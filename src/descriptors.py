# MolMod is a collection of molecular modelling tools for python.
# Copyright (C) 2005 Toon Verstraelen
#
# This file is part of MolMod.
#
# MolMod is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 2
# of the License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
#
# --


from scipy.special import lpmn as lpml, gammaln

import numpy

import math


__all__ = ["MolecularDescriptorTV1"]


def sph_harmonics(L, phi, theta):
    """Returns all spherical harmonics up to order L for theta and phi
       in the form [
            numpy.array([y^0_0]),
            numpy.array([y^(-1)_1, y^0_1, y^1_1]),
            ...,
            numpy.array([y^(-N)_N, ..., y^N_N])
       ].
    """
    x = math.cos(theta)
    val, foo = lpml(L,L,x)
    val = val.astype(complex)
    for l in xrange(0, L+1):
        val[:,l] *= math.sqrt((2*l+1)/4.0/math.pi)
        for m in xrange(0, l+1):
            val[m,l] *= math.exp(0.5*(gammaln(l-m+1)-gammaln(l+m+1)))
            val[m,l] *= numpy.exp(1j*m*phi)
    result = []
    for l in xrange(0, L+1):
        signs = ((numpy.array(range(l)) + l + 1)%2)*2-1
        result.append(numpy.concatenate((
            val[l:0:-1,l].conjugate()*signs,
            val[:l+1,l]
        )))
    return result


class FaseModulator(object):
    def __init__(self, order):
        self.factor = 0.5*math.pi*order
        self.scale = 1/self.factor

    def __call__(self, x):
        return numpy.exp(1j*self.factor*x)/self.scale


def total_mass(values):
    return sum(values["_masses"])
total_mass.label = "tm"
total_mass.internal = True
total_mass.shape = lambda parameters: (1,)
total_mass.tolerance = 1000.0


def mass_center(values):
    result = sum(
        coordinate * mass
        for coordinate, mass
        in zip(values["_coordinates"], values["_masses"])
    )/values["tm"]
    values["_deltas"] = values["_coordinates"] - result
    return result
mass_center.label = "mc"
mass_center.internal = False
mass_center.shape = lambda parameters: (3,)
mass_center.tolerance = 100.0


def inertia_tensor(values):
    tensor = numpy.zeros((3,3), float)
    for delta, mass in zip(values["_deltas"], values["_masses"]):
        tensor += mass*(
            numpy.dot(delta, delta)*numpy.identity(3, float)
           -numpy.outer(delta, delta)
        )
    return tensor
inertia_tensor.label = "it"
inertia_tensor.internal = False
inertia_tensor.shape = lambda parameters: (9,)
inertia_tensor.tolerance = 10.0


def rms_radius(values):
    distances = numpy.array([
        numpy.linalg.norm(delta)
        for delta
        in (values["_deltas"])
    ], float)
    values["_distances"] = distances
    epsilon = values["_distance_epsilon"]
    directions = numpy.array([
        delta/(distance + math.exp(-distance/epsilon))
        for delta, distance
        in zip(values["_deltas"], distances)
    ], float)
    values["_directions"] = directions
    result = math.sqrt(sum(distances**2*values["_masses"])/values["tm"])
    values["_reduced_distances"] = distances/result
    return result
rms_radius.label = "rmsr"
rms_radius.internal = True
rms_radius.shape = lambda parameters: (1,)
rms_radius.tolerance = 0.1


def spherical_harmonics_descriptor_internal(values):
    angles = numpy.array([
        [
            math.atan2(direction[1], direction[0]), # phi
            math.acos(direction[2]), # theta
        ] for direction
        in values["_directions"]
    ])
    values["_angles"] = angles
    L = values["_L"]
    result = []
    for modulator in values["_modulators"]:
        harmonics_sum = [
            numpy.zeros(2*l+1, complex)
            for l in xrange(L+1)
        ]
        for (phi, theta), reduced_distance, mass in zip(angles, values["_reduced_distances"], values["_masses"]):
            tmp = sph_harmonics(L, phi, theta)
            mod = modulator(reduced_distance) * mass
            for row, sum_row in zip(tmp, harmonics_sum):
                sum_row[:] += mod * row
        for sum_row in harmonics_sum:
            result.append(sum(abs(sum_row/values["tm"])**2))
            #result.append(sum(abs(sum_row)**2))
    result = numpy.array(result)
    result.shape = (-1,L+1)
    return result
spherical_harmonics_descriptor_internal.label = "shdi"
spherical_harmonics_descriptor_internal.internal = True
spherical_harmonics_descriptor_internal.shape = lambda parameters: (len(parameters["_modulators"]),parameters["_L"]+1)
spherical_harmonics_descriptor_internal.tolerance = 0.1


def spherical_harmonics_descriptor_external(values):
    L = values["_L"]
    result = []
    harmonics_sum = [
        numpy.zeros(2*l+1, complex)
        for l in xrange(L+1)
    ]
    for (phi, theta), reduced_distance, mass in zip(values["_angles"], values["_reduced_distances"], values["_masses"]):
        tmp = sph_harmonics(L, phi, theta)
        for row, sum_row in zip(tmp, harmonics_sum):
            sum_row[:] += row * reduced_distance * mass
    for sum_row in harmonics_sum:
        result.extend(abs(sum_row/values["tm"])**2)
    return numpy.array(result)
spherical_harmonics_descriptor_external.label = "shde"
spherical_harmonics_descriptor_external.internal = False
spherical_harmonics_descriptor_external.shape = lambda parameters: ((parameters["_L"]+1)**2,)
spherical_harmonics_descriptor_external.tolerance = 0.1


def static_fields(ordered_descriptors, parameters, internal, label, format):
    adresses = {}
    begin = 0
    for descriptor in ordered_descriptors:
        if descriptor.internal == internal:
            end = begin + reduce(lambda x,y: x*y, descriptor.shape(parameters), 1)
            adresses[descriptor.label] = (begin, end)
            end = begin
    size = end

    labels = []
    for descriptor in ordered_descriptors:
        if descriptor.internal == internal:
            size = reduce(lambda x,y: x*y, descriptor.shape(parameters), 1)
            l = len(str(size))
            labels.extend(["%s_%0*i" % (descriptor.label, l, index) for index in xrange(size)])

    table_name = "%s_%s" % (format, label)

    sql_create_descriptor_table = """
    CREATE TABLE %s (
        id BIGINT UNSIGNED NOT NULL AUTO_INCREMENT PRIMARY KEY, %s
    )
    """ % (
        table_name, ", ".join(
            "%s DOUBLE" % label
            for label in labels
        )
    )

    tolerance = []
    for descriptor in ordered_descriptors:
        if descriptor.internal == internal:
            size = reduce(lambda x,y: x*y, descriptor.shape(parameters), 1)
            tolerance.extend([descriptor.tolerance]*size)
    tolerance = numpy.array(tolerance, float)

    return labels, size, labels, table_name, sql_create_descriptor_table, tolerance


class MolecularDescriptorTV1(object):
    format = "mdtv1"

    ordered_descriptors = [
        total_mass,
        mass_center,
        inertia_tensor,
        rms_radius,
        spherical_harmonics_descriptor_internal,
        spherical_harmonics_descriptor_external,
    ]

    descriptors = dict(
        (descriptor.label, descriptor)
        for descriptor
        in ordered_descriptors
    )

    parameters = {
        "_distance_epsilon": 1e-6,
        "_L": 3,
        "_modulators": [FaseModulator(order) for order in xrange(1,6)]
    }


    (
        internal_adresses,
        internal_size,
        internal_labels,
        internal_table_name,
        sql_create_internal_descriptor_table,
        internal_tolerance
    ) = static_fields(ordered_descriptors, parameters, True, "internal", format)
    (
        external_adresses,
        external_size,
        external_labels,
        external_table_name,
        sql_create_external_descriptor_table,
        external_tolerance
    ) = static_fields(ordered_descriptors, parameters, False, "external", format)


    def __init__(self, coordinates, masses):
        assert len(coordinates.shape) == 2
        assert len(masses.shape) == 1
        assert coordinates.shape[1] == 3
        assert len(coordinates) > 0
        assert len(coordinates) == len(masses)

        self.values = self.parameters.copy()
        self.values["_coordinates"] = coordinates
        self.values["_masses"] = masses

        for descriptor in self.ordered_descriptors:
            self.values[descriptor.label] = descriptor(self.values)

        # only depends on internal coordinates
        self.internal_fingerprint = numpy.concatenate([
            numpy.array([self.values[descriptor.label]]).ravel()
            for descriptor
            in self.ordered_descriptors
            if descriptor.internal
        ])

        # also depends on external coordinates
        self.external_fingerprint = numpy.concatenate([
            numpy.array([self.values[descriptor.label]]).ravel()
            for descriptor
            in self.ordered_descriptors
            if not descriptor.internal
        ])

        for key in self.values.keys():
            if key[0] == '_':
                del self.values[key]

        self.derive_properties()

    def from_fingerprints(cls, internal_fingerprint, external_fingerprint):
        assert len(internal_fingerprint) == cls.internal_size
        assert len(external_fingerprint) == cls.external_size
        class _EmptyClass(object): pass
        result = _EmptyClass()
        result.__class__ = cls
        result.internal_fingerprint = internal_fingerprint
        result.external_fingerprint = external_fingerprint
        result.values = {}
        for label, (begin, end) in cls.internal_adresses.iteritems():
            value = internal_fingerprint[begin,end]
            value.shape = cls.descriptors[label].shape(cls.parameters)
            result.values[label] = value
        for label, (begin, end) in cls.external_adresses.iteritems():
            value = external_fingerprint[begin,end]
            value.shape = cls.descriptors[label].shape(cls.parameters)
            result.values[label] = value
        return result
    from_fingerprints = classmethod(from_fingerprints)

    def derive_properties(self):
        shdi = self.values["shdi"]
        ratio = sum(shdi[:,1::2].ravel()) / sum(shdi[:,0::2].ravel())
        self.inversion_symmetric = ratio < 1e-6

    def compare_structure(self, other):
        shdi1 = self.values["shdi"].ravel()
        shdi2 = other.values["shdi"].ravel()
        ratio = sum((shdi1 - shdi2)**2)  / (0.5*sum((shdi1 + shdi2)**2))
        return ratio < 1e-10

    def compare_global_rotation(self, other):
        shde1 = self.values["shde"]
        shde2 = other.values["shde"]
        ratio = sum((shde1 - shde2)**2) / (0.5*sum((shde1 + shde2)**2))
        return ratio < 1e-10

    def compare_global_translation(self, other):
        mc1 = self.values["mc"].ravel()
        mc2 = other.values["mc"].ravel()
        ratio = sum((mc1 - mc2)**2) / (0.5*sum((mc1 + mc2)**2))
        return ratio < 1e-10

    def look_up_internal(self, connection):
        sql = "SELECT id FROM %s WHERE (%s)" % (
            self.internal_table_name, " AND ".join(
                "%s < %f AND %s > %f" % (label, reference + tolerance, label, reference - tolerance)
                for label, reference, tolerance
                in zip(self.internal_labels, self.internal_fingerprint, self.internal_tolerance)
            )
        )
        cursor = connection.cursor()
        cursos.execute(sql)
        result = cursor.fetchall()
        cursor.close()

    def look_up_external(self, connection):
        sql = "SELECT id FROM %s WHERE (%s)" % (
            self.external_table_name, " AND ".join(
                "%s < %f AND %s > %f" % (label, reference + tolerance, label, reference - tolerance)
                for label, reference, tolerance
                in zip(self.external_labels, self.external_fingerprint, self.external_tolerance)
            )
        )
        cursor = connection.cursor()
        cursos.execute(sql)
        result = cursor.fetchall()
        cursor.close()

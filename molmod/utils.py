# -*- coding: utf-8 -*-
# MolMod is a collection of molecular modelling tools for python.
# Copyright (C) 2007 - 2012 Toon Verstraelen <Toon.Verstraelen@UGent.be>, Center
# for Molecular Modeling (CMM), Ghent University, Ghent, Belgium; all rights
# reserved unless otherwise stated.
#
# This file is part of MolMod.
#
# MolMod is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 3
# of the License, or (at your option) any later version.
#
# MolMod is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, see <http://www.gnu.org/licenses/>
#
# --
# -*- coding: utf-8 -*-
# MolMod is a collection of molecular modelling tools for python.
# Copyright (C) 2007 - 2010 Toon Verstraelen <Toon.Verstraelen@UGent.be>, Center
# for Molecular Modeling (CMM), Ghent University, Ghent, Belgium; all rights
# reserved unless otherwise stated.
#
# This file is part of MolMod.
#
# MolMod is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 3
# of the License, or (at your option) any later version.
#
# MolMod is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, see <http://www.gnu.org/licenses/>
#
# -- \
"""Utilities that are used in all parts of the MolMod library"""


from builtins import range
import numpy as np
from future.utils import with_metaclass


__all__ = ["cached", "ReadOnlyAttribute", "ReadOnly", "compute_rmsd"]


class cached(object):
    """A decorator that will turn a method into a caching descriptor

       When an attribute is requested for the first time, the original method
       will be called and its return value is cached. Subsequent access to the
       attribute will just return the cached value.

       Usage::

             class Foo(object):
                 @cached
                 def some_property(self):
                     return self.x*self.y

       There are a few limitations on the ``cached`` decorator. Once the result
       is computed and cached, it can not be erased. This means that the values
       on which the result depends have to be read-only parameters that can not
       be changed afterwards. This is facilitated by deriving from the
       :class:`ReadOnly` object. See :class:`molmod.molecules.Molecule` for an
       example.
    """
    def __init__(self, fn):
        self.fn = fn
        self.attribute_name = "_cache_%s" % fn.__name__
        fn_doc_lines = fn.__doc__.split("\n")
        self.__doc__ = "*Cached attribute:* %s.\n" % fn_doc_lines[0] + \
            "\n".join(fn_doc_lines[1:])

    def __get__(self, instance, cls=None):
        # make sure that the class attribute is simply this cached object:
        if instance is None:
            return self
        value = getattr(instance, self.attribute_name, self)
        if value is self:
            #print "COMPUTING %s" % self.attribute_name
            value = self.fn(instance)
            setattr(instance, self.attribute_name, value)
            #Disabled because not strictly enforcable and incompatible with Cython memoryviews.
            #if isinstance(value, np.ndarray):
            #    value.setflags(write=False)
        return value


class ReadOnlyAttribute(object):
    """A descriptor that becomes read-only after the first assignment.

       The initial value of the attribute is None. As long is this initial value
       is not overwritten, one can assign an immutable value or a numpy
       array. As soon as some value is assigned, no changes can be made.
       If the value is an array (an important exception to the restriction of
       mutable types), it is copied and set read-only. Warning: it is strictly
       speaking not guaranteed that the contents of a read-only numpy array are
       really fixed. We rely on your common sense to not write crappy code that
       changes read-only numpy arrays.

       Read-only attributes should be used in classes derived from the ReadOnly
       class.
    """

    def __init__(self, ptype=None, none=True, check=None, npdim=None,
                 npshape=None, npdtype=None, doc=None):
        """
           One can impose detailed type checking on the attribute through the
           following options.

           Optional arguments:
            | ``ptype`` -- The expected (Python) type of the attribute.
            | ``none`` -- When False, it is not possible to assign None.
            | ``npdim`` -- In case of numpy arrays: the number of dimensions.
            | ``npshape`` -- In case of numpy arrays: the expected shape.
            | ``npdtype`` -- In case of numpy arrays: the expected dtype.
            | ``check`` -- A method to check the validity of the attribute.
            | ``doc`` -- Short description.

           The npshape option must be a tuple with integer or None values. In
           case of a None value, the corresponding shape in the attribute is
           not tested. Trailing None values from the npshape argument can be
           omitted.
        """
        self.ptype = ptype
        self.none = none
        self.check = check
        if doc is None:
            doc = "no documentation available"
        self.doc = doc
        # process the type checks
        if ptype is not None:
            if issubclass(ptype, list):
                raise TypeError("ptype can not be a subclass of list because "
                    "lists are mutable.")
            if issubclass(ptype, np.ndarray):
                if not (npdim is None or isinstance(npdim, int)):
                    raise TypeError("npdim must be None or and integer.")
                if not (npshape is None or (isinstance(npshape, tuple) and
                        all(isinstance(s, int) or s is None for s
                        in npshape))):
                    raise TypeError("npshape must be None or a tuple of "
                        "integer and/or None values.")
                self.npdim = npdim
                self.npshape = npshape
                self.npdtype = npdtype
        elif not (npdim is None and npshape is None and npdtype is None):
            raise ValueError("The arguments npdim, npshape and npdtype are "
                "only allowed when ptype is a subclass of np.ndarray.")
        # make a nice docstring
        self.__doc__ = "*Read-only attribute:* %s.\n\n" % doc
        check_lines = []
        if not self.none:
            check_lines.append("* May not be None.")
        if self.ptype is not None:
            check_lines.append("* Must be an instance of ``%s``" % self.ptype)
            if issubclass(ptype, np.ndarray):
                if self.npdim is not None:
                    check_lines.append("* Must have dimension %i." % self.npdim)
                if self.npshape is not None:
                    tmp = []
                    for s in self.npshape:
                        if s is None:
                            tmp.append("?")
                        else:
                            tmp.append(str(s))
                    shape_str = "(%s)" % (", ".join(tmp))
                    check_lines.append("* Must have shape %s." % shape_str)
                if self.npdtype is not None:
                    check_lines.append("* Must have dtype ``%s``." % self.npdtype)
        if check is not None:
            check_lines.append("* Special conditions: %s." % check.__doc__)
        if len(check_lines) > 0:
            self.__doc__ += "The attribute must satisfy the following conditions:\n\n" + "\n\n".join(check_lines)
        # construct an annoying name
        self.attribute_name = "_read_only_%i" % id(self)

    def __get__(self, instance, cls=None):
        # make sure that the class attribute is simply this cached object:
        if instance is None:
            return self
        # just get the value, or return None if not present
        result = getattr(instance, self.attribute_name, None)
        if isinstance(result, np.ndarray):
            result = result.view()
        return result

    def __set__(self, instance, value, do_check=True):
        if value is None and not self.none:
            raise TypeError("This attribute may not be assigned None.")
        if self.ptype is not None and issubclass(self.ptype, tuple) and value is not None:
            # convert to a tuple of it is supposed one.
            value = tuple(value)
        if not (self.ptype is None or value is None):
            # Only type check if there is one and the value is not None.
            if issubclass(self.ptype, np.ndarray):
                if hasattr(value, "__len__"):
                    # try to turn non-arrays into arrays.
                    value = np.array(value, dtype=self.npdtype, copy=False)
                else:
                    raise TypeError("Single values are not automatically "
                        "converted into arrays.")
                if self.npdim is not None and len(value.shape) != self.npdim:
                    raise TypeError("Value does not have the right dimension. "
                        "Got %i. Expected %i" % (len(value.shape), self.npdim))
                if self.npshape is not None:
                    for i in range(len(self.npshape)):
                        s = self.npshape[i]
                        if s is not None:
                            if len(value.shape) < i+1:
                                raise TypeError("The array must have at least "
                                    "dimension %i." % (i+1))
                            if s != value.shape[i]:
                                raise TypeError("Element %i of the shape is "
                                    "not allowed. Got %i. Expected %i." %
                                    (i, value.shape[i], s))
                if self.npdtype is not None:
                    if not np.issubdtype(value.dtype, self.npdtype):
                        raise ValueError("The dtype must be a subdtype of %s. "
                            "Got %s." % (self.npdtype, value.dtype))
            if not isinstance(value, self.ptype):
                raise TypeError("Value (%s) does not have the expected type "
                    "(%s)." % (value, self.ptype))
        #Disabled because not strictly enforcable and incompatible with Cython memoryviews.
        #if isinstance(value, np.ndarray):
            #if value.flags.writeable:
            #    value = value.copy()
            #    value.setflags(write=False)
        #else:
        if not isinstance(value, np.ndarray):
            # Call the hash function. If that does not work, the object is
            # mutable, which we don't like in case the object is not a numpy
            # array.
            hash(value)
        if do_check:
            self.check_wrapper(instance, value)
        setattr(instance, self.attribute_name, value)

    def check_wrapper(self, instance, value):
        if not (self.check is None or value is None):
            # Only check if there is one and the value is not None.
            self.check(instance, value)


class ReadOnlyType(type):
    """A meta class for ReadOnly classes

       This meta class makes sure ReadOnlyAttribute descriptors are inherited.
    """

    def __init__(cls, name, bases, dct):
        super(ReadOnlyType, cls).__init__(name, bases, dct)
        for base in bases:
            for key, descriptor in base.__dict__.items():
                if isinstance(descriptor, ReadOnlyAttribute):
                    setattr(cls, key, descriptor)


class ReadOnly(with_metaclass(ReadOnlyType, object)):
    """A base class for read-only objects

       An object that has nothing but read-only attributes. If an attribute is
       not assigned yet, it is writable. Some attributes are cached, i.e. they
       are computed from other attributes upon request.

       If you want to modify a ReadOnly object, just create a modified one from
       scratch. This is greatly facilitated by the method :meth:`copy_with`.
    """

    def __copy__(self):
        return self

    def __deepcopy__(self, memo):
        return self

    def __getstate__(self):
        """Part of the pickle protocol"""
        result = {}
        for key, descriptor in self.__class__.__dict__.items():
            if isinstance(descriptor, ReadOnlyAttribute):
                result[key] = descriptor.__get__(self)
        return result

    def __setstate__(self, state):
        """Part of the pickle protocol"""
        for key, val in state.items():
            descriptor = self.__class__.__dict__.get(key)
            if not isinstance(descriptor, ReadOnlyAttribute):
                # Got wrong class attribute during unpickling. Just ignore.
                continue
            # Do not call custom check routines before all attributes
            # are assigned.
            descriptor.__set__(self, val, do_check=False)
        # Now do the custon checks
        for key, val in state.items():
            descriptor = self.__class__.__dict__.get(key)
            if not isinstance(descriptor, ReadOnlyAttribute):
                # Got wrong class attribute during unpickling. Just ignore.
                continue
            descriptor.check_wrapper(self, val)

    def copy_with(self, **kwargs):
        """Return a copy with (a few) changed attributes

           The keyword arguments are the attributes to be replaced by new
           values. All other attributes are copied (or referenced) from the
           original object. This only works if the constructor takes all
           (read-only) attributes as arguments.
        """
        attrs = {}
        for key, descriptor in self.__class__.__dict__.items():
            if isinstance(descriptor, ReadOnlyAttribute):
                attrs[key] = descriptor.__get__(self)
        for key in kwargs:
            if key not in attrs:
                raise TypeError("Unknown attribute: %s" % key)
        attrs.update(kwargs)
        return self.__class__(**attrs)


def compute_rmsd(a, b):
    """Compute the root-mean-square deviation between two arrays

       Arguments:
         a, b  --  Two numpy arrays with the same shape
    """
    return np.sqrt(((a-b)**2).mean())

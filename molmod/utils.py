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


import numpy, types


__all__ = ["cached", "ReadOnly", "compute_rmsd"]


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
        self.__doc__ = fn.__doc__

    def __get__(self, instance, owner):
        # make sure that the class attribute is simply this cached object:
        if instance is None:
            return self
        value = getattr(instance, self.attribute_name, self)
        if value is self:
            #print "COMPUTING %s" % self.attribute_name
            value = self.fn(instance)
            setattr(instance, self.attribute_name, value)
            if isinstance(value, numpy.ndarray):
                value.setflags(write=False)
        return value


class ReadOnly(object):
    """A base class for read-only objects

       An object that has mainly read-only attributes. If an attribute is not
       assigned yet, it is writable. Some attributes are cached, i.e. they are
       computed from other attributes upon request.

       If you want to modify a ReadOnly object, just create a modified one from
       scratch. This is greatly facilitated by the :meth:`copy_with`.

       The constructor of a derived class must call the :meth:`init_attributes`
       method to tell the ReadOnly object which attributes must be read-only.
       Some of them can be mandatory, and some can be optional.
    """

    def __init__(self):
        object.__setattr__(self, "__hidden", {})

    def __copy__(self):
        return self

    def __deepcopy__(self, memo):
        return self

    def init_attributes(self, mandatory, optional):
        """Prepare read-only attributes

           This method is called in the __init__ routines of derived classes to
           define read only attributes. Both arguments are dictionaries with
           (name, value) pairs. In case of mandatory, the value is not allowed to
           be None.

           Example::

               class Foo(ReadOnly):
                   def __init__(self, numbers, coordinates=None, title=None):
                       ReadOnly.__init__(self)
                       mandatory = {"numbers": numbers}
                       optional = {"coordinates": coordinates, "title": title}
                       self.init_attributes(mandatory, optional)

           Note that ``init_attributes`` could have been impelemented as a
           constructor of the ``ReadOnly`` class, but this would not allow
           second generation derived classes to introduce new read-only
           attributes.
        """
        for name, value in optional.iteritems():
            self._register_attribute(name)
            setattr(self, name, value)
        for name, value in mandatory.iteritems():
            self._register_attribute(name)
            if value is None:
                raise ValueError("'%s' is a mandatory argument, can not be None" % name)
            setattr(self, name, value)

    def _register_attribute(self, name):
        """Register an attribute name as a read-only attribute"""
        hidden = getattr(self, "__hidden")
        if name in hidden:
            raise ValueError("The name '%s' is already registered." % name)
        hidden[name] = None

    def __getattr__(self, name):
        value = getattr(self, "__hidden").get(name)
        #value = self.__hidden.get(name)
        if value is None:
            return object.__getattribute__(self, name)
            #raise AttributeError("Attribute '%s' is not defined for '%s'." % (name, self))
        return value

    def __setattr__(self, name, value):
        hidden = getattr(self, "__hidden")
        if name in hidden:
            current = hidden.get(name)
            if current is None:
                if isinstance(value, numpy.ndarray):
                    value.setflags(write=False)
                hidden[name] = value
            else:
                raise AttributeError("Attribute '%s' is read-only." % name)
        else:
            object.__setattr__(self, name, value)

    def __setstate__(self, state):
        """Part of the pickle protocol"""
        for key, val in state.iteritems():
            object.__setattr__(self, key, val)


    def copy_with(self, **kwargs):
        """Return a copy with (a few) changed attributes

           The keyword arguments are the attributes to be replaced by new
           values. All other attributes are copied (or referenced) from the
           original object.
        """
        attrs = getattr(self, "__hidden").copy()
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
    return numpy.sqrt(((a-b)**2).mean())

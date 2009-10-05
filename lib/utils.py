# MolMod is a collection of molecular modelling tools for python.
# Copyright (C) 2007 - 2008 Toon Verstraelen <Toon.Verstraelen@UGent.be>
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
"""Utilities that are used in all parts of the MolMod library"""


import numpy


__all__ = ["cached", "ReadOnly", "rmsd"]


class cached(object):
    """A decorator that will turn a method into a caching descriptor

    When an attribute is requested for the first time, the original method will
    be called and its return value is cached. Subsequent access to the attribute
    will just return the cached value.
    """
    def __init__(self, fn):
        self.fn = fn
        self.attribute_name = "_cache_%s" % fn.__name__
        self.__doc__ = fn.__doc__

    def __get__(self, instance, owner):
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
       scratch.
    """

    def __init__(self):
        """Initialize a read-only object"""
        object.__setattr__(self, "__hidden", {})

    def __copy__(self):
        return self

    def __deepcopy__(self, memo):
        return self

    def _init_attributes(self, mandatory, optional):
        """Prepare read-only attributes

           This method is called in the __init__ routines of derived classes to
           define read only attributes. Both arguments are dictionaries with
           (name, value) pairs. In case of mandatory, the value is not allowed to
           be None.

           Example:
           def __init__(self, numbers, coordinates=None, title=None):
               ReadOnly.__init__(self)
               mandatory = {"numbers": numbers}
               optional = {"coordinates": coordinates, "title": title}
               self._init_attributes(mandatory, optional)
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

    def copy_with(self, **kwargs):
        """Return a copy with (a few) changed attributes"""
        attrs = getattr(self, "__hidden").copy()
        for key in kwargs:
            if key not in attrs:
                raise TypeError("Unknown attribute: %s" % key)
        attrs.update(kwargs)
        return self.__class__(**attrs)



def rmsd(a, b):
    """Compute the root-mean-square deviation between two arrays

       Arguments:
         a, b  --  Two numpy arrays with the same shape
    """
    return numpy.sqrt(((a-b)**2).mean())


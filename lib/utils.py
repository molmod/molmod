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


import numpy


__all__ = ["cached", "cached_writable", "ReadOnly", "rmsd"]


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


class cached_writable(cached):
    """The writable variant of the normal caching descriptor

       This can be useful when constructing an object of which some results
       are at available without further computation. See transformations.py for
       an example.
    """
    def __set__(self, instance, value):
        value = getattr(instance, self.attribute_name, self)
        if value is not self:
            raise AttributeError("Value is already set.")
        else:
            setattr(instance, self.attribute_name, value)


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
                raise AttributeError("Attribute '%s' is read-only.")
        else:
            object.__setattr__(self, name, value)


def rmsd(a,b):
    return numpy.sqrt(((a-b)**2).mean())


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


from __future__ import print_function, division

import sys
import os
import platform
import datetime
import getpass
import time
import codecs
import locale
import functools
from contextlib import contextmanager

from molmod.units import kjmol, kcalmol, electronvolt, angstrom, nanometer, \
    femtosecond, picosecond, amu, deg, gram, centimeter


__all__ = ['ScreenLog', 'TimerGroup']


class Unit(object):
    def __init__(self, kind, conversion, notation, format):
        self.kind = kind
        self.conversion = conversion
        self.notation = notation
        self.format = format

    def __call__(self, value):
        return self.format % (value/self.conversion)


class UnitSystem(object):
    def __init__(self, *units):
        self.units = units
        # check for duplicates
        for i0, unit0 in enumerate(self.units):
            for unit1 in self.units[:i0]:
                if unit0.kind == unit1.kind:
                    raise ValueError('The unit of \'%s\' is encountered twice.' % unit0.kind)

    def log_info(self, log):
        if log.do_low:
            with log.section('UNITS'):
                log('The following units will be used below:')
                log.hline()
                log('Kind          Conversion               Format Notation')
                log.hline()
                for unit in self.units:
                    log('%13s %21.15e %9s %s' % (unit.kind, unit.conversion, unit.format, unit.notation))
                log.hline()
                log('The internal data is divided by the corresponding conversion factor before it gets printed on screen.')

    def apply(self, some):
        for unit in self.units:
            some.__dict__[unit.kind] = unit


class ScreenLog(object):
    # log levels
    silent = 0
    warning = 1
    low = 2
    medium = 3
    high = 4
    debug = 5

    # screen parameters
    margin = 8
    width = 72

    # unit systems
    # TODO: the formats may need some tuning
    joule = UnitSystem(
        Unit('energy', kjmol, 'kJ/mol', '%10.1f'),
        Unit('temperature', 1, 'K', '%10.1f'),
        Unit('length', angstrom, 'A', '%10.4f'),
        Unit('invlength', 1/angstrom, 'A^-1', '%10.5f'),
        Unit('area', angstrom**2, 'A^2', '%10.3f'),
        Unit('volume', angstrom**3, 'A^3', '%10.3f'),
        Unit('time', femtosecond, 'fs', '%10.1f'),
        Unit('mass', amu, 'amu', '%10.5f'),
        Unit('charge', 1, 'e', '%10.5f'),
        Unit('force', kjmol/angstrom, 'kJ/mol/A', '%10.1f'),
        Unit('forceconst', kjmol/angstrom**2, 'kJ/mol/A**2', '%10.1f'),
        Unit('velocity', angstrom/femtosecond, 'A/fs', '%10.5f'),
        Unit('acceleration', angstrom/femtosecond**2, 'A/fs**2', '%10.5f'),
        Unit('angle', deg, 'deg', '%10.5f'),
        Unit('c6', 1, 'E_h*a_0**6', '%10.5f'),
        Unit('c8', 1, 'E_h*a_0**6', '%10.5f'),
        Unit('c12', 1, 'E_h*a_0**12', '%10.5f'),
        Unit('diffconst', angstrom**2/picosecond, 'A**2/ps', '%10.5f'),
        Unit('density', gram/centimeter**3, 'g/cm^3', '%10.3f'),
    )
    cal = UnitSystem(
        Unit('energy', kcalmol, 'kcal/mol', '%10.2f'),
        Unit('temperature', 1, 'K', '%10.1f'),
        Unit('length', angstrom, 'A', '%10.4f'),
        Unit('invlength', 1/angstrom, 'A^-1', '%10.5f'),
        Unit('area', angstrom**2, 'A^2', '%10.3f'),
        Unit('volume', angstrom**3, 'A^3', '%10.3f'),
        Unit('time', femtosecond, 'fs', '%10.1f'),
        Unit('mass', amu, 'amu', '%10.5f'),
        Unit('charge', 1, 'e', '%10.5f'),
        Unit('force', kcalmol/angstrom, 'kcal/mol/A', '%10.1f'),
        Unit('forceconst', kcalmol/angstrom**2, 'kcal/mol/A**2', '%10.1f'),
        Unit('velocity', angstrom/femtosecond, 'A/fs', '%10.5f'),
        Unit('acceleration', angstrom/femtosecond**2, 'A/fs**2', '%10.5f'),
        Unit('angle', deg, 'deg', '%10.5f'),
        Unit('c6', 1, 'E_h*a_0**6', '%10.5f'),
        Unit('c8', 1, 'E_h*a_0**6', '%10.5f'),
        Unit('c12', 1, 'E_h*a_0**12', '%10.5f'),
        Unit('diffconst', angstrom**2/femtosecond, 'A**2/fs', '%10.5f'),
        Unit('density', gram/centimeter**3, 'g/cm^3', '%10.3f'),
    )
    solid = UnitSystem(
        Unit('energy', electronvolt, 'eV', '%10.4f'),
        Unit('temperature', 1, 'K', '%10.1f'),
        Unit('length', angstrom, 'A', '%10.4f'),
        Unit('invlength', 1/angstrom, 'A^-1', '%10.5f'),
        Unit('area', angstrom**2, 'A^2', '%10.3f'),
        Unit('volume', angstrom**3, 'A^3', '%10.3f'),
        Unit('time', femtosecond, 'fs', '%10.1f'),
        Unit('mass', amu, 'amu', '%10.5f'),
        Unit('charge', 1, 'e', '%10.5f'),
        Unit('force', electronvolt/angstrom, 'eV/A', '%10.1f'),
        Unit('forceconst', electronvolt/angstrom**2, 'eV/A**2', '%10.1f'),
        Unit('velocity', angstrom/femtosecond, 'A/fs', '%10.5f'),
        Unit('acceleration', angstrom/femtosecond**2, 'A/fs**2', '%10.5f'),
        Unit('angle', deg, 'deg', '%10.5f'),
        Unit('c6', 1, 'E_h*a_0**6', '%10.5f'),
        Unit('c8', 1, 'E_h*a_0**6', '%10.5f'),
        Unit('c12', 1, 'E_h*a_0**12', '%10.5f'),
        Unit('diffconst', angstrom**2/femtosecond, 'A**2/fs', '%10.5f'),
        Unit('density', gram/centimeter**3, 'g/cm^3', '%10.3f'),
    )
    bio = UnitSystem(
        Unit('energy', kcalmol, 'kcal/mol', '%10.2f'),
        Unit('temperature', 1, 'K', '%10.1f'),
        Unit('length', nanometer, 'nm', '%10.6f'),
        Unit('area', nanometer**2, 'nm^2', '%10.4f'),
        Unit('volume', nanometer**3, 'nanometer^3', '%10.1f'),
        Unit('invlength', 1/nanometer, 'nm^-1', '%10.8f'),
        Unit('time', picosecond, 'ps', '%10.4f'),
        Unit('mass', amu, 'amu', '%10.5f'),
        Unit('charge', 1, 'e', '%10.5f'),
        Unit('force', kcalmol/angstrom, 'kcal/mol/A', '%10.5f'),
        Unit('forceconst', kcalmol/angstrom**2, 'kcal/mol/A**2', '%10.5f'),
        Unit('velocity', angstrom/picosecond, 'A/ps', '%10.5f'),
        Unit('acceleration', angstrom/picosecond**2, 'A/ps**2', '%10.5f'),
        Unit('angle', deg, 'deg', '%10.5f'),
        Unit('c6', 1, 'E_h*a_0**6', '%10.5f'),
        Unit('c8', 1, 'E_h*a_0**6', '%10.5f'),
        Unit('c12', 1, 'E_h*a_0**12', '%10.5f'),
        Unit('diffconst', nanometer**2/picosecond, 'nm**2/ps', '%10.2f'),
        Unit('density', gram/centimeter**3, 'g/cm^3', '%10.3f'),
    )
    atomic = UnitSystem(
        Unit('energy', 1, 'E_h', '%10.6f'),
        Unit('temperature', 1, 'K', '%10.1f'),
        Unit('length', 1, 'a_0', '%10.5f'),
        Unit('invlength', 1, 'a_0^-1', '%10.5f'),
        Unit('area', 1, 'a_0^2', '%10.3f'),
        Unit('volume', 1, 'a_0^3', '%10.3f'),
        Unit('time', 1, 'aut', '%10.1f'),
        Unit('mass', 1, 'aum', '%10.1f'),
        Unit('charge', 1, 'e', '%10.5f'),
        Unit('force', 1, 'E_h/a_0', '%10.5f'),
        Unit('forceconst', 1, 'E_h/a_0**2', '%10.5f'),
        Unit('velocity', 1, 'a_0/aut', '%10.5f'),
        Unit('acceleration', 1, 'a_0/aut**2', '%10.5f'),
        Unit('angle', 1, 'rad', '%10.7f'),
        Unit('c6', 1, 'E_h*a_0**6', '%10.5f'),
        Unit('c8', 1, 'E_h*a_0**6', '%10.5f'),
        Unit('c12', 1, 'E_h*a_0**12', '%10.5f'),
        Unit('diffconst', 1, 'a_0**2/aut', '%10.2f'),
        Unit('density', 1, 'aum/a_0^3', '%10.3f'),
    )


    def __init__(self, name, version, head_banner, foot_banner, timer, f=None):
        self.name = name
        self.version = version
        self.head_banner = head_banner
        self.foot_banner = foot_banner
        self.timer = timer

        self._active = False
        self._level = self.medium
        self.unitsys = self.joule
        self.unitsys.apply(self)
        self.prefix = ' '*(self.margin-1)
        self._last_used_prefix = None
        self.stack = []
        self.add_newline = False
        if f is None:
            _file = sys.stdout
        else:
            _file = f
        self.set_file(_file)

    do_warning = property(lambda self: self._level >= self.warning)
    do_low = property(lambda self: self._level >= self.low)
    do_medium = property(lambda self: self._level >= self.medium)
    do_high = property(lambda self: self._level >= self.high)
    do_debug = property(lambda self: self._level >= self.debug)

    def _pass_color_code(self, code):
        if self._file.isatty():
            return code
        else:
            return ''

    reset =     property(functools.partial(_pass_color_code, code="\033[0m"))
    bold =      property(functools.partial(_pass_color_code, code="\033[01m"))
    teal =      property(functools.partial(_pass_color_code, code="\033[36;06m"))
    turquoise = property(functools.partial(_pass_color_code, code="\033[36;01m"))
    fuchsia =   property(functools.partial(_pass_color_code, code="\033[35;01m"))
    purple =    property(functools.partial(_pass_color_code, code="\033[35;06m"))
    blue =      property(functools.partial(_pass_color_code, code="\033[34;01m"))
    darkblue =  property(functools.partial(_pass_color_code, code="\033[34;06m"))
    green =     property(functools.partial(_pass_color_code, code="\033[32;01m"))
    darkgreen = property(functools.partial(_pass_color_code, code="\033[32;06m"))
    yellow =    property(functools.partial(_pass_color_code, code="\033[33;01m"))
    brown =     property(functools.partial(_pass_color_code, code="\033[33;06m"))
    red =       property(functools.partial(_pass_color_code, code="\033[31;01m"))

    def set_file(self, f):
        # Wrap sys.stdout into a StreamWriter to allow writing unicode.
        self._file = f
        self.add_newline = False

    def set_level(self, level):
        if level < self.silent or level > self.debug:
            raise ValueError('The level must be one of the ScreenLog attributes.')
        self._level = level

    def __call__(self, *words):
        s = u' '.join(str(w) for w in words)
        if not self.do_warning:
            raise RuntimeError('The runlevel should be at least warning when logging.')
        if not self._active:
            prefix = self.prefix
            self.print_header()
            self.prefix = prefix
        if self.add_newline and self.prefix != self._last_used_prefix:
            self._file.write(u'\n')
            self.add_newline = False
        # Check for alignment code '&'
        pos = s.find(u'&')
        if pos == -1:
            lead = u''
            rest = s
        else:
            lead = s[:pos] + ' '
            rest = s[pos+1:]
        width = self.width - len(lead)
        if width < self.width//2:
            raise ValueError('The lead may not exceed half the width of the terminal.')
        # break and print the line
        first = True
        while len(rest) > 0:
            if len(rest) > width:
                pos = rest.rfind(' ', 0, width)
                if pos == -1:
                    current = rest[:width]
                    rest = rest[width:]
                else:
                    current = rest[:pos]
                    rest = rest[pos:].lstrip()
            else:
                current = rest
                rest = u''
            self._file.write(u'%s %s%s\n' % (self.prefix, lead, current))
            if first:
                lead = u' '*len(lead)
                first = False
        self._last_used_prefix = self.prefix

    def warn(self, *words):
        self(u'WARNING!!&'+u' '.join(str(w) for w in words))

    def hline(self, char='~'):
        self(char*self.width)

    def center(self, *words, **kwargs):
        if len(kwargs) == 0:
            edge = ''
        elif len(kwargs) == 1:
            if 'edge' not in kwargs:
                raise TypeError('Only one keyword argument is allowed, that is edge')
            edge = kwargs['edge']
        else:
            raise TypeError('Too many keyword arguments. Should be at most one.')
        s = u' '.join(str(w) for w in words)
        if len(s) + 2*len(edge) > self.width:
            raise ValueError('Line too long. center method does not support wrapping.')
        self('%s%s%s' % (edge, s.center(self.width-2*len(edge)), edge))

    def blank(self):
        self._file.write(u'\n')

    def _enter(self, prefix):
        if len(prefix) > self.margin-1:
            raise ValueError('The prefix must be at most %s characters wide.' % (self.margin-1))
        self.stack.append(self.prefix)
        self.prefix = prefix.upper().rjust(self.margin-1, ' ')
        self.add_newline = True

    def _exit(self):
        self.prefix = self.stack.pop(-1)
        if self._active:
            self.add_newline = True

    @contextmanager
    def section(self, prefix):
        self._enter(prefix)
        try:
            yield
        finally:
            self._exit()

    def set_unitsys(self, unitsys):
        self.unitsys = unitsys
        self.unitsys.apply(self)
        if self._active:
            self.unitsys.log_info()

    def print_header(self):
        # Suppress any logging as soon as an exception is not caught.
        def excepthook_wrapper(type, value, traceback):
            self.set_level(self.silent)
            sys.__excepthook__(type, value, traceback)
        sys.excepthook = excepthook_wrapper

        if self.do_warning and not self._active:
            self._active = True
            print(self.head_banner, file=self._file)
            self._print_basic_info()
            self.unitsys.log_info(self)

    def print_footer(self):
        if self.do_warning and self._active:
            self._print_basic_info()
            self.timer._stop('Total')
            self.timer.report(self)
            print(self.foot_banner, file=self._file)

    def _print_basic_info(self):
        if self.do_low:
            with self.section('ENV'):
                self('User:          &' + getpass.getuser())
                self('Platform:      &' + platform.platform())
                self('Time:          &' + datetime.datetime.now().isoformat())
                self('Python version:&' + sys.version.replace('\n', ''))
                self('%s&%s' % (('%s version:' % self.name).ljust(15), self.version))
                self('Current Dir:   &' + os.getcwd())
                self('Command line:  &' + ' '.join(sys.argv))


class Timer(object):
    def __init__(self):
        self.cpu = 0.0
        self._start = None

    def start(self):
        assert self._start is None
        self._start = time.clock()

    def stop(self):
        assert self._start is not None
        self.cpu += time.clock() - self._start
        self._start = None


class SubTimer(object):
    def __init__(self, label):
        self.label = label
        self.total = Timer()
        self.own = Timer()

    def start(self):
        self.total.start()
        self.own.start()

    def start_sub(self):
        self.own.stop()

    def stop_sub(self):
        self.own.start()

    def stop(self):
        self.own.stop()
        self.total.stop()


class TimerGroup(object):
    def __init__(self):
        self.parts = {}
        self._stack = []
        self._start('Total')

    def reset(self):
        for timer in self.parts.values():
            timer.total.cpu = 0.0
            timer.own.cpu = 0.0

    @contextmanager
    def section(self, label):
        self._start(label)
        try:
            yield
        finally:
            self._stop(label)

    def _start(self, label):
        # get the right timer object
        timer = self.parts.get(label)
        if timer is None:
            timer = SubTimer(label)
            self.parts[label] = timer
        # start timing
        timer.start()
        if len(self._stack) > 0:
            self._stack[-1].start_sub()
        # put it on the stack
        self._stack.append(timer)

    def _stop(self, label):
        timer = self._stack.pop(-1)
        assert timer.label == label
        timer.stop()
        if len(self._stack) > 0:
            self._stack[-1].stop_sub()

    def get_max_own_cpu(self):
        result = None
        for part in self.parts.values():
            if result is None or result < part.own.cpu:
                result = part.own.cpu
        return result

    def report(self, log):
        max_own_cpu = self.get_max_own_cpu()
        #if max_own_cpu == 0.0:
        #    return
        with log.section('TIMER'):
            log('Overview of CPU time usage.')
            log.hline()
            log('Label             Total      Own')
            log.hline()
            bar_width = log.width-33
            for label, timer in sorted(self.parts.items()):
                #if timer.total.cpu == 0.0:
                #    continue
                if max_own_cpu > 0:
                    cpu_bar = "W"*int(timer.own.cpu/max_own_cpu*bar_width)
                else:
                    cpu_bar = ""
                log('%14s %8.1f %8.1f %s' % (
                    label.ljust(14),
                    timer.total.cpu, timer.own.cpu, cpu_bar.ljust(bar_width),
                ))
            log.hline()

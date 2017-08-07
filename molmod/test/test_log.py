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


from io import StringIO

from molmod import *


def test_line_wrapping():
    with StringIO() as f:
        timer = TimerGroup()
        log = ScreenLog('TEST', '0.0', '', '', timer, f)
        log._active = True
        with log.section('NVE'):
            log('This is just a long test message that should get splitted into two lines properly.')
        assert f.getvalue() == '\n    NVE This is just a long test message that should get splitted into two\n    NVE lines properly.\n'


def test_center():
    with StringIO() as f:
        timer = TimerGroup()
        log = ScreenLog('TEST', '0.0', '', '', timer, f)
        log._active = True
        with log.section('NVE'):
            log.center('Gets centered.', edge='***')
            log.center('Gets centered.')
        assert f.getvalue() == '\n    NVE ***                          Gets centered.                          ***\n    NVE                              Gets centered.                             \n'


def test_levels():
    timer = TimerGroup()
    log = ScreenLog('TEST', '0.0', '', '', timer)
    assert log.do_medium
    assert not log.do_high
    log.set_level(log.low)
    assert not log.do_medium


def test_header_footer():
    with StringIO() as f:
        timer = TimerGroup()
        log = ScreenLog(u'TEST', u'0.0', u'', u'', timer, f)
        log(u'test')
        assert f.tell() > 10
        log.print_footer()


def test_colors():
    timer = TimerGroup()
    log = ScreenLog('TEST', '0.0', '', '', timer)
    assert log.red == ''


def test_hline():
    with StringIO() as f:
        timer = TimerGroup()
        log = ScreenLog('TEST', '0.0', '', '', timer, f)
        log._active = True
        with log.section('FOOBAR'):
            log.hline()
        assert f.getvalue() == '\n FOOBAR ' + '~'*log.width + '\n'


def test_unitsys():
    from molmod.units import kjmol, kcalmol
    with StringIO() as f:
        timer = TimerGroup()
        log = ScreenLog('TEST', '0.0', '', '', timer, f)
        assert abs(log.energy.conversion - kjmol) < 1e-10
        log.set_unitsys(log.cal)
        assert abs(log.energy.conversion - kcalmol) < 1e-10


def test_lead():
    with StringIO() as f:
        timer = TimerGroup()
        log = ScreenLog('TEST', '0.0', '', '', timer, f)
        log._active = True
        with log.section('AAA'):
            log('Some prefix:&followed by a long text that needs to be wrapped over multiple lines.')
        assert f.getvalue() == '\n    AAA Some prefix: followed by a long text that needs to be wrapped over\n    AAA              multiple lines.\n'


def test_enter_leave():
    with StringIO() as f:
        timer = TimerGroup()
        log = ScreenLog('TEST', '0.0', '', '', timer, f)
        with log.section('FOO'):
            assert log.prefix == '    FOO'
            with log.section('BAR'):
                assert log.prefix == '    BAR'
            assert log.prefix == '    FOO'
        assert log.prefix == '       '

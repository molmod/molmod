# PyChem is a general chemistry oriented python package.
# Copyright (C) 2005 Toon Verstraelen
# 
# This file is part of PyChem.
# 
# PyChem is free software; you can redistribute it and/or
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

BEGIN {
    print "{";
    gradient=0;
    coordinate=0;
    converged=0;
    input_molecule=0;
    complete=0;
}

/Value of the MolecularEnergy:/ {
    printf("\t'energy': %s,\n", $5);
}

(gradient) {
    if (NF==0) {
        printf("\t],\n");
        gradient=0;
        converged=1;
    } else {
        if (coordinate==0) { printf("\t\t["); }
        printf("%s, ", $2);
        if (coordinate==2) { printf("],\n"); }
        coordinate = (coordinate + 1) % 3;
    }
}

/Gradient of the MolecularEnergy:/ {
    printf("\t'gradient': [\n");
    gradient=1;
}

(input_molecule) {
    if (NF==1 && $1="}") {
        input_molecule = 0;
        printf("\t],\n");
    } else {
        printf("\t\t[%s, %s, %s,\n", $4, $5, $6);
    }
}

(converged && ($4=="geometry")) {
    input_molecule = 1;
    printf("\t'input_coordinates': [\n");
}

/End Time/ {
    complete=1;
}

END {
    if (complete==1) {
        printf("\t'job_completed': True,\n");
    } else {
        printf("\t'job_completed': False,\n");
    }
    print "}";
}

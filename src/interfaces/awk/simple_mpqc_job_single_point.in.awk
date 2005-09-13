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
    printf("{\n");
    molecule=0;
    input_xyz=""
    numatoms = 0;
}

/method:/ {
    method = $0
    sub("method: ", "", method)
    printf("\t\"method\": \"%s\",\n", method);
}

/basis:/ {
    printf("\t\"basis\": \"%s\",\n", $2);
}

/charge:/ {
    printf("\t\"charge\": %s,\n", $2);
}

/multiplicity:/ {
    printf("\t\"multiplicity\": %s,\n", $2);
}

/gradient:/ {
    if ($2=="yes") {
        printf("\t\"do_gradient\": True,\n");
    } else {
        printf("\t\"do_gradient\": False,\n");
    }
}

(molecule==1) {
    if (match($1, ":")>0) {
        printf("\t\"input_xyz\": \"\"\"%i%s\"\"\",\n", numatoms, input_xyz);
        molecule=0
    } else {
        if (NF > 0) {
            numatoms += 1
        }
        input_xyz = sprintf("%s\n%s", input_xyz, $0);
    }
}

/molecule:/ {
    molecule=1;
}

END {
    printf("}\n");
}

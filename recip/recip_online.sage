"""
RECIP -- REpository of Complex multIPlication SageMath code.
See the file README.txt for version information and instructions.

#*****************************************************************************
# Copyright (C) 2010 -- 2020
# Marco Streng <marco.streng@gmail.com>
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License along
# with this program; if not, write to the Free Software Foundation, Inc.,
# 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
#*****************************************************************************

"""

load_attach_mode(load_debug=True)

abspath = "https://raw.githubusercontent.com/mstreng/recip/master/recip/"
#abspath = "https://bitbucket.org/mstreng/recip/raw/master/"

files = [   "basic.sage",
            "polynomials.sage",
            "symplectic_matrices.sage",
            "period_matrices.sage",
            "theta_trans.sage",
            "cm_types.sage",
            "find_class_invariant.sage",
            "igusa_invariants.sage",
            "list_fields.sage",
            "denominators.sage",
            "class_polynomials.sage",
            "orders.sage",
            "bissonstreng.sage",
            "bls.sage"]

for f in files:
    load(abspath + f)

recip_verbose = 0

def set_recip_verbose(level):
    r"""
    Sets the verbosity level for the Recip package.
    r"""
    global recip_verbose
    recip_verbose = level

def get_recip_verbose():
    r"""
    Returns the verbosity level for the Recip package.
    r"""
    try:
        return recip_verbose
    except NameError:
        recip_verbose = 0
        return 0
    

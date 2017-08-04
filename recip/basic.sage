"""
RECIP -- REpository of Complex multIPlication SageMath code.
See the file README.txt for version information and instructions.

#*****************************************************************************
# Copyright (C) 2010, 2011, 2012, 2013 Marco Streng <marco.streng@gmail.com>
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

To use this package, start Sage with the .sage files from this package in your
working directory. Then type::

    sage: from recip import *

This file implements some basic functions: base change for matrices and
uniformizers of ideals.

"""

from sage.matrix.constructor import Matrix


def mat_convert(M, ring_or_map):
    """
    Applies `ring_or_map` to the coefficients of M,
    i.e. given M = (m_ij)_ij, returns
    (ring_or_map(m_ij))_ij
    
    EXAMPLES::
    
        sage: from recip import *
        sage: M = Matrix([[1/2, 1/5], [7, 9/2]])
        sage: mat_convert(M, floor)
        [0 0]
        [7 4]
        sage: mat_convert(M, GF(3))
        [2 2]
        [1 0]
    """
    return Matrix([[ring_or_map(M[i,j]) for j in range(M.ncols())] \
                    for i in range(M.nrows())])


def uniformizer(p):
    """
    Given a prime ideal p, returns a uniformizer.

    EXAMPLES::

        sage: from recip import *
        sage: P = QuadraticField(-5,'a').ideal(5).factor()[0][0]
        sage: uniformizer(P)
        a
    """
    for pi in p.gens():
        if pi.valuation(p) == 1:
            return pi


def lift_small(a):
    """
    Given a in Z/muZ, output b in ZZ with |b| <= mu/2 and (b mod mu) = a.
    Given a in Z, output a.
    
    EXAMPLES::
    
        sage: from recip import *
        sage: lift_small(Zmod(8)(7))
        -1
        sage: lift_small(Zmod(7)(10))
        3
        sage: lift_small(Zmod(12)(20))
        -4
        sage: lift_small(-501)
        -501
        sage: lift_small(QQ(7))     
        7
    """
    if a.parent() is ZZ or (a.parent() is QQ and a in ZZ):
        return a
    mu = a.parent().order()
    if not mu in ZZ:
        raise ValueError, "a (=%s) must be in ZZ or ZZ/mu*ZZ for some mu, " \
                          "but is in %s of order %s" % (a, a.parent(), mu)
    b = ZZ(a) % mu
    if b > mu/2:
        b -= mu
    return b


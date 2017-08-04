"""
RECIP -- REpository of Complex multIPlication SageMath code.
See the file README.txt for version information and instructions.

#*****************************************************************************
# Copyright (C) 2010,2011,2012,2013 Marco Streng <marco.streng@gmail.com>
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

This file contains functions for listing quartic CM-fields.

"""

def iterate_DAB(stopD, stopA, startD=5, startA=1, prime_D=False):
    """
    Given a range of positive discriminants D and integers A, returns all
    minimal triples DAB corresponding to primitive quartic CM-fields with D and
    A in the range.
    
    The number of such triples is proportional to
    (maxA^2 - minA^2) * (sqrt(maxD) - sqrt(minD)).
    
    If prime_D is true, restrict to D of the form p for a prime p=1 mod 4
    and 4*p for p a prime p!=1 mod 4.
    
    EXAMPLES::
    
        sage: from recip import *
        sage: list(iterate_DAB(10, 10))
        [[5, 5, 5], [5, 7, 11], [5, 9, 19], [8, 4, 2], [8, 6, 7], [8, 8, 14]]
        sage: list(iterate_DAB(20, 20, prime_D=True))
        [[5, 5, 5], [5, 7, 11], [5, 9, 19], [5, 10, 20], [5, 11, 29], [5, 12, 31], [5, 13, 41], [5, 14, 44], [5, 15, 55], [5, 15, 45], [5, 16, 59], [5, 17, 71], [5, 17, 61], [5, 18, 76], [5, 19, 89], [5, 19, 79], [8, 4, 2], [8, 6, 7], [8, 8, 14], [8, 10, 23], [8, 10, 17], [8, 12, 34], [8, 12, 18], [8, 14, 47], [8, 14, 41], [8, 14, 31], [8, 16, 62], [8, 16, 46], [8, 18, 79], [8, 18, 73], [8, 18, 63], [12, 6, 6], [12, 8, 13], [12, 10, 22], [12, 10, 13], [12, 12, 33], [12, 14, 46], [12, 14, 37], [12, 14, 22], [12, 16, 61], [12, 16, 37], [12, 18, 78], [12, 18, 69], [12, 18, 33], [13, 5, 3], [13, 9, 17], [13, 10, 12], [13, 12, 23], [13, 13, 39], [13, 13, 13], [13, 15, 53], [13, 16, 51], [13, 17, 69], [13, 17, 43], [13, 18, 68], [13, 18, 29], [13, 19, 87], [13, 19, 61], [17, 5, 2], [17, 11, 26], [17, 12, 19], [17, 13, 38], [17, 15, 52], [17, 15, 18], [17, 16, 47], [17, 17, 68], [17, 17, 34], [17, 19, 86]]
        sage: l = list(iterate_DAB(30, 30))
        sage: l = [(DAB, CM_Field(DAB)) for DAB in l]
        sage: l = [(DAB, K, K.ideal(2).factor()) for (DAB, K) in l]
        sage: l = [(DAB, K) for (DAB, K, f) in l if len(f) == 1 and f[0][1] == 1]
        sage: l = [(DAB, K.class_number()/K.real_field().class_number()) for (DAB, K) in l]
        sage: [(DAB, h) for (DAB, h) in l if h > 3]
        [([21, 25, 109], 4), ([29, 25, 149], 5)]
    """
    # d is D//4
    for D in range(startD, stopD):
        if prime_D:
            if D % 4 == 1 and not ZZ(D).is_prime():
                continue
            elif D % 4 == 0 and not ZZ(D/4).is_prime():
                continue
        if not is_fundamental_discriminant(D):
            continue
        K0 = QuadraticField(D, 'a')
        for DAB in iterate_DAB_given_D(D, stopA, startA, K0):
            yield DAB


def iterate_DAB_given_D(D, stopA, startA=1, K0=None):
    """
    Given a real quadratic discriminant and a range of integers A,
    returns all minimal triples DAB corresponding to primitive quartic
    CM-fields with this D and in the A-range.
    
    The number of such triples is proportional to
    1/2 (maxA^2 - minA^2) / sqrt(D).
    
    K0, if given, should be the quadratic field with discriminant D.
    """
    # D is known, also m^2 D = A^2 - 4B,
    # so if D is even, then so is A
    if D % 2 == 0:
        e = 2
        if startA % 2 == 1:
            startA = startA + 1
    else:
        e = 1
    
    for A in range(startA, stopA, e):
        for DAB in iterate_DAB_given_DA(D, A, K0):
            yield DAB


def iterate_DAB_given_DA(D, A, K0):
    """
    Given a real quadratic discriminant and a positive integer A,
    returns all minimal triples DAB corresponding to primitive quartic
    CM-fields with this D and in the A-range.

    The number of such triples is proportional to A/sqrt(D).

    K0, if given, should be the quadratic field with discriminant D.
    """
    if K0 is None:
        K0 = QuadraticField(D,'a')
    # if D is even, then both D and A^2 are divisible by 4
    # if D is odd, then it is 1 mod 4, so m and A have the same parity
    # Take steps for m of length e, and start with start.
    if D % 2 == 0:
        e = 1
    else:
        e = 2
    start = A % 2
    if start == 0:
        start = e
    # Now m^2 D = A^2 - 4B < A^2, so m < A/sqrt(D), so m <= floor(A/sqrt(D))
    for m in range(start, A/sqrt(D)+1, e):
        DAB = DAm_to_DAB(D, A, m)
        if not DAB[2].is_square(): # equivalent to primitive CM-type
            min_dab = DAB_to_minimal(DAB, K0, m)
            if min_dab == DAB:
                yield DAB


def DAm_to_DAB(D, A, m):
    """
    Given D a fundamental discriminant, A a positive integer, and m
    a positive integer with m^2 D < A^2 and m^2 * D == A^2 mod 4,
    returns [D,A,B] with m^2*D = A^2 - 4B
    """
    return [D, A, ZZ((A**2-m**2*D)/4)]


def wamelen_dab_list():
    """
    Returns the list of triples DAB from Paul van Wamelen's list of CM curves.
    
    EXAMPLE::
    
        sage: from recip import *
        sage: lst = wamelen_dab_list(); lst
        [[5, 5, 5], [8, 4, 2], [13, 13, 13], [5, 10, 20], [5, 65, 845], [29, 29, 29], [5, 85, 1445], [37, 37, 333], [8, 20, 50], [13, 65, 325], [13, 26, 52], [53, 53, 53], [61, 61, 549]]
    
    We check that this list contains of 13 cyclic Galois CM-fields, 
    and that there are 19 curves corresponding to these 13
    fields.::
    
        sage: len(lst)
        13
        sage: flds = [CM_Field(DAB) for DAB in lst]
        sage: all([K.is_galois() for K in flds])
        True
        sage: sum([K.class_number()/K.real_field().class_number() for K in flds])
        19

    """
    # In van Wamelen's article, the fields are K=QQ(sqrt(-a+b*sqrt(d))),
    # so K is generated by a root of
    #   (x^2+a)^2-b^2*d = x^4 + 2*a*x^2 + a^2 - b^2*d.
    # The real quadratic subfield is QQ(sqrt(D)) for D=(2*a)^2 - 4*(a^2-b^2*d)
    #                                                 = 4*b^2*d
    # We start with the list of a,b,d's:
    lst = [[2,1,2], [13,2,13], [5,1,5], [65, 26, 5], [29, 2, 29], [85, 34, 5],
           [37, 6, 37], [10, 5, 2], [65, 10, 13], [13, 3, 13], [53, 2, 53],
           [61, 6, 61]]
    # And we convert that to DAB triples:
    lst = [[5,5,5]] + [[d if d%4 in [0,1] else 4*d, 2*a, a**2 - b**2*d] for [a,b,d] in lst]
    return [DAB_to_minimal(DAB) for DAB in lst]

def wamelen_curves():
    """
    Returns van Wamelen's CM curves
    (taken from https://www.math.lsu.edu/~wamelen/CMcurves.txt)
    in such a way that wamelen_curves()[i] corresponds to wamelen_dab_list()[i]
    
    """
    x = QQ['x'].gen()
    return [[x^5 - 1], [x^5 - 3*x^4 - 2*x^3 + 6*x^2 + 3*x - 1], [x^5 - 3*x^4 + 4*x^3 - 3*x^2 + 16/13*x - 11/52], [x^5 - 15/2*x^3 + 45/4*x - 11/2, x^5 - 10/11*x^4 - 430/121*x^3 + 380/121*x^2 - 100/121*x + 8/121], [x^5 + 32/11*x^4 - 1400/1573*x^3 + 152/1573*x^2 - 464/102245*x + 8/102245, x^5 - 2172/1271*x^4 + 23857440/21000733*x^3 - 7771776/21000733*x^2 + 80937216/1365047645*x - 5112832/1365047645], [x^5 - 8*x^4 + 22*x^3 - 25*x^2 + 373/29*x - 289/116], [x^5 - 870/71*x^4 + 454280/85697*x^3 - 73680/85697*x^2 + 89040/1456849*x - 2336/1456849, x^5 - 259000/302621*x^4 + 439853471360/1556850983897*x^3 - 68527385600/1556850983897*x^2 + 81861201920/26466466726249*x - 1820295168/26466466726249], [x^5 - 3496/121*x^4 - 7952768/131769*x^3 - 2415616/43923*x^2 - 1052672/43923*x - 557056/131769], [x^5 - 11540/23*x^4 - 64300/529*x^3 + 160/529*x^2 + 530/529*x - 8/529, x^5 - 26/23*x^4 - 14564/25921*x^3 + 14984/18515*x^2 - 1822/13225*x + 84/13225], [x^5 - 260/9*x^3 + 130/3*x^2 + 130/3*x + 1183/36, x^5 - 246/53*x^4 + 27784/2809*x^3 - 130608/14045*x^2 + 272592/70225*x - 42336/70225], [x^5 - 59/23*x^4 + 9332/4761*x^3 - 328/529*x^2 + 571/6877*x - 27/6877, x^5 - 9319/3013*x^4 - 24427743/18156338*x^3 + 131339171/18156338*x^2 + 36128207/36312676*x - 70399443/36312676], [x^5 + 19419/8381*x^4 + 195362586/70241161*x^3 + 121565502/70241161*x^2 + 36129645/70241161*x + 4057371/70241161], [x^5 + 8653/3075*x^4 - 333674198/85100625*x^3 + 73469002/28366875*x^2 - 5977848533/1730379375*x + 1730805683/1038227625]]

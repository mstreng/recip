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

this file:
denominators.sage

This file computes proven Igusa class polynomials.


"""


def class_polynomials(K, factor=False, prec=None, D=None,
                      verbose=False):
    """
    Returns the Igusa class polynomials of K
    (proven!).
    
    INPUT:
    
     - `K` -- a CM-field
     - `factor` -- whether the output consists of sets of irreducible factors
                   (currently not supported)
     - `prec` -- starting precision (will be increased until it is enough)
     - `D` -- known multiple of the denominator (if omitted, computes one)
     - `verbose` -- whether to print specific verbose output for this function.
     
    OUTPUT:
    
    Igusa class polynomials H_1, Hhat_2, Hhat3 of K with respect to the
    absolute Igusa invariants of the author's thesis.
    
    EXAMPLES:
    
    We compute proven class polynomials for some of the entries of
    van Wamelen's table. With these class polynomials, it is trivial to
    check correctness of the table entries. The full list, and more, is
    proven in [BouyerS].::
    
        sage: from recip import *
        sage: K = CM_Field([5,5,5])
        sage: class_polynomials(K) # long time: 6 seconds
        [x, 0, 0]

        sage: class_polynomials(CM_Field((x^2+2)^2-2), verbose=True) # long time: 10 seconds
        starting with precision 50
        output has 40.8562309773420 too few bits of precision
        increasing precision to 95
        output has 3.47594808128218 more bits of precision than needed
        Denominator: 1 out of 2^46
        [x + 7290, 437400, 2952450000]

        sage: class_polynomials(CM_Field((x^2+13)^2-2^2*13), verbose=True) # long time: 10 seconds
        starting with precision 50
        output has 8.62029143144967 too few bits of precision
        increasing precision to 63
        output has 4.68096300785080 more bits of precision than needed
        Denominator: 1 out of 2^14
        [x + 7840, 102400, -204800000]

        sage: class_polynomials(CM_Field((x^2+5)^2-5)) # long time: 22 seconds
        [x^2 - 183708000*x, 37826743837500/14641*x - 601817074425000000/14641, 1994141034144140625000/14641*x - 423741159843750000000/14641]

        sage: class_polynomials(CM_Field((x^2+65)^2-26^2*5)) # long time: 30 seconds
        [x^2 - 209024611260948000/1615441*x + 739519963620480000000/1615441, 1923741818956270886400000/2609649624481*x - 3630281807913088204800000000/2609649624481, 174715258193111891291126400000000000/2609649624481*x - 618134989860658268345733120000000000000/2609649624481]

        sage: bruinier_yang_applies(CM_Field((x^2+29)^2-4*29))
        True
        sage: class_polynomials(CM_Field((x^2+29)^2-2^2*29), verbose=True) # long time: 10 seconds
        starting with precision 50
        output has 35.4684649508012 too few bits of precision
        increasing precision to 90
        output has 4.50158705782075 more bits of precision than needed
        Denominator: 1 out of 2^14 * 5^4
        [x + 2589408, 131383296, -60466176000000]
        
        sage: bruinier_yang_applies(CM_Field((x^2+61)^2-6^2*61))
        True
        sage: class_polynomials(CM_Field((x^2+61)^2-6^2*61), verbose=True) # long time: 11 seconds
        starting with precision 50
        output has 54.2686150981159 too few bits of precision
        increasing precision to 109
        output has 4.58030311169402 more bits of precision than needed
        Denominator: 3^3 * 41^4 out of 2^14 * 3^8 * 5^4 * 41^4
        [x - 88833024/1681, -14055214415872/76295547, 9663676416000000/2825761]

    """
    if factor:
        raise NotImplementedError, "factor=True not yet implemented"
    Zs = list(K.period_matrices_iter())
    h = igusa_modular_forms()
    if D is None:
        D = denominator_bound(K)
    elif type(D) is str:
        D = denominator_bound(K, bound=D)
    elif not D in ZZ:
        raise ValueError
    else:
        D = ZZ(D)
        
    # compute the period matrices and denominators here
    if prec is None:
        # prec = ZZ(floor(40 + log(D)/log(2)))
        # or maybe some more sensible guess?
        prec = 50
        if verbose:
            print "starting with precision %s" % prec
            sys.stdout.flush()

    while True:
        hom_values = [[f(Z, prec=prec, interval=True) for f in h] for Z in Zs] 
        abs_values = [igusa_modular_forms_to_absolute(Is) for Is in hom_values]
        abs_values = [[z[k] for z in abs_values] for k in range(3)]
        x = hom_values[0][0].parent()['x'].gen()
        pols = [prod([x-i1 for i1 in abs_values[0]])]
        pols = pols + [short_interpolation(abs_values[0], abs_values[k]) for k in [1,2]]
        pols = [[x.real() for x in (D*p).list()] for p in pols]
        # These polynomials are approximations of integer polynomials.
        err = max(flatten([[x.absolute_diameter() for x in p] for p in pols]))
        err = (log(err)/log(2)).n()
        if err < 0:
            if verbose:
                print "output has %s more bits of precision than needed" % -err
            # we found the polynomial!
            try:
                pols = [[c.unique_integer() for c in p] for p in pols]
            except ValueError as e:
                raise RuntimeError, "Incorrect denominator or bug, if you did "\
                                    "not specify a denominator, please report "\
                                    "%s" % e
            ret = [QQ['x'](p)/D for p in pols]
            if verbose:
                print "Denominator: %s out of %s" % \
                      (lcm([p.denominator() for p in ret]).factor(), D.factor())
            return ret
        else:
            if verbose:
                print "output has %s too few bits of precision" % err
            # intervals are still pretty large, so we increase precision
            prec = prec + ZZ(floor(err)) + 5
            if verbose:
                print "increasing precision to %s" % prec
                sys.stdout.flush()


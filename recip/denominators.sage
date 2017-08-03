"""
RECIP -- REpository of Complex multIPlication SageMath code.

this file:
denominators.sage

This file contains implementations of the denominator bounds of [BY], [GL], and
[LV]. The implementations of [LV] were written for [BouyerS] and contain
additional original work.

See the file README.txt for version information, instructions, and references.

#*****************************************************************************
# Copyright (C) 2010, 2011, 2012, 2013, 2014, 2015, 2016,2017 Marco Streng
#                                                  <marco.streng@gmail.com>
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


def denominator_bound(K, c=2**14, k=2, d=None, Phi=None, bound='default', check=True, safe=True):
    """
    If ``Phi`` is None, returns a positive integer `D`
    such that `D N(f h_{10}^{-k}(Z))`
    is integral whenever `Z` has CM by `O_K`,
    `c*f` is a modular
    form of weight `k` and level `1` with integral
    Fourier coefficients, and
    `d` (if K/QQ is Galois) and `2d` (if K/QQ is non-Galois)
    is the degree of the algebraic number 
    `f h_{10}^{-k}(Z)`.
    
    For the choice of invariants in the author's thesis, take k=2 and c=2^14.
    These are also the defaults.
        
    If ``Phi`` is a CM-type, returns a number `D`
    in the real quadratic subfield of the reflex field
    of ``Phi`` such that
    that `D f h_{10}^{-k}(Z)`
    is integral whenever `Z` has CM by `O_K`
    of type ``Phi`` and `c*f` is a modular
    form of weight `k` and level `1` with integral
    Fourier coefficients.
    
    EXAMPLE:
    
    The following example shows how horribly unusable the Goren-Lauter bounds
    are in practice. The additional number of digits (about 2000) needed to
    compute the class polynomial with certainty using the algorithm of the
    author's thesis is huge even though it is also easy to see that actually
    the class polynomials are all equal to X, so the algorithm as given is not
    practical. The importance of the algorithm is in proving the asymptotic
    running time::
    
        sage: from recip import *
        sage: gl = denominator_bound(CM_Field([5,5,5]), bound='gl')
        sage: gl.factor()
        2^14 * 5^36 * 19^20 * 29^18 * 59^14 * 79^14 * 89^14 * 109^12 * 139^12 * 149^12 * 179^12 * 199^12 * 229^12 * 239^12 * 269^10 * 349^10 * 359^10 * 379^10 * 389^10 * 409^10 * 419^10 * 439^10 * 449^10 * 479^10 * 499^10 * 509^10 * 569^10 * 599^10 * 619^10 * 659^10 * 709^10 * 719^10 * 739^10 * 769^10 * 809^10 * 829^10 * 839^10 * 859^10 * 919^10 * 929^10 * 1009^8 * 1019^8 * 1039^8 * 1049^8 * 1069^8 * 1109^8 * 1129^8 * 1229^8 * 1249^8 * 1259^8 * 1279^8 * 1289^8 * 1319^8 * 1399^8 * 1409^8 * 1429^8 * 1439^8 * 1459^8 * 1489^8 * 1499^8 * 1549^8 * 1559^8 * 1579^8 * 1609^8 * 1619^8 * 1669^8 * 1699^8 * 1709^8 * 1759^8 * 1789^8 * 1879^8 * 1889^8 * 1949^8 * 1979^8 * 1999^8
        sage: (log(gl)/log(10)).n()
        1947.25745024934
        sage: (log(gl)/log(2)).n()
        6468.64923196203
    
    Fortunately the Bruinier-Yang bounds apply::
    
        sage: denominator_bound(CM_Field([5,5,5]), bound='by').factor()
        2^14
    
    Next, we try a different field::
                        
        sage: K = CM_Field((x^2+2)^2-2,'a')
        sage: bruinier_yang_applies(K)
        False
        sage: gl = denominator_bound(K, bound='gl')
        sage: (log(gl)/log(2)).n()
        313.633570483095
        sage: Z = K.one_period_matrix(); Z
        Period Matrix
        [ 0.?e-15 + 1.08239220029240?*I 0.414213562373095? + 0.?e-14*I]
        [             0.41421356237310?            1.53073372946036?*I]
        sage: h = igusa_modular_forms()
        sage: v = [f(Z, prec=100, interval=True) for f in h]; v
        [5.14511662180777441738853320? + 0.?e-28*I, 1.5657740361761366552190052? + 0.?e-27*I, -0.00110508779417351004090108656? + 0.?e-30*I, 0.0201781755750066641900463261? + 0.?e-29*I]
        sage: w = igusa_modular_forms_to_absolute(v); w
        [-7290.0000000000000000000000? + 0.?e-23*I, 437400.00000000000000000000? + 0.?e-21*I, 2.9524500000000000000000000?e9 + 0.?e-17*I]
        sage: [(log((gl*y).real().absolute_diameter())/log(2)).n() for y in w]
        [240.728133100306, 246.473959365695, 259.180686194823]

    These numbers are very recognizable, but we need 346 more bits of output
    precision to prove they are correct.::
 
        sage: v = [f(Z, prec=450, interval=True) for f in h]
        sage: w = igusa_modular_forms_to_absolute(v)
        sage: x = [gl*y.real() for y in w]
        sage: [(log(y.absolute_diameter())/log(2)).n() for y in x]
        [-107.565143763258, -101.889231246203, -89.1476167961862]

    So now we have more than enough precision and can simply round::
    
        sage: z = [y.real().unique_integer() for y in x]
        sage: invs = [y/gl for y in z]; invs
        [-7290, 437400, 2952450000]
    
    These numbers are proven correct, so let's prove the corresponding entry
    of van Wamelen's table with it. First, we double-check that the curve
    in the table corresponds with the invariants given in the table::
    
        sage: x = QQ['x'].gen()
        sage: f = -x^5+3*x^4+2*x^3-6*x^2-3*x+1
        sage: [I2,I4,I6,I10] = igusa_clebsch_invariants(f)
        sage: (I2^5/I10) == 2^7*3^15
        True
        sage: I2^3*I4/I10 == 2^5*3^11*5
        True
        sage: I2^2*I6/I10 == 2^4*3^9*31
        True
    
    Next, we show that these invariants coincide with the value we just proved.
    This seems to be the first complete proof of correctness for any entry
    of van Wamelen's table, other than the classical and easy QQ(zeta_5).
    It was already proven by van Wamelen that the curve `y^2 = f` has CM
    by an order of `K`, and a sketch was given of how to extend these results
    to show the order is the maximal order. Now we have a complete proof (by
    different means).::
    
        sage: I6pr = (I2*I4-3*I6)/2
        sage: invs == [I4*I6pr/I10, I4^2*I2/I10, I4^5/I10^2]
        True

    Check that a bug with incorrect embeddings is fixed::
    
        sage: K = CM_Field([5,26,149])
        sage: denominator_bound(K, bound='by').factor()
        2^14 * 5^4 * 7^4

    Our implementation computes Lauter-Viray bounds faster than Bruinier-Yang
    bounds, but Bruinier-Yang bounds are significantly sharper, so we
    use them by default::
    
        sage: K = CM_Field([389,37,245])
        sage: D1 = denominator_bound(K, bound='lv') # less than a second
        sage: D2 = denominator_bound(K, bound='by') # long time: 6 seconds
        sage: class_polynomials(K, D=D1) # long time: 50 seconds
        [x^2 - 502951680/303601*x - 215457915617280/5768419,
         102098407935836160/11153001631321*x + 1497833452550013478502400/4026233588906881,
         2213545771008000000/92173567201*x - 36646804488714190848000000/33274657759561]
        sage: class_polynomials(K, D=D2) # still a long time, but only 26 seconds
        [x^2 - 502951680/303601*x - 215457915617280/5768419,
         102098407935836160/11153001631321*x + 1497833452550013478502400/4026233588906881,
         2213545771008000000/92173567201*x - 36646804488714190848000000/33274657759561]

    """
    if K.degree() != 4:
        raise ValueError, "Input field %s is not quartic." % K
    if len(K.subfields()) > 3:
        raise ValueError, "Input field %s is biquadratic." % K
    if d is None:
        d = K.class_number()
        if Phi is None and not K.is_galois():
            d = d * 2
    if not Phi is None:
        raise NotImplementedError
    if bound == 'default':
        if bruinier_yang_applies(K):
            bound = 'by'
            check = False # as the check has already been done
        else:
            try:
                return denominator_bound(K=K, c=c, k=k, d=d, Phi=Phi, bound='lv', safe=safe)
            except NotImplementedError:
                pass
            if get_recip_verbose():
                print "[LV] bounds not implemented, falling back to GL"
            bound = 'gl'
    if bound == 'lv':
        b = lauter_viray_bound(K, safe=safe)    
    elif bound == 'gl':
        b = goren_lauter_bound(K, Phi=Phi)**d
    elif bound == 'by':
        b = bruinier_yang_bound(K, check=check)
    else:
        raise ValueError, "Unkown bound: %s" % bound
    gal = K.is_galois()
    if not gal:
        c = sqrt(c)
    ret = c**d * b**k
    if ret in QQ:
        return QQ(ret)
    else:
        # There may be some cases where the output is the square root of an
        # integer. We take the largest integer "dividing" that square root:
        return QQ([p**ZZ(floor(e/2)) for (p,e) in QQ(ret**2).factor()])


def bruinier_yang_applies(K, proof=True, reason=False):
    """
    Return ``True`` if and only if the Bruinier-Yang results are proven for
    the input quartic CM field in Yang's preprint.
    
    With ``proof=False`` returns ``True`` only if the Bruinier-Yang results are
    formulated and conjectured for the input quartic CM field
    as in Yang's preprint.
    
    With ``reason=True`` also return a string with an explanation.
    
    EXAMPLES::
    
        sage: from recip import *
        sage: bruinier_yang_applies(CM_Field([5,5,5]))
        True
    """
    #1Assumption from page 2: D is a prime that is 1 mod 4
    # Notation: F = QQ(sqrt(D))) is real quadratic subfield of K
    #           Delta is such that K=F(sqrt(Delta))
    #           Delta' is complex conjugate of Delta
    #           Dtilde is relative norm of Delta
    #           Ktilde is reflex field
    #           Ftilde is real quadratic subfield of reflex field
    #2Assumption from page 2: Dtilde is the discriminant of Ftilde
    #           (this is presented as notation, but does a Delta such that this
    #            holds always exist? also, Delta appears later on, and
    #            if we let Dtilde be the discriminant, then
    #            Dtilde * m^2 = Delta*Delta' for a positive integer m,
    #            but m does not have to be square in F. Also, it is part of the
    #            assumption a few lines lower.)
    #3Assumption from page 3: O_K is free over O_F
    #4Assumption: O_K can actually be written as O_F + O_F((w+sqrt(Delta))/2)
    #           (Yang's phrasing suggeststhat this is not stronger than
    #            "O_K is free over O_F", still, this seems to really be an
    #            assumption. Moreover, this Delta is the same as the one in
    #            Dtilde = Delta*Delta')
    #5Assumption from Theorem 1.2: Dtilde = Delta*Delta' is 1 modulo 4 and prime
    #
    if K.degree() != 4:
        if reason:
            return False, "Bruinier-Yang formula not formulated for fields of degree different from 4."
        return False
    if len(K.subfields()) > 3:
        if reason:
            return False, "Bruinier-Yang formula not formulated for biquadratic fields."
        return False
    if K.minimal_DAB() == [5,5,5]:
        if reason:
            return True, "Bruinier-Yang's conjecture has been proven in this case."
        return True
        # Proof: an old version of the code found the correct kind of generator

    # Now we may assume that K is quartic non-biquadratic with only 2 roots of unity.

    delta0 = K.real_field().discriminant()
    if not is_prime(delta0):
        # 1 is not satisfied
        if reason:
            return False, "Bruinier-Yang formula not formulated as K0=QQ(sqrt(%s)) does not have prime discriminant" % delta0
        return False
    # now we know 1 is satisfied

    # Let A = {x-xbar : x in OK}, which generates the relative different diff.
    # If the hypothesis is satisfied, then 
    # A = sqrt(Delta)*OF, so A*OK=diff.
#    Krel = K.relative_field()
    diff = K.relative_different()
    if not diff.is_principal():
        if reason:
            return False, "Bruinier-Yang formula not formulated as OK is not monogenic (in fact, the relative different is non-principal)."
        return False
    eta = diff.gens_reduced()[0]
    eps = K.unit_group().gens()[1]
    assert eps.multiplicative_order() == oo
#    eps = Krel.structure()[1](eps) # move it to the relative field

    if K.is_totally_real_element(eta^2):
        pass
    elif K.is_totally_real_element((eps*eta)^2):
        eta = eta*eps
    else:
        # 4 is not satisfied
        if reason:
            return False, "Bruinier-Yang formula not formulated as OK is not monogenic of the right form over OF"
        return False
    Delta = K.self_to_real(eta^2)
    
    # If we reach this point, then we have found Delta in F with K=F(sqrt(Delta))
    # and Delta*Delta' is prime, hence is the discriminant of Ftilde.
    # So 1, 2, 5 are satisfied. If 4 is satisfied, then so is 3.
    # Next, considering the minimal polynomial of the algebraic number
    # (w+sqrt(Delta))/2,
    # we get
    # ((w+sqrt(Delta))/2)^2 = w^2/4 + w*sqrt(Delta)/2 + Delta/4
    # = w * (w+sqrt(Delta))/2 - w^2/4 + Delta/4, so 
    # Delta is w^2 modulo 4 if and only if the number is an algebraic integer.
    OF = K.real_field().maximal_order()
    for w in (4*OF).residues():
        if (4*OF).divides(Delta-w^2):
            # 4 (and hence 3) is satisfied
            break
    else: # this "else" belongs to the "for"
        # 4 is not satisfied. Indeed, suppose 4 is true.
        # Delta in the Sage code above is uniquely defined, up to squares of units
        # in K. So up to squares of units in K, it is the Delta in 4.
        # Now squares of units in K, if they are in F, are also squares of units
        # in F (Lemma II.3.3 of my thesis).
        # So Delta is correct up to squares of units.
        if reason:
            return False, "Bruinier-Yang formula not formulated as OK is not monogenic of the right form over OF"
    
    # Now we know 1, 3, 4 are satisfied.
    
    Dtilde = QuadraticField(K.DAB()[2]).discriminant()
    if Delta.norm() != Dtilde:
        # We know 1, 3, 4, so Delta generates the relative discriminant
        # and is totally negative, 
        # hence is defined up to totally positive units. Its norm is therefore
        # defined by 1, 3, 4. So we now know that there is no Delta
        # at all satisfying 1 -- 4.
         
        # If Delta.norm() is prime, then we never reach this line.
        assert not ZZ(Delta.norm()).is_prime()
        assert (Delta.norm() / Dtilde).sqrt() in QQ
        
        if reason:
            return False, "Bruinier-Yang formula not formulated, at least not if both Dtilde's are equal (Delta*Delta' and disc(Ftilde)) and both Delta's are equal (in the generator of OK/OF and in Delta*Delta'=Dtilde)"
        return False
    
    # Now we know 1 -- 4 are satisfied.
    
    if not proof:
        if reason:
            return True, "Bruinier-Yang's conjecture has been stated as a conjecture in this case."
        return True
    
    if not is_prime(ZZ(K.discriminant()/delta0**2)):
        if reason:
            return False, "Bruinier-Yang formula conjectured, but not proven as the relative discriminant (%s) is not prime" % ZZ(K.discriminant()/delta0**2)
        return False
        
    if reason:
        return True, "Bruinier-Yang's conjecture has been proven in this case."
    return True


def bruinier_yang_bound(K, check=True, proof=True):
    """
    Returns the Bruinier-Yang bounds. Assumes
    that the input satisfies the hypotheses of the
    Bruinier-Yang conjecture. Bounds are proven
    only if ``bruinier_yang_applies`` returns ``True``.
    
    EXAMPLES::
    
        sage: from recip import *
        sage: bruinier_yang_bound(CM_Field([5,5,5]))
        1
    
    Let's check for a few fields how sharp these denominators are. First,
    all elements of van Wamelen's table where the conjecture has been
    proven.::
    
        sage: lst = wamelen_dab_list()
        sage: lst = [(DAB, CM_Field(DAB)) for DAB in lst]
        sage: lst = [(DAB, K, bruinier_yang_applies(K, proof=True)) for (DAB, K) in lst if bruinier_yang_applies(K, proof=False)]
        sage: len(lst)
        6
        sage: [(DAB, b) for (DAB, K, b) in lst]
        [([5, 5, 5], True),
         ([13, 13, 13], True),
         ([29, 29, 29], True),
         ([37, 37, 333], True),
         ([53, 53, 53], True),
         ([61, 61, 549], True)]
        sage: lst = [(K, denominator_bound(K, bound='by'), class_polynomials(K)) for (DAB, K, b) in lst if b] # long time
        sage: len(lst) # long time
        6
        sage: [(K.DAB(), den.factor(), lcm([c.denominator() for c in flatten([list(p) for p in pols])]).factor()) for (K,den,pols) in lst] # long time
        [([5, 5, 5], 2^14, 1), ([13, 13, 13], 2^14, 1), ([29, 29, 29], 2^14 * 5^4, 1), ([37, 37, 333], 2^14 * 3^4 * 11^4, 11^2), ([53, 53, 53], 2^14 * 17^4 * 29^4, 17^2 * 29^4), ([61, 61, 549], 2^14 * 3^8 * 5^4 * 41^4, 3^3 * 41^4)]

    So we see that den contains too many 2's, 5's, 11's in the first cases, but
    den is optimal wrt 29 and 41 in the last two cases. The 2, 3, 5, 11, 17
    could be cancellation, but is still a bit fishy. Anyway, here is the list
    of 6 triples of class polynomials that we just proved::
    
        sage: [(K.DAB(), pols) for (K,den,pols) in lst] # long time
        [([5, 5, 5], [x, 0, 0]), ([13, 13, 13], [x + 7840, 102400, -204800000]), ([29, 29, 29], [x + 2589408, 131383296, -60466176000000]), ([37, 37, 333], [x + 505440, -810393600/121, 1642291200000]), ([53, 53, 53], [x + 16982274516480/14297, 30036725005359513600/204404209, -16680511371779820748800000/707281]), ([61, 61, 549], [x - 88833024/1681, -14055214415872/76295547, 9663676416000000/2825761])]

    Here is the final curve in the author's preprint with Florian Bouyer::
    
        sage: K = CM_Field([389, 37, 245])
        sage: class_polynomials(K, verbose=True)
        starting with precision 50
        output has 120.772235075431 too few bits of precision
        increasing precision to 175
        output has 4.22059275336305 more bits of precision than needed
        Denominator: 11^2 * 19^6 * 29^4 out of 2^14 * 5^4 * 11^8 * 19^8 * 29^4
        [x^2 - 502951680/303601*x - 215457915617280/5768419,
         102098407935836160/11153001631321*x + 1497833452550013478502400/4026233588906881,
         2213545771008000000/92173567201*x - 36646804488714190848000000/33274657759561]
    """
    if check:
        b, s = bruinier_yang_applies(K, proof=proof, reason=True)
        if not b:
            raise ValueError, s
        if get_recip_verbose():
            print s
    if get_recip_verbose():
        print "Computing Bruinier-Yang's bound"

    # We need to take the number of roots of unity into account.
    # The logarithm of what we have below is a factor #OK^*^tors / 2 too small.
    # This factor is 5 if K=QQ(zeta_5) and 1 otherwise.
    # We can ignore the 1, and in the case of QQ(zeta_5), our output is 1
    # anyway, so raising to the power 5 does not change that.
    # Therefore, we can safely ignore the factor W_K/2 in [Yang].
    
    # We start with some standard fields and discriminants.
    Phi = K.CM_types()[0]
    Kr = Phi.reflex_field()                 # denoted Ktilde in Yang
    Kr_rel = Kr.relative_field()            # denoted Ktilde in Yang
#    K0r = Phi.reflex_field().real_field()   # denoted Ftilde in Yang
    K0 = K.real_field()                     # denoted F in Yang
    D = K0.disc()                           #         D
    delta1 = K.disc()/D**2                  # this equals Dtilde by assumptions
    d = Kr_rel.relative_discriminant()      # denoted d_{Ktilde/Ftilde} in Yang
    K0r = Kr_rel.base_field()               # denoted Ftilde in Yang
    Dtilde = K0r.disc()                     #         Dtilde
    B = sqrt(K0r(Dtilde))                   # in Yang: sqrt(Dtilde) (gen. of F)
    O0r = K0r.maximal_order()
    Kr_rel_to_Kr = Kr_rel.structure()[0]
    def cu(u):
        """
        Given u, which is one of the numbers t in F in the sum in (1.3),
        return exp(sum_{P prime of Ftilde | p} B_t(P) log(p)).
        In other words, if we sum over S, we get exp(2*RHS of (1.11)).
        """
        if get_recip_verbose()>1:
            print "Computing contribution of\nt = %s\n  = d^-1 * (%s)\n  = (k+(D-n^2)/4*sqrt(Dtilde))/(2D),\nwhere k=%s and (D-n^2)/4=%s and n=%s" % (u, u*d.gens_reduced()[0], (t*D).trace(), (D*t*2-(D*t).trace())/B, sqrt(D-(D*t*2-(D*t).trace())/B*4))
        f = (u*d).factor()
        def splittype(p):
            """
            Given a prime P of Ftilde, returns
            0 if P splits in Ktilde/Ftilde
            1 if P is inert in Ktilde/Ftilde
            2 if P ramifies in Ktilde/Ftilde
            """
            t = Kr_rel_to_Kr(p).factor()
            if len(t) == 2:
                return 0
            return t[0][1]
            
        # For each prime P of Ftilde, take (P, val_P(u*d), splittype)
        # if val_P(u*d) is non-zero:
        f2 = [(p[0], p[1], splittype(p[0])) for p in f]
        
        if get_recip_verbose()>1:
            print "Primes over t*d: (given by [p,e,f] in Ftilde/QQ)"
            for p in f2:
                if p[2] == 0:
                    s = "  Split in Ktilde/Ftilde:    "
                elif p[1] == 1:
                    s = "  Inert in Ktilde/Ftilde:    "
                else:
                    s = "  Ramified in Ktilde/Ftilde: "
                print s + "[%s,%s,%s] to power %s" % (p[0].gens_two()[0], p[0].ramification_index(), p[0].residue_class_degree(), p[1])
        
        def prerho(k, t):
            """
            Given an integer k and the output of splittingtype t for a prime
            P of Ftilde, returns the number of ideals of Ktilde of relative
            norm P^k, i.e., returns rho(P^k) as below (1.4).
            """
            if k < 0:  # not an integral ideal
                return 0
            if k == 0: # the unit ideal
                return 1
            if t == 2: # ramified
                return 1
            if t == 0: # split
                return k+1
            # inert:
            if is_even(k):
                return 1
            return 0

        def rho(k):
            """
            Given the index k of an entry (P,v,t)=f2[k] of f2, returns
            the number of ideals of Ktilde of relative norm u*d*P^-1,
            i.e., returns rho(u*d*P^-1) as below (1.4).
            
            Actually, it returns this number only if P is not split. If P
            is split, then because of (1.4) the number 0 is returned instead.
            """
            if f2[k][2] == 0:
                return 0
            return prod((prerho(f2[l][1], f2[l][2]) for l in range(len(f2)) if k != l))*prerho(f2[k][1]-1,f2[k][2])

        # ord_P(t) = f2[k][1] - d.valuation(f2[k][0])
        ret = prod(
                   f2[k][0].norm() ** ((f2[k][1] - d.valuation(f2[k][0]) + 1) * (rho(k))) for k in range(len(f2))
                  )
        if get_recip_verbose() > 1:
            print "Contribution of t is %s" % ret
        return ret
        
    # Some notation changes within Yang:
    # n of page 5 is denoted n by us
    # n of page 3 is denoted k by us
    # m of page 5 is 1 in our case
    # m of pages 2 and 3 is denoted (D-n^2)/4 by us
    
    # a is the set of pairs (k,n) where n is an odd integer and k is an
    #                             integer with k^2 <= (D-n^2)^2 Dtilde /4    
    a = [[(k,n) for k in range(-floor((D-n^2)*sqrt(Dtilde)/4),floor((D-n^2)*sqrt(Dtilde)/4)+1)] for n in range(1,floor(sqrt(D))+1,2)]

    # b is the set of elements t = k+m*sqrt(Dtilde))/(2D) in F, where
    # (k,n) is in a, and m = (D-n^2)/4 is the subscript of "b".
    b = [(4*k+(D-n^2)*B)/(8*D) for (k,n) in flatten(a, max_level=1)]

    # Equation (1.3) says we should restrict to t in d_{Ktilde/Ftilde}^-1 = d.
    # So we define S to be the subset of b satisfying this condition.
    # The further condition |n|<m sqrt(Dtilde) in (1.3) is in our notation
    # |k| < (D-n^2)/4 * Dtilde, which is exactly the restriction we already
    # have.
    # So S is the set of t in equation (1.3), but where m also ranges
    # over all the numbers m = (D-n^2)/4 of page 5.
    S = [u for u in b if u in d^-1]
    
    #maxnrofls = max([len([elt for elt in S if elt[1] == n]) for n in range(1,floor(sqrt(D))+1,2)]+[0])
    # sizeofS = len(S)
    # maxsizeforoneeltofS = max([log((cu(u))+0.0) for u in S]+[0])
    ret = prod([cu(t) for t in S])
    # By Remark 7.1 of [GJLSVW], if K/QQ is Galois, then we have counted the
    # denominator twice, so:
    if K.is_galois():
        if ret.is_square():
            ret = ret.sqrt()
        else:
            print "K is Galois, but the Bruinier-Yang arithmetic intersection number is not the logarithm of a square. Interesting..."
            ret = ret.sqrt()
    return ret


def goren_lauter_bound(K, Phi=None):
    """
    If ``Phi`` is None, returns a positive integer `D`
    such that `c^l D^{dk} N(f h_{10}^{-k}(Z))`
    is integral whenever `Z` has CM by `O_K`,
    `f` is a modular
    form of weight `k` and level `1` with integral
    Fourier coefficients, and
    `d` (if K/QQ is Galois) and `2d` (if K/QQ is non-Galois)
    is the degree of the algebraic number 
    `f h_{10}^{-k}(Z)`.
    
    Here l is the degree and is l=d if K/QQ is Galois and l=2*d otherwise.
    
    For the choice of invariants in the author's thesis, take k=2 and note that
    2^14*i_n satisfies the hypotheses.
        
    If ``Phi`` is a CM-type, returns a number `D`
    in the real quadratic subfield of the reflex field
    of ``Phi`` such that
    that `D^k f h_{10}^{-k}(Z)`
    is integral whenever `Z` has CM by `O_K`
    of type ``Phi`` and `f` is a modular
    form of weight `k` and level `1` with integral
    Fourier coefficients.
    """
    if not Phi is None:
        raise NotImplementedError
    gal = K.is_galois()
    prime_powers = []
    [D,A,B] = DAB_to_minimal(K.DAB())
    a = A/2
    d = ZZ(D)
    if d%4 == 0:
        d = ZZ(d/4)
    b = sqrt((A**2-4*B)/d)/2
    k = ZZ(floor(min([a.valuation(2), b.valuation(2)])/2))
    a = ZZ(a/2**(2*k))
    b = ZZ(b/2**(2*k))
    primes = prime_range(4*d*a**2)
    prime_powers = []
    disc = K.discriminant()
    for p in primes:
        # this is the place to test whether the splitting of p is right
        if gal: # for some splitting types we can take 1 even if gal is False
            if not (disc % p):
                galfactor = 1
            elif len(K.ideal(p).factor()) == 2:
                galfactor = 1
            else:
                galfactor = 0
        else:
            if Phi is None:
                Kr = K.CM_types()[0].reflex_field()
            else:
                Kr = Phi.reflex_field()
            if len(K.ideal(p).factor()) == 2:
                f = Kr.ideal(p).factor()
                if len(f) == 2:
                    galfactor = 2
                elif len(f) == 3:
                    galfactor = 1
                else:
                    galfactor = 0
            else:
                galfactor = 0
            galfactor = 2
        base = p**galfactor # for the case Phi is not None, we can do much better than this
        
        if p == 2 and K.discriminant() % 2 == 0:
            prime_powers.append(base**ZZ(floor(4*8*log(2*d*a**2)/log(2)+2)))
        elif p == 3 and K.discriminant() % 3 == 0: # actually, we only need to do this if e(3) == 4 in the normal closure
            prime_powers.append(base**ZZ(floor(4*8*log(2*d*a**2)/log(3)+2)))
        else:
            prime_powers.append(base**ZZ(floor(4*log(2*d*a**2)/log(p)+1)))
    return prod(prime_powers)


def lauter_viray_bound(K, safe=True, num_etas=None, implementation="bound", bound_calJ=1, bound_two=2):
    """
    Returns a denominator bound for Igusa class polynomials based on [LV].
    See below for the output format.
        
    INPUT:
    
     - `K`   -- a quartic CM-field
     - `safe` -- if True, do a few redundant computations, to detect possible
                 bugs
     - `num_etas` -- if O_K is not of the form O_{K_0}[eta], how many etas
                 with O_K containing O_{K_0}[eta] to use. If None, uses
                 the default of :func:`find_eta`
     - `implementation` -- "bound" or "gjlsvw", if "bound", then output a
                 bound, if "gjlsvw" and the code of [GJLSVW] is loaded
                 into a Magma session in a global variable `m`, then
                 compute an exact denominator with the aid of that code
                 (very slow).
     - `bound_calJ` - where we use upper bounds for calJ (see [LV]) instead
                 of exact numbers, multiply them by this number.
     - `bound_two` - replace the number 2 in Theorems 2.1 and 2.3 of [LV]
                 by this number.
    
    OUTPUT:
    
     - If O_K = O_{K_0}[eta] and implementation="bound" and bound_calJ is an
       integer >= 1 and bound_two is an integer >= 2, then returns a multiple
       of the denominator of the Igusa class polynomials.
     - If O_K is not of the form O_{K_0}[eta], takes a few values of eta
       such that O_K contains O_{K_0}[eta] and returns the output for each
       of those.
     - If implementation="gjlsvw", does not return bounds, but returns
       actual values (if it returns anything at all), hardly tested!
     - If bound_calJ and bound_two are polynomial variables, returns primes p
       and the orders e_p with which they appear in the denominators. Each e_p
       is of the form given, where 0<=bound_calJ<=1 and 0<=bound_two<=2 depend
       on p.

    EXAMPLES::
    
        sage: from recip import *

    The following gave output 4147360000 for an older version of Sage. I don't
    know whether that was correct too (it was a factor 5^4 lower).::

        sage: lauter_viray_bound(CM_Field([5,11,29]), safe=True)
        2592100000000

        sage: P.<a,b> = QQ[]
        sage: lauter_viray_bound(CM_Field([5,11,29]), safe=True, bound_two=b, bound_calJ=a)
        [[2, 4], [5, 3*a + 1], [7, 1], [23, 1]]

        sage: lauter_viray_bound(CM_Field([389,37,245]), safe=True)
        1913578584632940061725625
        sage: lauter_viray_bound(CM_Field([389,37,245]), safe=True, bound_two=b, bound_calJ=a)
        [[5, b], [11, 2*b], [19, 2*b], [29, 1]]
        
    The actual denominators in the class polynomials for the last one are:
    11^2 * 19^6 * 29^4
    A factor 4 and something related to the primes 2, 3, 5 can be expected,
    so really this is
    [[2, *], [3, *], [5, *], [11, 8*b], [19, 8*b], [29, 4]]
    The 29^4 is sharp, the 19^16 is quite far off, and the 11^16 even more.
    This can be due to cancellation, but probably isn't.
    
        sage: K = CM_Field([17, 255, 15300])
        sage: lauter_viray_bound(K, safe=True, bound_two=b, bound_calJ=a)
        [[2, 20*a*b + 15*b],
         [3, 16*a + 4],
         [5, 0],
         [13, 23*a + 2],
         [19, 4*a + 3],
         [67, 2],
         [83, 3*a + 3],
         [157, 8*a + 2],
         [883, 2]]

    Observed denominator in class polynomial:
    13^4 * 19^6 * 67^4 * 83^6 * 157^4 * 883^4
    Seems to suggest cancellation at 3 and non-sharpness or
    cancellation at 2, 13, 19, 83, 157.
    When getting rid of all cancellation, but without caring too much for the
    small primes, we get something like:
    2^(20/3) * 3^(2/3) * 13^2 * 19^4 * 67^2 * 83^4 * 157^2 * 883^2
    for I10 (product over all curves)
    So we have non-sharpness for 13, 19, 83, 157 and maybe 2 and 3.
    We get for the prime 13 an exponent
    23*a+2 = 2, so out of the 23 take 0.
    For the prime 83, we get 3*a+3 = 4, so out of 3, take only 1
    For the prime 157, we get 8*a+2=2, so out of 8, take 0.
    The big discrepancy is in the prime 13.
    
    In Magma, we get the following, which agrees with these conclusions.
    ListAllEmb(-1, 0, 319, -30, 17, 13);
    ...
    2
    ListAllEmb(-1, 0, 319, -30, 17, 19);
    ...
    4
    

    """
    if implementation == "bound":
        impl = lauter_viray_bound_given_eta
    elif implementation == "gjlsvw":
        impl = lauter_viray_win
    else:
        print "Unknown implementation %s" % (impl,)

    L = find_eta(K, how_many=num_etas)
    if len(L) == 1:
        ret = impl(K, L[0], safe=safe, factor_two=False, bound_calJ=bound_calJ, bound_two=bound_two)
        if get_recip_verbose():
            print "[LV] results: %s" % (ret,)
        return ret[0]
    return [(impl(K, safe=safe, eta=eta[0], bound_calJ=bound_calJ, bound_two=bound_two), eta[1]) for eta in L]


def lauter_viray_bound_given_eta(K, eta, safe=True, verbose=False, factor_two=True, bound_calJ=1, bound_two=2):
    """
    Returns a list of pairs (l_i, e_) such that sum of e_i with l_i=l is
    at least twice the right hand side of Theorem 2.1 or 2.3 of [LV].
    
    If `factor_two` is True, uses Theorem 2.3, otherwise, uses Theorem 2.1.
    
    INPUT:
    
     - `K`   -- a quartic CM-field
     - `eta` -- an integral element of K, not in `K_0`
     - `safe` -- if True, do a few redundant computations, to detect possible
                 bugs
     - `factor_two` -- If True, multiply by the factor two at the beginning
       of the right hand side as in Theorem 2.3
       of [LV]. If False, use factors two as in Theorem 2.1.
                        
    
    """
    if get_recip_verbose():
        print "Computing a denominator bound using [LV]"
    ret = []
    cc = K.complex_conjugation()
    F = K.real_field()
    D = ZZ(F.discriminant())
    sqrtD = F(D).sqrt()
    omega = 1/2*(D + sqrtD)
    
    factors_two = []
    # factors_two is a list of primes for which a factor 2 appears at the
    # beginning of the right hand side of [LV, Thm 2.1]
    
    our_basis_to_sage_basis = Matrix([list(F(1)), list(omega)]).transpose()
    sage_basis_to_our_basis = our_basis_to_sage_basis.inverse()

    (alpha0, alpha1) = sage_basis_to_our_basis*vector(list(K.relative_trace(eta)))
    (beta0, beta1) =   sage_basis_to_our_basis*vector(list(K.relative_norm(eta)))
    assert K.real_to_self()(alpha0+alpha1*omega) == eta + cc(eta)
    assert K.real_to_self()(beta0+beta1*omega) == eta * cc(eta)

    if get_recip_verbose():
        print "[LV] alpha0=%s; alpha1=%s; beta0=%s; beta1=%s" % (alpha0, alpha1, beta0, beta1)

    c_K = (alpha0^2 + alpha0*alpha1*D + 1/4*alpha1^2*(D^2-D)-4*beta0-2*beta1*D)
    
    OF = F.maximal_order()
    
    D_rel_of_eta = K.self_to_real((eta-cc(eta))^2)
    Dtilde = ZZ(D_rel_of_eta.norm())

    divisor_of_deltasq = ((Dtilde-c_K^2)/(4*D)).denominator()
    # We restrict to delta such that delta^2 is divisible by this number, as
    # explained below.
    
    for S in xsrange(ZZ(D%2), ZZ(ceil(sqrt(D))), 2):
        # S := the non-negative square root of D-4*delta.
        # Relation to [Yang]: S = (n in [Yang, Cor. 1.6]).
        #                     delta = (the subscript m in [Yang's] b_m)
        delta = ZZ((D-S^2)/4)
        
        if delta^2 % divisor_of_deltasq == 0:
            factors_two = factors_two + delta.prime_divisors()

            if get_recip_verbose():
                print "[LV] delta = %s" % delta

            C_delta = 2 if S == 0 else 1

            t_u = ZZ(alpha1*(D-S^2)/4)     # [LV, page 3]
            t_wxplus = 2*alpha0 + D*alpha1 # Correction to [LV, page 3], by email
            t_wxmin = alpha1*S             # [LV, page 3]

            t_w = ZZ((t_wxplus+t_wxmin)/2)
            t_x = ZZ((t_wxplus-t_wxmin)/2)

            a = ZZ((D-S)/2)                     # [LV, page 10]
            assert t_x == alpha0+a*alpha1       # [LV, page 10]
            assert t_w == t_x+(D-2*a)*t_u/delta # [LV, proof of Lemma 3.4, p.12]
            assert t_u == alpha1*delta          # [Lv, page 10]

            # Next, we list n. Conditions:
            # 4D | (delta^2*Dtilde-n^2) > 0
            # 2D | (n + c_K*delta)
            # Assuming the second condition, the first is equivalent to
            #    0 = delta^2*Dtilde-n^2 = delta^2*Dtilde-c_K^2*delta^2
            #      = delta^2*(Dtilde-c_K^2) mod 4D and
            #    n <= delta*sqrt(Dtilde)
            # So we can require
            #    4D | delta^2*(Dtilde-c_K^2)
            # and list n with n = -c_K*delta mod 2D and n <= delta*sqrt(Dtilde)
            # 
            # Write n = -c_K*delta + 2*D*k
            # then k = (n+c_K*delta) / (2*D) is in the interval
            #          (\pm delta*sqrt(Dtilde) + c_K*delta) / (2*D)
            # Relation to [Yang]: |n| there also is <= delta*sqrt(Dtilde)
            #                     Moreover, the condition
            #                     (n+delta*sqrt(Dtilde))/(2D) is in d_{Ktilde/Ftilde}^-1
            #                     is a congruence condition on n. 
            kmin = ZZ(ceil((-delta*sqrt(Dtilde) + c_K*delta) / (2*D)))
            kmax = ZZ(floor((delta*sqrt(Dtilde) + c_K*delta) / (2*D)))
            for k in srange(kmin, kmax+1):
                n = -c_K*delta + 2*D*k
                if get_recip_verbose():
                    print "[LV] n = %s" % n
                quo = (delta**2*Dtilde - n**2)/(4*D)
                
                assert quo in ZZ
                # This is equivalent to delta^2 % divisor_of_deltasq == 0
                
                if get_recip_verbose() > 2:
                    print "(delta**2*Dtilde - n**2)/(4*D) = %s" % quo
                if not quo in [-1, 1]:
                    # We next list prime divisors ell of quo, so +/- 1 is useless
                    n_u = ZZ(-delta*(n+c_K*delta)/(2*D)) # [LV, page 3]
                    t_xuvee = ZZ(beta1*delta + S*(n+c_K*delta)/(2*D)) # [LV, page 3]
                    n_wxplus = 2*beta0+D*beta1-2*n_u/delta # correction to [LV, page 3], by email
                    n_wxmin = beta1*S # [LV, page 3]
                    n_x = ZZ((n_wxplus-n_wxmin)/2)
                    n_w = ZZ((n_wxplus+n_wxmin)/2)
                    
                    assert n_x == beta0+a*beta1-n_u/delta # [LV, bottom of page 11]
                    assert t_xuvee == beta1*delta - (D-2*a)*n_u/delta # [LV, bottom of page 11]
   
                    d_u = t_u**2 - 4*n_u
                    # d_w = t_w**2 - 4*n_w # never used
                    d_x = t_x**2 - 4*n_x
                    assert t_u == alpha1*delta                 # [ABLPV]
                    assert t_x == alpha0 + 1/2 * (D-S)*alpha1  # [ABLPV]
                    assert   d_u == (alpha1*delta)**2 + 4*(n+c_K*delta)*delta/(2*D) # corrected version of [ABLPV], based on [LV]
                    assert   d_x == (alpha0+1/2*(D-S)*alpha1)**2 - 4*(beta0+1/2*(D-S)*beta1+(n+c_K*delta)/(2*D)) # corrected version of [ABLPV], based on [LV]
                    assert   t_xuvee == beta1*delta + S*(n+c_K*delta)/(2*D) # corrected version of [ABLPV], based on [LV]
                    assert d_u == (alpha1*delta)**2 - 4*n_u    # [LV], top of page 12
                    assert d_x == (alpha0+a*alpha1)**2 - 4*n_x # [LV], top of page 12

                    for l in ZZ(quo).prime_divisors():
                        if get_recip_verbose()>2:
                            print "[LV] l = %s" % l
                        mu_l = quo.valuation(l)
                        if d_u % l != 0 or d_x % l != 0:
                            mu_l = (mu_l + 1)/2
                        for f_u in conductor(d_u).divisors():
                            assert f_u > 0
                            if d_u > 0:
                                f_u = -f_u
                            should_be_discr = ZZ(d_u / f_u**2)
                            assert should_be_discr < 0
                            assert should_be_discr % 4 in [0, 1]
                            
                            if not ((should_be_discr/l**2) in ZZ and ZZ(should_be_discr/l**2) % 4 in [0, 1]):
                                # now we know should_be_discr is a discriminant of an imaginary quadratic order that is maximal at l
                                if get_recip_verbose() > 2:
                                    print "[LV] f_u = %s" % f_u
                                t_nfu = (d_x*d_u-f_u*(t_x*t_u-2*t_xuvee))/(2*f_u**2)
                                computed_calJ_bound = calJ_bound(n=n, delta=delta, Dtilde=Dtilde, d_u=d_u, f_u=f_u, d_x=d_x, t=t_nfu, l=l, D=D, safe=safe, bound_calJ=bound_calJ)
                                if computed_calJ_bound != 0:
                                    ret.append((l, C_delta * mu_l * calI(delta=delta, f_u=f_u, t_w=t_w, l=l, d_u=d_u, n_w=n_w, t_u=t_u, safe=safe) * computed_calJ_bound, delta, n, f_u))
                                else:
                                    ret.append((l, 0, delta, n, f_u))
    fact = []
    for e in ret:
        for i in range(len(fact)):
            if fact[i][0] == e[0]:
                fact[i][1] = fact[i][1] + e[1]
                break
        else:
            fact.append([e[0], e[1]])
    fact.sort()
    if factor_two:
        # Use a factor two everywhere, as in Theorem 2.3
        for i in range(len(fact)):
            fact[i][1] = bound_two*fact[i][1]
    else:
        # Only use a factor two as in Theorem 2.1, specified by factors_two.
        for i in range(len(fact)):
            if fact[i][0] in factors_two:
                fact[i][1] = bound_two*fact[i][1]
    # The "**2" in the following formula is in [Yang, Section 9], converting
    # intersection numbers to valuations of I_10. It is actually a
    # "**(number of roots of unity of K)", but we know that if that number
    # is > 2, then the denominator is trivial as K=QQ(zeta_5).
    if not (bound_two in ZZ and bound_calJ in ZZ):
        return fact, fact, ret, factors_two
    return prod([p**e for [p, e] in fact])**2, fact, ret, factors_two


def calJ_conjecture(d1, d2, t, l):
    """
    Return the right hand side of Conjecture 2.6,
    or None if the hypothesis of the conjecture is not satisfied.
    
    Also returns None if one of the Legendre symbols is ambiguous (i.e.,
    if l>2 and 2 divides m).
    
    Since this is a conjecture, it is only used as a check: counterexamples
    to the conjecture will be reported when safe=True in lauter_viray_bound,
    but otherwise nothing is done with this conjecture.
    If counterexamples are found, it is probably a good idea to check the code
    of this function, as not much care was given to it (it is only used as a
    redundant check).
    """
    if not t in ZZ:
        return None
    d1 = ZZ(d1)
    d2 = ZZ(d2)
    f1 = conductor(d1)
    f2 = conductor(d2)
    if not (d1 < 0 and d2 < 0):
        raise ValueError
    m = ZZ((d1*d2-(d1*d2-2*t)^2)/4)
    if m%2 == 0 and l > 2:
        return None
    if gcd([f1, f2, m]) != 0:
        return None
    ret = 1
    for p in m.prime_divisors():
        if l != p:
            dps = [d for d in [d1,d2] if not ((d/p^2) in ZZ) and ZZ(d/p^2) % 4 in [0,1]]
            assert dps == [d for d in [d1,d2] if conductor(d)%p != 0]
            assert len(dps) > 0
            g = []
            for dp in dps:
                f = []
                if legendre_symbol(dp, p) == 1 and f1%p != 0:
                    f.append(1+m.valuation(p))
                if legendre_symbol(dp, p) == 1 and f1%p == 0:
                    f.append(2)
                if dp%p == 0 and hilbert_symbol(dp, -m, p) == 1 and f1%p != 0:
                    f.append(2)
                if legendre_symbol(dp, p) == -1 and f1%p != 0 and m.valuation(p)%2 == 0:
                    f.append(1)
                if dp%p == 0 and hilbert_symbol(dp, -m, p) == 1 and f1%p == 0 and m.valuation(p) == 2:
                    f.append(1)
                if len(f) == 0:
                    f = [0]
                g = g + f
            assert all([h==g[0] for h in g])
            if g[0] == 0:
                return 0
            ret = ret * g[0]
    return ret
            

def conductor(d):
    """
    Returns the largest integer f such that d/f^2 is a discriminant.
    """
    ret = (sqrt(d/fundamental_discriminant(d)))
    if not ret in ZZ:
        raise ValueError, "Non-discriminant %s pretending to be a discriminant" % d
    return ZZ(ret)


def calJ_bound(n, delta, Dtilde, d_u, f_u, d_x, t, l, D, safe=True, bound_calJ=1):
    """
    Returns the upper bound on calJ of [LV] given in Theorem 2.4 of [LV].
    
    If safe is True, then also checks whether Conjecture 2.6 holds in this case
    and raises an error if it does not.
    """
    thm24 = calJ_bound_thm24(n=n, delta=delta, Dtilde=Dtilde, d_u=d_u, f_u=f_u, d_x=d_x, t=t, l=l, D=D, safe=safe)
    if safe:
        conj = calJ_conjecture(d_u/f_u**2, d_x, t, l)
        if get_recip_verbose() > 2:
            print "calJ <= %s" % thm24
            print "calJ conjectured to be %s" % conj
        if not (conj is None):
            assert conj > thm24
            if gcd(ZZ((delta**2*Dtilde-n**2)/(4*D*f_u^2)), conductor(d_u/f_u**2)) == 1:
                assert conj == thm24
    if thm24 == 0:
        return 0
    if gcd(ZZ((delta**2*Dtilde-n**2)/(4*D*f_u^2)), conductor(d_u/f_u**2)) == 1:
        return thm24
    return bound_calJ*thm24


def calJ_bound_thm24(n, delta, Dtilde, d_u, f_u, d_x, t, l, D, safe=True):
    """
    Returns the upper bound on calJ from Theorem 2.4.
    """
    if hilbert_symbol_trivial_outside(d_u, D*(n**2-delta**2*Dtilde), l):
        if safe:
            # Just to be sure, we'll use the other hilbert symbol too, which is supposed to be the same.
            assert hilbert_symbol_trivial_outside(d_u, (d_u*f_u**-2*d_x-2*t)**2 - d_u*f_u**-2*d_x, l)
        ret = count_ideals(d_u/f_u**2, (delta**2*Dtilde-n**2)/(4*D*l*f_u**2), safe=safe)
        if ret == 0:
            return 0
        extra_factor1 = 2**len([p for p in ZZ(d_u*f_u**-2).prime_divisors() if p != 2 and p != l and t.valuation(p) >= (d_u*f_u**-2).valuation(p)])
        extra_factor2 = rhotilde(d_u*f_u**-2, t, d_x, l)
        if safe and gcd(ZZ((delta**2*Dtilde-n**2)/(4*D*f_u**2)), ZZ(conductor(d_u/f_u**2))) == 1:
            # [LV, Remark 2.5]
            N = (delta**2*Dtilde-n**2)/(4*D)
            if not extra_factor1*extra_factor2 == 2**len([p for p in gcd(ZZ(N*f_u**-2), ZZ(d_u*f_u**-2)).prime_divisors() if p != l]):
                print "Problem with [LV, Remark 2.5], l=%s, factor1=%s, factor2=%s*%s, product should be 2^%s, sets are %s and %s, d_u=%s, f_u=%s, s_0=%s, s_1=%s, N=%s" % (l, extra_factor1.factor(), rhotilde(d_u*f_u**-2,t,d_x,l,part=1), rhotilde(d_u*f_u**-2,t,d_x,l,part=2), len([p for p in gcd(ZZ(N*f_u**-2), ZZ(d_u*f_u**-2)).prime_divisors() if p != l]) , [p for p in ZZ(d_u*f_u**-2).prime_divisors() if p != 2 and p != l and t.valuation(p) >= (d_u*f_u**-2).valuation(p)], [p for p in gcd(ZZ(N*f_u**-2), ZZ(d_u*f_u**-2)).prime_divisors() if p != l], d_u, f_u, t, d_x, N)
                if l != 2:
                    raise RuntimeError
                else:
                    pass
            else:
                if get_recip_verbose() > 1:
                    print "NO prob with [LV, Remark 2.5], factor1=%s, factor2=%s*%s, product should be 2^%s, sets are %s and %s" % (extra_factor1.factor(), rhotilde(d_u*f_u**-2,t,d_x,l,part=1), rhotilde(d_u*f_u**-2,t,d_x,l,part=2), len([p for p in gcd(ZZ(N*f_u**-2), ZZ(d_u*f_u**-2)).prime_divisors() if p != l]) , [p for p in ZZ(d_u*f_u**-2).prime_divisors() if p != 2 and p != l and t.valuation(p) >= (d_u*f_u**-2).valuation(p)], [p for p in gcd(ZZ(N*f_u**-2), ZZ(d_u*f_u**-2)).prime_divisors() if p != l])
        return ret * extra_factor1 * extra_factor2
    else:
        if safe:
            assert not hilbert_symbol_trivial_outside(d_u, (d_u*f_u**-2*d_x-2*t)**2 - d_u*f_u**-2*d_x, l)
        return 0


def count_ideals(disc, norm, safe=True):
    """
    Returns the number of invertible (equivalently: proper) ideals of the given
    norm in the quadratic order of the given discriminant.
    
    If safe=True, then calculates it in a few (possibly slow) ways and checks
    whether the outputs agree.
    """
    if not norm in ZZ:
        return 0
    norm = ZZ(norm)
    c1 = count_ideals1(disc, norm)
    if not safe:
        return c1
    c2 = count_ideals2(disc, norm)
    c3 = count_ideals3(disc, norm)
    if c1 != c2 or c2 != c3:
        raise RuntimeError, "counted ideals incorrectly for disc=%s and norm=%s: %s, %s, %s" % (disc, norm, c1, c2, c3)
    return c1


def count_ideals1(disc, norm):
    """
    Returns the number of invertible ideals of the given norm inside the
    quadratic order of the given discriminant.
    
    Does this quickly by checking some congruences if norm and the conductor
    of the order are coprime. Otherwise, falls back to a slower algorithm:
    :func:`count_ideals3`.
    
    EXAMPLE::

        sage: from recip import *
        sage: count_ideals1(-183, 3)
        1
        sage: count_ideals1(-183, 3*11*13)
        4
        sage: count_ideals1(-183, 7)
        0

    """
    ret = 1
    for p in norm.prime_divisors():
        if p > 2:
            leg = legendre_symbol(disc, p)
        elif disc % 8 == 1: # x^2 + x - c with 1+4*c, so c is even, so 2 roots mod 2, split
            leg = 1
        elif disc % 4 == 0: # ramified
            leg = 0
        elif disc % 8 == 5: # x^2 + x - c with 1+4*c, so c is odd, so no roots mod 2, inert
            leg = -1
        else:
            raise RuntimeError
        if leg == -1:
            if norm.valuation(p) % 2 == 1:
                return 0
            # inert prime with even valuation, does not change ret
        elif leg == 1:
            ret = ret * (norm.valuation(p) + 1)
        else:
            assert disc % p == 0
            if disc % p**2 == 0:
                # singular prime, the hardest case
                return count_ideals3(disc, norm)
                raise NotImplementedError
            # ramified prime, does not change ret
    return ret


def count_ideals2(disc, norm):
    """
    Returns the number of invertible ideals of the given norm inside the
    quadratic order of the given discriminant.
    
    This is a slow implementation, which counts binary quadratic forms. It is
    used only as a redundant safety check in :func:`count_ideals`.
    
    EXAMPLES:
    
    Here is an example of an ideal of norm 3::

        sage: QuadraticField(-183,'a').ideal(3).factor()[0][0].basis()
        [3, 1/2*a + 3/2]
    
    Let n=3, then this ideal is n*(tau*ZZ+ZZ), where tau = (a+3)/6 is a root
    of the primitive polynomial 3*x^2-3*x+16. So in the code of this function,
    this is given by A=3, B=-3, C=16. This is not in the form we are listing,
    but by changing tau by -1, we find tau=(a-3)/6, so B=3, with the same
    A and C. Also N=n/A=1.
    This is the only ideal of norm 3, since 3 is ramified. And indeed, we get::

        sage: from recip import *
        sage: count_ideals2(-183, 3)
        1

    """
    # The fractional ideals are n*(tau*ZZ+ZZ), A*tau^2 + B*tau + C = 0, A>0,
    # gcd(A, B, C) = 1. Note that aa = 
    # A*(tau*ZZ+ZZ) is an invertible ideal of the order O = A*tau*ZZ+ZZ of
    # norm A, not divisible by any integer>1 in ZZ.
    # Without loss of generality take Im(tau) > 0 (embed K into CC).
    # Proof: (A*tau)^2 = A*(-B*tau-C), so indeed an order,
    # (A*tau)*tau = -B*tau-C, so an ideal, index A is trivial, now
    # aa*aabar = A^2*(tau*taubar, tau, taubar, 1)
    #          = A*(C, A*tau, A*(tau+taubar), A)
    #          = A*(A*tau, C, B, A)
    #          = A*(A*tau, 1) = A*O, so invertible.
    # of norm A, so for the whole thing to be integral, take n = N*A.
    # The norm is then N^2*A, so A <= norm.
    # Here B^2 - 4*A*C = disc(O)
    # Every ideal has a unique way of being written as n*aa, where aa is
    # not divisible by an integer. For every aa, the norm AA is fixed.
    # For every aa, tau is defined up to addition by ZZ:
    # B |--> B+2*A*ZZ. So B in [0,2*A) from the translations
    # and then B >= 0 from the sign changes, so B in [0,2*A) is uniquely
    # determined by aa. Conversely, B and A together determine C.
    # Now N, A, B, C together determine n*(tau*ZZ+ZZ) up to complex
    # conjugation
    ret = 0
    for B in srange(disc%2, 2*norm, 2): # we may restrict to B of the same parity as disc
        AC = ZZ((B**2 - disc)/4) # A times C
        for A in gcd(norm, AC).divisors(): # A has to divide the norm and A*C.
            # Now we check if B and A are valid: n has to exist, B has to be in
            # the interval, and there is a gcd condition.
            if (norm/A).is_square() and B < 2*A and gcd([B, A, ZZ(AC/A)]) == 1:
                ret = ret + 1
    return ret


def count_ideals3(disc, norm):
    """
    Returns the number of invertible ideals of the given norm inside the
    quadratic order of the given discriminant.
    
    This is a slow implementation, which counts binary quadratic forms. It is
    used as a fall-back in difficult cases in :func:`count_ideals1`.
    
    EXAMPLES::
    
        sage: from recip import *
        sage: count_ideals3(-400*27, 5)
        0
        sage: count_ideals3(-400*27, 7)
        2
        sage: count_ideals3(-400*27, 7*11)
        0
        sage: count_ideals3(-400*27, 16)
        6
        sage: count_ideals3(-400*27, 32)
        0
        sage: count_ideals2(-400*27, 32)
        0

    """
    ret = 0
    # The fractional ideals are n*(tau*ZZ+ZZ), A*tau^2 + B*tau + C = 0, A>0
    # (tau*ZZ+ZZ) is the inverse of a proper integral ideal for A*tau*ZZ+ZZ
    # of norm A, so for the whole thing to be integral, take n = N*A.
    # The norm is then n^2/A = N^2*A
    # Here B^2 - 4*A*C = disc
    # tau can be replaced by +/- tau + ZZ, so wlog 0<=B<=A
    # A <= norm
    for N in prod([p**floor(e/2) for (p,e) in norm.factor()]).divisors():
        A = ZZ(norm/N**2)
        for B in srange(disc%2, 2*A, 2):
            C = (B**2-disc)/(4*A)
            if C in ZZ and gcd([A, B, ZZ(C)])==1:
                ret = ret + 1
    return ret


def rhotilde(d, s0, s1, l, part=0):
    """
    Returns rhotilde_d^(2)(s0, s1) from Theorem 2.4
    
    If part=0, really returns rhotilde, if part equals 1 or 2, only returns
    the first or second factor.
    """
    if not s0 in ZZ:
        raise RuntimeError
        # This s0 is t(n,f_u), which is (according to an email by authors of
        # [LV]) an integer whenever (delta^2Dtilde-n^2)/(4Dellf_u^2) is.
        # So by calling the function rhotilde only when count_ideals returns
        # something non-zero, this RuntimeError should not occur.
    
    if not part in [0,1,2]:
        raise ValueError
    
    if l == 2:
        return 1 # Addition by email from authors of [LV] (4th May 2013):
    
    # a is the first factor 1 or 2
    if part == 2:
        a = 1
    elif d % 16 == 12 and (s0-s1) % 2 == 0:
        a = 2
    elif d % 8 == 0 and s0.valuation(2) >= d.valuation(2)-2:
        a = 2
    else:
        a = 1
    
    if part == 1:
        return a
    elif d % 32 == 0 and (s0-2*s1) % 4 == 0:
        return 2*a
    return a


def hilbert_symbol_trivial_outside(a, b, l, negative=True):
    """
    Returns True if and only if (a,b)_p != -1 for all primes p not
    equal to l.
    
    If negative=True, then only accepts inputs with a and b negative,
    and performs some additional checks.
    """
    if negative:
        assert a<0 and b<0
    supp_a = [p for (p,e) in a.factor()]
    supp_b = [p for (p,e) in b.factor()]
    for p in union([2], union(supp_a, supp_b)):
        if p != l and hilbert_symbol(a, b, p) == -1:
            return False
    if negative:
        assert hilbert_symbol(a, b, l) == -1 # because the number of local obstructions is even
    return True


def calI(delta, f_u=None, t_w=None, l=None, d_u=None, n_w=None, t_u=None, safe=True):
    """
    Returns calI(n, f_u) as in Theorem 2.1 of [LV].
    """
    ret = 1
    for p in delta.prime_divisors():
        if p != l:
            v_p = delta.valuation(p)
            r_p = max([v_p-min([f_u.valuation(p), ((d_u-t_u*f_u)/(2*f_u)).valuation(p)]), 0])
            assert r_p <= v_p # email by authors of [LV]
            ret = ret * sum([calIp(C=j-r_p, p=p, a1=t_w, a0=n_w, safe=safe) for j in srange(v_p % 2, v_p+1, 2)])
    return ret


def calIp(C, p, a1, a0, safe=True):
    """
    Returns calI_C^{(p)}(a1, a2) as in [LV, Thm 2.1]
    
    With ``safe=True``, an additional redundant slow step is done to verify
    the result.
    """
    if not safe:
        return calIp_fast(C, p, a1, a0)
    I1 = calIp_fast(C, p, a1, a0)
    I2 = calIp_enumerate(C, p, a1, a0)
    if I1 != I2:
        raise RuntimeError, "%s, %s" % (I1, I2)
    return I1


def calIp_fast(C, p, a1, a0):
    """
    Returns calI_C^{(p)}(a1, a2) as in [LV, Thm 2.1], complicated but fast
    implementation.
    
    TESTS::

        sage: from recip import *
        sage: S = [(C, p, a1, a0) for C in srange(-1, 5) for p in [2,3,5] for a0 in srange(15) for a1 in srange(15)]
        sage: [(C, p, a1, a0) for (C, p, a1, a0) in S if calIp_enumerate(C, p, a1, a0) != calIp_fast(C, p, a1, a0)]
        []
    """
    if C < 0:
        return 0
    if C == 0:
        return 1
    disc = a1**2 - 4*a0
    v = disc.valuation(p)
    if p > 2:
        if v >= C:
            return p**floor(C/2)
        if v % 2 != 0:
            return 0
        u = ZZ(disc/p**v)
        if legendre_symbol(u, p) == -1:
            return 0
        return 2*p**ZZ(v/2)
    if v >= C+2:
        if a1 % 2 == 1:
            return 0
        return 2**floor(C/2)
    if v % 2 != 0:
        return 0
    u = ZZ(disc/p**v)
    if u % min(8, 2**(C+2-v)) != 1:
        return 0
    k = ZZ(v/2)
    if (a1 % 2 == 0 and k == 0) or (a1 % 2 == 1 and k > 0):
        return 0
    return 2**(k-1)*min(4, 2**(C+1-2*k))


def calIp_enumerate(C, p, a1, a0):
    """
    Returns calI_C^{(p)}(a1, a2) as in [LV, Thm 2.1], simple but slow
    implementation.
    """
    if C < 0:
        return 0
    if C == 0:
        return 1
    R = Zmod(p**C)
    return len([t for t in R if t**2 - a1*t + a0 == 0])
    

def find_eta(K, how_many=None, proof=True):
    """
    Given a primitive quartic CM-field K, returns a list [eta] of one element
    eta such that O_F[eta] = O_K, if it exists.
    
    If not, then returns a list L of pairs (eta, S), where S is thet set of
    primes dividing [O_K : O_f[eta]].
    The sets S don't contain any primes <= D/4 (D = Disc(F)), and the
    intersection of all sets S is empty.
    
    If how_many is None, use as few as possible. Otherwise, return that number.
    
    EXAMPLE::
    
        sage: from recip import *
        sage: K = CM_Field((x^2+2)^2-3) # a biquadratic field
        sage: find_eta(K)
        [alpha]
    """
    ret = None
    if K.minimal_DAB() == [5,5,5]:
        # return zeta, as OK = ZZ[zeta] = OF[zeta]
        return [K(K.unit_group(proof=proof).gens()[0])]
    Krel = K.relative_field()
    OK = K.maximal_order()
    F = Krel.base_field()
    D = F.discriminant()
    rel_diff = Krel.relative_different()
    if rel_diff.is_principal():
        eta = find_eta_pid(K, Krel, F, rel_diff)
        if not eta is None:
            if how_many is None:
                return [eta]
            ret = [eta]
    if F.class_number() == 1 and ret is None:
        raise RuntimeError, "Is it possible that F has class number one, but eta does not exist? %s" % K

    if ret is None:
        ret = []
    prime_bound = floor(max(D/4, 2))
    if how_many is None:
        how_many = 2
    omega = (K(D).sqrt()-D)/2
    for i in range(how_many - len(ret)):
        eta = find_eta_coprime(K, Krel, F, D, rel_diff, prime_bound)
        new_prime_bound = K.order([eta, omega]).index_in(OK)
        assert new_prime_bound > prime_bound
        prime_bound = new_prime_bound
        if get_recip_verbose() > 1:
            print "Found eta with [OK:OF[eta]] = %s" % prime_bound.factor()
        ret.append((eta, prime_bound))
    
    return ret


def find_eta_coprime(K, Krel, F, D, rel_diff, prime_bound, proof=True):
    eta = K.gen()
    eta = eta*eta.denominator()
    cc = K.complex_conjugation()
    x = Krel.structure()[1](eta-cc(eta))/rel_diff
    xsqr = x.relative_norm()
    assert xsqr == x^2
    F = Krel.base_field()
    x_F = prod([(p**ZZ(e/2)) for (p,e) in xsqr.factor()])
    assert x_F == x
    omega = (K(D).sqrt()-D)/2
    
    if get_recip_verbose() > 1:
        print "index of OF[eta] must be a prime larger than %s" % prime_bound
    
    P, y = next_prime_in_class(x_F, prime_bound, True, True)
    
    # now y*(eta-etabar)/rel_diff is coprime to the bound
    eta = Krel.structure()[0](y)*eta
    
    # Next, to make eta integral...
    
    x = (eta - cc(eta))
    
    # This x now satisfies x in D
    #                  and x/D coprime to prime_range(D/4+1)
    #                  and xbar == -x
    # Next, we need x in OF + 2*OK
    
    if get_recip_verbose() > 1:
        print "Trying (eta-etabar) that would yield [OK:OF[eta]] = %s" % P.norm()
    
    for a in range(2):
        for b in range(2):
            y = a+omega*b
            if x - y in K.ideal(2):
                return (x-y)/2
    
    # So apparently x is not in OF + 2*OK. We need to do something about this.
    raise NotImplementedError, "TODO: Sorry, I never needed this case so far, it should still be implemented"


def next_prime_in_class(a, n, gen=False, allow_one=False):
    """
    Returns a prime P of prime norm p>n in the same ideal class as a.
    
    If gen is True, also return an element y such that y*a=P.
    
    If allow_one is True, then instead of a prime P, also allow P=1, regardless
    of n, if a is principal.
    """
    if get_recip_verbose():
        print "Taking the next prime in a class for a=%s, n=%s" % (a, n)
    K = a.number_field()
    if allow_one and a.is_principal():
        P = K.ideal(1)
        if gen:
            return P, a.gens_reduced()[0]**-1
        return P
    p = ZZ(n)
    while True:
        p = p.next_prime()
        for P in K.ideal(p).prime_factors():
            if P.norm() == p:
                if get_recip_verbose() > 5:
                    print p
                y = P/a
                if y.is_principal():
                    if gen:
                        return P, y.gens_reduced()[0]
                    return P
    

def find_eta_pid_old(K, Krel, F, D, rel_diff):
    """
    Given a primitive quartic CM-field K, returns an one-element list
    [eta] such that O_F[eta] = O_K, if it exists, and return None otherwise.
    Assumes F is a PID. See the source code of find_eta for what the
    input is.
    
    THIS version had bugs and has been deprecated.
    
    """
    assert rel_diff.is_principal()
    rel_diff0 = rel_diff.gens_reduced()[0]
    # Now for every valid eta, we have
    # (eta-etabar) = rel_diff0 up to units in K
    # Iterate over all possible units
    # The unit group is <-1> x <eps>, and -1 is irrelevant
    eps = K(K.unit_group().gens()[1])
    k = 1
    twoK = K.ideal(2)
    twoKrel = Krel.ideal(2)
    while not eps**k - 1 in twoK:
        k = k + 1
        if k > (2^4-1)^4:
            raise RuntimeError
    # The relevant units are 1, eps, ..., eps^(k-1)
    for l in range(k):
        rel_diff = Krel.structure()[1](eps)**k * rel_diff0
        # eta-etabar is totally imaginary
        if rel_diff**2 in F and F(-rel_diff**2).is_totally_positive():
            # Suppose eta-etabar == rel_diff.
            # We have eta+etabar = A for some A in O_F,
            # so 2*eta = rel_diff + A.
            # Goal: find A in O_F, which is rel_diff mod 2*O_K.
            OKrel = Krel.maximal_order()
            qK = OKrel.quotient_ring(twoKrel, 'a')
            OF = F.maximal_order()
            bas = OF.basis()
            for b in cartesian_product_iterator([[0,1], [0,1]]):
                A = sum([bas[i]*b[i] for i in range(2)])
                if Krel(A) - rel_diff in twoKrel:
                    eta = (rel_diff - Krel(A)) / 2
                    # Now eta is in OK, and has the correct discriminant
                    assert eta in OKrel
                    assert eta.trace(F)^2 - 4*eta.norm(F) == Krel.relative_discriminant()
                    return [Krel.structure()[0](eta)]


def find_eta_pid(K, Krel, F, rel_diff):
    """
    Given a primitive quartic CM-field K with no roots of unity
    other than +/- 1, returns an element
    eta of K as a relative field such that O_F[eta] = O_K, if it exists, and
    return None otherwise.
    Assumes F is a PID. See the source code of find_eta for what the
    input is.
    
    EXAMPLE::
    
        sage: from recip import *
        sage: K = CM_Field((x^2+2)^2-3) # a biquadratic field
        sage: Krel = K.relative_field()
        sage: find_eta_pid(K, Krel, K.real_field(), Krel.relative_different())
        alpha
    """
    # Goal: find eta such that O_F[eta] = O_K
    # Suppose such eta exists.
    # Then the relative different is generated by eta - etabar.
    assert rel_diff.is_principal()
    OKrel = Krel.maximal_order()
    rel_diff = rel_diff.gens_reduced()[0]
    # Now eta - etabar = u*rel_diff for a u in O_K^*
    # Note that the left hand side is totally imaginary, hence so is the
    # right.
    #
    # Claim: the group O_F^* of totally real elements inside O_K^*
    # has index at most 2 unless
    # Proof: O_K^*/O_F^* is generated only by 1 element eps by Dirichlet and
    # the assumption that #tors((O_K)^*)=2. The relative norm of eps is
    # eps*epsbar = +/- eps^2, hence eps^2 is in O_F^*.
    #
    # Conclusion: as u*rel_diff is totally imaginary, we get that
    # either rel_diff or eps*rel_diff is totally imaginary.
    eps = K(K.unit_group().gens()[1])
    if K.is_totally_imaginary_element(K(rel_diff)):
        pass
    elif K.is_totally_imaginary_element(eps*K(rel_diff)):
        rel_diff = Krel(eps)*rel_diff
        assert K.is_totally_imaginary_element(K(rel_diff))
    else:
        # contradiction, so eta does not exist
        return None
    # At this stage, we have eta - etabar = u*rel_diff with u totally real
    # So u^-1 eta = (rel_diff + v) / 2 with totally real v.
    #
    # We are free to change eta by a unit and by adding elements of O_F,
    # so we now take u = 1 and v in a fixed set of residues for O_F / 2O_F.
    twoKrel = Krel.ideal(2)
    for v in F.ideal(2).residues():
        eta = (rel_diff - Krel(v))/2
        if eta in OKrel:
            assert eta.trace(F)^2 - 4*eta.norm(F) == Krel.relative_discriminant()
            return K(eta)
    return None


def win_implementation(alpha0, alpha1, beta0, beta1, D, p, m):
    """
    Assumes the code from [GJLSVW] has been loaded into the Magma session m
    (requires magma).
    Returns the output of length of the output of ListAllEmb.
    """
    if get_recip_verbose():
        print "calling [GJLSVW] function with p=%s" % p
    lst = m.ListAllEmb(m(ZZ(alpha0)), m(ZZ(alpha1)), m(ZZ(beta0)), m(ZZ(beta1)), m(ZZ(D)), m(ZZ(p)), nvals=3)
    return lst[2].sage()


def lauter_viray_win(K, eta, safe=True, verbose=False, factor_two=True, m=magma, bound_calJ=1, bound_two=2):
    """
    Returns a list of pairs (l_i, e_) such that sum of e_i with l_i=l is
    at least twice the right hand side of Theorem 2.1 or 2.3 of [LV].
    
    If `factor_two` is True, uses Theorem 2.3, otherwise, uses Theorem 2.1.
    
    INPUT:
    
     - `K`   -- a quartic CM-field
     - `eta` -- an integral element of K, not in `K_0`
     - `safe` -- if True, do a few redundant computations, to detect possible
                 bugs
     - `factor_two` -- If True, multiply by the factor two at the beginning
       of the right hand side as in Theorem 2.3
       of [LV]. If False, use factors two as in Theorem 2.1.
     - bound_calJ and bound_two are not used.
                        
    
    TESTS::
    
        sage: from recip import *
        sage: magma.load("gjlsvw.magma") # optional - magma gjlsvw
        'Loading "gjlsvw.magma"'
        sage: set_recip_verbose(1)
        sage: lauter_viray_bound(CM_Field([8,4,2]), safe=True, implementation="gjlsvw") # indirect doctest, long time 23 minutes, optional - magma gjlsvw
        16
    
    """
    primes = []
    
    if get_recip_verbose():
        print "Computing a denominator bound using [LV]"
    cc = K.complex_conjugation()
    F = K.real_field()
    D = ZZ(F.discriminant())
    sqrtD = F(D).sqrt()
    omega = 1/2*(D + sqrtD)
    
    factors_two = []
    # factors_two is a list of primes for which a factor 2 appears at the
    # beginning of the right hand side of [LV, Thm 2.1]
    
    our_basis_to_sage_basis = Matrix([list(F(1)), list(omega)]).transpose()
    sage_basis_to_our_basis = our_basis_to_sage_basis.inverse()

    (alpha0, alpha1) = sage_basis_to_our_basis*vector(list(K.relative_trace(eta)))
    (beta0, beta1) = sage_basis_to_our_basis*vector(list(K.relative_norm(eta)))
    assert K.real_to_self()(alpha0+alpha1*omega) == eta + cc(eta)
    assert K.real_to_self()(beta0+beta1*omega) == eta * cc(eta)

    if get_recip_verbose():
        print "[LV] alpha0=%s; alpha1=%s; beta0=%s; beta1=%s" % (alpha0, alpha1, beta0, beta1)

    c_K = (alpha0^2 + alpha0*alpha1*D + 1/4*alpha1^2*(D^2-D)-4*beta0-2*beta1*D)
    
    OF = F.maximal_order()
    
    D_rel_of_eta = K.self_to_real((eta-cc(eta))^2)
    Dtilde = ZZ(D_rel_of_eta.norm())

    divisor_of_deltasq = ((Dtilde-c_K^2)/(4*D)).denominator()
    # We restrict to delta such that delta^2 is divisible by this number, as
    # explained below.
    
    for S in xsrange(ZZ(D%2), ZZ(ceil(sqrt(D))), 2):
        # S := the non-negative square root of D-4*delta.
        # Relation to [Yang]: S = (n in [Yang, Cor. 1.6]).
        #                     delta = (the subscript m in [Yang's] b_m)
        delta = ZZ((D-S^2)/4)
        
        if delta^2 % divisor_of_deltasq == 0:
            factors_two = factors_two + delta.prime_divisors()

            if get_recip_verbose():
                print "[LV] delta = %s" % delta

            C_delta = 2 if S == 0 else 1

            t_u = ZZ(alpha1*(D-S^2)/4)     # [LV, page 3]
#            t_wxplus = 2*alpha0 + D*alpha1 # Correction to [LV, page 3], by email
#            t_wxmin = alpha1*S             # [LV, page 3]

#            t_w = ZZ((t_wxplus+t_wxmin)/2)
#            t_x = ZZ((t_wxplus-t_wxmin)/2)

            a = ZZ((D-S)/2)                     # [LV, page 10]
#            assert t_x == alpha0+a*alpha1       # [LV, page 10]
#            assert t_w == t_x+(D-2*a)*t_u/delta # [LV, proof of Lemma 3.4, p.12]
            assert t_u == alpha1*delta          # [Lv, page 10]

            # Next, we list n. Conditions:
            # 4D | (delta^2*Dtilde-n^2) > 0
            # 2D | (n + c_K*delta)
            # Assuming the second condition, the first is equivalent to
            #    0 = delta^2*Dtilde-n^2 = delta^2*Dtilde-c_K^2*delta^2
            #      = delta^2*(Dtilde-c_K^2) mod 4D and
            #    n <= delta*sqrt(Dtilde)
            # So we can require
            #    4D | delta^2*(Dtilde-c_K^2)
            # and list n with n = -c_K*delta mod 2D and n <= delta*sqrt(Dtilde)
            # 
            # Write n = -c_K*delta + 2*D*k
            # then k = (n+c_K*delta) / (2*D) is in the interval
            #          (\pm delta*sqrt(Dtilde) + c_K*delta) / (2*D)
            # Relation to [Yang]: |n| there also is <= delta*sqrt(Dtilde)
            #                     Moreover, the condition
            #                     (n+delta*sqrt(Dtilde))/(2D) is in d_{Ktilde/Ftilde}^-1
            #                     is a congruence condition on n. 
            kmin = ZZ(ceil((-delta*sqrt(Dtilde) + c_K*delta) / (2*D)))
            kmax = ZZ(floor((delta*sqrt(Dtilde) + c_K*delta) / (2*D)))
            for k in srange(kmin, kmax+1):
                n = -c_K*delta + 2*D*k
                if get_recip_verbose():
                    print "[LV] n = %s" % n
                quo = (delta**2*Dtilde - n**2)/(4*D)
                
                assert quo in ZZ
                # This is equivalent to delta^2 % divisor_of_deltasq == 0
                
                if get_recip_verbose() > 2:
                    print "(delta**2*Dtilde - n**2)/(4*D) = %s" % quo
                if not quo in [-1, 1]:
                    # We next list prime divisors ell of quo, so +/- 1 is useless
                    n_u = ZZ(-delta*(n+c_K*delta)/(2*D)) # [LV, page 3]
                    d_u = t_u**2 - 4*n_u
                    # d_w = t_w**2 - 4*n_w # never used
                    assert   d_u == (alpha1*delta)**2 + 4*(n+c_K*delta)*delta/(2*D) # based on [LV]
                    assert d_u == (alpha1*delta)**2 - 4*n_u    # [LV], top of page 12
                    for l in ZZ(quo).prime_divisors():
                        minus_N = (n**2 - delta**2*Dtilde)/(4*D)
                        if hilbert_symbol_trivial_outside(d_u, minus_N, l) and hilbert_symbol(d_u, minus_N, l) == -1:
                            if factor_two or delta % l == 0:
                                mult = 2
                            else:
                                mult = 1
                            primes.append((l, mult, delta, n))
    sorted_primes = []
    
    for e in primes:
        for i in range(len(sorted_primes)):
            if sorted_primes[i][0] == e[0]:
                if e[1] == 2:
                    sorted_primes[i] = ((e[0], e[1]))
                break
        else:
            sorted_primes.append((e[0], e[1]))
    sorted_primes.sort()
    
    if get_recip_verbose():
        print "[LV] relevant primes ell: %s" % (sorted_primes,)
    fact = []
    
    for (p, mult) in sorted_primes:
        e = mult*win_implementation(alpha0, alpha1, beta0, beta1, D, p, m)
        fact.append((p, e))

    # The "**2" in the following formula is in [Yang, Section 9], converting
    # intersection numbers to valuations of I_10. It is actually a
    # "**(number of roots of unity of K)", but we know that if that number
    # is > 2, then the denominator is trivial as K=QQ(zeta_5).
    return prod([p**e for [p, e] in fact])**2, fact, primes

    

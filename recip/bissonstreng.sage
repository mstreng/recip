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

This file gives the details of computations mentioned and used in
Bisson-Streng -- On polarised class groups of orders in quartic CM fields
http://arxiv.org/abs/1302.3756

#----------------------------------------------------------------------------

It may be very Sage version dependent. The tests take a long time, so use
"sage -t --long --timeout 3600 bissonstreng.sage" for long-testing this file.

The computations for Section 6.1.4 "The case of maximal orders" and Table 1
are in the examples section of the documentation of
:func:`orders_proven_with_wamelen`.

Here are the computations for the field `K=\QQ(\zeta_5)`, that is, the second
halve of Section 6.1.5::

        sage: from recip import *
        sage: lst = wamelen_dab_list()
        sage: assert(lst[0]) == [5,5,5]
        sage: K = CM_Field(lst[0])
        sage: OK = K.maximal_order()
        sage: zeta = K(K.unit_group().gens()[0])
        sage: assert zeta.multiplicative_order() == 10
        sage: Phi = K.CM_types()[0]
        
We start by computing the orders `O = \ZZ[\mu_i^5 : i] + p^k O_K`, where ``mu*mubar``
is in `\QQ`, `(mu) = N_{\Phi^r}(A)` and A ranges over the ray class group of
`K^r` mod `p^k`. We let `k=1,2,3,...` and notice that the sequences stabilize::

        sage: minorders = minimal_orders_of_class_number_one(Phi); minorders
        [(2,
          2,
          16,
          Order in CM Number Field in alpha with defining polynomial x^4 + 5*x^2 + 5),
         (3,
          1,
          9,
          Order in CM Number Field in alpha with defining polynomial x^4 + 5*x^2 + 5),
         (5,
          3,
          15625,
          Order in CM Number Field in alpha with defining polynomial x^4 + 5*x^2 + 5)]

In particular, the sequences stabilize at `2^2` with `[O_K:O_2]=16`, at `3^1`
with `[O_K:O_3]=9`, at `5^3` with `[O_K:O_5]=15625`, and at `p^0` with `O=O_K`
for all `p>5`. The orders `\ZZ[\zeta_5^{e_i}\mu_i : i]` contain
`O_2\cap O_3\cap O_5`.

We will not try to reduce the powers of 2 and 3 (as we know from scratch
computations that they are optimal). But we will do something about the
power of 5, because the following takes too long for F = 2^2*3*5^3 or even
F = 5^3::

        sage: mu = CM_type_to_mus(Phi, F) # not tested
        sage: from recip import _order_mod
        sage: [_order_mod(m, K.ideal(F)) for m in mu] # not tested

So let's reduce the power of 5::

        sage: F = 25
        sage: mu = CM_type_to_mus(Phi, F)
        sage: [_order_mod(m, K.ideal(F)) for m in mu] # exact output can change when different generators are used
        [50, 5, 5, 5, 5]

There are at most `5^5=3125` orders `\ZZ[\zeta^{e_i}\mu_i]` for `F=25` and we
compute them all. At the time of original writing, this took about 5 minutes.
Now it is faster.::

        sage: gens = [F*b for b in OK.basis()]
        sage: orders = [K.order(additive_gens_to_basis(gens + [zeta**i * mu[0], zeta**j * mu[1], zeta**k * mu[2], zeta**l * mu[3], zeta**m * mu[4]], K)) for i in range(5) for j in range(5) for k in range(5) for l in range(5) for m in range(5)] # long time: 5 minutes
        sage: indices = [A.index_in(OK) for A in orders] # long time
        sage: [a for a in indices if a != 1] # long time
        [25]
        sage: [Omin5] = [orders[k] for k in range(len(orders)) if indices[k]>1] # long time
        sage: Omin5.basis() # long time
        [1/2*alpha^3 + 1/2, 1/2*alpha^3 + 5/2*alpha^2 + 5/2*alpha, 5*alpha^2, alpha^3]
        sage: OK.basis()
        [1/2*alpha^3 + 1/2, 1/2*alpha^3 + 1/2*alpha^2 + 1/2*alpha, alpha^2, alpha^3]

So there is a unique minimal order `Omin5=O_{min,25}` with
`S_{Omin5} = I_{K^r}` and it has index 25 and satisfies `5O_K\subset Omin5`.
So if O has `S_O=I_{K^r}`, then `O + 25 O_K \supset \ZZ + 5 O_K`, so that is
where this sequence stabilizes and `O + 5^l O_K \supset Omin5`.
We get that `O` contains `O_2\cap Omin5 \cap O_5`. We can therefore use
`F = 2^2*3*5`, which we do::

        sage: F = 4*3*5
        sage: mu = CM_type_to_mus(Phi, F)
        sage: mu = [(_order_mod(m, K.ideal(F)), m) for m in mu]
        sage: mu.sort(reverse=True)
        sage: [m[0] for m in mu]
        [20, 10, 10, 10, 2, 1]
        sage: mu = [m[1] for m in mu]
        sage: gens = [F*b for b in OK.basis()] + mu[4:]
        sage: orders = [K.order(additive_gens_to_basis(gens + [zeta**i * mu[0], zeta**j * mu[1], zeta**k * mu[2], zeta**l * mu[3]], K)) for i in range(5) for j in range(5) for k in range(5) for l in range(5)] # long time
        sage: orders = [(A.index_in(OK), A) for A in orders] # long time
        sage: orders = [a for a in orders if a[0] != 1] # depends on a line with long time
        sage: orders.sort(reverse=True) # depends on a line with long time
        sage: [a[0].factor() for a in orders] # depends on a line with long time
        [5^2, 2^4, 3^2]
        sage: orders = [a[1] for a in orders] # depends on a line with long time

This shows that every order with CM class number one
contains one of three orders that we found, of index 2^4=16, 5^2=25,
and 3^2=9 respectively. In particular, there are no mixed cases,
only prime power indices.::
        
        sage: [len(superorders_stable_under_complex_conjugation(A))-1 for A in orders] # long time
        [2, 3, 1]
    
So besides the maximal order, there are 2 orders of index a power of 5
(one has index 5, one has index 25),
there are 3 of index a power of 2, and there is one of index 3^2=9.
In total, 2+3+1=6 non-maximal orders.
        
Now let's compute class polynomials for all six orders. The actual coefficients
of these polynomials are only numerically correct, but the degrees for all
orders are exactly correct.::
    
        sage: s5 = superorders_stable_under_complex_conjugation(orders[0]) # long time
        sage: pols5 = [class_polynomials_order(O, 5, Phi, 300) for O in s5] # long time, about a minute
        sage: pols5 # long time
        [[y, 0, 0], [y - 816480, 0, 2418647040000], [1, 0, 0]]
        sage: s2 = superorders_stable_under_complex_conjugation(orders[1]) # long time
        sage: pols2 = [class_polynomials_order(O, 4, Phi, 300) for O in s2] # long time, 3 minutes
        sage: pols2 # long time
        [[y, 0, 0], [y - 637875/4, 4556250, 576650390625/4], [y^2 - 4005855/4*y + 47003299245/16, 8901090*y - 67147570350, 15699911214375/4*y - 183904605019209375/16], [y^2 + 21381570*y - 21353299380, 715703040*y - 428554022400, 5986604920320000*y - 5978328612480000000]]
        sage: [len(p[0].factor()) for p in pols2] # long time
        [1, 1, 1, 1]
        sage: s3 = superorders_stable_under_complex_conjugation(orders[2]) # long time
        sage: pols3 = [class_polynomials_order(O, 3, Phi, 300) for O in s3] # long time, half a minute
        sage: pols3 # long time
        [[y, 0, 0], [y^2 - 238369068000/361*y + 67885776000000/361, 335606284800000/130321*y - 58147891200000000/130321, 226410725508480000000000/130321*y - 64480155909120000000000000/130321]]
        sage: [len(p[0].factor()) for p in pols3] # long time
        [1, 1]

If we were to believe these numerical coefficients, then the quadratic
polynomials are irreducible, and the linear polynomials correspond exactly
to the curves from Theorem 6.2 and Conjecture 6.3. The curve in Theorem 6.2
is proven using Magma and AVIsogenies, which is not part of this file. So we
now only need to prove that the three quadratic polynomials are really
quadratic, and that the orders s5[1] and s3[1] are as claimed in Theorem 6.2.

We first prove that the period matrix for the order of index a power of 5
is (5,5)-isogenous to the one for the maximal order::

        sage: [Zmax] = period_matrices(s5[0], 1, Phi) # long time
        sage: [Zother] = period_matrices(s5[1], 5, Phi) # long time
        sage: are_nn_isogenous(Zmax, Zother, 5, 1, 5, transformation=True) # long time
        (True, 1/10*alpha^3 - 1/2)

Next, we prove that the period matrices for the order of index a power of 3
is (3,3)-isogenous to the one for the maximal order::

        sage: [Zmax] = period_matrices(s3[0], 1, Phi) # long time
        sage: Zothers = period_matrices(s3[1], 3, Phi) # long time
        sage: are_nn_isogenous(Zmax, Zothers[0], 3, 1, 3, transformation=True) # long time
        (True, 1)
        sage: are_nn_isogenous(Zmax, Zothers[1], 3, 1, 3, transformation=True) # long time
        (True, 1)

A calculation modulo 23 using AVIsogenies as in the article shows that there
exist no (3,3)-isogenous curves to `y^2 = x^5+1` over `\QQ`, which shows that
there exist no curves with endomorphism ring of index a power of 3 in
`\ZZ[\zeta_5]`. A calculation using AVIsogenies directly over `\QQbar` shows
that the curve C given in Theorem 6.2 really is (2,2)-isogenous to
`y^2 = x^5+1`. So there are now only two things to do for this field:
show that the endomorphism rings of Theorem 6.2 are correct, and show that
the other curves with 2-power index endomorphism ring are not defined over
`\QQ`. The second we do for all fields simultaneously later. The first, we do
now::

        sage: O = s5[1]; # long time
        sage: C.<zeta> = CyclotomicField(5)
        sage: OinC = [C.order([phi(b) for b in O.basis()]) for phi in K.embeddings(C)] # long time
        sage: len(OinC) == 4 and all([A == OinC[0] for A in OinC[1:]]) # long time
        True
        sage: OinC[0].basis() # long time
        [1, 3*zeta^3 + zeta, zeta^3 + zeta^2, 5*zeta^3]
        
        sage: O = s2[1]; # long time
        sage: C.<zeta> = CyclotomicField(5)
        sage: OinC = [C.order([phi(b) for b in O.basis()]) for phi in K.embeddings(C)] # long time
        sage: len(OinC) == 4 and all([A == OinC[0] for A in OinC[1:]]) # long time
        True
        sage: OinC[0].basis() # long time
        [1, 2*zeta, zeta^3 + zeta^2, 2*zeta^3]

The documentation of :func:`all_period_matrices_two` proves Lemma 6.6 and
a little more. It proves that every period matrix for a non-maximal order
with S_O=I_Kr
in every field other than `\QQ[\zeta_5]` is
(2,2)-isogenous to a unique period matrix for a maximal order, and it proves
that every period matrix for an order of 2-power index 
with S_O=I_Kr in `\ZZ[\zeta_5]` is
(2,2)-isogeny connected to the maximal one. The actual (2,2)-isogenous curves
are then computed using AVIsogenies (not in this file), which proves that they
are correct and not defined over `\QQ`.

"""

def orders_proven_with_wamelen(latex_output=False):
    """
    Compute the orders O''' in the article, and prove that the non-maximal
    O''' have index 2 in the maximal order and no proper polarized ideal
    classes. Then output Table 1 of the article.
    
    EXAMPLES::
    
        sage: from recip import *
        sage: orders_proven_with_wamelen()
        [([5, 5, 5], 1, 1, 1, 1, 1), ([8, 4, 2], 1, x^4 + 4*x^2 + 2, 1, 1, 1), ([13, 13, 13], 1, x^4 - x^3 + 2*x^2 + 4*x + 3, 3, 1, 1), ([5, 10, 20], 2, x^4 + 10*x^2 + 20, 4, 4, 2), ([5, 65, 845], 2, x^4 - x^3 + 16*x^2 - 16*x + 61, 19, 1, 1), ([29, 29, 29], 1, x^4 - x^3 + 4*x^2 - 20*x + 23, 7, 1, 1), ([5, 85, 1445], 2, x^4 - x^3 + 21*x^2 - 21*x + 101, 29, 1, 1), ([37, 37, 333], 1, x^4 - x^3 + 5*x^2 - 7*x + 49, 21, 3, 1), ([8, 20, 50], 2, x^4 + 20*x^2 + 50, 25, 25, 1), ([13, 65, 325], 2, x^4 - x^3 + 15*x^2 + 17*x + 29, 23, 1, 1), ([13, 26, 52], 2, x^4 + 26*x^2 + 52, 36, 36, 2), ([53, 53, 53], 1, x^4 - x^3 + 7*x^2 + 43*x + 47, 13, 1, 1), ([61, 61, 549], 1, x^4 - x^3 + 8*x^2 - 42*x + 117, 39, 3, 1)]
    """
    M = [[] for i in range(19)]
    M[1] = [[0, 0, 0, -1], [0, 0, -1, 0], [2, 2, 0, 0], [2, 1, 0, 0]]
    M[2] = [[-1, 1, -2, -1], [-1, -1, -2, 0], [3, -1, 3, 2], [0, -1, -1, 0]]
    M[3] = [[0, 0, 0, 1], [0, 0, 1, -1], [-2, -4, 0, 0], [-4, 2, 0, 0]]
    M[4] = [[0, -2, 2, 0], [3, 0, 0, -4], [1, 0, 0, -3], [0, 0, 2, 0]]
    M[5] = [[1, 1, -1, 1], [-1, -2, 1, -2], [11, 6, 0, 2], [4, 9, 0, 2]]
    M[6] = [[1, 7, 4, 0], [6, 6, -1, 4], [-13, -13, 1, -7], [-13, -26, -6, -7]]
    M[7] = [[0, 0, 1, 1], [0, 0, 1, 2], [-13, 4, -3, -5], [4, -3, 1, 4]]
    M[8] = [[-8, 5, -9, 5], [4, -13, 4, -13], [9, -4, 9, -4], [-3, 14, -3, 13]]
    M[9] = [[2, 0, -3, 1], [-1, -1, 1, -2], [5, 4, -1, 0], [5, 9, -1, 1]]
    M[10] = [[-2, -3, 1, -1], [2, 1, -2, 1], [-1, 0, 3, -4], [5, 4, 1, -1]]
    M[11] = [[-1, 0, -1, -1], [-1, 0, -1, -2], [11, -10, 1, 1], [-10, 15, 0, 0]]
    M[12] = [[-1, 2, 1, -1], [-1, 0, -1, -1], [-9, 5, 1, 1], [5, -2, -2, 0]]
    M[13] = [[-2, 2, 0, -1], [-7, 8, -1, -5], [1, -9, 0, 4], [-8, 12, -1, -5]]
    M[14] = [[-1, -1, -3, -3], [-2, 1, -2, -5], [3, -1, 0, 3], [0, 4, 3, 1]]
    M[15] = [[-1, -1, -1, 0], [-1, -1, 0, -1], [6, 8, 1, 1], [8, 24, 1, 1]]
    M[16] = [[-2, -3, -2, 0], [0, 0, 0, -2], [4, 3, 2, 0], [3, 11, 3, 0]]
    M[17] = [[-5, 3, 1, 5], [-2, 1, 4, 7], [1, -1, 2, 3], [-4, 0, 1, 3]]
    M[18] = [[3, -3, 2, 3], [0, 0, 3, 3], [-7, 4, -5, -9], [-5, 2, 2, 3]]
    pols = [(Matrix(QQ, m).characteristic_polynomial()) for m in M]
    duplicates = [1 if k>1 and pols[k-1] == pols[k] else 0 for k in range(19)]
    pols = [pols[k] for k in range(19) if not duplicates[k]]

    ncurves = []
    for k in range(19):
        if duplicates[k]:
            ncurves[-1] = ncurves[-1]+1
        else:
            ncurves.append(1)
    
    Ks = [CM_Field([5,5,5])] + [CM_Field(p) for p in pols[1:]]

    wamelen_orders = [K.order(K.gen()) for K in Ks]
    wamelen_orders[0] = Ks[0].maximal_order()
    wamelen_indices = [O.index_in(O.number_field().maximal_order()) for O in wamelen_orders]

    orders_with_conjugate = [K.order([K.gen(), K.complex_conjugation()(K.gen())]) for K in Ks]
    orders_with_conjugate[0] = wamelen_orders[0]
    indices_with_conjugate = [O.index_in(O.number_field().maximal_order()) for O in orders_with_conjugate]

    orders_with_lemma = [minimal_order_cl_nr_one_F(O.number_field().CM_types()[0], Ostart=O) for O in orders_with_conjugate]    
    indices_with_lemma = [O.index_in(O.number_field().maximal_order()) for O in orders_with_lemma]

    proper_pol_ideals = [polarized_ideal_classes(orders_with_lemma[k], 2)  for k in range(13) if indices_with_lemma[k] != 1]

    assert proper_pol_ideals == [[],[]]
    assert all([i in [1,2] for i in indices_with_lemma])

    DAB = [[5,5,5]] + [CM_Field(p).minimal_DAB() for p in pols[1:]]

    out = zip(DAB, ncurves, map(latex, pols) if latex_output else pols, wamelen_indices, indices_with_conjugate, indices_with_lemma)

    if not latex_output:
        return out

    t = ""
    for o in out[1:]:
        t = t + "%s & %s & %s & %s & %s & %s \\\\  \\hline \n" % o

    return t


def minimal_orders_of_class_number_one(Phi, output_type='order'):
    r"""
    Let K be the domain of the input CM-type Phi. If K is not `\QQ(\zeta_5)`,
    returns the list of tuples (p, k, i, O), where p ranges over all primes
    such that there is an order that is not p-maximal, but does have
    Shimura class
    number one. Here O is the smallest such order of index a power i of p, and
    `p^k` is the exponent of the additive group OK/O.
    
    If K is `\QQ(\zeta_5)`, then the output is as above, except that the list
    may contain too many primes p, and the orders O may be too small
    (in other words, for a classification, the output will be on the safe
    side).
        
    EXAMPLES::
    
        sage: from recip import *
        sage: lst = wamelen_dab_list()

        sage: K = CM_Field([5,5,5])
        sage: Phi = K.CM_types()[0]
        sage: minorders = minimal_orders_of_class_number_one(Phi); minorders # long time 1 second
        [(2, 2, 16, Order in CM Number Field in alpha with defining polynomial x^4 + 5*x^2 + 5), (3, 1, 9, Order in CM Number Field in alpha with defining polynomial x^4 + 5*x^2 + 5), (5, 3, 15625, Order in CM Number Field in alpha with defining polynomial x^4 + 5*x^2 + 5)]

        sage: Phi = CM_Field(lst[1]).CM_types()[0]
        sage: m = minimal_orders_of_class_number_one(Phi); m
        [(2, 1, 4, Order in CM Number Field in alpha with defining polynomial x^4 + 4*x^2 + 2)]
    """
    K = Phi.domain()
    OK = K.maximal_order()
    ret = []
    for p in _primes_of_interest(K):
        if get_recip_verbose():
            print "p = %s" % p
        k = 1
        previous = OK
        while True:
            if get_recip_verbose():
                print "k = %s" % k
            O = minimal_order_cl_nr_one_F(Phi, p**k, output_type=output_type)
            if O.index_in(previous) == 1:
                if k != 1:
                    ret.append((p,k-1,previous.index_in(OK),previous))
                break
            k = k + 1
            previous = O
    return ret


def minimal_order_cl_nr_one_F(Phi, F=None, output_type='order', Ostart=None):
    """
    Returns an order O in the domain K of Phi such that for any order O' in K
    containing ZZ+F*OK,
    if the image of I_Kr(F) in CCC(O) is trivial,
    then O contains O'.
    
    Returns None every O' satisfies this property, i.e., if O does not exist.
    
    If K is not QQ(zeta_5), then actually the "if" in the first paragraph
    is an "if and only if".
    
    Here F is an integer.
    """
    K = Phi.domain()

    OK = K.maximal_order()

    if F is None:
        if Ostart is None:
            raise ValueError
        F = Ostart.index_in(OK)

    gens = CM_type_to_mus(Phi, F)


    if K.unit_group().torsion_generator().multiplicative_order() != 2:
        if K.unit_group().torsion_generator().multiplicative_order() != 10:
            raise RuntimeError, K
        gens = [mu**5 for mu in gens]
    
    # explanation of why it is going to work: up to roots of unity,
    # these elements are in every O', and since "-a in O' <=> a in O'",
    # we know that the element to the power e is in O'.
    
    gens = [F*b for b in OK.basis()] + gens
    if not Ostart is None:
        gens = Ostart.basis() + gens

    # Sage's order construction is slow for these long lists of elements, while
    # we already give a reasonable suborder to start with. This can be improved
    # by giving a basis of the subgroup generated by our input elements.
    
    if output_type == 'gens':
        return gens
    if output_type == 'direct_from_gens':
        return K.order(gens)

    gens = additive_gens_to_basis(gens, K)

    if output_type == 'order':
        return K.order(gens)
    elif output_type == 'magma':
        return magma(gens).Order()
    else:
        raise ValueError, "Unkown output_type: '%s'" % output_type


def CM_type_to_mus(Phi, F):
    """
    Returns mu generating N_Phir(aaa) with mu*mubar in QQ for aaa ranging
    over generators of the ray class group mod F of the reflex field of Phi.
    """
    K = Phi.domain()
    Psi = Phi.reflex()
    Kr = Psi.domain()
    Fidl = Kr.ideal(F)
    if not F in ZZ:
        raise TypeError
    
    gens = []
    for C in Kr.class_group().gens():
        mu = a_to_mu(Psi, lift_ray_class_group_element(C.ideal(), Kr.ideal(1), Fidl))
        if mu is None:
            return None
        gens.append(mu)
    
    for b in Fidl.idealstar(2).gens():
        gens.append(Psi.type_norm(Kr(b)))

    return gens


def _order_mod(a, I):
    
    for i in range(1, I.norm()):
        if a^i-1 in I:
            return i
    raise RuntimeError

    
def _primes_of_interest(K):
    r"""
    Returns a set of primes, containing all primes dividing the index of the
    minimal order of Shimura class number one.
    
    Assumes that the input is a primitive quartic CM-field and that the maximal
    order has Shimura class number one. Uses results from
    Bisson-Streng.
    
    More precisely, if K is not `\QQ(\zeta_5)`, then Theorem 5.1 can be used
    directly, as it does not depend on the computations that we did with this
    Sage file. It is equivalent to combining Theorem 4.3 with Corollary 5.3.
    
    If K is `\QQ(zeta_5)`, then part of the purpose of this Sage file is to
    prove Theorem 5.1, so we cannot use that theorem. Instead, we use
    Proposition 4.2, and Corollary 5.3, together with the fact that the
    discriminant of K is a power of 5.
    """
    if K.minimal_DAB() == [5,5,5]:
        return prime_range(19)
    return prime_divisors(2*3*
                ZZ(K.discriminant()/K.real_field().discriminant()**2))
                

def all_period_matrices_two(lst):
    """
    Find all period matrices with CM by an order of 2-power index in OK such
    that the field of moduli is in Kr, where K ranges over the fields given by
    lst.
    
    And tests whether for all but [5,5,5] that there are no orders with
    non-2-power index such that the field of moduli is in Kr
    (by raising an error if there is such an order).
    
    OUTPUT:
    
    - List of triples (K, Phi, F, Zs, orders), where
      Zs[k] is a period matrix has CM by orders[k] for each k, and
      F*OK is contained in each element of orders.
      
    
    EXAMPLES::
    
        sage: from recip import *
        sage: lst = wamelen_dab_list() + [[5, 15, 45], [5, 30, 180], [5, 35, 245], [5, 105, 2205], [8, 12, 18], [17, 119, 3332], [17, 255, 15300]]
        sage: Zs = all_period_matrices_two(lst) # long time, 2 minutes
        sage: [len(l) for (K,Phi,F,l,orders) in Zs] # long time
        [6, 3, 1, 6, 2, 1, 2, 1, 6, 2, 6, 1, 1, 4, 12, 4, 4, 6, 2, 4]
        
        sage: matr = [two_two_matrix(l, F) for (K, Phi, F, l, orders) in Zs] # long time
        sage: matr # long time, output is random
        [
        [0 1 1 1 0 0]                [0 1 1 0 1 0]
        [1 0 0 0 1 1]                [1 0 0 1 0 1]
        [1 0 0 0 1 1]                [1 0 0 1 0 1]
        [1 0 0 0 1 1]  [1 1 1]       [0 1 1 0 1 0]
        [0 1 1 1 0 0]  [1 1 1]       [1 0 0 1 0 1]  [0 0]       [0 0]
        [0 1 1 1 0 0], [1 1 1], [0], [0 1 1 0 1 0], [0 0], [0], [0 0], [0],
        <BLANKLINE>
        [1 0 0 1 0 1]         [0 1 1 0 1 0]
        [0 1 1 0 1 0]         [1 0 0 1 0 1]
        [0 1 1 0 1 0]         [1 0 0 1 0 1]            [0 1 0 0]
        [1 0 0 1 0 1]         [0 1 1 0 1 0]            [1 0 0 0]
        [0 1 1 0 1 0]  [0 0]  [1 0 0 1 0 1]            [0 0 0 1]
        [1 0 0 1 0 1], [0 0], [0 1 1 0 1 0], [0], [0], [0 0 1 0],
        <BLANKLINE>
        [0 0 1 0 1 0 1 0 0 0 0 0]
        [0 0 0 1 0 1 0 1 0 0 0 0]
        [1 0 0 0 0 0 0 0 1 0 1 0]
        [0 1 0 0 0 0 0 0 0 1 0 1]
        [1 0 0 0 0 0 0 0 1 0 1 0]
        [0 1 0 0 0 0 0 0 0 1 0 1]
        [1 0 0 0 0 0 0 0 1 0 1 0]                        [1 0 0 1 0 1]
        [0 1 0 0 0 0 0 0 0 1 0 1]                        [0 1 1 0 1 0]
        [0 0 1 0 1 0 1 0 0 0 0 0]  [0 0 0 1]  [0 0 0 0]  [0 1 1 0 1 0]
        [0 0 0 1 0 1 0 1 0 0 0 0]  [0 0 1 0]  [0 0 0 0]  [1 0 0 1 0 1]
        [0 0 1 0 1 0 1 0 0 0 0 0]  [0 1 0 0]  [0 0 0 0]  [0 1 1 0 1 0]  [0 0]
        [0 0 0 1 0 1 0 1 0 0 0 0], [1 0 0 0], [0 0 0 0], [1 0 0 1 0 1], [0 0],
        <BLANKLINE>
        [0 0 0 0]
        [0 0 0 0]
        [0 0 0 0]
        [0 0 0 0]
        ]

    The output above was with an older version of Sage. We now check
    that the actual output agrees up to graph isomorphism, and then we apply
    a graph isomorphism and continue with the result after isomorphism.::

        sage: L = [[[0, 1, 1, 1, 0, 0],
        ....:   [1, 0, 0, 0, 1, 1],
        ....:   [1, 0, 0, 0, 1, 1],
        ....:   [1, 0, 0, 0, 1, 1],
        ....:   [0, 1, 1, 1, 0, 0],
        ....:   [0, 1, 1, 1, 0, 0]],
        ....:  [[1, 1, 1], [1, 1, 1], [1, 1, 1]],
        ....:  [[0]],
        ....:  [[0, 1, 1, 0, 1, 0],
        ....:   [1, 0, 0, 1, 0, 1],
        ....:   [1, 0, 0, 1, 0, 1],
        ....:   [0, 1, 1, 0, 1, 0],
        ....:   [1, 0, 0, 1, 0, 1],
        ....:   [0, 1, 1, 0, 1, 0]],
        ....:  [[0, 0], [0, 0]],
        ....:  [[0]],
        ....:  [[0, 0], [0, 0]],
        ....:  [[0]],
        ....:  [[1, 0, 0, 1, 0, 1],
        ....:   [0, 1, 1, 0, 1, 0],
        ....:   [0, 1, 1, 0, 1, 0],
        ....:   [1, 0, 0, 1, 0, 1],
        ....:   [0, 1, 1, 0, 1, 0],
        ....:   [1, 0, 0, 1, 0, 1]],
        ....:  [[0, 0], [0, 0]],
        ....:  [[0, 1, 1, 1, 0, 0],
        ....:   [1, 0, 0, 0, 1, 1],
        ....:   [1, 0, 0, 0, 1, 1],
        ....:   [1, 0, 0, 0, 1, 1],
        ....:   [0, 1, 1, 1, 0, 0],
        ....:   [0, 1, 1, 1, 0, 0]],
        ....:  [[0]],
        ....:  [[0]],
        ....:  [[0, 1, 0, 0], [1, 0, 0, 0], [0, 0, 0, 1], [0, 0, 1, 0]],
        ....:  [[0, 0, 1, 0, 1, 0, 1, 0, 0, 0, 0, 0],
        ....:   [0, 0, 0, 1, 0, 1, 0, 1, 0, 0, 0, 0],
        ....:   [1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0],
        ....:   [0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1],
        ....:   [1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0],
        ....:   [0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1],
        ....:   [1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0],
        ....:   [0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1],
        ....:   [0, 0, 1, 0, 1, 0, 1, 0, 0, 0, 0, 0],
        ....:   [0, 0, 0, 1, 0, 1, 0, 1, 0, 0, 0, 0],
        ....:   [0, 0, 1, 0, 1, 0, 1, 0, 0, 0, 0, 0],
        ....:   [0, 0, 0, 1, 0, 1, 0, 1, 0, 0, 0, 0]],
        ....:  [[0, 0, 0, 1], [0, 0, 1, 0], [0, 1, 0, 0], [1, 0, 0, 0]],
        ....:  [[0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0]],
        ....:  [[1, 0, 0, 1, 0, 1],
        ....:   [0, 1, 1, 0, 1, 0],
        ....:   [0, 1, 1, 0, 1, 0],
        ....:   [1, 0, 0, 1, 0, 1],
        ....:   [0, 1, 1, 0, 1, 0],
        ....:   [1, 0, 0, 1, 0, 1]],
        ....:  [[0, 0], [0, 0]],
        ....:  [[0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0]]]
        sage: matrL = [Matrix(c) for c in L]
        sage: matr == matrL # exact equality, was False at latest Sage version, long time
        False
        sage: t = [(Graph(matrL[c]).is_isomorphic(Graph(matr[c]), certificate=True), matr[c]) for c in range(len(L))] # long time
        sage: all([a[0] for (a,b) in t]) # here is what matters, long time
        True
        sage: matr = [Matrix([[c[a[1][i],a[1][j]] for j in range(c.ncols())] for i in range(c.nrows())]) for (a,c) in t] # long time
        sage: matr == matrL # now we do have the same graph as before (when forgetting the orders, so still doctests can fail at later stages in later versions of Sage). long time
        True
        sage: matr # long time
        [
        [0 1 1 1 0 0]                [0 1 1 0 1 0]
        [1 0 0 0 1 1]                [1 0 0 1 0 1]
        [1 0 0 0 1 1]                [1 0 0 1 0 1]
        [1 0 0 0 1 1]  [1 1 1]       [0 1 1 0 1 0]
        [0 1 1 1 0 0]  [1 1 1]       [1 0 0 1 0 1]  [0 0]       [0 0]
        [0 1 1 1 0 0], [1 1 1], [0], [0 1 1 0 1 0], [0 0], [0], [0 0], [0],
        <BLANKLINE>
        [1 0 0 1 0 1]         [0 1 1 1 0 0]
        [0 1 1 0 1 0]         [1 0 0 0 1 1]
        [0 1 1 0 1 0]         [1 0 0 0 1 1]            [0 1 0 0]
        [1 0 0 1 0 1]         [1 0 0 0 1 1]            [1 0 0 0]
        [0 1 1 0 1 0]  [0 0]  [0 1 1 1 0 0]            [0 0 0 1]
        [1 0 0 1 0 1], [0 0], [0 1 1 1 0 0], [0], [0], [0 0 1 0],
        <BLANKLINE>
        [0 0 1 0 1 0 1 0 0 0 0 0]
        [0 0 0 1 0 1 0 1 0 0 0 0]
        [1 0 0 0 0 0 0 0 1 0 1 0]
        [0 1 0 0 0 0 0 0 0 1 0 1]
        [1 0 0 0 0 0 0 0 1 0 1 0]
        [0 1 0 0 0 0 0 0 0 1 0 1]
        [1 0 0 0 0 0 0 0 1 0 1 0]                        [1 0 0 1 0 1]
        [0 1 0 0 0 0 0 0 0 1 0 1]                        [0 1 1 0 1 0]
        [0 0 1 0 1 0 1 0 0 0 0 0]  [0 0 0 1]  [0 0 0 0]  [0 1 1 0 1 0]
        [0 0 0 1 0 1 0 1 0 0 0 0]  [0 0 1 0]  [0 0 0 0]  [1 0 0 1 0 1]
        [0 0 1 0 1 0 1 0 0 0 0 0]  [0 1 0 0]  [0 0 0 0]  [0 1 1 0 1 0]  [0 0]
        [0 0 0 1 0 1 0 1 0 0 0 0], [1 0 0 0], [0 0 0 0], [1 0 0 1 0 1], [0 0],
        <BLANKLINE>
        [0 0 0 0]
        [0 0 0 0]
        [0 0 0 0]
        [0 0 0 0]
        ]

    The matrices above correspond to the (2,2)-isogeny graphs, one matrix for
    every order, one row/column for every period matrix. Next, we look at the
    orders corresponding to these period matrices. An entry a in row i and
    column j means that a is non-zero if and only if period matrices i and j
    have the same endomorphism ring, and a is the index of that endomorphism
    ring in the maximal order.::

        sage: vis = [visualize_orders(orders) for (K, Phi, F, l, orders) in Zs] # long time
        sage: flatten([[Matrix(lst[k]),Matrix([[vis[k][t[k][0][1][i],t[k][0][1][j]] for j in range(vis[k].ncols())] for i in range(vis[k].nrows())]), matr[k]] for k in range(len(matr))]) # long time
        [
                 [ 1  0  0  0  0  0]  [0 1 1 1 0 0]
                 [ 0  4  0  0  0  0]  [1 0 0 0 1 1]
                 [ 0  0  8  8  0  0]  [1 0 0 0 1 1]
                 [ 0  0  8  8  0  0]  [1 0 0 0 1 1]           [1 0 0]  [1 1 1]
                 [ 0  0  0  0 16 16]  [0 1 1 1 0 0]           [0 2 2]  [1 1 1]
        [5 5 5], [ 0  0  0  0 16 16], [0 1 1 1 0 0], [8 4 2], [0 2 2], [1 1 1],
        <BLANKLINE>
                                          [1 1 0 0 0 0]  [0 1 1 0 1 0]
                                          [1 1 0 0 0 0]  [1 0 0 1 0 1]
                                          [0 0 4 4 4 4]  [1 0 0 1 0 1]
                                          [0 0 4 4 4 4]  [0 1 1 0 1 0]
                                          [0 0 4 4 4 4]  [1 0 0 1 0 1]
        [13 13 13], [1], [0], [ 5 10 20], [0 0 4 4 4 4], [0 1 1 0 1 0],
        <BLANKLINE>
                       [1 1]  [0 0]
        [  5  65 845], [1 1], [0 0], [29 29 29], [1], [0], [   5   85 1445],
        <BLANKLINE>
                                                           [1 1 0 0 0 0]
                                                           [1 1 0 0 0 0]
                                                           [0 0 2 2 2 2]
                                                           [0 0 2 2 2 2]
        [1 1]  [0 0]                                       [0 0 2 2 2 2]
        [1 1], [0 0], [ 37  37 333], [1], [0], [ 8 20 50], [0 0 2 2 2 2],
        <BLANKLINE>
        [1 0 0 1 0 1]                                           [1 1 0 0 0 0]
        [0 1 1 0 1 0]                                           [1 1 0 0 0 0]
        [0 1 1 0 1 0]                                           [0 0 4 4 4 4]
        [1 0 0 1 0 1]                                           [0 0 4 4 4 4]
        [0 1 1 0 1 0]                 [1 1]  [0 0]              [0 0 4 4 4 4]
        [1 0 0 1 0 1], [ 13  65 325], [1 1], [0 0], [13 26 52], [0 0 4 4 4 4],
        <BLANKLINE>
        [0 1 1 1 0 0]
        [1 0 0 0 1 1]
        [1 0 0 0 1 1]
        [1 0 0 0 1 1]
        [0 1 1 1 0 0]
        [0 1 1 1 0 0], [53 53 53], [1], [0], [ 61  61 549], [1], [0],
        <BLANKLINE>
                    [1 1 1 1]  [0 1 0 0]
                    [1 1 1 1]  [1 0 0 0]
                    [1 1 1 1]  [0 0 0 1]
        [ 5 15 45], [1 1 1 1], [0 0 1 0], [  5  30 180],
        <BLANKLINE>
        [1 1 1 1 0 0 0 0 0 0 0 0]  [0 0 1 0 1 0 1 0 0 0 0 0]
        [1 1 1 1 0 0 0 0 0 0 0 0]  [0 0 0 1 0 1 0 1 0 0 0 0]
        [1 1 1 1 0 0 0 0 0 0 0 0]  [1 0 0 0 0 0 0 0 1 0 1 0]
        [1 1 1 1 0 0 0 0 0 0 0 0]  [0 1 0 0 0 0 0 0 0 1 0 1]
        [0 0 0 0 4 4 4 4 4 4 4 4]  [1 0 0 0 0 0 0 0 1 0 1 0]
        [0 0 0 0 4 4 4 4 4 4 4 4]  [0 1 0 0 0 0 0 0 0 1 0 1]
        [0 0 0 0 4 4 4 4 4 4 4 4]  [1 0 0 0 0 0 0 0 1 0 1 0]
        [0 0 0 0 4 4 4 4 4 4 4 4]  [0 1 0 0 0 0 0 0 0 1 0 1]
        [0 0 0 0 4 4 4 4 4 4 4 4]  [0 0 1 0 1 0 1 0 0 0 0 0]
        [0 0 0 0 4 4 4 4 4 4 4 4]  [0 0 0 1 0 1 0 1 0 0 0 0]
        [0 0 0 0 4 4 4 4 4 4 4 4]  [0 0 1 0 1 0 1 0 0 0 0 0]
        [0 0 0 0 4 4 4 4 4 4 4 4], [0 0 0 1 0 1 0 1 0 0 0 0], [  5  35 245],
        <BLANKLINE>
        [1 1 1 1]  [0 0 0 1]                    [1 1 1 1]  [0 0 0 0]
        [1 1 1 1]  [0 0 1 0]                    [1 1 1 1]  [0 0 0 0]
        [1 1 1 1]  [0 1 0 0]                    [1 1 1 1]  [0 0 0 0]
        [1 1 1 1], [1 0 0 0], [   5  105 2205], [1 1 1 1], [0 0 0 0],
        <BLANKLINE>
                    [1 1 0 0 0 0]  [1 0 0 1 0 1]
                    [1 1 0 0 0 0]  [0 1 1 0 1 0]
                    [0 0 2 2 2 2]  [0 1 1 0 1 0]
                    [0 0 2 2 2 2]  [1 0 0 1 0 1]
                    [0 0 2 2 2 2]  [0 1 1 0 1 0]                    [1 1]
        [ 8 12 18], [0 0 2 2 2 2], [1 0 0 1 0 1], [  17  119 3332], [1 1],
        <BLANKLINE>
                                    [1 1 1 1]  [0 0 0 0]
                                    [1 1 1 1]  [0 0 0 0]
        [0 0]                       [1 1 1 1]  [0 0 0 0]
        [0 0], [   17   255 15300], [1 1 1 1], [0 0 0 0]
        ]
        
    The first matrix is for `\QQ[\zeta_5]`, and from the others, we conclude
    Lemma 6.6. For `\QQ[\zeta_5]`, we see that the graph is connected.
    """    
    # the things below are a copy of another similar function
    ret = []
    for DAB in lst:
#        print DAB
        K = CM_Field(DAB)
        Phi = K.CM_types()[0]
        m = minimal_orders_of_class_number_one(Phi)
        if len(m):
            Omin = m[0][3]
            if m[0][0] != 2 or (DAB != [5,5,5] and len(m) > 1):
                raise RuntimeError, "tests fail!"
            F = m[0][0]**m[0][1]
        else:
            Omin = K.maximal_order()
            F = 1
        Zs = []
        orders = []
        sup = superorders_stable_under_complex_conjugation(Omin)
#        print "%s orders" % len(sup)
        for O in sup:
#            print "Computing period matrices for an order"
            Zs_this_order = period_matrices(O, F, Phi)
            Zs = Zs + Zs_this_order
            orders = orders + [O for k in range(len(Zs_this_order))]
        ret.append((K, Phi, F, Zs, orders))
    return ret
            
    
def two_two_matrix(Zs, F):
    """
    Returns an n x n matrix A (where n=len(Zs)) with
    A[k,l] = 1 if Zs[k] and Zs[l] are (2,2)-isogenous, and A[k,l] = 0 otherwise
    """
    n = len(Zs)
    A = zero_matrix(n)
    for k in range(n):
#        print "testing k = %s and l=0,...,%s" % (k,k)
        for l in range(k+1):
            entry = 1 if are_nn_isogenous(Zs[k], Zs[l], 2, F, F) else 0
            A[k,l] = entry
            A[l,k] = entry
#        print A[k][0:k+1]
    return A


def visualize_orders(orders):
    """
    Returns an n x n matrix A (where n=len(Zs)) with
    A[k,l] = [OK:orders[k]] if orders[k]=orders[l], and A[k,l] = 0 otherwise
    """
    n = len(orders)
    A = zero_matrix(n)
    for k in range(n):
        for l in range(k+1):
            if orders[k] == orders[l]:
                OK = orders[k].number_field().maximal_order()
                entry = orders[k].index_in(OK)
                A[k,l] = entry
                A[l,k] = entry
    return A
    

def is_S_O_equal_S_OK(O):
    """
    Given an order O in a non-biquadratic quartic CM-field O, returns True if
    and only if S_O = S_{O_K}, in the notation of Bisson-Streng.
    
    EXAMPLES:
    
    The same example as in is_trivial_in_shimura_group, but much faster::
    
        sage: from recip import *
        sage: K = CM_Field([149,13,5])
        sage: P.<x> = QQ[]
        sage: [alpha1,alpha2]=(x^4+2*x^3+16*x^2+15*x+19).roots(K, multiplicities=False)
        sage: O = K.order(alpha1)
        sage: is_S_O_equal_S_OK(O)
        False
        
    """
    K = O.number_field()
    OK = K.maximal_order()
    F = minimal_F(O)
    Phi = K.CM_types()[0]
    Kr = Phi.reflex_field()
    Psi = Phi.reflex()
    for g in Kr.ideal(F).idealstar(flag=2).gens():
        g = Kr.ideal(Kr(g))
        A = Psi.type_norm(g)
        alpha = g.norm()
        if not is_trivial_in_shimura_group(A, alpha, O):
            # These are all in S_OK, since they are obviously principal
            return False
    # Now we know that both S_O and S_OK contain the principal ideals,
    # so we can look at S_O/P_Kr \subset S_OK/P_Kr \subset Kr.class_group()
    # TODO: with linear algebra on the group Kr.class_group(), the following
    # can be made faster.
    idF = Kr.ideal(F)
    for c in Kr.class_group():
        I = c.ideal()
        I = I * I.idealcoprime(idF)
        A = Psi.type_norm(I)
        alpha = I.norm()
        if (is_trivial_in_shimura_group(A, alpha, OK) and not
                is_trivial_in_shimura_group(A, alpha, O)
                ):
            return False
    # Now S_O/P_Kr = S_OK/P_Kr.
    return True
                
                
"""

Some additional tests of orders.sage, made possible by the code in bissonstreng.sage::

        sage: from recip import *
        sage: lst = wamelen_dab_list()
        sage: Phi = CM_Field(lst[3]).CM_types()[0]
        sage: m = minimal_orders_of_class_number_one(Phi); m
        [(2, 1, 4, Order in CM Number Field in alpha with defining polynomial x^4 + 10*x^2 + 20)]
        sage: O = m[0][3]
        sage: polarized_ideal_classes(O, 2)
        [([1/2*alpha^2 + 1, alpha, alpha^2, alpha^3], 1/40*alpha), ([1/2*alpha^2 + 1, alpha, alpha^2, alpha^3], 1/40*alpha^3 + 3/40*alpha), ([1/2*alpha^2 + 5, alpha, alpha^2, alpha^3], 1/200*alpha^3 + 1/40*alpha), ([1/2*alpha^2 + 5, alpha, alpha^2, alpha^3], -1/100*alpha^3 - 1/40*alpha), ([1/2*alpha^3 + 1/2*alpha^2 + 1, alpha, alpha^2, alpha^3], 1/40*alpha), ([1/2*alpha^3 + 1/2*alpha^2 + 1, alpha, alpha^2, alpha^3], 1/40*alpha^3 + 3/40*alpha), ([1/2*alpha^3 + 1/2*alpha^2 + 5, alpha, alpha^2, alpha^3], 1/200*alpha^3 + 1/40*alpha), ([1/2*alpha^3 + 1/2*alpha^2 + 5, alpha, alpha^2, alpha^3], -1/100*alpha^3 - 1/40*alpha)]

"""



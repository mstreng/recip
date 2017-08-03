"""
RECIP -- REpository of Complex multIPlication SageMath code.

this file:
bls.sage

This file gives the details of computations for Table 2 of [BLS].

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

::

    sage: from recip import *
    sage: lst = recognize_all_from_article(3, print_results=True)
    0 -7
    (1, [(x^2 - x + 2, -7), (x^2 - x + 2, -7)])
    0 -6
    (1, [(x^2 - x + 1, -3), (x^2 + 3, -12)])
    (2, [(x^2 + 2, -8), (x^2 + 2, -8)])
    0 -5
    (1, [(x^2 - x + 1, -3), (x^2 - x + 1, -3)])
    (1, [(2*x^2 - x + 2, -15), (x^2 - x + 4, -15)])
    (2, [(x^2 - x + 1, -3), (x^2 - x + 1, -3)])
    (2, [(2*x^2 - x + 2, -15), (x^2 - x + 4, -15)])
    1 -5
    (1, [(x^2 - x + 1, -3), (x^2 - x + 1, -3)])
    0 -3
    (1, [(2*x^2 - x + 2, -15), (x^2 - x + 4, -15)])
    (2, [(2*x^2 - 2*x + 3, -20), (x^2 + 5, -20)])
    (2, [(2*x^2 - 2*x + 3, -20), (x^2 + 5, -20)])
    (2, [(2*x^2 - x + 2, -15), (x^2 - x + 4, -15)])
    0 -2
    (1, [(x^2 - x + 1, -3), (x^2 - x + 1, -3)])
    (1, [(x^2 - x + 1, -3), (x^2 + 3, -12)])
    (2, [(x^2 + 6, -24), (x^2 + 6, -24)])
    (2, [(2*x^2 + 3, -24), (x^2 + 6, -24)])
    2 -2
    (1, [(x^2 - x + 1, -3), (x^2 - x + 1, -3)])
    (1, [(x^2 + 1, -4), (x^2 + 1, -4)])
    3 1
    (1, [(x^2 - x + 1, -3), (x^2 - x + 1, -3)])
    (2, [(2*x^2 - x + 2, -15), (x^2 - x + 4, -15)])

So all of them are polarized products that we can write down easily
or are (2,2)-isogenous to such. In the latter case, we can list the
possibilities and do numerical things. And maybe sometimes better than
numerical things: how about degrees of fields of moduli?

    sage: Z = lst[1][1][1][0]
    sage: j = igusa_invariants_absolute()[2]
    sage: j_val = j(Z.complex_matrix().base_extend(CIF), interval=True); j_val
    -50000.000000? + 0.?e-6*I
    sage: K.<sqrt2>=QuadraticField(2)
    sage: a = (j_from_a()-8000).numerator().roots(K)[0][0]
    sage: a
    2*sqrt2 + 3
    sage: P.<y>=K[]
    sage: Q = P.fraction_field()
    sage: p = (Q(j_from_a())-8000).numerator()
    sage: l = [b for b in p.roots(multiplicities=False) if b != a]
    sage: [out] = [(b, ic_from_Cab(a,b)) for b in l if (j_val-CIF(j_from_Cab(a, b))).contains_zero()]
    sage: out
    (-2*sqrt2 + 3, (716800*sqrt2 - 1013760, 72666316800*sqrt2 - 102765690880, -14732865254195200*sqrt2 + 20835417855098880, 30280661510083082557849600*sqrt2 - 42823322185188459983929344))
    sage: mult_ic(reduce_ic(out[1]),sqrt2)
    [20, -20, -40, 8]
    sage: ic_on_humb8(_)
    True

So the first table entry that is not a polarized product would be something
like "C_{-8}, C(3+2*sqrt2, 3-2*sqrt2), (20: -20: -40: 8)".

The next is a bit more complicated, and appears for two fields::

    sage: Z = lst[2][1][3][0]
    sage: j_val = j(Z.complex_matrix().base_extend(CIF), interval=True); j_val
    3.9106607...?e6 + 0.00000?*I
    sage: Z2 = lst[4][1][3][0]
    sage: Z.complex_matrix(QQbar) == Z2.complex_matrix(QQbar)
    True
    sage: Z3 = lst[7][1][1][0]
    sage: Z.complex_matrix(QQbar) == Z3.complex_matrix(QQbar)
    True
    sage: K.<sqrt5>=QuadraticField(5)
    sage: [j1, j2] = hilbert_class_polynomial(-15).roots(K, multiplicities=False)
    sage: P.<x> = K[]
    sage: L.<sqrtminus3> = NumberField(x^2+3)
    sage: Q = L['x'].fraction_field()
    sage: a = (Q(j_from_a())-j1).numerator().roots()[0][0]
    sage: a
    (-7/2*sqrt5 + 8)*sqrtminus3 + 1/2
    sage: p = (Q(j_from_a())-j2).numerator()
    sage: l = [b for b in p.roots(multiplicities=False) if b != a]
    sage: [out] = [(b, ic_from_Cab(a,b)) for b in l if j_from_Cab(a,b) in QQ and (j_val-CIF(j_from_Cab(a,b))).contains_zero()]
    sage: out
    ((-7/2*sqrt5 - 8)*sqrtminus3 + 1/2, (-3790080*sqrt5 - 8474880, 36135532339200*sqrt5 + 80801506713600, -161288631717325701120*sqrt5 - 360652344517868912640, 37539075121066454177959882260480*sqrt5 + 83939923783175739234512919330816))
    sage: out = mult_ic(out[1], 20/out[1][0])
    sage: out
    [20, 225, 1185, -384]
    sage: ic_on_humb8(_)
    False

So the next such entry is
"C_{-15}, C((-7/2*sqrt5 + 8)*sqrtminus3 + 1/2, (-7/2*sqrt5 - 8)*sqrtminus3 + 1/2),
[20, 225, 1185, -384]", and it appears not only for this field (0, -5), but also
for the fields (0, -3) and (3, 1).

There is another one for this field.::

    sage: Z = lst[2][1][2][0]
    sage: j = igusa_invariants_absolute()[2]
    sage: j_val = j(Z.complex_matrix().base_extend(CIF), interval=True); j_val
    5.12578125000?e6 + 0.?e-7*I
    sage: K.<zeta> = CyclotomicField(6)
    sage: j_from_a()(zeta) == 0
    True
    sage: a = zeta
    sage: [b] = [b for b in j_from_a().numerator().roots(K, multiplicities=False) if b != a]
    sage: b == zeta**-1
    True
    sage: reduce_ic(ic_from_Cab(a, b))
    [40, 45, 555, 6]
    sage: ic_on_humb8(_)
    True

So the next such entry is
"C_{-3}, C(zeta_6, zeta_6^{-1}), [40, 45, 555, 6]", and that finishes this
field (0, -5).

For the next field (0, -3), we have already done one curve. The other
is a pair of curves, which may be the most difficult.

    sage: Z1 = lst[4][1][1][0]
    sage: j_val1 = j(Z1.complex_matrix().base_extend(CIF), interval=True); j_val1 # the following was + in an older version of Sage, and should be different from the one after, so something is wrong.
    8.1894000000?e7 + 1.50920000000?e7*I
    sage: Z2 = lst[4][1][2][0]
    sage: j_val2 = j(Z2.complex_matrix().base_extend(CIF), interval=True); j_val2
    8.1894000000?e7 + 1.50920000000?e7*I
    sage: K.<sqrt5>=QuadraticField(5)
    sage: [j1, j2] = hilbert_class_polynomial(-20).roots(K, multiplicities=False)
    sage: P.<x> = K[]
    sage: L.<twoa> = NumberField(x^2 - 2*x - 4*sqrt5 + 9)
    sage: a = twoa/2
    sage: Q.<y> = L[]
    sage: M.<c> = NumberField(y^2 - 2*y + 4*sqrt5 + 9)
    sage: j_from_a()(a) == j1
    True
    sage: p = (P.fraction_field()(j_from_a())-j2).numerator()
    sage: l = [b for b in p.roots(M, multiplicities=False) if b != a]
    sage: len(l) == 6 # so we have found all candidate b's
    True
    sage: out = [b for b in l if j_from_Cab(a,b).absolute_minpoly().degree() <= 2]
    sage: len(out) # swapping j1 and j2 does nothing, so here are exactly our two curves C:
    2
    sage: Mabs.<d> = M.absolute_field()
    sage: N.<i> = QuadraticField(-1)
    sage: phi = N.embeddings(Mabs)[1]
    sage: Mrel = Mabs.relativize(phi, 'alpha')
    sage: phi = Mrel.structure()[1] * Mabs.structure()[1]
    sage: phi(j_from_Cab(a, out[0]))
    15092000*i + 81894000
    sage: ic = ic_from_Cab(a,out[0])
    sage: js = [ic[0]^5/ic[3], ic[0]^3*ic[1]/ic[3], ic[0]^2*ic[2]/ic[3]]
    sage: js = [phi(ji) for ji in js]
    sage: all([ji == ji[0] for ji in js])
    True
    sage: js = [ji[0] for ji in js]
    sage: ic = [N(1),0,0,0]
    sage: ic[3] = js[0]**-1
    sage: ic[1] = ic[3]*js[1]
    sage: ic[2] = ic[3]*js[2]
    sage: mult_ic(reduce_ic(ic), i)
    [46*i + 42, 210*i - 220, 50*i - 4350, -73*i + 161]
    sage: ic_on_humb8(_)
    True
    sage: out[0].minpoly()
    x^2 - x + sqrt5 + 9/4
    sage: out[1].minpoly() == out[0].minpoly()
    True
    sage: out[1] == out[0]
    False
    sage: a.minpoly()
    x^2 - x - sqrt5 + 9/4

So this entry is "C_{-20}^i, C((1+2sqrt(-sqrt5-2))/2, (1+2sqrt(sqrt5-2))/2),
[46*i + 42, 210*i - 220, 50*i - 4350, -73*i + 161]", where the
signs of the square root of 5 must be picked consistent, and the other
square roots can be varied, giving two distinct curves.
This finishes (0, -3), except that the Igusa-Clebsch invariants are different
from what is written in the preprint. That was due to a sign error in
an absolute invariant that fortunately was saved in the comments, and
can now be fixed.

Now only for the field (0, -2), something needs to be computed.

    sage: Z1 = lst[5][1][3][0]
    sage: j_val1 = j(Z1.complex_matrix().base_extend(CIF), interval=True); j_val1
    1.76433163200?e9 + 0.00000?*I
    sage: Z2 = lst[5][1][2][0]
    sage: j_val2 = j(Z2.complex_matrix().base_extend(CIF), interval=True); j_val2
    2.55091680000?e7 + 0.0000?*I
    sage: K.<sqrt2>=QuadraticField(2)
    sage: [j1, j2] = hilbert_class_polynomial(-24).roots(K, multiplicities=False)
    sage: P.<x> = K[]
    sage: L.<twoa> = NumberField(x^2 - 2*x + 12*sqrt2 -17)
    sage: a = twoa/2
    sage: j_from_a()(a) == j1
    True
    sage: p = (P.fraction_field()(j_from_a())-j2).numerator()
    sage: l = [b for b in p.roots(L, multiplicities=False) if b != a]
    sage: len(l) == 6 # so we have found all candidate b's
    True
    sage: [b] = [b for b in l if j_from_Cab(a,b).absolute_minpoly()(j_val1).contains_zero()]
    sage: b
    (sqrt2 + 3/2)*twoa - sqrt2 - 1
    sage: all([c == c[0] for c in ic_from_Cab(a, b)])
    True
    sage: ic = [c[0] for c in ic_from_Cab(a, b)]
    sage: mult_ic(reduce_ic(ic), sqrt2)
    [76, 252, 5160, 24]
    sage: ic_on_humb8(_)
    True
    sage: b.minpoly()
    x^2 - x - 3*sqrt2 - 17/4
    sage: a.minpoly()
    x^2 - x + 3*sqrt2 - 17/4

So we found the first one, C_{-6}^2 in the preprint, also
C(a, 2*sqrt2*a+3*a-sqrt2-1) where a^2 - 2*a + 3*sqrt2 - 17 = 0,
also [76, 252, 5160, 24]. These a and be can also be denoted
in notation very similar to the previous example, except that we have
to be very careful about signs.

As for the other one::

    sage: p = (P.fraction_field()(j_from_a())-j1).numerator()
    sage: l = [b for b in p.roots(L, multiplicities=False) if b != a]
    sage: len(l) == 5 # so we have found all candidate b's
    True
    sage: [b] = [b for b in l if j_from_Cab(a,b).absolute_minpoly()(j_val2).contains_zero()]
    sage: b
    -1/2*twoa + 1
    sage: all([c == c[0] for c in ic_from_Cab(a, b)])
    True
    sage: ic = [c[0] for c in ic_from_Cab(a, b)]
    sage: mult_ic(reduce_ic(ic), sqrt2)
    [-92, 108, -4104, -24]
    sage: ic_on_humb8(_)
    True
    sage: b.minpoly()
    x^2 - x + 3*sqrt2 - 17/4
    sage: a.minpoly()
    x^2 - x + 3*sqrt2 - 17/4

This time a and b are conjugates over K. We get the final curve,
C_{-6}^1 with [92, 108, 4104, 24].

For all the curves in Table 2, we now have proofs and models.
In particular, for all the curves in Table 1, we now have proofs,
up to a correction for one of the entries, and we had
models at some stage. For some curves, we forgot the models:
C_{-7}, C_{-4, -7}, C_{-4, -8}, but the proof went via computing such a model,
and we can recompute it easily using the methods above.

    sage: K.<sqrtmin7> = QuadraticField(-7)
    sage: l = (j_from_a()+3375).numerator().roots(K, multiplicities=False)
    sage: l.sort()
    sage: [CC(b) for b in l]
    [0.0312500000000000 - 0.248039185412305*I, 0.0312500000000000 + 0.248039185412305*I, 0.500000000000000 - 3.96862696659689*I, 0.500000000000000 + 3.96862696659689*I, 0.968750000000000 - 0.248039185412305*I, 0.968750000000000 + 0.248039185412305*I]
    sage: len(l)
    6

We use the following pari/gp code to select the a from the article (the
one for which we proved b=1/a)
? g(z) = ellwp([1, (sqrt(-7)+1)/2], z)
? a = (g((sqrt(-7)+1)/4 + 1/2) - g(1/2)) / (g((sqrt(-7)+1)/4) - g(1/2))
%2 = 0.96875000000000000000000000000000000003 - 0.24803918541230536785952647690368066482*I
So here they are::

    sage: a = l[4]
    sage: b = 1/a
    sage: b in l
    True
    sage: reduce_ic(ic_from_Cab(a, b))
    [10840, 2004345, 7846230105, 131736761856]
    sage: ic_on_humb8(_)
    True
    sage: a
    -3/32*sqrtmin7 + 31/32
    sage: b
    3/32*sqrtmin7 + 31/32
    sage: a.minpoly()
    x^2 - 31/16*x + 1
    sage: b.minpoly()
    x^2 - 31/16*x + 1

Similarly for Z[i] and Z[sqrt-2], we do
? beta1 = -1 + sqrt(-1)
%4 = -1 + 1.0000000000000000000000000000000000000*I
? beta2 = sqrt(-2)
? a(tau) = (ellwp([1, tau], (tau+1)/2) - ellwp([1,tau],1/2)) / (ellwp([1,tau],tau/2) - ellwp([1,tau],1/2))
%15 = (tau)->(ellwp([1,tau],(tau+1)/2)-ellwp([1,tau],1/2))/(ellwp([1,tau],tau/2)-ellwp([1,tau],1/2))
? a(beta1)
%16 = 2.0000000000000000000000000000000000000 + 0.E-37*I
? a(beta2)
%17 = 0.82842712474619009760337744841939615860 + 0.E-38*I
? algdep(a(beta2),2)
%19 = x^2 + 4*x - 4

    sage: (j_from_a()-1728).numerator().factor()
    (4096) * (x - 2)^2 * (x - 1/2)^2 * (x + 1)^2
    sage: hilbert_class_polynomial(-8)
    x - 8000
    sage: (j_from_a()-8000).numerator().factor()
    (4096) * (x^2 - 6*x + 1) * (x^2 - x - 1/4) * (x^2 + 4*x - 4)

So we take a=2 and take b to be a root of x^2+4x-4 and get

    sage: x = polygen(QQ)
    sage: K.<b> = NumberField(x^2+4*x-4)
    sage: reduce_ic(ic_from_Cab(2, b))
    [24, 30, 366, 2]
    sage: ic_on_humb8(_)
    False

This proves C_{-4, -8}, and the final curve C_{-4, -7} comes from a
(3,3)-isogeny, hence is not of this form.

As for whether the points are on the Humbert surface::

    sage: ic_on_humb8([40,45,555,6])
    True
    sage: ic_on_humb8([92,108,4104,24])
    True
    sage: ic_on_humb8([76,252,5160,24])
    True
    sage: ic_on_humb8([10840,2004345,7846230105,131736761856])
    True
    sage: ic_on_humb8([20,-20,-40,8])
    True
    sage: ic_on_humb8([20,225,1185,-384])
    False
    sage: K.<i> = QuadraticField(-1)
    sage: ic_on_humb8([46*i+42,210*i-220,50*i-4350,-73*i+161])
    True
    sage: K.<a> = QuadraticField(7)
    sage: ic_on_humb8([8+20*a,-1035,-450*a,87246+33606*a,25164+9504*a])
    False
    sage: ic_on_humb8([24,30,366,2])
    False

For the last two, we did not need a computation, the endomorphism ring is
contained in QQ(sqrt(-1)) x QQ(sqrt(D)) with D = -7, -8, which does not
contain a square root of 2.

"""




def cosets_of_GammaN_iter(p):
    """
    Returns exactly one representative of every left coset of Gamma_0(N)
    in Sp_4(ZZ).
    """
    if not p.is_prime():
        raise ValueError
    F = srange(p)
    for a in F:
        for b in F:
            for c in F:
                yield GSp_element([[1,0,0,0],[0,1,0,0],[a,b,1,0],[b,c,0,1]])
    yield GSp_element([[0,0,1,0],[0,0,0,1],[-1,0,0,0],[0,-1,0,0]])
    for c in F:
        yield GSp_element([[0,0,0,1],[1,0,0,0],[0,-1,0,0],[c,0,1,0]])
        for a in F:
            yield GSp_element([[1,0,0,0],[a,0,0,1],[c,a,1,0],[0,-1,0,0]])        


def nn_isogenous_matrices_iter(Z, p):
    """
    Returns all Sp_4(ZZ)-classes of period matrices that are (p,p)-isogenous to
    Z.
    """
    for M in cosets_of_GammaN_iter(p):
#        print "."
#        sys.stdout.flush()
        yield (Z.Sp_action(M)*p).reduce()


def nn_isogenous_matrices_iter2(Z, p):
    """
    Returns all Sp_4(ZZ)-classes of period matrices that are (p,p)-isogenous to
    Z.
    """
    for M in cosets_of_GammaN_iter(p):
#        print "."
#        sys.stdout.flush()
        yield _reduce(Sp_action(M.matrix(), Z)*p)[0]

def count_nn_endomorphisms(Z, p):
    """
    Gives a lower bound (hopefully tight) on the number of (p,p)-endomorphisms
    for the input prime p.
    """
    return len([V for V in nn_isogenous_matrices_iter2(Z, p) if V == Z])

def is_polarized_product(Z, info=False, reduce=True):
    """
    Returns True or False. If True, then we know Z is a polarized product
    of CM elliptic curves. If False, then likely it is not (no proof!).
    If info is True, then returns information on the elliptic curves.
    If reduce is True, then reduce first (use reduce=False if already reduced).
    """
    if reduce:
        Z = Z.reduce()
    if Z[0][1] != 0 or Z[1][0] != 0:
        if info:
            return False, None
        return False
    if info:
        pols = [Z[i][i].minpoly() for i in [0,1]]
        pols = [p*p.denominator() for p in pols]
        discs = [p.discriminant() for p in pols]
        return True, zip(pols, discs)
    return True


def recognize_matrix(Z, bound=2):
    """
    Tries to write Z as a polarized product of CM elliptic curves.
    If that fails, tries this for all (p,p)-isogenous varieties
    for p below the given bound.
    """
    for p in [1] + prime_range(bound):
        for V in [Z] if p == 1 else nn_isogenous_matrices_iter(Z, p):
            b, i = is_polarized_product(V, info=True, reduce=False)
            if b:
                return p, i
    return False, None


def recognize_all_from_article(bound=2, print_results=False):
    """
    Returns a list of pairs ((a1,a2),l).
    Here a1, a2 runs over the numbers a1 and a2 from the article, and each
    l is itself a list.
    The elements of l are triples (Z, o0, o1), where the Z are period matrices.
    """
    X = polygen(QQ)
    a1a2list = [(0,-7), (0,-6), (0,-5), (1,-5), (0,-3), (0,-2), (2,-2), (3,1)]
    ret = []
    for (a1,a2) in a1a2list:
        if print_results:
            print a1, a2
            sys.stdout.flush()
        f = X^4-a1*X^3+(4+a2)*X^2-2*a1*X+4
        K = CM_Field(f)
        alpha = K.gen()
        alphabar = K.complex_conjugation()(alpha)
        p =  (alpha+alphabar).minpoly()
        assert p.degree() == 2 and p.discriminant() > 0
        O = K.order([alpha, alphabar])
        s = superorders_stable_under_complex_conjugation(O)
        if not len(s) <= 2:
            raise NotImplementedError
        if len(s) == 2:
            if not s[1].index_in(s[0]) in [2,4,8]:
                raise NotImplementedError
        l = []
        for A in s:
            for Z in period_matrices(A, 2, reduced=False):
                Z = Z.reduce()
                o = recognize_matrix(Z, bound)
                l.append((o[0], o[1], Z))
        l.sort()
        if print_results:
            for (o0, o1, Z) in l:
                print (o0, o1)
                sys.stdout.flush()
        l = [(Z,o0,o1) for (o0,o1,Z) in l] # restructure for backwards compatibility
        ret.append(((a1,a2),l))
    return ret
        

def j_from_a():
    a = polygen(QQ)
    return (4096*a^6 - 12288*a^5 + 24576*a^4 - 28672*a^3 + 24576*a^2 - 12288*a + 4096)/(16*a^4 - 32*a^3 + 16*a^2)


def ic_from_Cab(a, b):
    x = polygen(Sequence([a,b]).universe())
    H = HyperellipticCurve((x^2-1)*(x^2-b/a)*(x^2-(b-1)/(a-1)))
    return H.igusa_clebsch_invariants()


def j_from_Cab(a, b):
    ic = ic_from_Cab(a, b)
    return ic[1]^5/ic[3]^2


def ideal_power(I, n):
    """
    Returns I^n if I is principal.
    (because I**n is too slow in Sage 5.4.1)
    """
    return I.number_field().ideal(I.gens_reduced()[0]**n)


def lower_root(ideal, n):
    return prod([ideal_power(p, floor(e/n)) for (p,e) in ideal.factor()])

def mult_ic(ic, d):
    return [ic[0]*d, ic[1]*d^2, ic[2]*d^3, ic[3]*d^5]

def reduce_ic(ic):
    K = ic[0].parent()
    O = K.maximal_order()
    id = lower_root(ic[0]*O, 1)+lower_root(ic[1]*O,2)+lower_root(ic[2]*O,3)+lower_root(ic[3]*O,5)
    if id.is_principal():
        d = id.gens_reduced()[0]
        ic = [ic[0]/d, ic[1]/d^2, ic[2]/d^3, ic[3]/d^5]
        d = ic[0].factor().unit()
        ic = [ic[0]/d, ic[1]/d^2, ic[2]/d^3, ic[3]/d^5]
        return ic
    raise NotImplementedError


def humbert8():
    """
    Gives an equation for the Humbert surface of discriminant 8, due to
    David Gruenewald.
    
    See also http://echidna.maths.usyd.edu.au/~davidg/thesis.html
    which contains this surface with respect to other invariants.
    """
    P.<i1,i2,i3> = QQ[]
    return (-12008187649867604184*i1^10 + 46966009581837720*i1^9*i2^2 -
    281796057491026320*i1^9*i2*i3 + 74323144140993953726400*i1^9*i2 +
    422694086236539480*i1^9*i3^2 - 188308985155317427104000*i1^9*i3 +
    6525527376491882166412800000*i1^9 - 40922652296790*i1^8*i2^4 +
    491071827561480*i1^8*i2^3*i3 + 80797331265418968300*i1^8*i2^3 -
    2209823224026660*i1^8*i2^2*i3^2 - 600817177103562872352*i1^8*i2^2*i3 -
    5646435963973808150448000*i1^8*i2^2 + 4419646448053320*i1^8*i2*i3^3 +
    1464035661952137007200*i1^8*i2*i3^2 + 11816513613940876191360000*i1^8*i2*i3

    - 3314734836039990*i1^8*i3^4 - 1164270951989800848000*i1^8*i3^3 -
    9692302860*i1^7*i2^6 + 174461451480*i1^7*i2^5*i3 +
    13070519254692720*i1^7*i2^5 - 1308460886100*i1^7*i2^4*i3^2 -
    166191670052927424*i1^7*i2^4*i3 - 151905897476921833920*i1^7*i2^4 +
    5233843544400*i1^7*i2^3*i3^3 + 828167627209721904*i1^7*i2^3*i3^2 +
    935454732407866725120*i1^7*i2^3*i3 - 11776147974900*i1^7*i2^2*i3^4 -
    2008024036434920736*i1^7*i2^2*i3^3 - 1442710323979397280000*i1^7*i2^2*i3^2 +

    14131377569880*i1^7*i2*i3^5 + 2342758351512643200*i1^7*i2*i3^4 -
    7065688784940*i1^7*i3^6 - 1031195565889344000*i1^7*i3^5 - 690120*i1^6*i2^8 +

    16562880*i1^6*i2^7*i3 + 768285653265*i1^6*i2^7 - 173910240*i1^6*i2^6*i3^2 -

    13581227720916*i1^6*i2^6*i3 + 54907566009968574*i1^6*i2^6 +
    1043461440*i1^6*i2^5*i3^3 + 101234335327542*i1^6*i2^5*i3^2 -
    673276229026228512*i1^6*i2^5*i3 + 3037457591848026771264*i1^6*i2^5 -
    3912980400*i1^6*i2^4*i3^4 - 407375781811164*i1^6*i2^4*i3^3 +
    2473553295898899120*i1^6*i2^4*i3^2 + 9391152960*i1^6*i2^3*i3^5 +
    933192344543553*i1^6*i2^3*i3^4 - 2847845695533497280*i1^6*i2^3*i3^3 -
    14086729440*i1^6*i2^2*i3^6 - 1153079429223816*i1^6*i2^2*i3^5 +
    12074339520*i1^6*i2*i3^7 + 599830229102400*i1^6*i2*i3^6 -
    4527877320*i1^6*i3^8 - 16*i1^5*i2^10 + 480*i1^5*i2^9*i3 + 16407360*i1^5*i2^9

    - 6480*i1^5*i2^8*i3^2 - 370821888*i1^5*i2^8*i3 + 3370305406482*i1^5*i2^8 +
    51840*i1^5*i2^7*i3^3 + 3573802944*i1^5*i2^7*i3^2 -
    37884332178108*i1^5*i2^7*i3 - 219088920070495152*i1^5*i2^7 -
    272160*i1^5*i2^6*i3^4 - 19051044480*i1^5*i2^6*i3^3 +
    157074803277714*i1^5*i2^6*i3^2 + 513572895493723584*i1^5*i2^6*i3 +
    979776*i1^5*i2^5*i3^5 + 60699222720*i1^5*i2^5*i3^4 -
    283564834757424*i1^5*i2^5*i3^3 - 2449440*i1^5*i2^4*i3^6 -
    115641561600*i1^5*i2^4*i3^5 + 186903505656720*i1^5*i2^4*i3^4 +
    4199040*i1^5*i2^3*i3^7 + 122024522304*i1^5*i2^3*i3^6 -
    4723920*i1^5*i2^2*i3^8 - 55031778432*i1^5*i2^2*i3^7 + 3149280*i1^5*i2*i3^9 -

    944784*i1^5*i3^10 + 80*i1^4*i2^11 - 1920*i1^4*i2^10*i3 - 4770576*i1^4*i2^10

    + 20160*i1^4*i2^9*i3^2 + 122472000*i1^4*i2^9*i3 + 7836835007961*i1^4*i2^9 -

    120960*i1^4*i2^8*i3^3 - 1040312160*i1^4*i2^8*i3^2 -
    46478808076104*i1^4*i2^8*i3 + 453600*i1^4*i2^7*i3^4 +
    4037376960*i1^4*i2^7*i3^3 + 69880444045344*i1^4*i2^7*i3^2 -
    1088640*i1^4*i2^6*i3^5 - 7430726160*i1^4*i2^6*i3^4 + 1632960*i1^4*i2^5*i3^6

    + 5283232128*i1^4*i2^5*i3^5 - 1399680*i1^4*i2^4*i3^7 + 524880*i1^4*i2^3*i3^8

    - 160*i1^3*i2^12 + 2880*i1^3*i2^11*i3 - 36920448*i1^3*i2^11 -
    21600*i1^3*i2^10*i3^2 + 317540736*i1^3*i2^10*i3 - 446201642736*i1^3*i2^10 +

    86400*i1^3*i2^9*i3^3 - 908392320*i1^3*i2^9*i3^2 - 194400*i1^3*i2^8*i3^4 +
    864162432*i1^3*i2^8*i3^3 + 233280*i1^3*i2^7*i3^5 - 116640*i1^3*i2^6*i3^6 +
    160*i1^2*i2^13 - 1920*i1^2*i2^12*i3 + 25973784*i1^2*i2^12 +
    8640*i1^2*i2^11*i3^2 - 85753728*i1^2*i2^11*i3 - 17280*i1^2*i2^10*i3^3 +
    12960*i1^2*i2^9*i3^4 - 80*i1*i2^14 + 480*i1*i2^13*i3 - 720*i1*i2^12*i3^2 +
    16*i2^15)

def ic_on_humb8(ic):
    return humbert8()([ic[0]^5/ic[3], ic[1]*ic[0]^3/ic[3], ic[2]*ic[0]^2/ic[3]]) == 0

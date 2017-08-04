"""
RECIP -- REpository of Complex multIPlication SageMath code.

See the file README.txt for version information and instructions.

#*****************************************************************************
# Copyright (C) 2010,2011,2012,2013,2016 Marco Streng <marco.streng@gmail.com>
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

This file implements the theta transformation formula and a few different
approaches to a modular-form-in-terms-of-theta class.

For the theta transformation formula, we take
[BL, Formula 8.6.1 and Exercise 8.11(9)],
but restrict to the case of theta constants (v=0) and principal
polarizations (delta=1). We use delta when [BL] use D
and A, B, C, D when [BL] use greek letters.

For v = 0, M = ((A B) (C D)) in Sp_2g(Z), z in H_g and c in R^(2g),
it reads
  theta[M[c]](0,M(Z)) = kappa(M) det(CZ+D)^(1/2) e(k(M,c)/2) theta[c](0,Z).
  
Here kappa(M) is defined up the the choice of square root of det(CZ+D).
We will define kappa(M), k(M,c), and M[c] below.

[BL] -- Birkenhake-Lange, Complex abelian varieties
"""

from sage.structure.sage_object import SageObject
from sage.rings.number_field.number_field import CyclotomicField
from sage.structure.unique_representation import UniqueRepresentation 
from sage.structure.element import CommutativeRingElement
from sage.groups.group import Group
from sage.structure.element import Element
from sage.rings.fraction_field import is_FractionField
from sage.structure.element import is_RingElement

def subscript_zero(A):
    """
    Given A in Mat_g(Z), returns A_0 in Z^g as defined in [BL],
    which is the diagonal of A.
    
    EXAMPLES::
    
        sage: from recip import *
        sage: A = Matrix([[1,2],[3,4]])
        sage: subscript_zero(A)
        (1, 4)
    """
    return vector(diag(A))


def is_den_even(den):
    """
    All code assumes den to be even. This function
    raises an error if den is odd, and returns nothing
    otherwise.
    
    EXAMPLES::
    
        sage: from recip import *
        sage: is_den_even(5)
        Traceback (most recent call last):
        ...
        ValueError: The integer den (=5) must be even.
        sage: is_den_even(8)
        sage: is_den_even(4)
        
    """
    if den % 2 == 1:
        raise ValueError, "The integer den (=%s) must be even." % den


def sq_brackets_inverse(M, nu_inv, c):
    """
    Given M in GSp_2g(Z/NZ) and c in (1/den)ZZ^{2g}, returns d as in my Shimura
    reciprocity article. If M is in Sp_2g, then M[d^1,d^2]^i == c^i in the
    notation of Birkenhake-Lange.
    
    INPUT:
    
    - M -- matrix in GSp_2g(Z/NZ), lifted to ZZ
    - nu_inv -- inverse of nu(M), lifted to ZZ
    - c -- vector in (1/den)Z^{2g}
    
    ..NOTE::
    
    The input is not checked, make sure that the congruences are satisfied when
    using this function.
    
    EXAMPLES::
    
        sage: from recip import *
        sage: M = Matrix([[1,0,2,0],[0,1,0,0],[-2,0,-3,0],[0,0,0,1]])
        sage: c = vector([1/2, 1/2, 0, 0])
        sage: sq_brackets_inverse(M, 1, c)
        (-1/2, 1/2, -2, 0)
    """
    A,B,C,D = ABCD(M)
    g = A.ncols()
    d = M.transpose() * (vector(c) - 1/2 * nu_inv * \
            vector(Sequence(subscript_zero(C*D.transpose())) + \
                   Sequence(subscript_zero(A*B.transpose()))))
    return d


def kappa(M, k = 1):
    """
    Given M in Sp_2g(Z) and k in ZZ, returns kappa(M)^k in <zeta_8^{gk}>.
    Note that the output is defined only up to sign if k is odd.

    The output is defined by
      M --> kappa(M)^2 is a group homomorphism
      kappa((0 -1) (1 0))^2 = (-i)^g
      kappa((1 B) (0 1))^2 = 1
      kappa((A 0) (0 D))^2 = det(A) in {-1,1}
      
    EXAMPLE::
    
        sage: from recip import *
        sage: kappa(Matrix([[0,5,0,1],[-1,0,-2,0],[1,4,2,1],[1,-5,1,-1]]))
        Traceback (most recent call last):
        ...
        NotImplementedError: sorry, kappa(M, k) only implemented in trivial cases
    """
    g = ZZ(M.ncols()/2)
    if k*g % 8 == 0:
        return 1
    if M.base_ring() != ZZ:
        mu = M.base_ring().order()
        if not (mu % 8 == 0 and M.base_ring() == Zmod(mu)):
            raise ValueError, "M (=%s) must be well-defined over Zmod(8), but is defined over %s" % (M, M.base_ring())
    A,B,C,D = ABCD(M)
    if C == 0 and k % 2 == 0:
        return A.determinant()^(ZZ(k/2))
    raise NotImplementedError, "sorry, kappa(M, k) only implemented in trivial cases"


def theta_trans_k(M, nu_inv, c):
    """
    Returns k(M,c) in R as in Formula 8.6.1 of BL on page 227.
    Note: if c is in (1/n)ZZ, then k(M,c) is in (1/n^2)ZZ.
    
    Note: if M is 1 modulo 2*den^2, then the output is an even integer.
    
    EXAMPLE::
    
        sage: from recip import *
        sage: M = Matrix([[1, 0, 2, 0], [0, 1, 0, 0], [-2, 0, -3, 0], [0, 0, 0, 1]])
        sage: d = vector([1/2, 1/2, 0, 0])
        sage: theta_trans_k(M, 1, d)
        -3/2
    """
    A,B,C,D = ABCD(M)
    g = A.nrows()
    c1 = vector(c[0:g])
    c2 = vector(c[g:2*g])
    ret = nu_inv * (D*c1-C*c2) * (-B*c1+A*c2+subscript_zero(A*B.transpose())) - c1*c2
    if all([a in QQ for a in c]):
        n = LCM([QQ(a).denominator() for a in c])
        if not n**2*ret in ZZ:
            raise RuntimeError, "n^2 * ret not in ZZ for n=%s and ret=%s" % (n, ret)
    return ret


def theta_c_mod_Z(c):
    """
    Given c in RR^2g, returns d in [0,1)^2g and
    s in RR such that theta[c](0,Z) = e(s)*theta[d](0,Z).
    
    Note: depends only on c modulo den*ZZ.
    
    EXAMPLES::
    
        sage: from recip import *
        sage: theta_c_mod_Z([7/2, 7/3, 7/4, 7/5, 7/6, 7/7, 7/8, 7/9])
        ([1/2, 1/3, 3/4, 2/5, 1/6, 0, 7/8, 7/9], 35/6)
    """
    m = [floor(a) for a in c]
    g = ZZ(len(c)/2)
#    den = LCM([QQ(a).denominator() for a in c])
#    C.<zeta> = CyclotomicField(den)
    s = sum([m[i+g]*c[i] for i in range(g)])
    d = vector(c) - vector(m)
    return list(d), s


def _gottschling_kappa_from_list():
    zeta8 = CyclotomicField(8).gen()
    return [-zeta8^3 for i in range(11)] + [-zeta8^2 for i in range(27)]

gottschling_kappa = _gottschling_kappa_from_list()


def _determine_kappa_for_gottschling():
    """
    Returns the sequence of values of kappa(M) in QQ(zeta_{1/8})
    where M ranges over the Gottschling matrices.
    
    EXAMPLES::
    
        sage: from recip import *
        sage: from recip import _determine_kappa_for_gottschling
        sage: gottschling_kappa == _determine_kappa_for_gottschling()
        True
    """
    return [determine_kappa(M) for M in gottschling_matrices()]


def determine_kappa(M, prec=200):
    """
    Returns the element kappa(M) of QQ(zeta_{1/8}) such that for all c and Z:

      theta[M[c]](0, M(Z)) = kappa(M) det(CZ+D)^(1/2) e(k(M,c)/2) theta[c](0,Z).

    where the square root is chosen such that det(C*I+D)^(1/2) has argument
    in [0, 2*pi).
    """
    CIF = ComplexIntervalField(prec)
    i = CIF.gen()
    assert i.imag() > 0
    Z = i*identity_matrix(2)
    assert nu(M) == 1

    d = [0,0,0,0]
    LHS = evaluate_theta_interval(d, Sp_action(M, Z))

    (A,B,C,D) = ABCD(M)
    det = (C*Z+D).determinant()
    sqrtdet = det.sqrt()
    im = sqrtdet.imag()
    if C == 0:
        assert im.contains_zero()
        re = sqrtdet.real()
        assert not re.contains_zero()
        assert re > 0
    if C != 0:
        if  im.contains_zero():
#            print sqrtdet
#            print im.upper()
#            print im.lower()
            raise ValueError, "M too complicated in determine_kappa"
        if im < 0:
            sqrtdet = -sqrtdet

    (c, s) = theta_action_without_kappa(M, 1, d)
    e = exp(CIF(2*pi*I*s))

    th = evaluate_theta_interval(c, Z)

    kappa_numerical = LHS / sqrtdet / e / th
    assert kappa_numerical.parent() is CIF
    l = log(kappa_numerical)
    assert l.parent() is CIF
    k = (l.imag()/ RIF(pi/4)).unique_integer()

    g = CyclotomicField(8).gen()
    Cg = CIF(CC(g)) + RIF(-0.1,0.1) + i*RIF(-0.1,0.1)

    assert (log(Cg).imag() / RIF(pi/4)).unique_integer() == 1

    return g**k




def theta_action_without_kappa(M, nu_inv, d):
    """
    Given d in RR^2g and M in Sp_2g(ZZ), returns
    c in [0,1)^2g and s in RR such that
      theta[d](0, M(Z)) = kappa(M)*det(C*Z+D)^(1/2)*e(s)*theta[c](0, Z)
    
    This uses the theta transformation formula
      theta[M[c]](0, M(Z)) = kappa(M) det(CZ+D)^(1/2) e(k(M,c)/2) theta[c](0,Z).
      
    In fact, this is useful for slightly more, and also accepts elements M of
    GSp_2g(ZZ/NZZ) when given with the inverse of nu(M).    
      
    Note: depends only on M modulo LCM(2*den^2, 8):
    s1/2 mod ZZ depends only on M mod 2*den^2
    cpre mod 2*den depends only on M mod 2*den^2,
    hence s2 and c depend only on M mod 2*den^2.
    
    EXAMPLES::
    
        sage: from recip import *
        sage: M = Matrix([[0, 1], [-1, 0]])
        sage: theta_action_without_kappa(M, 1, [1/6, 1/2])
        ([1/2, 1/6], 1/12)
        sage: M = Matrix([[1, 1], [0, 1]])
        sage: theta_action_without_kappa(M, 1, [1/6, 1/6])
        ([1/6, 5/6], -7/72)
        sage: M = matrix_from_blocks(zero_matrix(2), identity_matrix(2), -identity_matrix(2), zero_matrix(2))
        sage: theta_action_without_kappa(M, 1, [1/2, 1/2, 1/2, 1/2])
        ([1/2, 1/2, 1/2, 1/2], 1/2)
        sage: theta_action_without_kappa(M, 1, [0, 1/2, 1/2, 0])
        ([1/2, 0, 0, 1/2], 0)
    """
    cpre = sq_brackets_inverse(M, nu_inv, d)
    # we have theta[d](0,Z) = kappa(M) det(CZ+D)^(1/2) e(k(M,cpre)/2) theta[cpre](0,Z).
    s1 = theta_trans_k(M, nu_inv, cpre)
    c, s2 = theta_c_mod_Z(cpre)
    # ..................... = kappa(M) det(CZ+D)^(1/2) e(s1/2) e(s2) theta[c](0,Z)
    return c, s1/2+s2


def theta_leading_term_complex(c, C=CC):
    """
    Returns a pair (a, o), where o is rational and a is an element of C such
    that theta[c] = a*q^o + h.o.t. The output a is in C.
    
    Here c=(c1,c2) is a pair of numbers, and this only works in
    the classical genus-one case.
    """
    # terms are exp(pi*i*(n+c1)^2*z + 2*pi*i*(n+c1)*c2) for integers n,
    # so we may as well translate c1 to the interval (-1/2, 1/2].
    (c1,c2) = c
    c1 = c1 + (1/2-c1).floor()
    assert -1/2 < c1 and c1 <= 1/2
    i = C.gen()
    assert i^2 == -1
    if c1 == 1/2:
        if c2 + (1/2-c2).floor() == 1/2:
            return (0, oo)
        # two terms: n+c1 = 1/2 and n+c1 = -1/2
        # exp(pi*i*(n+c1)^2*z + 2*pi*i*(n+c1)*c2)
        # = exp(2*pi*i*z/8) * exp(\pm 2*pi*i*c2/2)
        # so their sum is 2*Re(exp(2*pi*i*c2))*q^(1/8)
        return (C(2*exp(C(pi)*i*c2).real()), 1/8)
    # one term: exp(pi*i*c1^2*z + 2*pi*i*c1*c2)
    return (exp(2*C(pi)*i*c1*c2), c1^2/2)


def theta_leading_term_cyclotomic(c, n=None):
    """
    Returns a pair (a, o), where o is rational and a is an element of C such
    that theta[c] = a*q^o + h.o.t. The output a is in CyclotomicField(n)
    or ZZ. If n is None, then n is 2 times the product of the denominators
    of c1 and c2.
    
    Here c=(c1,c2) is a pair of numbers, and this only works in
    the classical genus-one case.
    """
    (c1,c2) = c
    if n is None:
        n = c1.denominator()*c2.denominator()*2
        C.<zeta> = CyclotomicField(n)
    # terms are exp(pi*i*(n+c1)^2*z + 2*pi*i*(n+c1)*c2) for integers n,
    # so we may as well translate c1 to the interval (-1/2, 1/2].
    c1 = c1 + (1/2-c1).floor()
    assert -1/2 < c1 and c1 <= 1/2
    if c1 == 1/2:
        if c2 + (1/2-c2).floor() == 1/2:
            return (0, oo)
        # two terms: n+c1 = 1/2 and n+c1 = -1/2
        # exp(pi*i*(n+c1)^2*z + 2*pi*i*(n+c1)*c2)
        # = exp(2*pi*i*z/8) * exp(\pm 2*pi*i*c2/2)
        # so their sum is 2*Re(exp(2*pi*i*c2))*q^(1/8)
        return (zeta**(ZZ(c2*n/2))+zeta**(-ZZ(c2*n/2)), 1/8)
    # one term: exp(pi*i*c1^2*z + 2*pi*i*c1*c2)
    if c1 == 0 or c2 == 0:
        return (1, c1^2/2)
    return (zeta**(ZZ(n*c1*c2)), c1^2/2)


def theta_lc_complex(c, C=CC):
    """
    Returns the number a from :func:`theta_leading_term_complex`.
    """
    return theta_leading_term_complex(c, C)[0]

def theta_lc_cyclotomic(c, n=None):
    """
    Returns the number a from :func:`theta_leading_term_cyclotomic`.
    """
    return theta_leading_term_cyclotomic(c, n)[0]
    

try:
    theta_rings = theta_rings
except NameError:
    theta_rings = []


def make_theta_ring(g, den):
    """
    See theta_ring.
    
    EXAMPLE::
    
        sage: from recip import *
        sage: make_theta_ring(1,6)[0]
        Multivariate Polynomial Ring in t0d6, t1d6, t2d6, t3d6, t4d6, t5d6, t6d6, t7d6, t8d6, t9d6, t10d6, t11d6, t12d6, t13d6, t14d6, t15d6, t16d6, t17d6, t18d6, t19d6, t20d6, t21d6, t22d6, t23d6, t24d6, t25d6, t26d6, t27d6, t28d6, t29d6, t30d6, t31d6, t32d6, t33d6, t34d6, t35d6 over Cyclotomic Field of order 72 and degree 24

    """
    is_den_even(den)
    n = den**(2*g)
    C = CyclotomicField(2*den**2)
    z = C.gen()
#    i = 1
#    while 10^i < n:
#        i += 1
#    names = ["t" + str(10^i + j)[1:] + "d" + str(den) for j in range(n)]
    if den == 2:
        names = ["t" + str(j) for j in range(n)]
    else:
        names = ["t" + str(j) + "d" + str(den) for j in range(n)]
    P = PolynomialRing(C, n, names = names)
    #return P
    ideal_gens = []
    eval = Sequence(P.gens())
    for i in range(n):
        c = num_to_c(i, g = g, den = den)
        d, k = theta_c_mod_Z(-vector(c))
        e = ZZ(2*den**2*k)
        j = c_to_num(d, den = den)
        if i > j:
            ideal_gens = ideal_gens + [P.gens()[i] - z**e * P.gens()[j]]
            eval[i] = z**e * P.gens()[j]
        elif i == j and e % (2*den**2) != 0:
            ideal_gens = ideal_gens + [P.gens()[i]]
            eval[i] = 0
        #print i,j,c,d,k,e,ideal_gens
    I = P.ideal(ideal_gens)
    return P, eval, P.quotient(I)


def theta_ring(g, den):
    """
    Returns a triple (P, eval, R),
    where P is a polynomial ring representing polynomials
    in theta[c](0,Z) with c in (1/den)*ZZ^2g
    and coefficients in CyclotomicField(LCM(8, 2*den^2)),
    R is a quotient of P by an ideal of zero modular forms,
    and eval is such that R(x(eval)) = R(x) for all x in P.
    
    EXAMPLES::
    
        sage: from recip import *
        sage: tr = theta_ring(3, 2); len(tr)
        3
        sage: tr[0]
        Multivariate Polynomial Ring in t0, t1, t2, t3, t4, t5, t6, t7, t8, t9, t10, t11, t12, t13, t14, t15, t16, t17, t18, t19, t20, t21, t22, t23, t24, t25, t26, t27, t28, t29, t30, t31, t32, t33, t34, t35, t36, t37, t38, t39, t40, t41, t42, t43, t44, t45, t46, t47, t48, t49, t50, t51, t52, t53, t54, t55, t56, t57, t58, t59, t60, t61, t62, t63 over Cyclotomic Field of order 8 and degree 4
    """
    global theta_rings
    for r in theta_rings:
        if r[0] == g and r[1] == den:
            return r[2]
    rng = make_theta_ring(g, den)
    theta_rings = theta_rings + [(g, den, rng)]
    return rng


def theta_ring_inclusion(g, den1, den2):
    """
    Returns the natural map from theta_ring(g, den1) to theta_ring(g, den2)
    as a sequence to evaluate elements of theta_ring(g, den1) in.
    
    EXAMPLE::
    
        sage: from recip import *
        sage: theta_ring_inclusion(2,2,4)
        [t0d4, t2d4, t8d4, t10d4, t32d4, t34d4, t40d4, t42d4, t128d4, t130d4, t136d4, t138d4, t160d4, t162d4, t168d4, t170d4]

    """
    P = theta_ring(g, den2)[0]
    gens = P.gens()
    return [gens[c_to_num(num_to_c(i, g, den1), den2)] for i in range(den1**(2*g))]


def c_to_num(c, den):
    """
    Returns an integer from 0 to den^2g (exclusive), one unique
    number for each c in [0,1)^2g with den*c integral.
    
    For g = 2 and den = 2, this is Dupont's notation.
    
    EXAMPLE::
    
        sage: from recip import *
        sage: c_to_num(num_to_c(1234, g = 5, den = 7), den = 7)
        1234
    """
    
    g = ZZ(len(c)/2)
    d = den*vector(c)
    for a in d:
        if not (a in ZZ and a >= 0 and a < den):
            raise ValueError, "coordinates of c (=%s) not all in [0, 1) or not all in (1/den)*ZZ for den=%s" % (c,den) 
    ret = 0
    r = c[0:g]
    s = c[g:2*g]
    t = den
    for a in Sequence(s) + Sequence(r):
        ret = ret + a * t
        t = t * den
    return ret


def num_to_c(n, g, den):
    """
    Inverse of c_to_num.
    
    For g = 2 and den = 2, this is Dupont's notation.
    
    EXAMPLE::
    
        sage: from recip import *
        sage: num_to_c(c_to_num([0,1/2,2/3,3/4,4/5,5/6],den=60), g=3, den=60)
        [0, 1/2, 2/3, 3/4, 4/5, 5/6]
    """
    if not (n in ZZ and n >= 0 and n < den**(2*g)):
        raise ValueError, "n (=%s) is not an integer between 0 (inclusive) and den^(2g) for den=%s and g=%s" % (n, den, g)
    n = ZZ(n)
    g = ZZ(g)
    den = ZZ(den)
    sr = [0 for i in range(2*g)]
    for i in range(2*g):
        sr[i] = (n % den) / den
        n = ZZ(n/den - sr[i])
    c = sr[g:2*g] + sr[0:g]
    return c       


def name_to_den(name):
    """
    Given a string ending in "d?????", returns
    the integer with decimal expansion "?????".
    
    EXAMPLE::
    
        sage: from recip import *
        sage: name_to_den('t12d034')
        34

    """
    i = name.find('d')
    if i < 0:
        return ZZ(2)
    return Integer(name[i+1:], base=10)


def name_to_den_c(name, g):
    """
    Given a string name = "t.....d?????", returns
    (den, num_to_c(.....), g, den))
    for den = ????? (both interpreted as decimal expansions.
    
    EXAMPLE::
    
        sage: from recip import *
        sage: name_to_den_c('t0123d45', 2)
        (45, [0, 0, 11/15, 2/45])

    """
    if not name[0] == 't':
        raise ValueError, "string name (=%s) does not start with 't'" % name
    i = name.find('d')
    if i < 0:
        return ZZ(2), num_to_c(Integer(name[1:], base=10), g, 2)
    den = name_to_den(name)
    return den, num_to_c(Integer(name[1:i], base=10), g, den)


def cycl_galois_action(alpha, amodb):
    """
    For alpha in CyclotomicField(b), return alpha with zeta_b raised to the power a
    
    EXAMPLES::
    
        sage: from recip import *
        sage: C = CyclotomicField(10)
        sage: amodb = Zmod(15)(7)
        sage: z = C.gen()
        sage: cycl_galois_action(z^5+z+1, amodb)
        -zeta15^3
        
    The modulus must match up::
        
        sage: f = Zmod(6)(5)
        sage: cycl_galois_action(z, f)
        Traceback (most recent call last):
        ...
        TypeError: Cannot coerce zeta10 into Cyclotomic Field of order 6 and degree 2

    """
    a = ZZ(amodb)
    b = amodb.modulus()
    C = CyclotomicField(b)
    pol = C(alpha).polynomial()
    x = pol.parent().gen()
    return C(pol(x^a))


def cycl_galois_action_on_polynomials(x, amodb):
    """
    On input a polynomial x and an element amodb of Zmod(b),
    return x where cycl_galois_action is aplied to each coefficient.
    
    
    EXAMPLES::
    
        sage: from recip import *
        sage: P = theta_ring(2,6)[0]
        sage: z = P.base_ring().gen()
        sage: x = P.gens()[1]*(z+z^2+z^3)+P.gens()[3]*z^-1; x
        (zeta72^3 + zeta72^2 + zeta72)*t1d6 + (-zeta72^23 + zeta72^11)*t3d6
        sage: cycl_galois_action_on_polynomials(x, Zmod(72)(73)) == x
        True
        sage: cycl_galois_action_on_polynomials(x, Zmod(72)(3))
        (zeta72^9 + zeta72^6 + zeta72^3)*t1d6 + (2*zeta72^21 - 2*zeta72^9)*t3d6
    """
    P = x.parent()
    return P(sum([cycl_galois_action(t[0], amodb)*t[1] for t in Sequence(x)]))


def my_sum(s):
    if len(s) == 0:
        return 0
    ret = s[0]
    for i in range(1,len(s)):
        ret = ret.__add__(s[i])
    return ret


def ThetaModForm(x, g = None, den = None):
    """
    Constructs a modular form given as a polynomial in theta series.

    EXAMPLES:
      
        sage: from recip import *
        sage: P = theta_ring(3, 2)[0]
        sage: x = sum(P.gens()[2:7])
        sage: t = ThetaModForm(x, g = 3); t
        t2 + t3 + t4 + t5 + t6
        sage: t.long_name()
        'th[0,0,0;0,1/2,0] + th[0,0,0;1/2,1/2,0] + th[0,0,0;0,0,1/2] + th[0,0,0;1/2,0,1/2] + th[0,0,0;0,1/2,1/2]'
    """
    if type(x) == Theta_element_polynomial_ring:
        return x
    if x in ZZ and not (g is None or den is None):
        c = num_to_c(x, g=g, den=den)
        return ThetaProduct(c)
    if is_RingElement(x):
        B = x.parent()
        if is_FractionField(B):
            y = x.denominator()
            x = x.numerator()
            B = B.ring()
        else:
            y = B(1)
        n = B.ngens()
        if den == None:
            den = name_to_den(B.variable_name())
        if g == None:
            if den == 1:
                raise ValueError, "g not supplied an cannot be read off from x (=%s) in %s" % (x, x.parent())
            g = 1
            while den**(2*g) < n:
                g += 1
        P, eval, R = theta_ring(g, den)
        try:
            x = P(x)
            y = P(y)
        except TypeError:
            raise TypeError, "x (=%s) or y (=%s) is not an element of a polynomial ring supplied by theta_ring for g=%s and den=%s" % (x, y, g, den)
        if R(y) == 0:
            if R(x) == 0:
                raise NotImplementedError, "cannot simplify x/y to something that is not 0/0 for x=%s and y=%s" % (x,y)
            raise ValueError, "denominator must be non-zero in ThetaModForm, but is %s" % y
        if not x.is_homogeneous():
            raise ValueError, "numerator %s not homogeneous" % num_pol
        if not y.is_homogeneous():
            raise ValueError, "denominator %s not homogeneous" % den_pol
        x = P(x(eval))
        y = P(y(eval))
        return Theta_element_polynomial_ring(x, y, g, den)
    if type(x) == list or is_Vector(x):
        if len(x) % 2 == 1:
            raise ValueError, "Length of x (=%s) is odd" % x
        g = ZZ(len(x)/2)
        den = LCM([2]+[a.denominator() for a in x])
        n = c_to_num(x, den)
        P, eval, R = theta_ring(g, den)
        return ThetaModForm(P.gens()[n], g, den)
    try:
        return ThetaModForm(x.rational_function())
    except AttributeError:
        raise TypeError, "type %s of x (=%s) not allowed in ThetaModForm" % (type(x), x)


        
class Theta_element(Element):
    def orbit(self, M, group_elements=False):
        """
        Returns the orbit of self under M.
        Note: see :meth:`__eq__`
        
        INPUT:
        
        - group_elements can be either False, True, or "list".
          If False, return only the orbit
          If "list", also return a list of lists `l` of elements of M such that
          the elements of the orbit are given by l*self
          If True, then give an element of the group generated by M instead
          of a list of elements of M.
        
        EXAMPLES::
        
            sage: from recip import *
            sage: gens = symplectic_generators(2)
            sage: th0 = ThetaModForm(theta_ring(2,2)[0].gens()[0])
            sage: h4 = my_sum((th0^8).orbit(gens)); h4 # long time
            t0^8 + t1^8 + t2^8 + t3^8 + t4^8 + t6^8 + t8^8 + t9^8 + t12^8 + t15^8
        """
        grp_elements_ret = [[]]
        grp_elements_nw = [[]]
        if not type(M) == list:
            M = [M]
        ret = [self]
        nw = [r for r in ret]
        while not len(nw) == 0:
            nw2 = []
            grp_elements_nw2 = []
            for A in M:
                for y_ind in range(len(nw)):
                    y = nw[y_ind]
                    g_y = grp_elements_nw[y_ind]
                    z = y^A
                    if not z in ret:
                        ret = ret + [z]
                        grp_elements_ret = grp_elements_ret + [g_y + [A]]
                        nw2 = nw2 + [z]
                        grp_elements_nw2 = grp_elements_nw2 + [g_y + [A]]
            nw = nw2
            grp_elements_nw = grp_elements_nw2
        if group_elements == False:
            return ret
        if len(grp_elements_ret) != len(ret):
            raise RuntimeError
        if group_elements == "list":
            return ret, grp_elements_ret
        raise NotIMplementedError

        
        
class ThetaSum(CommutativeRingElement, Theta_element):
    _sequence = None # consists of tuples, (s,sc) thetaproducts and coeffs
    _g = None
    _level = None
    _weight = None

    def __init__(self, sequence):
        if isinstance(sequence, int):
            sequence = ZZ(sequence)
        if isinstance(sequence, Integer):
            if sequence == 0:
                sequence = []
            else:
                sequence = [(ThetaProduct([]), sequence)]
        if type(sequence) != list:
            sequence = [(sequence,1)]
        sequence.sort()
        newsequence = []
        for (s, sc) in sequence:
            if s._unity != 0:
                N = s._unity.denominator()
                zeta = CyclotomicField(N).gen()
                x = ZZ(N*s._unity)
                sc = sc * zeta**x
                if sc in QQ:
                    sc = QQ(sc)
#                else:
#                    N = sc.parent().zeta().multiplicative_order()
#                    divs = N.divisors()
#                    divs.sort()
#                    for d in divs:
#                        if sc in CyclotomicField(d):
#                            sc = CyclotomicField(d)(sc)
#                            break

                s = ThetaProduct(s._sequence)
            if newsequence and newsequence[-1][0]._sequence == s._sequence:
                newsequence[-1][1] += sc
            else:
                newsequence.append([s, sc])
        sequence = [(s, sc) for [s, sc] in newsequence if sc != 0]
        self._sequence = sequence
        
        level = 1
        weight = None
        g = None
        for (s,cs) in sequence:
            if weight is None:
                weight = s.weight()
            elif weight != s.weight():
                raise ValueError
            if g is None:
                g = s.g()
            elif g != s.g():
                raise ValueError
            level = lcm(level, s.level())
        CommutativeRingElement.__init__(self, ThetaRing())

    def __call__(self, z, use_magma=False, prec=None):
        """
        For z in H_g, returns self evaluated at z.
        
        WARNING: unreliable if z is not reduced
        """
        return sum([sc*s(z, use_magma, prec) for (s, sc) in self._sequence])

    def _repr_(self, long=None):
        if not self._sequence:
            return "0"
        ret = ""
        for (t,c) in self._sequence:
            ret += " + "
            if c != 1:
                strc = str(c)
                if " + " in strc or " - " in strc:
                    strc = "("+strc+")"
                ret += strc + "*"
            ret += t._repr_(long)
        ret.replace('+ -', '- ')
        if ret[:3] == " + ":
            ret = ret[3:]
        return ret

    def __pow__(self, M):
        r"""
        Let the symplectic matrix M act on self from the right.
        Or if M is an integer, raise to that power.
        """
        if type(M) is Integer:
            return CommutativeRingElement.__pow__(self, M)
#        for (s, cs) in self._sequence:
#            if not cs in QQ:
#                raise NotImplementedError
        return ThetaSum([(s**M,cycl_galois_action(cs, M.nu())) for (s, cs) in self._sequence])
    
    def _eq_(self, other):
        r"""
        Check whether self and other are the same rational function
        in theta constants.
        This is not equivalent to being the same modular forms.
        """
        return self._sequence == other._sequence
               
    def _mul_(self, other):
        r"""
        Return the product of self and other.
        """
        s1 = self._sequence
        s2 = other._sequence
        ret = [(a*b,ca*cb) for (a,ca) in s1 for (b,cb) in s2]
        return ThetaSum(ret)
        
    def _neg_(self):
        return ThetaSum([(a,-ca) for (a,ca) in self._sequence])

    def _add_(self, other):
        return ThetaSum(self._sequence + other._sequence)

    def _invert_(self):
        raise NotImplementedError

    def _div_(self, other):
        return self._mul_(other._invert_())
        
    def _sub_(self, other):
        return self._add_(other._neg_())
        
    def __cmp__(self, other):
        l = self._repr_()
        r = other._repr_()
        if r > l:
            return -1
        if r < l:
            return 1
        return 0


class ThetaRing(CommutativeRing, UniqueRepresentation):
    def __init__(self):
        pass
    
    Element = ThetaSum

    def _coerce_map_from_(self, other):
        if other is ZZ or other is QQ or other is int:
            return True
        if other is ThetaGroup():
            return True
            

    def _convert_map_from_(self, other):
        if other is ZZ or other is QQ or other is int:
            return True
        if other is ThetaGroup():
            return True
    
    def __call__(self, other):
        """
        I'm probably not supposed to override this, but it wasn't working
        before
        """
        return ThetaSum(other)

   
class ThetaProduct(MultiplicativeGroupElement, Theta_element):

    # _sequence is a list of triples (c,d,e), where c is a theta characteristic
    # d*c is integral, and e is an exponent (non-zero integer)
    _sequence = None
    _unity = None
    _g = None
    _d = None
    
    def __init__(self, sequence, unity=0):
        if isinstance(sequence, Integer):
            if sequence < 16 and sequence >=0:
                sequence = num_to_c(sequence, 2,2)
        if type(sequence) == list:
            if not sequence:
                self._g = None
                self._d = None
                self._sequence = sequence
                self._unity = 0
                MultiplicativeGroupElement.__init__(self, parent=ThetaRing())
                return
            if len(sequence) > 0 and type(sequence[0]) != list and \
                                     type(sequence[0]) != tuple:
                d = lcm([x.denominator() for x in sequence])
                sequence = [(sequence, d, 1)]
        self._sequence = sequence
        self._unity = unity
        if len(sequence) == 0:
            self._d = 1
        else:
            self._d = lcm([d for (c,d,e) in sequence])
            self._g = ZZ(len(sequence[0][0])/2)
        MultiplicativeGroupElement.__init__(self, parent=ThetaRing())
    
    def __call__(self, z, use_magma=False, prec=None):
        """
        For z in H_g, returns self evaluated at z.
        
        WARNING: unreliable if z is not reduced
        """
        g = self._g
        den = self._d
        if _is_numerical_complex_field(z.base_ring()):
            if prec != None:
                raise ValueError, "prec!=None cannot be combined with a " \
                                  "non-exact period matrix"
            C = z.base_ring()
        else:
            C = ComplexField(prec)
            z = mat_convert(z, C)
        ret = 1
        for (c,d,e) in self._sequence:
            ret = ret * evaluate_theta(c, z, use_magma=use_magma)**e
        ret = ret * C(exp(2*pi*C.gen()*self._unity))
        return ret


    @cached_method
    def _repr_(self, long=None):
        if (not self._d is None) and long is None:
            long = (2 % self._d == 1) 
        unity = self._unity
        d = unity.denominator()
        if 4 % d == 0:
            if unity == 0:
                ret = ""
            elif unity == 1/2:
                ret = "-"
            elif unity == 1/4:
                ret = "i*"
            elif unity == 3/4:
                ret = "-i*"
            else:
                raise RuntimeError
        else:
            n = unity.numerator()
            if 2*n < d:
                if n == 1:
                    ret = "zeta%s*" % d
                else:
                    ret = "zeta%s^%s*" % (d, n)
            else:
                if n == d+1:
                    ret = "-zeta%s*" % d
                else:
                    ret = "-zeta%s^%s*" % (d, n-d)
        if self._d is None:
            if ret == "":
                return "1"
            return ret
        factors = []
        for (c,d,e) in self._sequence:
            if 2 % d == 0 and not long:
                num = c_to_num(c,den=d)
                r = "t%s" % num
            else:
                num = c
                r = "t[" + ",".join(["%s" % (x) for x in c[:self._g]]) + \
                   ";" + ",".join(["%s" % (x) for x in c[self._g:]]) + "]"
            eabs = abs(e)
            if eabs != 1:
                r = r + "^%s" % eabs
            factors.append((e,num,r))
        numfactors = [(num,r) for (e,num,r) in factors if e > 0]
        denfactors = [(num,r) for (e,num,r) in factors if e < 0]
        numfactors.sort()
        denfactors.sort()
        ret = ret + "*".join([f[1] for f in numfactors])
        if len(denfactors) > 0:
            if len(denfactors) > 1:
                ret = ret + "/(" + "*".join([f[1] for f in denfactors]) + ")"
            else:
                ret = ret + "/" + denfactors[0][1] 
        if ret[-1] == "*":
            ret = ret[:-1]
        return ret

    def __pow__(self, M):
        r"""
        Let the symplectic matrix M act on self from the right.
        Or if M is an integer, raise to that power.
        """
        if type(M) is Integer:
            return MonoidElement.__pow__(self, M)
        ret = []
        nu_inv = lift_small(M.nu()^-1)
        unity_ret = lift_small(M.nu())*self._unity
        M = mat_convert(M.matrix(), lift_small)
        for (c,d,e) in self._sequence:
            cret, s = theta_action_without_kappa(M, nu_inv, c)
            unity_ret = unity_ret + s*e
            ret.append((cret,d,e))
        n = unity_ret.numerator() % unity_ret.denominator()
        ret.sort()
        return ThetaProduct(ret, n/unity_ret.denominator())
    
    def _cmp_(self, other):
        if self._sequence < other._sequence:
            return -1
        if self._sequence > other._sequence:
            return 1
        return self._unity._cmp_(other._unity)

    __cmp__ = _cmp_
    

    
    def _eq_(self, other):
        r"""
        Check whether self and other are the same quotient of theta constants.
        I don't know whether this is equivalent to being the same modular forms.
        """
        return self._sequence == other._sequence and \
               self._unity    == other._unity
               
    def __eq__(self, other):
        r"""
        Check whether self and other are the same quotient of theta constants.
        I don't know whether this is equivalent to being the same modular forms.
        """
        return self._sequence == other._sequence and \
               self._unity    == other._unity

    def _mul_(self, other):
        r"""
        Return the product of self and other.
        """
        ret = self._sequence + other._sequence
        ret.sort()
        for i in range(len(ret)-1):
            if ret[i][0] == ret[i+1][0]:
                d = gcd(ret[i][1],ret[i+1][1])
                ret[i] = (ret[i][0], d, ret[i][2]+ret[i+1][2])
                ret[i+1] = (0,0,0)
        unity = self._unity+other._unity
        if unity >= 1:
            unity = unity - 1
        return ThetaProduct([r for r in ret if r[2]!=0], unity)
        
    def _neg_(self):
        if self._unity >= 1/2: 
            return ThetaProduct(self._sequence, self._unity-1/2)
        return ThetaProduct(self._sequence, self._unity+1/2)

    def _add_(self, other):
        return ThetaSum(self)+ThetaSum(other)
        
    def __add__(self, other):
        return ThetaSum(self)+ThetaSum(other)

    def _invert_(self):
        ret = [(c,d,-e) for (c,d,e) in self._sequence]
        if self._unity == 0:
            return ThetaProduct(ret, 0)
        return ThetaProduct(ret, 1-self._unity)

    def _div_(self, other):
        return self._mul_(other._invert_())
        
    def _sub_(self, other):
        return self._add_(other._neg_())
        
    def _coerce_map_from_(self, other):
        if other is ZZ or other is QQ:
            return True
        if is_CyclotomicField(other):
            return True
    
    def weight(self):
        return sum([s[2] for s in self._sequence])
        
    def g(self):
        return self._g
        
    def level(self):
        """
        Returns a multiple of the level of this function
        """
        if self._d is None:
            return ZZ(1)
        return self._d**2 * 2

            
class ThetaGroup(Group, UniqueRepresentation):
    def __init__(self):
        pass
    
    Element = ThetaProduct
    
    def is_abelian(self):
        return True
    
    def is_finite(self):
        return False
    
    def order(self):
        return PlusInfinity
    
    def _repr_(self):
        return "Group of formal products of theta constants"
    
    def _element_constructor_(self, *args, **kwds):
        return self.element_class(*args, **kwds)
            
            
class Theta_element_polynomial_ring(Theta_element):
    _num_pol = None
    _den_pol = None
    _g = None
    _den = None
    
    def __init__(self, num_pol, den_pol, g, den):
        self._num_pol = num_pol
        self._den_pol = den_pol
        self._g = g
        self._den = den
        CommutativeRingElement.__init__(self, ThetaRing())
    
    @cached_method
    def _cached_call(self, z, use_magma=False, prec=None, interval=False):
        return self._call_code(z, use_magma=use_magma, prec=prec, interval=interval)
    
    def __call__(self, z, use_magma=False, prec=None, interval=False):
        """
        For z in H_g, returns self evaluated at z.
        
        WARNING: unreliable if z is not reduced
        """
        if z.is_mutable():
            return self._call_code(z, use_magma=use_magma, prec=prec, interval=interval)
        return self._cached_call(z, use_magma=use_magma, prec=prec, interval=interval)
    
    def _call_code(self, z, use_magma=False, prec=None, interval=False):
        g = self._g
        den = self._den
        if _is_numerical_complex_field(z.base_ring()):
            if prec != None:
                raise ValueError, "prec!=None cannot be combined with a " \
                                  "non-exact period matrix"
            C = z.base_ring()
            if get_verbose() > 1:
                print "Working with a period matrix over CC"
            if interval:
                vals = [evaluate_theta_interval(num_to_c(i, g, den), z)
                             for i in range(den**(2*g))]
            else:
                vals = [evaluate_theta(num_to_c(i, g, den), z,
                             use_magma=use_magma) for i in range(den**(2*g))]
        elif hasattr(z, '_theta_vals'):
            if interval:
                C = ComplexIntervalField(prec)
            else:
                C = ComplexField(prec)
            if get_verbose() > 1:
                print "Working with a period matrix over a number field"
            vals = z._theta_vals(den, prec=prec, use_magma=use_magma, interval=interval)
        else:
            raise TypeError, "Period matrix z (=%s) has incorrect type " \
                             "(%s) or base ring" % (z, type(z))
        num_ret = sum([C(a[0])*a[1](vals) for a in Sequence(self._num_pol)])
        den_ret = sum([C(a[0])*a[1](vals) for a in Sequence(self._den_pol)])
        return num_ret / den_ret
        
    def __pow__(self, M):
        """
        For M an integer, simply raise self to the power M.
        
        For M in Sp_2g(ZZ), returns Z |--> kappa(M)^(-2k) * det(CZ+D)^-k self(M(Z)),
        where k is the weight of the numerator divided by the weight of the
        denominator.
        
        This works also for M in Sp_2g(Zmod(m)) with m divisible by 8 and 2*den^2.
        
        TODO: Use the Sage category framework properly.
        
        EXAMPLES:
        
            sage: from recip import *
            sage: P = theta_ring(2, 2)[0]
            sage: x = sum(P.gens()[2:6])
            sage: t = ThetaModForm(x, g = 2, den = 2)
            sage: M = Matrix(Zmod(8), [(5, 0, 7, 1), (6, 6, 0, 7), (1, 2, 2, 3), (2, 3, 5, 4)])
            sage: t^M
            (zeta8)*t1 + (zeta8^2)*t9 + t15
            sage: P = theta_ring(2,2)[0]
            sage: x = ThetaModForm(P.gens()[0] / P.gens()[1])
            sage: x.long_name()
            'th[0,0;0,0]/th[0,0;1/2,0]'
            sage: M = Matrix([(-10, 3, 6, -9), (-1, 0, 0, -2), (5, -2, -4, 3), (-13, 5, 9, -11)])
            sage: Z = Matrix(CC, [(9.2629056612 + 1.18488500061*I, 7.56582781329 + 0.330070006085*I), (7.56582781329 + 0.330070006085*I, 7.28011701001 + 9.0055158146*I)])
            sage: y = x^M; y
            (zeta8^2)*t6/t0
            sage: x(Sp_action(M,Z), use_magma=True) # optional - magma
            -0.69816125872... + 0.42001387598...*I
            sage: Z.set_immutable()
            sage: y(Z)
            -0.698161258723300 + 0.420013875983267*I
        """
        den = self._den
        g = self._g
        if isinstance(M, sage.rings.finite_rings.integer_mod.IntegerMod_int):
            if not M in Zmod(LCM(2*den**2, 8)):
                raise ValueError
            M = Zmod(LCM(2*den**2, 8))(M)
            return ThetaModForm(cycl_galois_action_on_polynomials(self._num_pol, M)/
                                cycl_galois_action_on_polynomials(self._den_pol, M), g, den)
        if M in ZZ:
            return ThetaModForm(self._num_pol**M / self._den_pol**M, g, den)

        try:
        # TODO: be absolutely sure the following line should have nu^-1
        #       and not nu^-1?
        #       Check this in paper, note: no problem for den=2,
        #       as Zmod(8)^* has exponent 2
            nu_inv = Zmod(2*den**2)(M.nu())**-1
            action = M.action_on_theta_generators(den)
            return ThetaModForm(cycl_galois_action_on_polynomials( \
                                    self._num_pol,nu_inv)(action) / \
                                cycl_galois_action_on_polynomials( \
                                    self._den_pol, nu_inv)(action), g, den)
        #  The following three lines should be an incorrect older version of
        #  the three lines above. I'm keeping them here just in case.
        #            return ThetaModForm(cycl_galois_action_on_theta(
        #            self._num_pol(Sp_action), nu_inv)/
        #            cycl_galois_action_on_theta(self._den_pol(
        #            Sp_action), nu_inv), g, den)
        except AttributeError:
            pass
        if is_Matrix(M):
            return self.__pow__(GSp_element(M))
        
    def __eq__(self, right):
        """
        Returns True if and only if self is equal to right
        as a rational function in the theta functions.
        
        Note that this may return False even if self and right
        are equal as functions on H_g.
        """
        g = self._g
        den = self._den
        P, eval, R = theta_ring(g, den)
        right = ThetaModForm(right, g=g)
        return R(self._den_pol*right._num_pol) == R(self._num_pol * right._den_pol)
        
    def weight(self):
        return self.num_weight() - self.den_weight()
        
    def num_weight(self):
        return self._num_pol.degree() / 2
    
    def den_weight(self):
        return self._den_pol.degree() / 2

    def long_name(self):
        ret = str(self.rational_function())
        g = self._g
        den = self._den
        for i in range(den**(2*g)-1,-1,-1):
            c = num_to_c(i, g, den)
            r = c[0:g]
            s = c[g:2*g]
            rstr = str(r).replace(" ", "").replace("]",";")
            sstr = str(s).replace(" ", "").replace("[","")
            longstr = "th" + rstr + sstr
            if den == 2:
                shortstr = "t" + str(i)
            else:
                shortstr = "t" + str(i) + "d" + str(den)
            ret = ret.replace(shortstr, longstr)
        return ret
        
    def _repr_(self):
        return str(self._num_pol / self._den_pol).replace("d", "_")
        
    def change_den(self, den):
        """
        Returns a representation of self with denominator den
        """
        if den == self._den:
            return self
        if get_verbose() > 1:
            print "Changing denominator of theta characteristic from %s to %s" % (self._den, den)
        if not (den % self._den) == 0:
            raise NotImplementedError, "Lowering den is not implemented: trying to change den from %s to %s" % (self._den, den)
        g = self._g
        change = change_theta_ring(g, self._den, den)
        return ThetaModForm(self.rational_function()(change), g, den)

    def __neg__(self):
        return ThetaModForm(-self.rational_function())

    def __sub__(self, right):
        return self.__add__(right.__neg__())

    def __add__(self, right):
        g = self._g
        den1 = self._den
        right = ThetaModForm(right)
        den2 = right._den
        if not g == right._g:
            raise ValueError, "self (=%s) and right (=%s) do not have same genus" % (self, right)
        den = LCM(den1, den2)
        left = self.change_den(den).rational_function()
        right = right.change_den(den).rational_function()
        return ThetaModForm(left+right, g, den)

    def __mul__(self, right):
        g = self._g
        den1 = self._den
        right = ThetaModForm(right)
        den2 = right._den
        if not g == right._g:
            raise ValueError, "self (=%s) and right (=%s) do not have same genus" % (self, right)
        den = LCM(den1, den2)
        left = self.change_den(den).rational_function()
        right = right.change_den(den).rational_function()
        return ThetaModForm(left*right, g, den)
                            
    def __inv__(self):
        g = self._g
        den = self._den
        P, eval, R = theta_ring(g, den)
        if R(self._den_pol) == 0:
            raise ValueError, "trying to invert zero element %s" % self
        return ThetaModForm(self.rational_function()**-1, g, den)
        
    def __div__(self, right):
        return self * right.__inv__()
        
    __invert__ = __inv__
    
    def is_fixed_by(self, M):
        """
        Returns the orbit of self under M.
        Note: see __eq__
        """
        if not type(M) == list:
            M = [M]
        for A in M:
            if not self^A == self:
                return False
        return True

    def rational_function(self):
        return self._num_pol / self._den_pol
    
    def is_power_of_theta(self, data=False):
        """
        Returns true if and only if self is a constant times
        a positive power of a non-zero theta constant.
        
        If data is True, then also output a triple (a, n, e)
        such that self = a*theta_n^e.
        """
        l = list(self._num_pol)
        if not (self._den_pol.is_constant() and len(l) == 1):
            if data:
                return False, None
            return False
        m = l[0][1]
        degs = m.degrees()
        e = 0
        for m in range(len(degs)):
            d = degs[n]
            if d != 0:
                if e != 0:
                    if data:
                        return False, None
                    return False
                e = d
                n = m
        if e == 0:
            if data:
                return False, None
            return False
        if data:
            return True, (l[0][0]/list(self._den_pol)[0][0], n, e)
        return True
        
    def multiplication_formula(self, n):
        """
        Gives Z |-> self(nZ)
        
        EXAMPLES::
        
            sage: from recip import *
            sage: M = Matrix([[I, 1/2+I/7], [1/2+I/7, 3/2*I-1/5]], ring=CC)
            sage: M.set_immutable()
            sage: P = theta_ring(2,2)[0]
            sage: N = 1/2 * M
            sage: N.set_immutable()
            sage: u = ThetaModForm(P.gens()[0]^2/P.gens()[1]^2)
            sage: u(M)
            1.39698675882520 + 0.0121273780920932*I
            sage: u.multiplication_formula(2)(N)
            1.39698675882520 + 0.0121273780920933*I
            
        Here is a Delta-quotient, but it will need to be simplified before it is useful::
        
            sage: i = igusa_modular_forms()[2]; i
            t0^2*t1^2*t2^2*t3^2*t4^2*t6^2*t8^2*t9^2*t12^2*t15^2
            sage: q = i.multiplication_formula(2)/i

        """
        if n != 2:
            raise NotImplementedError, "Sorry, multiplication formula only implemented for n=2 (i.e. duplication formula)"
        if self._den != 2:
            raise NotImplementedError, "Sorry, multiplication formula only implemented for den=2"
        g = self._g
        return ThetaModForm(dup_formula(self._num_pol, g)/dup_formula(self._den_pol, g))


_dup_data_cache = []


def dup_data(g):
    """
    Writes out the duplication formula on pages 129 and 130 of Dupont's thesis.
    
    EXAMPLES::
    
        sage: from recip import *
        sage: M = Matrix([[I, 1/2+I/7], [1/2+I/7, 3/2*I-1/5]], ring=CC)
        sage: M.set_immutable()
        sage: P = theta_ring(2,2)[0]
        sage: N = 1/2 * M
        sage: N.set_immutable()
        sage: ThetaModForm(dup_data(2)[1])(N)
        0.864595401455589 - 0.0220695738100844*I
        sage: ThetaModForm(P.gens()[1])(M)^2
        0.864595401455589 - 0.0220695738100845*I
        sage: [abs(1-ThetaModForm(dup_data(2)[i])(N) / ThetaModForm(P.gens()[i])(M)^2)<10^-10 for i in [0,1,2,3,4,6,8,9,12,15]] # long time
        [True, True, True, True, True, True, True, True, True, True]
    """
    global _dup_data_cache
    for a in _dup_data_cache:
        if a[0] == g:
            return a[1]
            
    P = theta_ring(g=g, den=2)[0]
    ret = []
    for n in range(2**(2*g)):
        current_formula = 0
        c = num_to_c(n=n, g=g, den=2)
        a = c[0:g]
        b = c[g:2*g]
        for b1 in cartesian_product_iterator([[0,1/2] for i in range(g)]):
            sgn = 1
            b2 = []
            for i in range(g):
                if b1[i] == b[i]:
                    b2.append(0)
                else:
                    b2.append(1/2)
                if b1[i] == 1/2 and a[i] == 1/2:
                    sgn = sgn * -1
            current_formula = current_formula + \
                sgn*P.gens()[c_to_num([0 for i in range(g)]+list(b1), den=2)] \
                * P.gens()[c_to_num([0 for i in range(g)] + b2, den=2)]            
        ret.append(current_formula / 2**g)
    _dup_data_cache.append([g, ret])
    return ret

        
def dup_formula(pol, g):
    """
    Applies the duplication formula on pages 129 and 130 of Dupont's thesis.
    
    EXAMPLES::
    
        sage: from recip import *
        sage: M = Matrix([[I, 1/2+I/7], [1/2+I/7, 3/2*I-1/5]], ring=CC)
        sage: M.set_immutable()
        sage: P = theta_ring(2,2)[0]
        sage: N = 1/2 * M
        sage: N.set_immutable()
        sage: t = P.gens()[1]^2*P.gens()[2]^4
        sage: ThetaModForm(dup_formula(t, 2))(N)
        1.14954909189867 + 0.0102975011948174*I
        sage: ThetaModForm(t)(M)
        1.14954909189867 + 0.0102975011948173*I
    """
    pol = list(pol)
    ret = 0
    dup_dat = dup_data(g)
    for a in pol:
        e = a[1].exponents()
        if not len(e) == 1:
            raise RuntimeError, "Bug in dup_formula, monomial is not a " \
                                "monomial? %s" % a
        e = e[0]
        for i in range(len(e)):
            if not (e[i]/2 in ZZ):
                raise ValueError, "Duplication formula can only handle even " \
                                  "powers of theta's"
        t = prod([dup_dat[i]**(e[i]/2) for i in range(len(e))])
        ret = ret + a[0]*t
    return ret
    
    
Theta = ThetaProduct

def even_theta_characteristics(dupont=False):
    """
    Returns the even theta characteristics.
    If ``dupont`` is True returns numbers in Dupont's
    notation instead of theta characteristics.
    
    EXAMPLE::
    
        sage: from recip import *
        sage: even_theta_characteristics()
        [[0, 0, 0, 0], [0, 0, 1/2, 0], [0, 0, 0, 1/2], [0, 0, 1/2, 1/2], [1/2, 0, 0, 0], [1/2, 0, 0, 1/2], [0, 1/2, 0, 0], [0, 1/2, 1/2, 0], [1/2, 1/2, 0, 0], [1/2, 1/2, 1/2, 1/2]]

    """
    ns = [0, 1, 2, 3, 4, 6, 8, 9, 12, 15]
    if dupont:
        return ns
    return [num_to_c(i, 2, 2) for i in ns]

"""
RECIP -- REpository of Complex multIPlication SageMath code.
See the file README.txt for version information and instructions.

When using this package in a publication, the author would appreciate citation
of the following, or preferably a future published version. This package is
based on formulas from that paper, is motivated by it, and contains examples
from it.

- Marco Streng -- An explicit version of Shimura's reciprocity law for Siegel
  modular functions. arXiv:1201.0020

#*****************************************************************************
# Copyright (C) 2010 -- 2020 Marco Streng <marco.streng@gmail.com>
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

See below for some examples.

The current version of this package is suitable when working with theta
characteristics in (1/den)ZZ and g-dimensional abelian varieties where
den^(2*g) is small and den is even (e.g. g=den=2). This already allows for
class invariants that yield a significant improvement in size relative to Igusa
invariants.

Main functions of this file:

 * reciprocity_map_image - gives generators of the image of H(1) cap I(N)
   in GSp under the reciprocity map (small N only)

 * visualize - visualizes the action of a set of elements of GSp on theta
   constants (small N only)

"""
    
def a_to_mu(Psi, A):
    r"""
    Given an ideal A and CM-type Psi, returns mu with mu*mubar in QQ and
    N_Psi(A) = (mu) if mu exists, and None otherwise
    r"""
    C = Psi.type_norm(A)
    if not C.is_principal():
        return None
    mu = C.gens_reduced()[0]
    # This may not be the correct mu yet. But every candidate mu is epsilon*mu
    # for some unit epsilon, abbreviated eps below.
    K = mu.parent()
    conj = K.complex_conjugation()
    epsepsbar = A.norm() / mu / conj(mu)  # eps*mu is correct iff eps*epsbar = epsepsbar
    # So we need to find a unit eps with eps*epsbar = epsepsbar (if it exists)
    # and multiply mu by it.
    # Then eps*mu is the output.
    # The unit group of K0 has finite index e in that of K.
    # If eps is any unit in K and epsbar is its complex conjugate, then
    # eps^e is in K0, so epsbar^e = eps^e, hence epsbar = eps * u for
    # a root of unity u
    # We get eps*epsbar = eps^2 * u.
    # So if the correct eps exists, then there is a root of unity
    # u with epsepsbar / u = eps^2 and epsbar = eps*u.
    # For any u, eps is well-defined up to sign, but sign is irrelevant.
    for u in K.roots_of_unity():
        if (epsepsbar/u).is_square():
            eps = (epsepsbar/u).sqrt()
            if conj(eps) == u*eps:
                return eps*mu
    return None
    
    
def a_to_mus(Psi, A):
    r"""
    Given an ideal A and CM-type Psi, returns a list containing
    all mu with mu*mubar in QQ and N_Psi(A) = (mu).

    EXAMPLE::

        sage: from recip import *
        sage: K = CM_Field([5,5,5])
        sage: Phi = K.Phi()
        sage: A = K.ideal(K.gen())
        sage: l = a_to_mus(Phi, A)
        sage: l.sort(key=str)
        sage: l
        [-1/2*alphar^2 - 5/2,
         -1/8*alphar^3 + 1/8*alphar^2 - 15/8*alphar + 15/8,
         -1/8*alphar^3 - 1/8*alphar^2 - 15/8*alphar - 15/8,
         -3/8*alphar^3 + 1/8*alphar^2 - 25/8*alphar - 5/8,
         -3/8*alphar^3 - 1/8*alphar^2 - 25/8*alphar + 5/8,
         1/2*alphar^2 + 5/2,
         1/8*alphar^3 + 1/8*alphar^2 + 15/8*alphar + 15/8,
         1/8*alphar^3 - 1/8*alphar^2 + 15/8*alphar - 15/8,
         3/8*alphar^3 + 1/8*alphar^2 + 25/8*alphar - 5/8,
         3/8*alphar^3 - 1/8*alphar^2 + 25/8*alphar + 5/8]

    EXAMPLE::

        sage: from recip import *
        sage: K = CM_Field([764, 28, 5])
        sage: Phi = K.Phi()
        sage: A = K.ideal([19, K.gen()-2])
        sage: a_to_mus(Phi, A^4)
        []
        sage: a_to_mus(Phi, A^8)
        [-687/2*alphar^3 + 19177*alphar^2 - 28345*alphar + 536769,
         687/2*alphar^3 - 19177*alphar^2 + 28345*alphar - 536769]
    r"""
    mu = a_to_mu(Psi, A)
    if mu is None:
        return []
    K = mu.parent()
    return [mu*u for u in K.roots_of_unity()]


def principal_type_norms(Psi, modulus):
    r"""
    Returns a set of generators of the image of
    (H(1) cap I(modulus)) / H(modulus)
    in
    (Or / modulus)^* / O^*
    under the type norm of Psi.

    Here I(N) is the group of ideals of the domain of Psi, coprime to N,
    and H(N) is the subgroup of ideals A such that the type norm of A
    is principal and generated by mu with mu*mubar in QQ.
    
    Or is the maximal order of the reflex field.
    r"""
    K = Psi.domain()
    O = K.maximal_order()
    M = modulus * O
    
    gens = [Psi.type_norm(K(g)) for g in M.idealstar(2).gens()]
    
    if K.degree() == 2:
        gens_of_candidate_group = []
        orders_of_gens = []
    elif K.degree() > 4:
        gens_of_candidate_group = K.class_group().gens()
        orders_of_gens = [g.order() for g in gens_of_candidate_group]
        # TODO: find a smaller set of generators for g > 2.
    else:
        gens_of_class_group = K.class_group().gens()
        orders_of_gens = []
        gens_of_candidate_group = []
        #g.order() for g in gens_of_class_group]
        e = K.real_field().class_group().exponent()
        for k in range(len(gens_of_class_group)):
            g = gens_of_class_group[k]
            o = gcd(g.order(), 2*e)
            if o != 1:
                orders_of_gens.append(o)
                gens_of_candidate_group.append(g**ZZ(g.order() / o))
        
        # If an ideal A has a principal type norm,
        # then the type norm of its type norm is also principal,
        # but that is A^2 * relative_norm(A)^(real conjugation),
        # hence A^(2e) is principal, where e is as above.
        # Any such A is in the group generated by gens_of_candidate_group.
    
    for k in cartesian_product_iterator([range(o) for o in orders_of_gens]):
        if any(k):
            for c in [gens_of_candidate_group[i]**k[i] for i in range(len(k))]:
                A = c.ideal()
                # Now we need A to be coprime to "modulus"
                if not A.is_coprime(M):
                    A = A * A.idealcoprime(M)
                mu = a_to_mu(Psi, A)
                if mu != None:
                    gens.append(mu)
    return gens


def reciprocity_map_image(Z, level, modulus=None):
    r"""
    Given a CM period matrix Z, returns Z.epsilon(b) mod level
    for a set of generators A of the group (H(1) cap I(modulus)) / H(modulus).

    Here I(N) is the group of ideals of the reflex field coprime to N,
    and H(N) is the subgroup of ideals A such that the reflex type norm of A
    is principal and generated by mu with mu*mubar in QQ.
    
    The b in Z.epsilon(b) is b=mu with mu as in the definition of H(1).
    
    If modulus is None, then take modulus=level
    
    EXAMPLE::
    
        sage: from recip import *
        sage: k = CM_Field((x^2+5)^2-4*5)
        sage: Z = k.one_period_matrix(); Z
        Period Matrix
        [-0.30901699437495? + 0.95105651629515?*I -0.50000000000000? + 0.36327126400268?*I]
        [-0.50000000000000? + 0.36327126400268?*I  0.30901699437495? + 0.95105651629515?*I]
        sage: reciprocity_map_image(Z, 6) # not tested, is this answer correct?
        [[1 1 4 1]
        [2 4 1 5]
        [2 2 0 3]
        [1 2 4 2], [3 2 2 0]
        [2 3 0 2]
        [4 0 3 0]
        [2 4 0 5]]
        sage: gammas = reciprocity_map_image(Z, 8)
        sage: len(gammas)
        4
        sage: gammas[0] # not tested, is this answer correct?
        [0 2 1 6]
        [0 1 6 7]
        [7 1 2 1]
        [2 7 7 2]
    r"""
    Phi = Z.CM_type()
    Psi = Phi.reflex()
    if modulus == None:
        modulus = level
    elif level % modulus != 0:
        raise ValueError( "modulus (=%s) must divide level (=%s)" %                (modulus, level))
    elts = principal_type_norms(Psi, modulus)
    
    gammas = [GSp_element(mat_convert(Z.epsilon(b), Zmod(level))) for b in elts]
    gammas = [g for g in gammas if not g == identity_matrix(2*g.g())]
    return gammas


def list_group(gens):
    r"""
    Given a set of elements of a finite group, lists all elements of the group
    generated by the given elements.
    
    NOTE:
    Use :func:`group_generators_to_list` instead, as it is better. The
    function list_group will be removed in a later version.
    
    EXAMPLES:
    
    The matrices M and N were output by an earlier version of
    reciprocity_map_image. The outputs are equivalent, as the following test
    shows.::
    
        sage: from recip import *
        sage: k = CM_Field((x^2+5)^2-4*5)
        sage: Z = next(k.period_matrices_iter())
        sage: gens = reciprocity_map_image(Z, 6)
        sage: lst = list_group(gens) # long time
        Warning: the function list_group will be removed in a later version, use group_generators_to_list instead
        sage: M = Matrix(Zmod(6), [[3,1,2,5],[0,2,5,1],[4,2,4,5],[1,4,4,4]])
        sage: M in lst # long time, not tested, apparently wrong??
        True
        sage: N = Matrix(Zmod(6), [[3,0,4,0],[0,5,0,4],[2,4,3,2],[0,2,2,3]])
        sage: N in lst # long time, not tested, apparently wrong??
        True
        sage: lst = list_group([M,N]) # long time
        Warning: the function list_group will be removed in a later version, use group_generators_to_list instead
        sage: all([g.matrix() in lst for g in gens]) # long time, not tested, apparently wrong??
        True
    r"""
    print("Warning: the function list_group will be removed in a later "            "version, use group_generators_to_list instead")
    try:
        gens = [g.matrix() for g in gens]
    except AttributeError:
        pass
    S = [gens[0]^0]
    fin = False
    while not fin:
        fin = True
        for s in S:
            for g in gens:
                if not g*s in S:
                    fin = False
                    S.append(g*s)
    return S

    
def table(gens, den, select_rows=None):
    r"""
    Returns a table showing the action of gens on the theta's
    with characteristics in 1/den.

    EXAMPLE::

        sage: from recip import *
        sage: gens = [GSp_element([[5,0,7,1],[6,6,0,7],[3,6,6,1],[6,1,7,4]],ring=Zmod(8))]
        sage: table(gens, 2, select_rows=[0,1,4,6,8])
        [            0|            1]
        [-------------+-------------]
        [           t0|          -t6]
        [           t1|           t0]
        [           t4|   (zeta8)*t1]
        [           t6|(-zeta8^3)*t8]
        [           t8|(-zeta8^2)*t4]
        [      (zeta8)|    (zeta8^3)]
    r"""
    g = gens[0].g()
    if select_rows == None:
        if g == 2:
            select_rows = [0,1,2,3,4,6,8,9,12,15]
        else:
            select_rows = range(den**(2*g))
    P = theta_ring(gens[0].g(), den)
    zeta = P[0].base_ring().gen()
    l = []
    gens = [gens[0]**0] + gens
    for k in range(len(gens)):
        g = gens[k]
        m = g.action_on_theta_generators(den=den)
        m = [m[i] for i in select_rows]
        m = [k]+m+[zeta**ZZ(Zmod(2*den**2)(g.nu())**-1)]
        l = l + [m]
    M = Matrix(l).transpose()
    M.subdivide(1,1)
    return M

def make_zero_big(n, big):
    r"""
    Returns `big` if `n==0` and returns `n` otherwise.
    We need this because permutation groups don't handle 0 well.
    r"""
    if n == 0:
        return big
    return n


def make_big_zero(n, big):
    r"""
    Returns `0` if `n==big` and returns `n` otherwise.
    We need this because permutation groups don't handle 0 well.
    r"""
    if n == big:
        return ZZ(0)
    return n


def _permutation(M, den, L=None):
    r"""
    Given a symplectic matrix M in GSp_2g(Zmod(2*den^2)), returns
    a sequence p such that M*theta_k = theta_{p[k]} for k in L
    and p[k]=k for k not in L.
    
    If L is None, then take L = {0,...,den**(2*g)}, where 2*g is the
    number of columns of M.
    
    Note that the output only makes sense if L is stable under M.
    r"""
    B = M.base_ring()
    #g = ZZ(M.ncols()/2)
    g = M.g()
    if not (B == ZZ or (B == Zmod(B.order()) and B.order() % 2*den**2 == 0 and B.order() % 8 == 0)):
        raise ValueError( "Invalid den (=%s) for element of Sp_(%s)" % (den, B))
#    if not nu(mat_convert(M, Zmod(2*den**2))) == 1:
#        raise ValueError( "Not a symplectic matrix")
    nu_inv = lift_small(M.nu()**-1)
    M = mat_convert(M.matrix(), ZZ)
    ret = [c_to_num(theta_action_without_kappa(M,nu_inv,num_to_c(k,g,den))[0],den)             for k in range(den**(2*g))]
    if not L == None:
        for a in range(den**(2*g)):
            if not a in L:
                ret[a] = a
    return ret


def _permutations(gens, den, L=None):
    r"""
    Returns a permutation group showing the action
    of gens on the 2*den^2-th powers of the theta's
    with characteristic in 1/den.
    
    If L is not None, restrict to those theta's
    corresponding to numbers in L.
    
    Note that this works only if L is stable under
    all elements of gens.
    
    EXAMPLES::
    
        sage: from recip import *
        sage: M = []
        sage: M.append([[5,0,7,1],[6,6,0,7],[3,6,6,1],[6,1,7,4]])
        sage: M.append([[1,6,4,2],[2,3,6,4],[0,4,5,2],[4,4,6,7]])
        sage: M.append([[7,0,6,4],[2,1,0,6],[0,2,1,2],[2,2,4,3]])
        sage: M.append([[5,0,4,0],[4,1,0,4],[0,4,1,4],[4,4,0,5]])
        sage: gens = [GSp_element(m, ring=Zmod(8)) for m in M]
        sage: from recip import _permutations
        sage: _permutations(gens, 2)      
        ([[16, 15, 9, 1, 7, 8, 14, 4, 2, 13, 11, 3, 5, 10, 12, 6], [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16], [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16], [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16]], Permutation Group with generators [(), (1,16,6,8,4)(2,15,12,3,9)(5,7,14,10,13)])
    r"""
    g = gens[0].g()
    n = den ** (2*g)
    big = n
    l = []
    for g in gens:
#        M = g.Sp_part()
        p = _permutation(g, den, L)
        p = [make_zero_big(k, big) for k in p]
        p = p[1:] + [0 for k in range(big-n)] + [p[0]]
        l = l + [Permutation(p)]
    H = PermutationGroup([p.cycle_tuples() for p in l])
    return l, H


def visualize(gens, den, select_rows = None):
    r"""
    Gives a visual representation of the action of gens
    on the theta's of characteristic in ZZ/den.

    EXAMPLE::

        sage: from recip import *
        sage: M = []
        sage: M.append([[5,0,7,1],[6,6,0,7],[3,6,6,1],[6,1,7,4]])
        sage: M.append([[1,6,4,2],[2,3,6,4],[0,4,5,2],[4,4,6,7]])
        sage: M.append([[7,0,6,4],[2,1,0,6],[0,2,1,2],[2,2,4,3]])
        sage: M.append([[5,0,4,0],[4,1,0,4],[0,4,1,4],[4,4,0,5]])
        sage: gens = [GSp_element(m, ring=Zmod(8)) for m in M]
        sage: visualize(gens, 2)
        On the 8-th powers, the action is (note that 16 means 0):
        1: (1,16,6,8,4)(2,15,12,3,9)
        2: ()
        3: ()
        4: ()
        Permutation Group with generators [(), (1,16,6,8,4)(2,15,12,3,9)] of order 5
        The action on the orbit [0, 1, 4, 6, 8] is as follows
        [            0|            1             2             3             4]
        [-------------+-------------------------------------------------------]
        [           t0|          -t6            t0            t0            t0]
        [           t1|           t0            t1           -t1            t1]
        [           t4|   (zeta8)*t1            t4  (zeta8^2)*t4           -t4]
        [           t6|(-zeta8^3)*t8            t6           -t6            t6]
        [           t8|(-zeta8^2)*t4            t8 (-zeta8^2)*t8           -t8]
        [      (zeta8)|    (zeta8^3)       (zeta8)    (-zeta8^3)      (-zeta8)]
        The action on the orbit [2, 3, 9, 12, 15] is as follows
        [            0|            1             2             3             4]
        [-------------+-------------------------------------------------------]
        [           t2|          t15           -t2  (zeta8^2)*t2           -t2]
        [           t3| (zeta8^2)*t9           -t3  (zeta8^2)*t3           -t3]
        [           t9|(-zeta8^2)*t2           -t9 (-zeta8^2)*t9           -t9]
        [          t12| (zeta8^3)*t3          -t12           t12           t12]
        [          t15|(zeta8^3)*t12          -t15 (zeta8^2)*t15          -t15]
        [      (zeta8)|    (zeta8^3)       (zeta8)    (-zeta8^3)      (-zeta8)]
    r"""
    g = gens[0].g()
    if select_rows == None:
        if g == 2 and den == 2:
            select_rows = [0,1,2,3,4,6,8,9,12,15]
        else:
            select_rows = range(den**(2*g))
    big = den**(2*g)
    print("On the %s-th powers, the action is (note that %s means 0):" % (2*den**2, big))
    l, H = _permutations(gens, den, L=select_rows)
    for k in range(len(gens)):
        print(str(k+1) + ": " + l[k].cycle_string())
    print("%s of order %s" % (H, H.order()))
    for h in H.orbits():
        i = [make_big_zero(k, big) for k in h]
        j = [k for k in i if k in select_rows]
        if len(j) > 0:
            j.sort()
            print("The action on the orbit %s is as follows" % j)
            print(table(gens, den, select_rows=j))






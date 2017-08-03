"""
RECIP -- REpository of Complex multIPlication SageMath code.

this file:
orders.sage

This file implements polarized ideal classes of orders in CM-fields.
Examples of how to use this code are in bissonstreng.sage

See the file README.txt for version information and instructions.

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


def class_polynomials_order(O, F, Phi, prec, include_period_matrices=False):
    """
    Returns the class polynomials of an order in a cyclic quartic CM-field K.
    
    TODO: more general fields, examples
    
    """
    Zs = period_matrices(O, F, Phi)
    if len(Zs) == 0:
        P = ZZ['X']
        ret = [P(1),P(0),P(0)]
    else:
        i = igusa_invariants_absolute()
        numerical_values = [[j(Z, prec=prec) for Z in Zs] for j in i]
        X = numerical_values[0][0].parent()['X'].gen()
        pols_numerical = ([prod([X-v for v in numerical_values[0]])] +
                        [short_interpolation(numerical_values[0], n) for n in numerical_values[1:]])
        
        ret = [recognize_polynomial(pol_num, QQ) for pol_num in pols_numerical]
    if include_period_matrices:
        return ret, Zs
    else:
        return ret


def superorders(O):
    """
    Given an order O in a number field K, returns all orders in K containing O.
    
    NOTE: this implementation is probably not optimal, as it uses all subgroups
    of OK/O and tests whether they yield rings.
    
    EXAMPLES::
    
        sage: from recip import *
        sage: K = CM_Field([5,5,5])
        sage: OK = K.maximal_order()
        sage: O = K.order([5*b for b in OK.basis()])
        sage: s = superorders(O)
        sage: [A.index_in(OK) for A in s]
        [1, 5, 25, 25, 25, 25, 25, 25, 125]

    """
    K = O.number_field()
    t = K.degree()
    OK = K.maximal_order()
    Omatrix =  Matrix([K(g).list() for g in  O.basis()])
    # O --> K (with bases) is v --> v*Omatrix 
    OKmatrix = Matrix([K(g).list() for g in OK.basis()]) 
    # OK --> K (with bases) is v --> v*OKmatrix
    mat = Matrix(ZZ, Omatrix*OKmatrix.inverse())
    # O --> OK (with bases) is v --> v*Omatrix --> (v*Omatrix)*OKmatrix^-1
    S, U, V = mat.smith_form()
    # S = U*mat*V, S is diagonal and the others are invertible.
    # so let OS and OKS be O and OK with the Smith normal form basis
    # then OS --> OKS is v --> v*S
    # and OS --> O is v --> v*U
    # and OK --> OKS is v --> v*V
    Sdiag = S.diagonal()
    
    # This is a workaround for distinct bugs in Sage 5.6 and Sage 5.7
    non_ones = [k for k in range(len(Sdiag)) if Sdiag[k] != 1]
    if len(non_ones) == 0:
        return [O]
    omitted_ones = non_ones[0]
    A = AbelianGroup(t - omitted_ones, Sdiag[omitted_ones:])
    assert len(A.gens()) + omitted_ones == t
    # we check that these ones are the first entries of S
    assert all([a==1 for a in S.diagonal()[0:omitted_ones]])
    zeroes = [0 for i in range(omitted_ones)]

    Vinv_times_OKmatrix = V.inverse()*OKmatrix
    B0 = [g for g in S]
    ret = []
    for s in A.subgroups():
        B = B0 + [vector(zeroes + g.list()) for g in s.gens()] # these are in OKS and generate O'
        M = Matrix(B).hermite_form(include_zero_rows=False)
        # now O' --> OKS is v --> v*M
        if _is_ring(M*Vinv_times_OKmatrix, K):
            O1 = K.order([K(x*Vinv_times_OKmatrix) for x in M])
            ret.append(O1)
            assert M.determinant() == O1.index_in(OK)
    return ret


def superorders_stable_under_complex_conjugation(O, c=None):
    """
    Given an order O in a number field K, and an automorphism c of K, returns
    all orders of K containing O and fixed under c.
    
    If c is omitted and K has a method complex_conjugation, then use that for
    c.

    EXAMPLES::
    
        sage: from recip import *
        sage: K = CM_Field([5,5,5])
        sage: OK = K.maximal_order()
        sage: O = K.order([5*b for b in OK.basis()])
        sage: s = superorders_stable_under_complex_conjugation(O)
        sage: [A.index_in(OK) for A in s]
        [1, 5, 25, 25, 125]

    """
    # TODO: By using linear algebra to directly list only subgroups fixed under
    # complex conjugation,
    # we could save an order of magnitude in the complexity compared to what
    # the following does.
    if c is None:
        K = O.number_field()
        c = K.complex_conjugation()
    ret = []
    for A in superorders(O):
        if all([c(b) in A for b in A.basis()]):
            ret.append(A)
    return ret


def _is_ring(mat, K):
    """
    INPUT:
    
     - ``mat`` - a square matrix over `\QQ`, the rows of ``mat`` are
       elements of the number field `K`, expressed in terms of the standard
       basis of `K`.
    
     - `K` - a number field
     
    OUTPUT:
    
    True or False - whether the given elements of `K` are a basis of an order
    in `K`.
    """
    mat_inv = mat.inverse()
    for x in range(mat.ncols()):
        for y in range(x+1):
            prd = vector((K(mat[x])*K(mat[y])).list())
            for c in (prd*mat_inv):
                if not c in ZZ:
                    return False
    return True
          
    
def proper_ideal_classes(O, F):
    """
    Input: O an order, F non-zero in ZZ and a multiple of the conductor.
    
    Output: representatives of all classes of proper ideals of O modulo
    principal ideals.
    
    Warning: not tested much except for primitive quartic CM-fields.
    
        sage: from recip import *
        sage: K = CM_Field([5,5,5])
        sage: Phi = K.CM_types()[0]
        sage: alpha = K.gen()
        sage: O = K.order([1/2*alpha^3 + 1/2, 1/2*alpha^3 + 5/2*alpha^2 + 5/2*alpha, 5*alpha^2, alpha^3])
        sage: OK = K.maximal_order()
        sage: O.index_in(OK)
        25
        sage: s = superorders_stable_under_complex_conjugation(O)
        sage: len(s)
        3
        sage: F = 5
        sage: [ideal_contains(A.basis(), [F*b for b in OK.basis()]) for A in s]
        [True, True, True]
        sage: proper_ideal_classes(s[0], 1)
        ([[1/2*alpha^3 + 1/2, 1/2*alpha^3 + 1/2*alpha^2 + 1/2*alpha, alpha^2, alpha^3]], [[1/2*alpha^3 + 1/2, 1/2*alpha^3 + 1/2*alpha^2 + 1/2*alpha, alpha^2, alpha^3]], [[1/2*alpha^3 + 1/2, 1/2*alpha^3 + 1/2*alpha^2 + 1/2*alpha, alpha^2, alpha^3]])
        sage: proper_ideal_classes(s[0], 5) == _ # long time
        True
        sage: proper_ideal_classes(s[1], 5) # long time
        ([[1/2*alpha^3 + 1/2, 1/2*alpha^3 + 1/2*alpha^2 + 5/2*alpha, alpha^2, alpha^3]], [[1/2*alpha^3 + 1/2, 5/2*alpha^3 + 5/2*alpha^2 + 5/2*alpha, alpha^2, alpha^3]], [[1/2*alpha^3 + 1/2, 1/2*alpha^3 + 1/2*alpha^2 + 5/2*alpha, alpha^2, alpha^3]])
        sage: proper_ideal_classes(s[2], 5) # long time
        ([[1/2*alpha^3 + 1/2, 1/2*alpha^3 + 1/2*alpha^2 + 1/2*alpha, 5*alpha^2, alpha^3], [1/2*alpha^3 + 1/2, 1/2*alpha^3 + 5/2*alpha^2 + 5/2*alpha, 5*alpha^2, alpha^3]], [[1/2*alpha^3 + 1/2, 1/2*alpha^3 + 1/2*alpha^2 + 1/2*alpha, 5*alpha^2, alpha^3], [1/2*alpha^3 + 1/2, 5/2*alpha^3 + 5/2*alpha^2 + 5/2*alpha, 5*alpha^2, alpha^3]], [[1/2*alpha^3 + 1/2, 1/2*alpha^3 + 5/2*alpha^2 + 5/2*alpha, 5*alpha^2, alpha^3]])
    """
    # We start with the class group of the maximal order:
    K = O.number_field()
    C = K.class_group()
    Fidl = K.ideal(F)
    g = [A.ideal() for A in C]
    g = [A*A.idealcoprime(Fidl) for A in g]
    # So far these are integral ideals of OK that are coprime to F,
    # we will make them into ideals of O.
    # Given AK = A*OK with A an ideal of O coprime to F, how to find A?
    # AK*F is contained in A, but not coprime to F
    # N := norm(AK), then N*O is coprime to F, and N*OK is contained in AK,
    # Now consider the integral ideal N*O + F*AK. Multiplication with OK yields
    # AK = N*AK + F*AK subset N*OK+F*AK subset AK, so we have found A.
    g = [additive_gens_to_basis([AK.norm()*b for b in O.basis()] + (F*AK).basis(), K) for AK in g]

    from_maximal_order = g
    
    # We now have (the basis of) one invertible ideal coprime to F, for each
    # class of invertible ideal
    # of the maximal order. Next we want to have every ideal I that tensors
    # to give OK itself back.

    OK = K.maximal_order()
    n = K.degree()
    
    # Note: F*OK = F*(OK*I) = (F*OK)*I \subset O*I = I

    # We compute all subgroups of FOK/OK, as our ideals I correspond to
    # certain of those subgroups.
    B = [g*F for g in OK.basis()]
    A = AbelianGroup(n, [F for i in range(n)])
    S = A.subgroups()

    # We now create matrices for multiplication by generators of O
    # with respect to bases of OK.
    # This will help us recognize ideals of O among the subgroups of A.
    IK_matrix = Matrix([K(g).list() for g in OK.basis()])
    # IK --> K is s --> s*IK_matrix
    IK_matrix_inverse = IK_matrix.inverse()
    mult_by_O = []
    for a in O.basis():
        a_matrix = Matrix([(K(a)*K(g)).list() for g in OK.basis()])
        # IK --> K through multiplication by a is s --> s*a_matrix
        a_matrix_IK = a_matrix*IK_matrix_inverse
        # IK --> IK through multiplication by a is s --> s*a_matrix_IK
        mult_by_O.append(a_matrix_IK)
    
    F_primary_bases = []
    
    # Next, we create matrices for multiplication by generators of O'
    # (with respect to bases of OK) for all superorders O' of O. This
    # will allow us to eliminate non-proper ideals.
    strict_superorders = [(O.index_in(a), a) for a in superorders(O) if a != O]
    strict_superorders.sort()
    minimal_superorders = []
    minimal_superorder_matrices = []
    for C in strict_superorders:
        c = C[1]
        if all([not b.is_suborder(c) for b in minimal_superorders]):
            minimal_superorders.append(c)
            mult_by_c = []
            for a in c.basis():
                a_matrix = Matrix([(K(a)*K(g)).list() for g in OK.basis()])
                # IK --> K through multiplication by a is s --> s*a_matrix
                a_matrix_IK = a_matrix*IK_matrix_inverse
                # IK --> IK through multiplication by a is s --> s*a_matrix_IK
                mult_by_c.append(a_matrix_IK)
            minimal_superorder_matrices.append(mult_by_c)
    if get_recip_verbose():
        print "Number of minimal superorders: %s" % len(minimal_superorder_matrices)
    
    for s in S:
        gns = ([g.list() for g in F*identity_matrix(4)] +
               [g.list() for g in s.gens()])
        M = Matrix(gns).echelon_form(include_zero_rows=False)
        if M.nrows() != M.ncols():
            raise RuntimeErorr,  "Matrix M=%s not square, it came from %s" % (M, gens)
        # I --> IK is s --> s*M
        gns = [K(r) for r in M*IK_matrix]
        actual_IK = K.ideal(gns)
        if (actual_IK == K.ideal(1) and
            _is_proper_ideal(mult_by_O, minimal_superorder_matrices, M)):
            F_primary_bases.append(gns)

    # Now F_primary_bases contains the bases of the proper ideals I with
    # I*OK = OK.

    # Next, take only one representative of each OK^*-orbit in F_primary_bases.
    # We can take OK^*/O^* here as O^* acts trivially.
    # We enumerate all elements of OK^*/O^*.
    unit_gens = K.unit_group().gens()
    unit_gens = [K(g) for g in unit_gens]
    unit_reps = [K(1)]
    for i in range(len(unit_gens)):
        # At this stage, unit_reps contains exactly one representative of every
        # element of the subgroup of OK^*/O^* generated by unit_gens[:i].
        u = unit_gens[i]
        # We will find the order e of u in OK^*/O^*/unit_reps and compute
        # the new unit_reps
        new_unit_reps = unit_reps
#        print "new_unit_reps" + str(new_unit_reps)
#        sys.stdout.flush()
        # e = 1
        v = u
        while all([not(u*t in O) for t in unit_reps]):
            new_unit_reps = new_unit_reps + [u*t for t in unit_reps]
#            print "new_unit_reps" + str(new_unit_reps)
#            sys.stdout.flush()
            u = u*v
            # e = e+1
        unit_reps = new_unit_reps
#        print "unit_reps" + str(unit_reps)
#        sys.stdout.flush()
    
    lst = []
    
    times_found = [0 for f in F_primary_bases]
    for k in range(len(F_primary_bases)):
        if not times_found[k]:
            A = F_primary_bases[k]
            lst.append(A)
            for e in unit_reps:
                B = [e*b for b in A]
                for l in range(len(F_primary_bases)):
                    if ideals_equal(B, F_primary_bases[l]):
                        times_found[l] = times_found[l] + 1

    if not times_found == [1 for f in F_primary_bases]:
        print times_found
    if not len(lst)*len(unit_reps) == len(F_primary_bases):
        print len(lst), len(unit_reps), len(F_primary_bases)

    F_primary_bases = lst


    out = [mult_ideals(a,b,K) for a in F_primary_bases for b in from_maximal_order]
    
    return (out, F_primary_bases, from_maximal_order)


def polarized_ideal_classes(O, F):
    """
    On input an order O of a quartic CM-field, returns all principally
    polarized ideal
    classes up to
      * changing (aaa, xi) into (aaa, -xi),
      * changing (aaa, xi) into (x*aaa, xi/x/xbar) for x in K.
    In the primitive quartic case, this is either the empty set or
    half of a
    torsor for the group CCC(O) in the article of Bisson and Streng.
    (Half of, because only one of (aaa, xi) and (aaa, -xi) is given.)
    If O is the maximal order of a primtive quartic CM-field, then the
    output is non-empty.
    
    INPUT:
    
     - `O` -- order in a quartic `CM_Field` object
     - `F` -- positive integer such that F*OK is contained in O (not checked)
    
    EXAMPLES::
    
        sage: from recip import *
        sage: K = CM_Field([5,5,5])
        sage: Phi = K.CM_types()[0]
        sage: alpha = K.gen()
        sage: O = K.order([1/2*alpha^3 + 1/2, 1/2*alpha^3 + 5/2*alpha^2 + 5/2*alpha, 5*alpha^2, alpha^3])
        sage: OK = K.maximal_order()
        sage: O.index_in(OK)
        25
        sage: s = superorders_stable_under_complex_conjugation(O)
        sage: len(s)
        3
        sage: F = 5
        sage: [ideal_contains(A.basis(), [F*b for b in OK.basis()]) for A in s]
        [True, True, True]
        sage: polarized_ideal_classes(s[0], 1)
        [([1/2*alpha^3 + 1/2, 1/2*alpha^3 + 1/2*alpha^2 + 1/2*alpha, alpha^2, alpha^3],
          1/5*alpha),
         ([1/2*alpha^3 + 1/2, 1/2*alpha^3 + 1/2*alpha^2 + 1/2*alpha, alpha^2, alpha^3],
          -1/5*alpha^3 - 3/5*alpha)]
        sage: polarized_ideal_classes(s[0], 5) == _ # long time
        True
        sage: polarized_ideal_classes(s[1], 5) # long time
        [([1/2*alpha^3 + 1/2, 1/2*alpha^3 + 1/2*alpha^2 + 5/2*alpha, alpha^2, alpha^3],
          2/25*alpha^3 + 1/5*alpha),
         ([1/2*alpha^3 + 1/2, 1/2*alpha^3 + 1/2*alpha^2 + 5/2*alpha, alpha^2, alpha^3],
          1/25*alpha^3 + 1/5*alpha)]
        sage: polarized_ideal_classes(s[2], 5) # long time
        []
    """
    a,_,_ = proper_ideal_classes(O,F)
    # we also need a generator of N_{K/K0}(O^*):
    K = O.number_field()
    cc = K.complex_conjugation()
    [z, u] = K.unit_group().gens()
    z = K(z)
    u = K(u)
    e = ZZ(z.multiplicative_order() / 2)
    k = 1
    while all([not (z**j * u**k in O) for j in range(e)]):
        k = k + 1
    # now z**j * u**k generates the free part of O^*, so the relative
    # norm of u**k generates N_{K/K0}(O^*).
    # That relative norm is z**l * u**(2*k) for some l, so the units of
    # OK^*/N_{K/K0}(O^*)/<-1>
    # are:
    units_mod_relative_norms = [z**j * u**m for j in range(e) for m in range(2*k)]
    if get_recip_verbose():
        print "units modulo relative norms: %s" % units_mod_relative_norms
    ret = []
    
    for A in a:
        # print "ideal %s" % A
        # we need xi*A to be the trace dual of Abar
        Abardual = trace_dual_complex_conjugate(A)
        # xi*A = Abardual, so xi*A*OK = Abardual*OK
        xi_ideal = K.ideal(Abardual) / K.ideal(A)
        if xi_ideal.is_principal():
            if get_recip_verbose() > 2:
                print "ideal is principal"
            xi = xi_ideal.gens_reduced()[0]
            for u in units_mod_relative_norms:
                if ideals_equal([xi*u*b for b in A], Abardual):
                    if get_recip_verbose() > 3:
                        print "found xi generating the correct ideal"
                    if cc(xi*u) == -xi*u:
                        if get_recip_verbose() > 3:
                            print "xi is totally imaginary"
                        ret.append((A, xi*u))
                        #print A, xi, u, xi*u
                    elif get_recip_verbose() > 2:
                        print "xi is not totally imaginary"
                elif get_recip_verbose() > 2:
                    print "This generator (%s) of AbardualOK / AOK (%s, %s) does not map A to Abardual (%s), but to (%s): " % (xi*u, xi*u*K.ideal(A) == K.ideal(Abardual), ideals_equal([xi*u*b for b in K.ideal(A).basis()], K.ideal(Abardual).basis()), Abardual, [xi*u*b for b in A]) 
        elif get_recip_verbose() > 2:
            print "ideal is not principal"
    return ret
                

def trace_dual_complex_conjugate(A, cc=None):
    """
    Given a basis A of a lattice, returns the trace dual of cc(A).
    If cc is not specified and K has a complex_conjugation method, use that.
    """
    K = A[0].parent()
    if cc is None:
        cc = K.complex_conjugation()
    M = Matrix([cc(c).list() for c in A]) # s |-> s*M  is  A --> K
    trace_pairing = Matrix([[(K.gen()**(i+j)).trace() for i in range(K.degree())] for j in range(K.degree())])
    # trace(s*y) for s in A and y in K is: s*M*trace_pairing*y^transpose
    # we want a matrix B such that M*trace_pairing*B = identity, so
    N = (M*trace_pairing).inverse()
    # now s^transpose |-> N*s^transpose  is  B --> K
    #     s |-> s*Ntranspose
    return [K(r) for r in N.transpose()]
    

def ideal_index(B, A, check_included=True):
    """
    Given two ideals B and A, assuming B contains A,
    returns the index of A in B.
    
    Here A and B are given by their bases, which are sequences of n elements
    of K, where n=[K:QQ] and A and B are lattices in K.
    
    If check_included is False, then B does not have to contain A, and
    the output is the quotient of the indices in B+A.
    """
    if check_included and not ideal_contains(B, A):
        raise ValueError
    return abs(Matrix([a.list() for a in A]).determinant() / Matrix([b.list() for b in B]).determinant())


def ideal_contains(B, A):
    """
    Given bases of two lattices B and A in K, returns true
    if and only if the lattice B contains the lattice A.
    """
    # s |-> s*Matrix(A) is the map A --> K
    M = Matrix([a.list() for a in A])*Matrix([b.list() for b in B]).inverse()
    for r in M:
        for c in r:
            if not c in ZZ:
                return False
    return True


def ideals_equal(A, B, K=None):
    """
    Given bases of two lattices A and B in K, returns true
    if and only if the lattices are equal.
    """
    # s |-> s*Matrix(A) is the map A --> K
    if K is None:
        K = Sequence(A+B).universe()

    M = Matrix([K(a).list() for a in A])*Matrix([K(b).list() for b in B]).inverse()
    for r in M:
        for c in r:
            if not c in ZZ:
                return False
    if not M.determinant() in [-1,1]:
        return False
    return True
    
    
def additive_gens_to_basis(gens, K):
    """
    Given generators of a lattice L in K, returns a basis of L.
    """
    M = Matrix([K(g).list() for g in gens])
    den = lcm([lcm([c.denominator() for c in r]) for r in M])
    N = Matrix(ZZ, den*M).hermite_form(include_zero_rows=False)
    if N.nrows() != N.ncols():
        print "Matrix %s is not square, so the lattice does not have full rank. The input was %s" % (N, gens)
    return [K(g)/den for g in N]


def mult_ideals(basis1, basis2, K):
    """
    Given bases of lattices in K, returns a basis of their product.

    EXAMPLES::

        sage: from recip import *
        sage: K.<a> = QuadraticField(15)
        sage: basis1 = [1, a]   # Maximal order OK
        sage: basis2 = [2, 2*a] # 2OK
        sage: basis3 = [1, 2*a] # Another order, O
        sage: basis4 = [2, a]   # Fractional ideal. Note: 4*basis4 is an ideal of basis3.
        sage: mult_ideals(basis1, basis2, K) # basis2 multiplied by its endomorphism ring, hence unchanged
        [2, 2*a]
        sage: mult_ideals(basis1, basis3, K) # Two orders multiplied, hence the largest is obtained
        [1, a]
        sage: mult_ideals(basis2, basis3, K) # Twice the above
        [2, 2*a]
        sage: mult_ideals(basis4, basis3, K) # basis4 multiplied by its endomorphism ring, hence unchanged
        [2, a]
        sage: mult_ideals(basis4, basis4, K) # the square of basis4, generated by 4, 2a and a^2 = 15.
        [1, 2*a]
    """
    return additive_gens_to_basis([b1*b2 for b1 in basis1 for b2 in basis2], K)


def intersect_ideals(basis1, basis2, K):
    """
    Given bases of lattices in K, returns a basis of their intersection.

    EXAMPLES::

        sage: from recip import *
        sage: K.<a> = QuadraticField(15)
        sage: basis1 = [1, a]   # Maximal order OK
        sage: basis2 = [2, 2*a] # 2OK
        sage: basis3 = [1, 2*a] # Another order, O
        sage: basis4 = [2, a]   # Fractional ideal. Note: 4*basis4 is an ideal of basis3.
        sage: intersect_ideals(basis1, basis2, K) # ideal2 intersected with an ambient ring, hence unchanged
        [2, 2*a]
        sage: intersect_ideals(basis1, basis3, K) # Two orders intersected, hence the smallest is obtained
        [1, 2*a]
        sage: intersect_ideals(basis2, basis3, K) # ideal2 intersected with an ambient ring, hence unchanged
        [2, 2*a]
        sage: intersect_ideals(basis4, basis3, K) # intersection
        [2, 2*a]
        sage: intersect_ideals(basis4, basis4, K) # intersection with itself
        [2, a]
    """
    basis1 = [K(b).list() for b in basis1]
    basis2 = [K(b).list() for b in basis2]
    n = K.degree()
    if len(basis1) != n or len(basis2) != n:
        raise ValueError
    ideal1 = span(basis1, ZZ)
    ideal2 = span(basis2, ZZ)
    ideal = ideal1.intersection(ideal2)
    basis = ideal.basis()
    assert len(basis) == n
    return [K(b) for b in basis]


def ideal_divide(basis1, basis2, K):
    """
    Given bases of ideals A1 and A2, returns a basis
    of (A1 : A2) = {a in K : a*A2 \subset A1}.

    EXAMPLES::

        sage: from recip import *
        sage: K.<a> = QuadraticField(15)
        sage: basis1 = [1, a]   # Maximal order OK
        sage: basis2 = [2, 2*a] # 2OK
        sage: basis3 = [1, 2*a] # Another order, O
        sage: basis4 = [2, a]   # Fractional ideal. Note: 4*basis4 is an ideal of basis3.
        sage: ideal_divide(basis1, basis2, K) # OK / 2OK = (1/2)OK
        [1/2, 1/2*a]
        sage: ideal_divide(basis1, basis3, K) # contains OK as OK*O = OK, and is contained in OK as it multiplies 1 into OK
        [1, a]
        sage: ideal_divide(basis2, basis3, K) # twice the above
        [2, 2*a]
        sage: ideal_divide(basis4, basis3, K) # ideal divided by its own order
        [2, a]
        sage: ideal_divide(basis4, basis4, K) # endomorphism ring of ideal 4
        [1, 2*a]

    """
    ideals = [[b^-1 * c for c in basis1] for b in basis2]
    ideal = ideals[0]
    for I in ideals[1:]:
        ideal = intersect_ideals(ideal, I, K)
    return ideal


def ideal_order(basis, K):
    """
    Given a basis of an ideal, returns the endomorphism ring inside K
    of this ideal, as an order.

    EXAMPLES::

        sage: from recip import *
        sage: K.<a> = QuadraticField(15)
        sage: basis1 = [1, a]   # Maximal order OK
        sage: basis2 = [2, 2*a] # 2OK
        sage: basis3 = [1, 2*a] # Another order, O
        sage: basis4 = [2, a]   # Fractional ideal. Note: 4*basis4 is an ideal of basis3.
        sage: O1 = ideal_order(basis1, K); O1
        Order in Number Field in a with defining polynomial x^2 - 15
        sage: O1.basis()
        [1, a]

        sage: ideal_order(basis2, K) == O1
        True

        sage: O3 = ideal_order(basis3, K)
        sage: O3.basis()
        [1, 2*a]
        sage: O1 == O3
        False

        sage: ideal_order(basis4, K) == O3
        True

        sage: is_ideal(O3, basis3)
        True
        sage: is_ideal(O1, basis3)
        False
        sage: is_ideal(O1, basis1)
        True
        sage: is_ideal(O3, basis1)
        True


    """
    return K.order(ideal_divide(basis, basis, K))


def ideal_inverse(basis, K, O=None):
    """
    Given a basis of an ideal A, returns
    a basis of (O : A) = {a in K : a*A \subset O}.

    If O is not specified, use ideal_order(basis, K).

    EXAMPLES::

        sage: from recip import *
        sage: K.<a> = QuadraticField(15)
        sage: basis1 = [1, a]   # Maximal order OK
        sage: basis2 = [2, 2*a] # 2OK
        sage: basis3 = [1, 2*a] # Another order, O
        sage: basis4 = [2, a]   # Fractional ideal. Note: 2*basis4 is an ideal of basis3.
        sage: ideal_inverse(basis1, K) # The order itself
        [1, a]
        sage: ideal_inverse(basis2, K) # (1/2)OK
        [1/2, 1/2*a]
        sage: ideal_inverse(basis3, K) # the order itself
        [1, 2*a]
        sage: basis5 = ideal_inverse(basis4, K); basis5 # happens to be self-inverse
        [2, a]
        sage: mult_ideals(basis4, basis5, K)
        [1, 2*a]
        sage: ideal_order(basis4, K).basis()
        [1, 2*a]

        sage: ideal_is_invertible(basis2, O = K.order(basis3))
        False
        sage: ideal_is_invertible(basis2, O = K.order(basis1))
        True
        sage: bases = [basis1, basis2, basis3, basis4]
        sage: all([ideal_is_invertible(b) for b in bases]) # quadratic ring, hence all invertible for some ring
        True

    Special case of Exercise 2.20 of
    http://websites.math.leidenuniv.nl/algebra/ant.pdf
    ::

        sage: from recip import *
        sage: K.<a> = NumberField(x^3+x+1)
        sage: A = K.order([1, a, a^2])
        sage: p = 7
        sage: R = K.order([1, p*a, p*a^2])
        sage: I = additive_gens_to_basis([1, p*a, p*a^2, a, p*a^2, p*a^3], K); I
        [1, a, 7*a^2]
        sage: ideal_order(I, K) == R
        True
        sage: Isq = mult_ideals(I, I, K); Isq
        [1, a, a^2]
        sage: ideal_order(Isq, K) == A
        True

    The final part of the exercise is to show that the above
    implies that I is a proper, but non-invertible R-ideal. In particular,
    it is for all orders non-invertible.
    We compute the ``inverse'' of this non-invertible ideal::

        sage: Iinv = ideal_inverse(I, K); Iinv
        [7, 7*a, 7*a^2]
        sage: mult_ideals(I, Iinv, K)
        [7, 7*a, 7*a^2]
        sage: ideal_is_invertible(I)
        False
        sage: ideal_is_invertible(Isq, O=A)
        True
        sage: ideal_is_invertible(Isq, O=R)
        False
        sage: ideal_is_invertible(Isq)
        True

    Something from work with Sorina Ionica, Chloe Martindale
    and Damien Robert::

        sage: from recip import *
        sage: K = CM_Field([40, 20, 90])
        sage: alpha = K.gen()
        sage: I1 = [2*alpha^3 + 2, 220/3*alpha^3 + 2/3*alpha, 3*alpha^3 + 3*alpha^2,    111*alpha^3]
        sage: xi1 = -1/5328*alpha^3 - 31/26640*alpha
        sage: I2 = [28*alpha^3 + 2*alpha^2 + 10, 61/3*alpha^3 + 2*alpha^2 + 2/3*alpha,   15*alpha^3 + 6*alpha^2,   39*alpha^3]
        sage: xi2 = -1/93600*alpha^3 + 1/9360*alpha
        sage: I3 = [2*alpha^3 + 2, 4/3*alpha^3 + 2*alpha^2 + 2/3*alpha, 6*alpha^2, 3*alpha^3]
        sage: xi3 = 1/1440*alpha^3 + 7/720*alpha

        sage: O = ideal_order(I1, K)
        sage: O == ideal_order(I2, K)
        True
        sage: O == ideal_order(I3, K)
        True
        sage: O.index_in(K.maximal_order())
        27

        sage: I1inv = ideal_inverse(I1, K)
        sage: P1 = mult_ideals(I1, I1inv, K); P1
        [3, 2*alpha^3 + alpha, 3*alpha^2, 3*alpha^3]
        sage: O1 = ideal_order(P1, K); O1 == O
        False
        sage: O1.index_in(K.maximal_order())
        1

        sage: I2inv = ideal_inverse(I2, K)
        sage: P2 = mult_ideals(I2, I2inv, K); P2
        [3, 2*alpha^3 + alpha, 3*alpha^2, 3*alpha^3]
        sage: O2 = ideal_order(P2, K)
        sage: O2 == O1
        True

        sage: I3inv = ideal_inverse(I3, K)
        sage: P3 = mult_ideals(I3, I3inv, K); P3
        [3, 2*alpha^3 + alpha, 3*alpha^2, 3*alpha^3]
        sage: O3 = ideal_order(P3, K)
        sage: O3 == O1
        True

    """
    if O is None:
        O = ideal_order(basis, K)
    try:
        O = O.basis()
    except AttributeError:
        pass
    O = [K(b) for b in O]
    return ideal_divide(O, basis, K)


def is_ideal(O, B):
    r"""
    Given a basis B of a lattice L in a field K, and an order O in K,
    return True if and only if `O L \subset L`.
    """
    K = O.number_field()
    M = Matrix([K(b).list() for b in B])
    
    mult_by_O = [Matrix([(K(b)*K.gen()**i).list() for i in range(K.degree())]) for b in O.basis()]
    # each element M of mult_by_O gives s |-> s*M as a map K --> K
    
    return _is_ideal(mult_by_O, M)
    

def _is_ideal(mult_by_O, M):
    """
    Given M: I --> IK and mult_by_O the multiplication by O matrices on IK,
    returns True if and only if I is an O-module.
    
    Here IK can be K itself or a sublattice, with any fixed basis.
    """
    for a in mult_by_O:
        # I --> IK --> IK --> I by multiplication by a is:
        # s --> s*M*a*M^-1
        for r in M*a*M.inverse():
            for c in r:
                if not c in ZZ:
                    return False
    return True


def _is_proper_ideal(mult_by_O, minimal_superorder_matrices, M):
    if not _is_ideal(mult_by_O, M):
        if get_recip_verbose() > 2:
            print "Is not an ideal"
        return False
    for c in minimal_superorder_matrices:
        if _is_ideal(c, M):
            if get_recip_verbose() > 2:
                print "Is an ideal of the over-order"
            return False
    return True


def ideal_is_invertible(basis, K=None, O=None):
    """
    Returns True if and only if the ZZ-lattice I spanned by basis is an
    invertible O-ideal.

    If O is not specified, returns True if and only if there is an order O
    such that I is an invertible O-ideal. Equivalently, it then returns True
    if and only if it is an invertible ideal of its endomorphism order.

    See :func:`ideal_inverse` for examples.
    """
    if K is None:
        K = Sequence(basis).universe()
    if O is None:
        O = ideal_order(basis, K)
    elif O.number_field() != K:
        raise ValueError
    Iinv = ideal_inverse(basis, K, O)
    A = mult_ideals(basis, Iinv, K)
    return ideals_equal(A, O.basis(), K)


def period_matrices(O, F, Phi=None, reduced=True):
    """
    Returns all isomorphism classes of principally polarized abelian surfaces
    with CM by O of type Phi. Here F is such that F*OK is contained in O.
    
    If Phi is None, use all CM-types and gives all elements
    up to complex conjugation of the CM type.
    
    EXAMPLES:
    
    Here is a bug, some numbers are close to zero and apparently CLF's > 0
    is not to be trusted.::
    
        sage: from recip import *
        sage: K = CM_Field([1837, 81, 1181])
        sage: Phi = K.CM_types()[0]
        sage: len(period_matrices(K.maximal_order(), 1, Phi)) # long time, half a minute
        66

    """
    pl = polarized_ideal_classes(O, F)
    if Phi is None:
        ret = [PeriodMatrix(None, a[0], a[1]) for a in pl]
    if not Phi is None:
        pl = pl + [(a[0], -a[1]) for a in pl]
        ret = [PeriodMatrix(Phi, aaa, xi) for (aaa, xi) in pl
                                          if Phi.is_totally_positive_imaginary(xi)]
    if reduced:
        ret = [a.reduce(CC) for a in ret]
    return ret


def are_nn_isogenous(Z1, Z2, n, F1, F2, transformation=False, double_check=False):
    """
    Given two CM period matrices with the same CM-type, returns True if and
    only if there is an (n,n)-isogeny between them. Assumes the orders Oi of Zi
    to satisfy Fi*OK \subset Oi.

    If double_check, also checks whether there is an isogeny the other way
    (which should be equivalent). Raises an error if there is an isogeny in only
    one direction.
    
    EXAMPLES::
    
        sage: from recip import *
        sage: K = CM_Field([5,5,5])
        sage: Phi = K.CM_types()[0]
        sage: alpha = K.gen()
        sage: orders = [K.order([1/2*alpha^3 + 1/2, 1/2*alpha^3 + 5/2*alpha^2 + 5/2*alpha, 5*alpha^2, alpha^3]), K.order([1, 2*alpha, 2*alpha^2, alpha^3]), K.order([3/2*alpha^3 + 1/2, 3/2*alpha^3 + 1/2*alpha^2 + 3/2*alpha, alpha^2, 3*alpha^3])]
        sage: OK = K.maximal_order()
        
        sage: O = orders[0]
        sage: O.index_in(OK)
        25
        sage: s = superorders_stable_under_complex_conjugation(O)
        sage: len(s)
        3
        sage: F = 5
        sage: [ideal_contains(A.basis(), [F*b for b in OK.basis()]) for A in s]
        [True, True, True]
        sage: Zs = [period_matrices(A, F, Phi) for A in s] # long time, half a minute
        sage: len(Zs[0]) == 1 and s[0] == OK # long time
        True
        sage: Zmax = Zs[0][0] # long time
        sage: len(Zs[2]) == 0 and len(Zs[1]) == 1 # long time
        True
        sage: Zother = Zs[1][0] # long time
        sage: are_nn_isogenous(Zmax, Zother, 5, 1, F, transformation=True) # long time
        (True, 1/10*alpha^3 - 1/2)
        sage: b, t = _ # long time
        sage: ideal_index([b*t for b in Zmax.basis()], Zother.basis()) # long time
        25
        sage: cc = K.complex_conjugation() # long time
        sage: Zother.xi()*t*cc(t)*5 == Zmax.xi() # long time
        True

    Some checks for the are_nn_isogenous function::
    
        sage: are_nn_isogenous(Zmax, Zmax, 5, 5, 5, transformation=True) # long time
        (True, -1/5*alpha^3 - 1/10*alpha^2 - 1/2*alpha - 1/2)
        sage: are_nn_isogenous(Zmax, Zmax, 1, 1, 1, transformation=True) # long time
        (True, 1)
        sage: are_nn_isogenous(Zother, Zother, 1, F, F, transformation=True) # long time
        (True, 1)
        sage: are_nn_isogenous(Zmax, Zother, 3, 1, F) # long time
        False
        sage: are_nn_isogenous(Zmax, 7*Zmax, 7, 1, 7) # todo, multiplication of period matrices by integers is not implemented
        True

    Another example::
    
        sage: O = orders[2]
        sage: O.index_in(OK)
        9
        sage: s = superorders_stable_under_complex_conjugation(O)
        sage: len(s)
        2
        sage: F = 3
        sage: [ideal_contains(A.basis(), [F*b for b in OK.basis()]) for A in s]
        [True, True]
        sage: Zs = [period_matrices(A, F, Phi) for A in s]
        sage: len(Zs[0]) == 1 and s[0] == OK
        True
        sage: Zmax = Zs[0][0]
        sage: len(Zs[1]) == 2
        True
        sage: Zothers = Zs[1]
        sage: are_nn_isogenous(Zmax, Zothers[0], 3, 1, F, transformation=True)
        (True, 1)
        sage: are_nn_isogenous(Zmax, Zothers[1], 3, 1, F, transformation=True)
        (True, 1)

    The 2-part will be dealt with for all fields together.

    """
    if Z1.CM_type() != Z2.CM_type():
        raise NotImplementedError, "are_nn_isogenous is only implemented for period matrices of equal CM type"
    if double_check:
        b1, alpha = are_nn_isogenous(Z1, Z2, n, F1, F2, transformation=True, double_check=False)
        b2 = are_nn_isogenous(Z2, Z1, n, F2, F1, transformation=False, double_check=False)
        if b1 != b2:
            raise RuntimeError, "Bug: contradiction with symmetry of being isogenous. Z1 = " + str(Z1) + " Z2 = " + str(Z2) + "(n, F1, F2) = " + str((n, F1, F2))
        if transformation:
            return b1, alpha
        return b1
    # Being (n,n)-isogenous is a symmetric relation. Write Zi = (Ai, xii).
    #
    # alpha is an (n,n)-isogeny if and only if both of the following hold:
    # * A2/alpha is contained in A1 with index n^g and quotient group (ZZ/nZZ)^2,
    # * xi1 = n*xi2*alpha*alphabar.
    #
    # We list all candidates alpha and check them. So suppose alpha is as above.
    if not is_squarefree(n):
        raise NotImplementedError
    rel_norm_alpha = Z1.xi()/Z2.xi()/n # alpha*alphabar as an element of K
    K = Z1.CM_field()
    cc = K.complex_conjugation()
    A1K = K.ideal(Z1.basis())
    A2K = K.ideal(Z2.basis())
    zeta = K(K.unit_group().gens()[0])
    zetaorder = zeta.multiplicative_order()
    roots_of_unity_mod_minus_one = [zeta**k for k in range(ZZ(zetaorder/2))]
    # Note that alpha satisfies the following:
    #   n*A1 \subset A2/alpha \subset A1
    # In particular, we get
    #   F1*n*alpha*A1K \subset n*alpha*A1 \subset A2 \subset A2K
    # AND
    #   F2*A2K \subset A2 \subset alpha*A1 \subset alpha*A1K
    # so
    # (NOT) F1*n*A1K/A2K \subset OK*alpha \subset F2^-1*A1K/A2K
    # (BUT) F2*A2K/A1K \subset alpha*OK\subset A2K/A1K/(F1*n)
    base = (F1*n)**-1 * A2K/A1K
    # In other words:
    #   alpha*OK     \subset base,     and
    #   F2*F1*n*base \subset alpha*OK
    # Let alpha*OK = mult*base, then
    #   mult in OK
    #   mult divides F2*F1*n
    # But also, don't forget alpha*alphabar = rel_norm_alpha,
    # which in particular implies (mult*base).norm()^2 == alpha.norm()^2
    #                                                  == rel_norm_alpha.norm()
    for mult in [a for a in divisors(K.ideal(F1*F2*n))
                   if (a.norm()*base.norm())^2 == rel_norm_alpha.norm()]:
        # We're starting a slow search, obviously it is much better to do
        # linear algebra on the powers of the primes dividing F1*F2*n,
        # but since these numbers are all small anyway in our application, we
        # don't bother in this implementation.
        alphaidl = mult*base
        if alphaidl.is_principal():
            # Take an arbitrary generator alpha0 of the ideal (alpha)
            alpha0 = alphaidl.gens_reduced()[0]
            assert K.ideal(alpha0) == alphaidl
            # Again, don't forget alpha*alphabar = rel_norm_alpha, which implies:
            if (cc(alphaidl)*alphaidl) == K.ideal(rel_norm_alpha):
                u = rel_norm_alpha / alpha0 / cc(alpha0) # element of K
                # We need a unit v of relative norm u, and take
                # alpha = alpha0 * v, because then
                # alpha generates the correct ideal and satisfies
                # n*xi2*alpha*alphabar = n*xi2*u*alpha0*alpha0bar
                #                      = n*xi2*rel_norm_alpha
                #                      = xi1
                # Note v = z*v0 with v0 in K0, and z a root of unity
                # (by a Lemma from my thesis), so u = v*vbar = z^2
                u_in_K0 = K.self_to_real(u)
                if u_in_K0.is_square():
                    v0 = K.real_to_self(u_in_K0.sqrt())
                    # Next, z is only relevant up to sign, hence only
                    # if K = QQ(zeta_5).
                    for z in roots_of_unity_mod_minus_one:
                        alpha = z*v0*alpha0
                        # Now alpha ranges over all numbers with
                        # alpha*alphabar*n*xi2 = xi1 for which *potentially*
                        # A2/alpha is contained in A1 with the 
                        # appropriate index.
                        if ideal_contains([alpha*b for b in Z1.basis()], Z2.basis()):
                            # This is the point where we assume n is prime:
                            if transformation:
                                return True, alpha
                            return True
    if transformation:
        return False, None
    return False


def minimal_F(O):
    """
    Given an order O in a number field K, returns the smallest integer F such
    that F*O_K \subset O.
    """
    K = O.number_field()
    OK = K.maximal_order()
    MO = Matrix([K(g).list() for g in O.basis()])
    # ZZ^d . MO is a set of row vectors giving the elements of O in terms of
    # a basis of K
    MOK = Matrix([K(g).list() for g in OK.basis()])
    # ZZ^d . MOK is a set of row vectors giving the elements of OK in terms of
    # a basis of K
    M = MOK*MO.inverse()
    # ZZ^d . M is a set of row vectors giving the elements
    # of OK in terms of a basis of O
    F = M.denominator()
    assert all([F*K(b) in O for b in OK.basis()])
    return F
    

def is_trivial_in_shimura_group(A, alpha, O, cc=None):
    """
    Given an integral ideal A of the maximal order O_K of O coprime to 
    F = [O_K : O] and given an element alpha of O with A*Abar = alpha*O_K,
    returns True if and only if
    there exists mu in K^* with mu*mubar = alpha and mu*O coprime to F*O as
    ideals of O and mu*O_K = A.
    
    cc is complex conjugation (if not given, assumes K is a CM-field object).
    
    Assumes K is a non-biquadratic quartic CM-field (TODO: does not warn for
    biquadratic fields).
    
    Assumes O is stable under complex conjugation.
    
    EXAMPLES:
    
    We first create a suitable order, and standard things like reflex
    fields etcetera::
    
        sage: from recip import *
        sage: P.<x> = QQ[]
        sage: CM_Field(x^4+2*x^3+16*x^2+15*x+19).minimal_DAB()
        [149, 13, 5]
        sage: K = CM_Field([149,13,5])
        sage: Phi = K.CM_types()[0]
        sage: Psi = Phi.reflex()
        sage: Kr = Phi.reflex_field()
        sage: [alpha1,alpha2]=(x^4+2*x^3+16*x^2+15*x+19).roots(K, multiplicities=False)
        sage: cc = K.complex_conjugation()
        sage: cc(alpha1) == alpha2
        True
        sage: O = K.order(alpha1)
        sage: alpha2 in O
        True
        sage: OK = K.maximal_order()
        sage: O.index_in(OK)
        7
    
    Now let's compute lots of ideal classes to try::
    
        sage: ideal_classes = [c.ideal() for c in Kr.class_group()]
        sage: idealstar = Kr.ideal(7).idealstar(flag=2).list()
        sage: ray_classes = [a*b for a in idealstar for b in ideal_classes]
        sage: len(ray_classes)
        2304
        sage: type_norms = [(Psi.type_norm(I), I.norm()) for I in ray_classes] # long time, 20 seconds
    
    Here we test the current function::
    
        sage: all([is_trivial_in_shimura_group(A, alpha, OK) for (A, alpha) in type_norms]) # long time, 30 seconds
          ***   Warning: precision too low for generators, not given.
          ...
          ***   Warning: precision too low for generators, not given.
        True
        sage: len([A for (A, alpha) in type_norms if is_trivial_in_shimura_group(A, alpha, O)]) # long time, 25 seconds
        384
        sage: 2304/384
        6
    
    1/6 of ideals maps to a principal ideal, so the image has order 6.

    """
    K = O.number_field()
    if K.degree() != 4:
        raise NotImplementedError
    OK = K.maximal_order()
    if gcd(O.index_in(OK), A.norm()) != 1:
        raise ValueError, "A is not coprime to F"
    if cc is None:
        cc = K.complex_conjugation()
    if not all([cc(b) in O for b in O.basis()]):
        raise ValueError, "Order is not stable under complex conjugation"
    if not A*cc(A) == K.ideal(alpha):
        raise ValueError
    if not cc(alpha) == alpha:
        raise ValueError, "alpha not totally positive"
    if not all([b>0 for b in alpha.minpoly(polygen(QQ)).roots(AA, multiplicities=False)]):
        raise ValueError, "alpha not totally positive"

    if not A.is_principal():
        return False
    mu0 = A.gens_reduced()[0]
    # If mu exists, then mu = mu0 times a unit u, so
    # alpha = mu*mubar = mu0*mu0bar*u*ubar 
    # assuming quartic non-biquadratic CM-field, we have
    # u = u0 * z for a root of unity z
    # and a totally real u0, so
    # alpha / ( mu0*mu0bar ) = u0^2
    u0sq = alpha / mu0 / cc(mu0)
    if not K.ideal(u0sq) == 1:
        return False
    if not u0sq.is_square():
        return False
    u0 = u0sq.sqrt()
    if len(K.roots_of_unity())>2:
        raise NotImplementedError
    else:
        # z = +/- 1, and is only relevant up to sign, so u=u0*z=u0
        u = u0
    mu = u*mu0
    # Now we know mu satisfies (and is unique with) mu*O_K = A and
    # mu*mubar = alpha. So we only need to check mu*O is coprime to F*O.
    if not mu in O:
        return False
    # Now we know mu*mubar = alpha is coprime to F and mubar is in O,
    # hence mu is coprime to F.
    return True
    

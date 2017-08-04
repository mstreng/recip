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

This file implements CM period matrices and evaluations of theta constants in
them. It contains Shimura's reciprocity law!

"""

from sage.rings.complex_field import ComplexField
from sage.misc.misc import get_verbose
from sage.matrix.constructor import (Matrix, identity_matrix)
from sage.misc.cachefunc import cached_method
from sage.rings.all import ZZ
from sage.structure.sequence import Sequence
from sage.modules.free_module_element import vector
from sage.functions.other import (imag, floor)
from sage.rings.number_field.number_field_ideal import \
                                                is_NumberFieldFractionalIdeal
from sage.matrix.matrix_generic_dense import Matrix_generic_dense



def evaluate_theta(c, z, u = None, use_magma=False):
    """
    Numerically evaluate `\theta[c](z,u) =
    \sum_{n \in \ZZ^2}\exp(\pi i (n+c')^t z (n+c') + 2\pi i (n+c')^t (z+c''))`,
    where `c = (c',c'')`.
    If u is unspecified, use u = (0,...,0).
    
    TODO:
    
        Uses bounds from my thesis, but according to errata, these should be
        changed slightly for actual correctness.
        
    NOTE:
    
        algorithm='magma' requires Sage 5.6 or higher
    
    EXAMPLES::
    
        sage: from recip import *
        sage: I = CC.gen()
        sage: evaluate_theta([0,0,0,0], Matrix(CC, [[I+1, 1/2+I/3], [1/2+I/3, 3/2*I+1/5]]), use_magma=True) # optional - magma
        0.933295835691982 + 0.0143249911726235*I
    """
    if use_magma:
        zmag = magma(z)
        g = z.nrows()
        if u == None:
            u = [0 for i in range(g)]
        umag = magma(Matrix([u]).transpose())
        cmag = magma(Matrix([c]).transpose())
        return cmag.Theta(umag, zmag).sage()
    # Multiply the required precision by the following number, just to be sure.
    extra_prec = 1.1
    # Multiply the required bound on n by the following number, just to be sure.
    extra_terms = 1.1
    prec = z.base_ring().precision()
    if get_verbose() == 2:
        print "Precision %s" % prec
    R = extra_terms*ceil((0.4*prec+2.2).sqrt())
    if get_verbose() == 2:
        print "%s terms up to %s" % ((2*R+1)**2, R)
    if len(c) != 4:
        raise NotImplementedError, "sorry, evaluate_theta is only " \
                                   "implemented for g=2"
    F_extra = ComplexField(prec*extra_prec)
    I = F_extra.gen()
    cpp = vector([c[2],c[3]])
    piI = F_extra(pi*I)
    s = F_extra(0)
    for n in range(-R, R+1):
        for m in range(-R, R+1):
            cp = vector([ZZ(m)+c[0],ZZ(n)+c[1]])
            s = s + exp(piI*(cp*z*cp + 2*cp*cpp))
    return ComplexField(prec)(s)


def evaluate_theta_interval(c, z, R=None, reduce_first=True):
    """
    Numerically evaluate `\theta[c](z) =
    \sum_{n \in \ZZ^2}\exp(\pi i (n+c')^t z (n+c') + 2\pi i (n+c')^t (z+c''))`,
    where `c = (c',c'')`.
    
    `R` specifies the square over which to sum for n. It is automatically
    determined if not specified.

    TODO: Make a more sensible choice for `R` in case g=2:
    one that is different for
    the two coordinates m and n in the code below, and which depends on Z.
    
    TODO: It turns out that sometimes `R` is the limiting factor, and
    increasing precision by 100 bits has no effect as `R` is not raised.
    Check that `R` is large enough, especially in the non-interval setting!

    EXAMPLES::
    
        sage: from recip import *
        sage: I = ComplexIntervalField(100).gen()
        sage: Z = Matrix([[I+1, 1/2+I/3], [1/2+I/3, 3/2*I+1/5]])
        sage: c = [0,0,0,0]
        sage: evaluate_theta_interval(c, Z, 1)
        0.?e1 + 0.?e1*I
        sage: evaluate_theta_interval(c, Z, 2)
        1.0? + 0.1?*I
        sage: evaluate_theta_interval(c, Z, 3)
        0.93330? + 0.0144?*I
        sage: evaluate_theta_interval(c, Z, 4)
        0.9332958357? + 0.014324992?*I
        sage: evaluate_theta_interval(c, Z, 5)
        0.933295835691983? + 0.0143249911726235?*I
        sage: evaluate_theta_interval(c, Z, 6)
        0.93329583569198215423213? + 0.01432499117262351644467?*I
        sage: evaluate_theta_interval(c, Z, 7)
        0.9332958356919821542321268818? + 0.01432499117262351644467161792?*I
        sage: v = [evaluate_theta_interval(c, Z, k) for k in range(1, 8)+[20]]
        sage: vru = [x.real().upper() for x in v]
        sage: all([vru[k] >= vru[k+1] for k in range(6)])
        True
        sage: vrl = [x.real().lower() for x in v]
        sage: all([vrl[k] <= vrl[k+1] for k in range(6)])
        True
        sage: [RR(x.diameter()) for x in v]
        [4.38675115849014, 1.97112678302085, 0.000734124840204993, 1.23050838590045e-8, 8.92050564190827e-15, 2.79460042829428e-22, 1.53160261406578e-28, 9.18961568439468e-28]

    TESTS::
    
        sage: cs = [(1/2, 0), (0, 1/2), (0,0), (5/3, 111/7), (-3/4, 8), (1/100, -3/2), (-1,-1)]
        sage: I = sqrt(-1)
        sage: zs = [I/20, I*0.9, 2*I+8.3, I+0.1, I/3+1/7, I/5-10/3]
        sage: pairs = [(c,z) for c in cs for z in zs]
        sage: CIF_values = [evaluate_theta_interval(c, Matrix([[CIF(z)]]), reduce_first=True) for c in cs for z in zs]
        sage: CIF_values_non_reduced_first = [evaluate_theta_interval(c, Matrix([[CIF(z)]]), reduce_first=False) for c in cs for z in zs]
        sage: magma_values = [evaluate_theta(c, Matrix([[CC(z)]]), use_magma=True) for c in cs for z in zs] # optional - magma
        sage: quotients = [CIF_values[i]/CIF_values_non_reduced_first[i] for i in range(len(CIF_values))]
        sage: magma_quotients = [CIF_values_non_reduced_first[i].center()/magma_values[i] for i in range(len(CIF_values))] # optional - magma
        sage: [i for i in range(len(pairs)) if (magma_quotients[i]-1).abs() > 10^-8] # optional - magma
        []
        sage: [i for i in range(len(pairs)) if (quotients[i]-1).abs().lower() != 0]
        []
        sage: [i for i in range(len(pairs)) if (quotients[i]-1).abs().upper() > 10^-7]
        []

    """
    C = z.base_ring()
    I = C.gen()
    piI = C(pi*I)
    s = C(0)

    prec = C.precision()
    if get_verbose() == 2:
        print "Precision %s" % prec
    
    if len(c) == 2:
    
        z = z[0,0]
        if z.imag().lower() < 0:
            raise ValueError, "z must be in the upper half plane, but is %s" % z
        prefactor = 1
        if reduce_first:
            zeta8 = C((I+1)/sqrt(2))
            reduced = False
            z_for_reduction = z.center()
            while not reduced:
                t = (z_for_reduction.real()+1/2).floor()
#                print evaluate_theta_interval(c, Matrix([[z]]), reduce_first=False)
                if t != 0:
                    z_new = z - t
                    N = Matrix([[1, t],[0,1]])
                    # f(z) = f(N*z_new) = f|N(z_new)*f(N*u)/f|N(u), u->ioo
                    # So what is that limit? Write q = exp(2*pi*i*u),
                    # f(u) = lc * q^{e} + h.o.t.
                    # f(N*u) = f(u+t) = lc * exp(2*pi*i*e*t) * q^{e} + h.o.t.
                    # f|N(u) = lcN * q^{eN} + h.o.t.
                    c_trans, e = theta_action_without_kappa(N, 1, c)
                    lc, e = theta_leading_term_complex(c, C)
                    lcN, eN = theta_leading_term_complex(c_trans, C)
                    factor_for_this_transformation =  lc * exp(2*piI*e*t) / lcN
#                    print "factor", factor_for_this_transformation
                    prefactor = prefactor * factor_for_this_transformation
                    z = z_new
                    z_for_reduction = z_for_reduction - t
#                    print "T", z_new, c_trans
                    c = c_trans
#                    print "T"
#                    print prefactor*evaluate_theta_interval(c, Matrix([[z]]), reduce_first=False)
                
                if z_for_reduction.abs() < 0.9:
                    z_new = -1/z
                    N = Matrix([[0, -1], [1, 0]])
                    # f(z) = f(N*z_new) = sqrt(z_new)*f|N(z_new)*f(N*u)/f|N(u)/sqrt(u), u=I
                    #      = sqrt(z_new)*f|N(z_new) * f(I)/f|N(I)/((1+I)/sqrt(2))
                    #      = sqrt(z_new)*f|N(z_new) * f(I) * sqrt(2) / (1+I) / f|N(I)
                    sqrt_z_new = z_new.sqrt()
                    c_trans, e = theta_action_without_kappa(N, 1, c)
                    if sqrt_z_new.real() < 0:
                        sqrt_z_new = -sqrt_z_new
                    assert sqrt_z_new.imag() > 0
#                    print c, c_trans, e
#                    print prefactor
                    factor_for_this_transformation =  sqrt_z_new * evaluate_theta_interval(c, Matrix([[I]])) / zeta8 / evaluate_theta_interval(c_trans, Matrix([[I]]))
#                    print factor_for_this_transformation
                    prefactor = prefactor * factor_for_this_transformation
#                    print prefactor
                    z = z_new
                    z_for_reduction = -1/z_for_reduction
                    c = c_trans
#                    print "S", z_new, c_trans
#                    print "S"
#                    print prefactor*evaluate_theta_interval(c, Matrix([[z]]), reduce_first=False)
                else:
                    reduced=True
        
        c1 = c[0]
        c2 = c[1]
        RF = RealIntervalField(prec)
        Q = exp(RF(-pi*z.imag()))
        if not Q.upper() < 1:
            print Q
            raise ValueError
        if R is None:
            R = ZZ(ceil(sqrt((prec+3)*log(2)/-log(Q.upper())).n()))
        for n in srange(-R, R+1):
            cp = n+c1
            s = s + exp(piI * (cp**2 * z + 2*cp*c2))
        # and now to add the error term:
        if not c1 >= 0 and c1 < 1 and c2 >= 0 and c2 < 1:
            raise ValueError
        # sum of absolute value of remaining terms:
        # cp is at most -R-1+c1 or at least R+1+c1,
        # |cp| is at least R-c1
        # cp^2 is at least (R-c1)^2
        # so sum is at most: 2 * sum_{n>R}(|exp(pi*I * (n-c1)**2 * z)|)
        # Recall Q = |exp(pi*I*z)|
        # sum <= 2 * sum(Q^{(n-c1)^2})
        #      = 2 * Q^{(R+1-c1)^2} * sum(Q^{(n-c1)^2 - (R+1-c1)^2})
        #      = 2 * Q^{R^2} * sum(Q^{[(n-c1) - (R+1-c1)] * [(n-c1) + (R+1-c1)]})
        #      = 2 * Q^{R^2} * sum(Q^{[(n - (R+1))] * [(n-c1) + (R+1-c1)]})
        #     <= 2 * Q^{R^2} * sum(Q^{2*R*k})     ( k = n-R-1 >= 0 )
        #      = 2 * Q^{R^2} / (1 - Q^{2*R}) 
        error_term_bound = (2*Q**(R**2) / (1-Q**(2*R))).upper()
        error_term_bound = RF(-error_term_bound, error_term_bound)
        s = s + C(error_term_bound, error_term_bound)
#        print s
        return prefactor * s
        

    if len(c) != 4:
        raise NotImplementedError, "sorry, evaluate_theta_interval is only " \
                                   "implemented for g<=2"
    if R is None:
        R = ceil((0.4*prec+2.2).sqrt())

    if get_verbose() == 2:
        print "%s terms up to %s" % ((2*R+1)**2, R)
        
    cpp = vector([c[2],c[3]])
    verb = get_recip_verbose()
    if verb:
        print "Evaluating a theta constant in interval arithmetic, R=%s" % R
    for n in range(-R, R+1):
        if verb and R > 40:
            print "Evaluating a theta constant, R=%s, n=%s, so %s percent" % (R, n, 100*((n+R)/2/R).n())
        for m in range(-R, R+1):
            cp = vector([ZZ(m)+c[0],ZZ(n)+c[1]])
            s = s + exp(piI*(cp*z*cp + 2*cp*cpp))

    # and now to add the error term:
    RF = RealIntervalField(prec)
    error_term_bound = (RF(7.247)*sum([exp(RF(-1/2*pi*R**2)*z[k,k].imag()) for k in [0,1]])).upper()
    error_term_bound = RF(-error_term_bound, error_term_bound)
    s = s + C(error_term_bound, error_term_bound)

    return s
    
    
def _riemann_form(basis, xi, bar):
    """
    Returns a matrix for the Riemann form (x,y) --> trace(xi*bar(x)*y)
    with respect to the given basis
    """
    return Matrix(ZZ, [[(bar(x)*xi*y).trace() for y in basis] for x in basis])


def _symplectic_basis(A, xi, bar, double_check=False):
    """
    Returns a symplectic basis of A for the Riemann form
    (x,y) --> trace(xi*bar(x)*y)
    """
    if type(A) is list:
        bas = A
    else:
        bas = A.basis()
    E = _riemann_form(bas, xi, bar)
    T = E.symplectic_form()[1]
    ret = Sequence(T*vector(bas))
    if double_check:
        g = len(bas)/2
        E = _riemann_form(ret, xi, bar)
        assert E == _matrix_omega(g)

    return ret






def _big_period_matrix(Phi, bas):
    """
    Returns the matrix with columns Phi(b) for b in bas 
    """
    return Matrix([[phi(b) for b in bas] for phi in Phi])


def _small_period_matrix(Phi, bas):
    """
    Returns a small period matrix Omega_2^-1 Omega_1 from a given
    period matrix (Omega_1; Omega_2)
    """
    bigone = _big_period_matrix(Phi, bas)
    g = bigone.nrows()
    assert bigone.ncols() == 2*g
    bigone.subdivide(g,g)
    Omega1 = bigone.subdivision(0,0)
    Omega2 = bigone.subdivision(0,1)
#    return Omega2**-1 * Omega1
    return Omega2.adjoint() * Omega1 / Omega2.det()


def my_ceil(a):
    ret = ceil(a)
    if ret.parent() is SR:
        # happens for a an integer in QQbar.
        return ZZ(a)
    return ret


def my_floor(a):
    ret = floor(a)
    if ret.parent() is SR:
        # happens for a an integer in QQbar.
        return ZZ(a)
    return ret


def _realred(Z):
    """
    Given a real matrix Z, returns (Z', T), where Z' = Z + T
    has coefficients in the interval [-1/2, 1/2) and T
    has integer coefficients.
    """
    V = Matrix([[my_ceil(t.real()-1/2) for t in u] for u in Z])
    return Z - V, -V


def _reduce_symm_matrix(Y):
    #assert Y.determinant() > 0
    #assert Y[0,1] == Y[1,0]
    #assert Y[0,0] > 0
    T = Matrix([[0,1],[-1,0]])
    finished = False
    U = Matrix([[1,0],[0,1]])
    while not finished:
        r = my_floor(-Y[0,1]/Y[0,0]+1/2)
        #R = Matrix([[1,0],[r,1]], ring = QQbar)
        R = Matrix([[1,0],[r,1]])
        U = R*U
        #Y = QQbar_to_AA(R*Y*R.transpose())
        Y = R*Y*R.transpose()
        if Y[0,0] > Y[1,1]:
            #T = Matrix([[0,1],[-1,0]], ring = QQbar)
            U = T*U
            #Y = QQbar_to_AA(T*Y*T.transpose())
            Y = T*Y*T.transpose()
        else:
            finished = True
    if Y[0,1] < 0:
        T = Matrix([[1,0],[0,-1]])
        U = T*U
        Y = T*Y*T.transpose()
    return Y, mat_convert(U, ZZ)


def _imagred(Z):
    Y, U = _reduce_symm_matrix(mat_convert(Z, imag))
    return U*Z*U.transpose(), U

    
# The following variable is the cache for the function gottschling_matrices.
_gottschling = None


def gottschling_matrices():
    """
    Returns a set of matrices in Sp(4,ZZ) containing
    gottschling's set.
    """
    global _gottschling
    if _gottschling == None:
        from sage.misc.mrange import cartesian_product_iterator
        _gottschling = [Matrix([[0,0,-1,0],[0,1,0,0],[1,0,e,0],[0,0,0,1]]) for e in [-1,0,1]] \
            + [Matrix([[1,0,0,0] ,[0,0,0,-1],[0,0,1,0],[0,1,0,e]]) for e in [-1,0,1]] \
            + [Matrix([[0,0,-1,0],[0,1,0,0],[1,-1,d,0],[0,0,1,1]]) for d in [-2,-1,0,1,2]] \
            + [Matrix([[0,0,-1,0],[0,0,0,-1],[1,0,e1,e3],[0,1,e3,e2]]) \
                 for [e1,e2,e3] in cartesian_product_iterator([[-1,0,1] for k in range(3)])]
    return _gottschling


def _matrix_omega(g):
    """

    EXAMPLE::

        sage: from recip import _matrix_omega
        sage: _matrix_omega(2) == Matrix([[0,0,1,0],[0,0,0,1],[-1,0,0,0],[0,-1,0,0]])
        True
    """
    g = ZZ(g)
    return matrix_from_blocks(zero_matrix(g), identity_matrix(g), -identity_matrix(g), zero_matrix(g))


def is_symplectic(gamma, g=None):
    """
    Tests whether a 4x4 matrix is an element of Sp_g(ZZ).
    
    EXAMPLES::

        sage: from recip import *
        sage: g = gottschling_matrices()
        sage: all([is_symplectic(gamma) for gamma in g]) and len(g) == 38
        True
        sage: is_symplectic(g[0]+1)
        False
    """
    if g is None:
        g = ZZ(gamma.ncols()/2)
    om = _matrix_omega(g)
    return gamma.transpose() * om * gamma == om


def _det_bottom_part(gamma, tau):
    gamma.subdivide(2,2)
    c = gamma.subdivision(1,0)
    d = gamma.subdivision(1,1)
    return (c*tau+d).determinant()

    
def Sp_action(gamma, Z):
    """
    Given a matrix gamma in GSp_2g(QQ)^+ and a gxg period matrix Z,
    returns gamma*Z as an immutable gxg Matrix.
    """ 
    gamma.subdivide(2,2)
    a = gamma.subdivision(0,0)
    b = gamma.subdivision(0,1)
    c = gamma.subdivision(1,0)
    d = gamma.subdivision(1,1)
    ret = (a*Z+b)*(c*Z+d)**-1
    ret.set_immutable()
    return ret

    
def _gottschling_reduce(Z):
    """
    Given Z, returns a pair (Z', g) with g an element of the set
    gottschling_matrices() and Z' = g*Z, such that det Im(g*Z) > det Im(Z)
    holds, if such a g exists. If such a g does not exist, return (Z, I) with I
    a 4x4 identity matrix.
    """
    improvement = 1
    for g in gottschling_matrices():
        d = abs(_det_bottom_part(g, Z))
        if d < improvement:
            if get_verbose() == 2:
                print "The matrix %s would improve the determinant of the " \
                      "imaginary part by a factor %s" % (g, 1/d)
            gamma = g
            improvement = d
    if improvement == 1:
        if get_verbose() == 2:
            print "Cannot increase the imaginary part of %s" % Z
        return (Z, identity_matrix(4))
    Z = Sp_action(gamma, Z)
    if get_verbose() == 2:
        print "Improving the imaginary part by a factor %s using %s to get " \
              "%s" % (1/improvement, gamma, Z)
    return (Z, gamma)


def _reduce(Z, reduction_sequence=False):
    """
    Given Z, returns (Z', M) such that Z' is reduced, M is in Sp(2g,ZZ) and
    MZ = Z'
    
    Only implemented for g=2 and maybe g=1.
    
    If reduction_sequence is True, then instead of M returns
    [M1,...,Mk] such that M = Mk*...*M1 and Mi is in a certain list
    of generators of Sp(2g,ZZ).
    """
    g = Z.nrows()
    if get_verbose() == 2:
        print "Beginning reduction of " + str(Z)

    if g > 2:
        raise NotImplementedError

    if reduction_sequence:
        M = []
    else:
        M = identity_matrix(2*g)
    while True:
        Z, U = _imagred(Z)
        if get_verbose() == 2:
            print "The imaginary part is made reduced by %s. This yields %s" % \
                  (U,Z)
        Uti = U.transpose().inverse()
        U = Matrix([[U[0,0],U[0,1],0,0],[U[1,0],U[1,1],0,0],[0,0,Uti[0,0],
                     Uti[0,1]],[0,0,Uti[1,0],Uti[1,1]]])
        assert is_symplectic(U)
        if reduction_sequence:
            M.append(U)
        else:
            M = U*M
            assert is_symplectic(M)
        Z, V = _realred(Z)
        # Now it may happen that V is not symmetric (because of rounding),
        # so we use V[0,1] twice and ignore V[1,0].
        if get_verbose() == 2:
            print "The real part is made reduced by %s. This yields %s" % (V, Z)
        V = Matrix([[1,0,V[0,0],V[0,1]],[0,1,V[0,1],V[1,1]],[0,0,1,0],[0,0,0,1]])
        assert is_symplectic(V)
        if reduction_sequence:
            M.append(V)
        else:
            M = V*M
            assert is_symplectic(M)
        Z, gamma = _gottschling_reduce(Z)
        if gamma == 1:
            return Z, M
        if reduction_sequence:
            M.append(V)
        else:
            M = gamma*M
            assert is_symplectic(M)


def is_positive_definite(m):
    """
    Returns `True` if and only if symmetric real matrix `m`
    is positive definite.
    
    EXAMPLES::
    
        sage: from recip import *
        sage: is_positive_definite(Matrix(AA,[[2,2],[2,4]]))
        True
        sage: is_positive_definite(Matrix(AA,[[1,2],[2,4]]))
        False
        sage: is_positive_definite(Matrix(AA,[[-3,2],[2,-4]]))
        False
    """
    try:
        return QuadraticForm(m).is_positive_definite()
    except NotImplementedError:
        pass
    B = m.base_ring()
    if m.ncols() == 2 and (B == AA or B == RDF or is_RealField(B)):
        if m.nrows() != 2 or m[1,0] != m[0,1]:
            raise ValueError, "Matrix m (=%s) is not symmetric"
        return m[0,0] > 0 and m[0,0]*m[1,1] - m[1,0]**2 > 0
    raise NotImplementedError, "is_positive_definite not implemented for matrices of dimension %s over %s" % (m.ncols(), B)

    
def is_period_matrix(m):
    """
    Returns `True` if and only if m is a period matrix.
    
    EXAMPLES::
    
        sage: from recip import *
        sage: i = CC.gen()
        sage: is_period_matrix(Matrix([[1+i, 2+i], [2+i, 5*i]]))
        True
        sage: is_period_matrix(Matrix([[1+i, 2+i], [2+i, i]]))
        False
    """
    height = m.nrows()
    width = m.ncols()
    if 2*height == width:
        raise NotImplementedError
    elif height != width:
        raise ValueError
    if not m.is_symmetric():
        if get_verbose() == 2 != 0:
            print "Matrix %s is not symmetric in is_period_matrix" % m
        return False
    if not is_positive_definite(mat_convert(m, imag)):
        if get_verbose() == 2 != 0:
            print "Imaginary part of %s is not positive definite" % m
        return False
    return True
    

def PeriodMatrix(arg1, arg2=None, arg3=None, check=True):
    r"""
    Create a period matrix.
    
    Currently, period matrices with CM by
    non-maximal orders may not completely work.
    
    INPUT:
    
    - `arg1` -- a CM-type `\Phi` of a CM-field `K`
    - `arg2` -- an ideal `A` of a maximal order of `K`
    - `arg3` -- an element `\xi` of `K` with totally positive imaginary part
      and real part 0, such that `\xi A` is the trace dual of `\bar{A}`.
      (If `A` is an ideal of the maximal order `O_K` and `D` is the different,
      then the final condition is equivalent to `\xi A\bar{A} D = O_K`.)
    - `check` (default:True) -- whether to check the condition on `arg3`.
    
    Instead of an ideal of a maximal order, it is possible to give a basis
    of any lattice in `K`. This basis does not have to be symplectic. For this,
    call directly the class (it is not implemented in this constructor yet).
    
    OUTPUT:
    
    A period matrix corresponding to some symplectic basis for the input data.
    """
    if arg2 == None and arg3 == None:
        if get_verbose() == 2:
            print "PeriodMatrix got a single input (%s), interpreting it as a PeriodMatrix_CM and creating a new PeriodMatrix_CM (copy) out of it" % arg1
        return PeriodMatrix_CM(arg1.CM_type(), arg1.ideal(), arg1.xi(), arg1.basis(), arg1, check=check)
    if get_verbose() == 2:
        print "PeriodMatrix got inputs CM_type = %s, ideal = %s, xi = %s" % (arg1, arg2, arg3)
    return PeriodMatrix_CM(arg1, arg2, arg3, check=check)


class PeriodMatrix_CM(Matrix_generic_dense):
    r"""
    Object representing a CM period matrix. See :func:`PeriodMatrix`.
    """
    def __init__(self, CM_type=None, ideal=None, xi=None, basis=None,
                 matrix=None, check=True):
        if xi is None:
            raise NotImplementedError, "xi must be supplied"
        
        if CM_type is None:
            CM_type = CM_Type(xi, check=False)
        
        g = CM_type.g()
        bar = CM_type.domain().complex_conjugation()
        if basis is None:
            if ideal is None:
                raise ValueError, "Either basis or ideal must be supplied"
            basis = _symplectic_basis(ideal, xi, bar)
        no_ideal_supplied = (ideal is None or type(ideal) is list)
        if no_ideal_supplied:
            ideal = CM_type.domain().ideal(basis)
            if not xi*bar(ideal)*ideal*CM_type.domain().different() == 1:
                ideal = None
        if check:
            if not (ideal is None or
                    xi*bar(ideal)*ideal*CM_type.domain().different() == 1):
                raise ValueError(
                      "We don't have xi*bar(ideal)*ideal*different = O_K")
            for i in range(len(basis)):
                if not (ideal is None or basis[i] in ideal):
                    raise ValueError, "(%s)th basis element %s not in " \
                                      "ideal %s" % (i, basis[i], ideal)
                for j in range(i):
                    b = (xi*bar(basis[i])*basis[j]).trace()
                    if j+g == i:
                        if b != -1:
                            raise ValueError, "Basis %s is not symplectic: " \
                                              "-1 expected in position " \
                                              "(%s, %s), but %s found. Matrix:\n %s" % \
                                              (basis,i,j,b,Matrix([[(xi*bar(basis[i])*basis[j]).trace() for j in range(len(basis))] for i in range(len(basis))]))
                    elif b != 0:
                        raise ValueError, "Basis %s is not symplectic: " \
                                          "0 expected in position (%s, %s), " \
                                          "but %s found.\n Matrix: %s" % (basis,i,j,b,Matrix([[(xi*bar(basis[i])*basis[j]).trace() for j in range(len(basis))] for i in range(len(basis))]))
            
                    
        if matrix == None:
            matrix = _small_period_matrix(CM_type, basis)
        elif check:
            if Sequence(matrix) != \
               Sequence(_small_period_matrix(CM_type, basis)):
                raise ValueError, "Period matrix belonging to basis %s is " \
                         "%s, but %s was supplied" % \
                         (basis, _small_period_matrix(CM_type, basis), matrix)
        from sage.matrix.matrix_space import MatrixSpace
        M = MatrixSpace(CM_type.codomain(), nrows=g, ncols=g, sparse=False)
        from sage.misc.flatten import flatten
        entries = flatten([[a for a in b] for b in matrix])
        Matrix_generic_dense.__init__(self, M, entries=entries, copy=True,
                                      coerce=True)
        self._CM_type = CM_type
        self._ideal = ideal
        self._xi = xi
        self._basis = basis
        self.set_immutable()
        # rows of Bt are bas[i] expressed in terms of standard basis
        self._Bt = Matrix([b.vector() for b in basis])

    def g(self):
        return self.CM_type().g()

    def __repr__(self):
        return "Period Matrix\n" + Matrix_generic_dense.__repr__(self)

    def complex_matrix(self, prec=None):
        """
        Returns self as a matrix over a complex field.
                
        INPUT:
        
          - prec (default=None) -- Either None or a positive integer.
          
        OUTPUT:
        
          - If prec is an integer, returns self as a matrix over ComplexField(prec)
          
          - If prec=None, returns self as a matrix over self.base_ring().embedding().codomain()
        """
        if prec == None:
            K = self.base_ring()
            if _is_accepted_complex_field(K):
                return self
            emb = K.embedding()
            if emb == None:
                raise RuntimeError, "Incorrect period matrix: field not embedded"
            return mat_convert(self, emb)
        if prec in ZZ:
            prec = ComplexField(prec)
        return mat_convert(self.complex_matrix(None), prec)

    def complex_conjugate(self, transformation=False):
        """
        Returns minus the complex conjugate of self.
        
        If transformation is True, then also returns
        a matrix in Sp_2g(QQ) (if possible in Sp_2g(ZZ)) that maps
        self to the output.
        
        EXAMPLES::
        
            sage: from recip import *
            sage: U = CM_Field(x^4+5*x^2+5).period_matrices_iter().next(); U
            Period Matrix
            [ 0.30901699437...? + 0.95105651629...?*I  0.50000000000...? + 0.36327126400268?*I]
            [ 0.50000000000...? + 0.36327126400268?*I -0.30901699437495? + 0.95105651629...?*I]
            sage: Ubar, M = U.complex_conjugate(transformation=True); Ubar
            Period Matrix
            [-0.30901699437...? + 0.95105651629...?*I -0.50000000000...? + 0.36327126400268?*I]
            [-0.50000000000...? + 0.36327126400268?*I  0.30901699437495? + 0.95105651629...?*I]
            sage: M.base_ring()
            Integer Ring
            sage: M
            [-1  0  0  0]
            [ 0 -1  0  0]
            [ 0  1 -1  0]
            [ 1 -1  0 -1]
            sage: U.complex_matrix(20)
            [ 0.30902 + 0.95106*I  0.50000 + 0.36327*I]
            [ 0.50000 + 0.36327*I -0.30902 + 0.95106*I]
            sage: Ubar.complex_matrix(20)
            [-0.30902 + 0.95106*I -0.50000 + 0.36327*I]
            [-0.50000 + 0.36327*I  0.30902 + 0.95106*I]
            sage: U.Sp_action(M) == Ubar
            True
        """
        # Minus the complex conjugate of self has:
        # Phi     : Phi
        # ideal   : idealbar
        # basis   : basisbar with the first g multiplied by -1
        # xi      : xi
        # Z       : -Zbar
        # So we calculate these quantities where they differ.
        Phi = self.CM_type()
        ideal = self.ideal()
        conjugation_k = ideal.number_field().complex_conjugation()
        idealbar = conjugation_k(ideal)
        basisbar = [conjugation_k(b) for b in self.basis()]
        for i in range(self.g()):
            basisbar[i] = -basisbar[i]
        xi = self.xi()
        # conjugation_l = self.base_ring().complex_conjugation()
        # Z = - mat_convert(self, conjugation_l)
        
        # If possible, then we replace idealbar by ideal so that we
        # can identify ideal with idealbar.
        # In other words, if ideal/idealbar is principal and generated by
        # x, then we replace idealbar by ideal, basisbar by x*basisbar,
        # and xi by xi/(x*xbar).
        if (ideal/idealbar).is_principal():
            gns = (ideal/idealbar).gens_reduced()
            if len(gns)>1:
                raise RuntimeError
            x = gns[0]
            idealbar = ideal
            basisbar = [x*b for b in basisbar]
            xi = xi / (x*conjugation_k(x))            
        
        Z = PeriodMatrix_CM(Phi, idealbar, xi, basisbar) #, Z)
        if not transformation:
            return Z
        M = (Z._Bt * self._Bt.inverse())    
        if not M*vector(self.basis()) == vector(Z.basis()):
            raise RuntimeError, "bug in complex_conjugate, wrong matrix"
        #if not Sp_action(M, self) == Z:
        #    raise RuntimeError, "bug in complex_conjugate, wrong action"
        try:
            M = mat_convert(M, ZZ)
        except TypeError:
            pass

        return Z, M

    def CM_field(self):
        return self.CM_type().domain()
    
    def reflex_field(self):
        return self.CM_type().reflex_field()
        
    @cached_method
    def complex_conjugation_symplectic_matrix(self, level, mu=None, A=None):
        r"""
        Returns the matrix U of Proposition 2.7.
        
        INPUT:
        
        - ``level`` -- positive integer
        - ``mu`` -- element of the CM-field such that mu*mubar is in QQ and
          such that an ideal A exists as below. If not specified, then in
          some cases one will be determined automatically.
        - ``A`` -- invertible fractional ideal of the reflex field satisfying
          N_{reflextype, O} * Bbar = mu*B (where B is the ideal of self).
          It must exist. If it is not provided, then it is assumed that it
          exists and is coprime to the level. If provided, then it is scaled
          so that it is coprime to the level, and the input is checked.
        
        The cases where A and mu are determined automatically are as follows:
        
        - B is stable under complex conjugation: A=1, mu=1
        - g=1 and B is an ideal of the maximal order: start with
          A=B/Bbar, mu=1, then scale if A is not coprime to the level
        - g=2 and B is an ideal of the maximal order: start with
          A=N_Phi(B), mu = N(B), then scale if A is not coprime to the level

        TODO: verify doctest also by hand.
        
        OUTPUT:
        
        Matrix U in GSp(ZZ/level*ZZ) such that, if f(self) is a class
        invariant, then f^U(self) = f(self) if and only if f(self) is in M_0.
        
        EXAMPLES::
        
            sage: from recip import *
            sage: k = CM_Field((x^2+5)^2-4*5)
            sage: it = k.period_matrices_iter()
            sage: Z = it.next(); Z
            Period Matrix
            [-0.30901699437495? + 0.95105651629515?*I -0.50000000000000? + 0.36327126400268?*I]
            [-0.50000000000000? + 0.36327126400268?*I  0.30901699437495? + 0.95105651629515?*I]
            sage: Z.complex_conjugation_symplectic_matrix(8)
            [1 0 0 0]
            [0 1 0 0]
            [0 7 7 0]
            [7 1 0 7]

        """
        B = self.ideal()
        rho = self.CM_field().complex_conjugation()
        if mu is None:
            if B == rho(B):
                A = self.reflex_field().ideal(1)
                mu = self.CM_field()(1)
            elif self.g() == 1:
                A = self.CM_type().type_norm(B/rho(B))
                mu = self.CM_field()(1)
            elif self.g() == 2:
                A = self.CM_type().type_norm(B)
                mu = self.CM_field()(B.norm())
            else:
                raise NotImplementedError, "Finding A and mu is only " \
                                           "implemented for g<=2 and for " \
                                           "ideals that are invariant under " \
                                           "complex conjugation."
        elif not mu in self.CM_field():
            raise ValueError, "mu not in CM-field of self"
        if not rho(mu)*mu in QQ:
            raise ValueError, "mu*mubar not in QQ for mu = %s" % mu
        if not A is None:
            if self.CM_type().reflex().type_norm(A)*rho(B) != mu*B:
                raise ValueError, "A (=%s) and mu (=%s) do not satisfy the " \
                                  "hypothesis"
            A, x = lift_ray_class_group_element(A, 1, level, generator=True)
            x = self.reflex_field()(x)
            mu = mu / self.CM_type().reflex().type_norm(x)
            if self.CM_type().reflex().type_norm(A)*rho(B) != mu*B:
                raise RuntimeError, "A (=%s) and mu (=%s) do not satisfy the "\
                                  "hypothesis"
        Ct = Matrix([(mu**-1*rho(b)).vector() for b in self.basis()])
        # C = B Mt, so Ct = M Bt, so M^-1 = Bt Ct^-1
        # U = M^-1 = Bt Ct^-1
        U = self._Bt * Ct.inverse()
        U = mat_convert(U, Zmod(level))
        return GSp_element(U)
            
        
        
#    @cached_method
    def reduce(self, prec=None, transformation=False):
        """
        Returns a reduced period matrix Z that is Sp_{2g}-equivalent to self.
        
        INPUT:
        
          - prec (default=None) -- the working precision when reducing.
            If None, use the embedding of self.base_ring() into the complex
            numbers.
            If a positive integer, use a complex field with that many bits of
            precision.
            
          - transformation (default=False) -- whether to also return a matrix M
            such that Z = M*self.
        """
        self_numerical = self.complex_matrix(prec)
        Z_numerical, M = _reduce(self_numerical)
        if not nu(M) == 1:
            raise RuntimeError, "M should be a symplectic matrix, but is %s over %s with nu=%s" % (M, M.base_ring(), nu(M))
        Z_basis = Sequence(M * vector(self.basis()))
        Z = PeriodMatrix_CM(self.CM_type(), self.ideal(), self.xi(), Z_basis) #, Sp_action(M, self))
        if not M*vector(self.basis()) == vector(Z.basis()):
            raise RuntimeError, "bug in reduce, wrong matrix"        
        if transformation:
            return (Z, M)
        return Z

    def CM_type(self):
        return self._CM_type

    def ideal(self):
        return self._ideal

    def xi(self):
        return self._xi

    def basis(self):
        return self._basis

    def evaluate_theta(self, c, prec, use_magma=False, interval=False):
        if get_verbose() == 2:
            print "Evaluating theta constant of characteristic %s at %s" % \
                  (c, self)
        if interval:
            if use_magma:
                raise NotImplementedError, "Cannot use both Magma and interval arithmetic"
            return evaluate_theta_interval(c,
                mat_convert(self, ComplexIntervalField(prec)))
        return evaluate_theta(c, self.complex_matrix(prec),
                use_magma=use_magma)

    @cached_method
    def _theta_vals(self, den, prec, use_magma=False, interval=False):
        g = self.g()
        return [self.evaluate_theta(num_to_c(i,g,den),prec,
                use_magma=use_magma, interval=interval) for i in range(den**(2*g))]

    def Sp_action(self, M):
        """
        Returns the period matrix M*self for any M in GSp_{2g}(QQ)^+.
        """
        try:
            M = M.matrix()
        except AttributeError:
            pass
        Z_basis = Sequence(M * vector(self.basis()))
        return PeriodMatrix_CM(self.CM_type(), self.ideal(), self.xi(), Z_basis)
        
    def negative_inverse(self):
        """
        Returns -1/self.
        """
        basis = self.basis()
        g = self.g()
        basis = basis[g:] + [-basis[i] for i in range(g)]
        return PeriodMatrix_CM(self.CM_type(), self.ideal(), self.xi(), basis)

    def galois_action(self, A, transformation=False, reduce=100):
        """
        Returns a period matrix Z with the same CM-type as self,
        and with Z.ideal() = N_Phi(A)^-1 * self.ideal()
        and Z.xi() = N(A)^-1 * self.xi().
        
        For a modular function f of level 1, we have
        f(self)^sigma(A) = f(Z), where sigma is the Artin map. 

        INPUT:
        
        - A -- an ideal of self.CM_type().reflex_field()
          
        - N (default=1) -- interpret x as an element of I(N)/R(N),
            where I(N) is the set of ideals coprime to N and R(N)
            is the set of ideals generated by elements that are 1 mod N.
            
        - transformation (default=False) -- whether to also
            return a matrix M in Mat(ZZ,2g) with M*Z.basis() = self.basis()
            (and M*Z = self)
            
        - reduce (default=100) -- if False, do not reduce the basis,
            if a positive integer, reduce the period matrix numerically
            with that precision.
            
        OUTPUT:
        
        With transformation=False, a matrix Z as above.
        With transformation=True, a pair (Z,M) as above.
        """
        Phi = self.CM_type()
        ideal = Phi.reflex().type_norm(A)**-1 * self.ideal()
        xi = A.norm() * self.xi()
        Z = PeriodMatrix(Phi, ideal, xi)
        M = (self._Bt * Z._Bt.inverse())
        if get_verbose()==2:
            print "galois action obtained by %s" % M
        if reduce != False:
            Z, M2 = Z.reduce(prec=reduce, transformation=True)
            if get_verbose()==2:
                print "reduction obtained by the inverse\n %s of\n %s" % (M2.inverse(), M2)
            M = M * M2.inverse()
        if not M*vector(Z.basis()) == vector(self.basis()):
            raise RuntimeError, "bug in galois_action, wrong matrix"
        if transformation:
            return (Z, M)
        return Z
        
    def epsilon(self, x):
        """
        The map epsilon of page 57 of Shimura's "on certain reciprocity laws..."
        Returns the transpose of the matrix of multiplication by x wrt the basis self.basis()
        """
        # columns of M are x*bas[i] expressed in terms of standard basis
        # M = Matrix([(x*b).vector() for b in bas]).transpose()
        Mt = Matrix([(x*b).vector() for b in self.basis()])
        # return (B.inverse()*M).transpose()
        return Mt*(self._Bt.inverse())

    def Shimura_reciprocity(self, A, n, m=None, reduce=100, period_matrix=False, transformation=True):
        """
        Returns matrices M in GSp(QQ)^+ and u in GSp(ZZ/n*ZZ) that
        give the Galois action of the ray class A mod m on a
        modular function f of level n.
        
        INPUT:
        
        - A -- an ideal of K = self.CM_type().reflex_field()
        
        - n -- a positive integer
        
        - m (default:None) -- a positive integer or None.
          If None, then use m=n.
          For every prime P of K dividing m, we must have
          ord_P(A)=0.
          
        - reduce (default:100) -- a positive integer or False.
          If a positive integer, then M*self is reduced
          (with precision reduce).
          
        - period_matrix (default:False) -- whether to
          also return Z = (M^-1)*self
          
        - transformation (default:True) -- whether to
          return M. period_matrix and transformation
          can not both be false.
                      
        OUTPUT:
        
        A pair (M, u) or (Z, u) or a triple (Z, M, u)
        depending on period_matrix and transformation.
        Here we have M in GSp(QQ)^+, u in GSp(ZZ/n*ZZ),
        and Z = M^-1*self a period matrix
        such that the following holds.
        
        Let f be a modular function of level n
        and assume that f(self) is in the ray class field
        of K of conductor m. (This assumption
        automatically holds in the case m=n.)
        
        Then f(self)^Artin(A) = f^u(Z).
        
        NOTE:
        
        If m=n, then M and Z are as in
        self.galois_action(x, transformation=True).
        
        EXAMPLES::
        
            sage: from recip import *
            sage: k = CM_Field([5,13,41])
            sage: Z = k.one_period_matrix()
            sage: a = Z.CM_type().reflex_field().gen()
            sage: i = igusa_invariants_absolute()[0]
            
        The following is supposed to be real, so the precision of 100 is not
        really 100. This should be corrected (TODO)::
        
            sage: i(Z, prec=100)
            6464.1219280862927825807751509 + 1.4177571165844968616946581863e-26*I
            sage: (U, M, u) = Z.Shimura_reciprocity(a.parent().ideal(a), m=1, n=8, period_matrix=True)
            sage: M
            [ 1 -2 -1 -1]
            [-2 -1 -1 -2]
            [ 0  0  1 -2]
            [ 0  0 -2 -1]
            sage: u
            [1 6 7 7]
            [6 7 7 6]
            [0 0 1 6]
            [0 0 6 7]
            sage: (i^u)(U, prec=100)
            6464.12192808629278258077515...

        The following sign error points out a bug, or not? a is not 1 mod 8,
        so it is allowed to have a non-trivial effect::

            sage: P = theta_ring(2, 2)[0]
            sage: t = ThetaModForm(P.gens()[4]/P.gens()[6])
            sage: t(Z, prec=100)
            0.99854006288205177918601876434 - 0.054016134807185886527908776...*I
            sage: t
            t4/t6
            
        Simplification of theta quotient expressions does not wcompletely work
        yet in all cases (TODO)::

            sage: t^u
            ((-zeta8)*t6)/((zeta8)*t4)
            sage: (t^u)(U, prec=100)
            -0.99854006288205177918601876434 + 0.054016134807185886527908776...*I
            sage: P = theta_ring(2,2)[0]

            
        """
        if not (n in ZZ and n > 0):
            raise TypeError, "n (=%s) must be a positive integer" % n
        if m == None:
            m = n
        elif not (m in ZZ and n > 0):
            raise TypeError, "m (=%s) must be a positive integer or None" % m
        
        A = lift_ray_class_group_element(A, m, n)
        (Z, M) = self.galois_action(A, transformation=True, reduce=reduce)
        
        
        u = GSp_element(mat_convert(M, Zmod(n)))
        
        if period_matrix and transformation:
            return (Z, M, u)
        if period_matrix:
            return (Z, u)
        if transformation:
            return (M, u)
        raise ValueError, "period_matrix and transformation are not " \
                          "allowed to be both False"
    
    def has_real_moduli(self):
        """
        Returns true if and only if self is isomorphic
        to its complex conjugate over CC.
        
        EXAMPLES::

            sage: from recip import *
            sage: K = CM_Field(x^4+x^3+x^2+x+1, embedding=QQbar) # the embedding is not added afterwards, and the reflex field is taken to be the field itself, so this is necessary for now
            sage: Z = K.one_period_matrix()
            sage: Z.has_real_moduli()
            True
            sage: K = CM_Field([521,27,52])
            sage: mat = list(K.period_matrices_iter())
            sage: ZZ(sum([Z.has_real_moduli() for Z in mat])) / ZZ(len(mat))
            1/7
        """
        M = (self.complex_conjugate(transformation=True)[1])
        return M.base_ring() == ZZ
    
    def __mul__(self, other):
        if not other in ZZ:
            raise ValueError
        basis = self._basis
        g = self.g()
        basis = ([basis[i] * other for i in range(g)] +
                 [basis[i] for i in range(g, 2*g)])
        return PeriodMatrix_CM(xi=self._xi/other, basis=basis,
                               CM_type=self._CM_type, check=True)
    
    def __div__(self, other):
        if not other in ZZ:
            raise ValueError
        basis = self._basis
        g = self.g()
        basis = ([basis[i] / other for i in range(g)] +
                 [basis[i] for i in range(g, 2*g)])
        return PeriodMatrix_CM(xi=self._xi*other, basis=basis,
                               CM_type=self._CM_type, check=True)

        
def lift_ray_class_group_element(A, M, N, generator=False):
    """
    Given a ray class mod M, lift it to a ray class mod N.
    
    INPUT:
    
    - A -- a fractional ideal of a field K
    
    - M -- an integral ideal of K, such that
      the valuation of A at every divisor of M
      is zero. Actually, the current implementation
      requires that M can be defined already over the
      rationals.
          
    - N -- an integral ideal of K divisible by M
    
    - generator -- whether to also output an element x
      with x = 1 mod* M and output = A / x.
    
    OUTPUT:
    
    A fractional ideal B such that every valuation
    of A at every divisor of N is zero, and such that
    A is in the same ray class modulo M as A.
    
    Or a pair (B, x) if generator is True
    
    
    EXAMPLES::
    
        sage: from recip import *
        sage: P.<x> = QQ[]
        sage: k = NumberField(x^2+5,'a')
        sage: A = k.class_group().gen().ideal()
        sage: B = lift_ray_class_group_element(A, M=1, N=2); B
        Fractional ideal (3, a + 1)
        sage: (B/A).is_principal()
        True
        
        sage: B, x = lift_ray_class_group_element(A, M=3, N=6, generator=True)
        sage: B
        Fractional ideal (23/2, 1/2*a + 15/2)
        sage: B*x == A
        True
        sage: (x - 1)/3
        2/23*a - 7/23

    A must be coprime to M to start with::

        sage: lift_ray_class_group_element(A, M=2, N=2)
        Traceback (most recent call last):
        ...
        ValueError: A (=Fractional ideal (2, a + 1)) and M (=Fractional ideal (2)) are not coprime
    """
    A_den = A.denominator()
    if A_den != 1:
        (B_den, x_den) = lift_ray_class_group_element(A_den, M=M, N=N,
                                                      generator=True)
        (B_num, x_num) = lift_ray_class_group_element(A.numerator(), M=M, N=N,
                                                      generator=True)
        B = B_num/B_den
        if generator:
            return (B, x_num/x_den)
        return B
    K = A.number_field()
    M = K.ideal(M)
    N = K.ideal(N)
    if not (M.is_integral() and (N/M).is_integral()):
        raise ValueError, "The ideals M (=%s) and N (=%s) must be integral and M must divide N" % (M, N)
    if A + M != 1:
        raise ValueError, "A (=%s) and M (=%s) are not coprime" % (A, M)
    if A + N == 1:
        if generator:
            return (A, 1)
        return A
    # The following makes sure that y*A is coprime to N.
    y = A.idealcoprime(N)
    if M == 1:
        x_inv = y
    else:
        # If y is integral, then we can simply do y.inverse_mod(M)
        # We do know that both y and y*A are coprime to M, hence so is y.
        # We need to write y as a quotient q/r of things that are both coprime to M.
        # This is easiest if we simply assume M to be defined over QQ.
        # Then we can take r to be a rational integer denominator,
        # which is automatically coprime to M.
        # If M is not defined over QQ, then the following may fail.
        r = y.denominator_ideal().gens_two()[0]
        q = y*r
        x_inv = y * q.inverse_mod(M) / r.inverse_mod(M)
    B = x_inv * A
    if generator:
        return (B, x_inv**-1)
    return B


def random_period_matrix(prec=53, g=2):
    """
    Outputs a pseudorandom Z in H_g. Only implemented for g=2.
    
    The period matrix is in the block with diagonal entries uniformly
    in [-1/2, 1/2] + [0.1, 2]*i and off-diagonal entries with real part
    in [-1/2, 1/2] and imaginary part in [0, min(Im(diagonal entries))].
    
    EXAMPLE::

        sage: from recip import *
        sage: random_period_matrix(200, 2).parent()
        Full MatrixSpace of 2 by 2 dense matrices over Complex Field with 200 bits of precision
        sage: random_period_matrix().parent()
        Full MatrixSpace of 2 by 2 dense matrices over Complex Field with 53 bits of precision
        sage: random_period_matrix(g=3)
        Traceback (most recent call last):
        ...
        NotImplementedError

    """
    C = ComplexField(prec)
    I = C.gen()
    R = RealField(prec)
    if g != 2:
        raise NotImplementedError
    x1 = R.random_element(-1/2, 1/2)
    y1 = R.random_element(0.1, 2)
    x2 = R.random_element(-1/2, 1/2)
    y2 = R.random_element(0.1, 2)
    x3 = R.random_element(-1/2, 1/2)
    y3 = R.random_element(0, min(y1,y2))
    return Matrix(C, [[x1+I*y1, x3+I*y3],[x3+I*y3,x2+I*y2]])



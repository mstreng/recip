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

This file has code for basic polynomial operations:

 * recognize_polynomial -- for recognizing polynomials over number fields
   from numerical data
   
 * short_interpolation -- gives the `Hecke representation'

"""

from sage.matrix.constructor import Matrix
from sage.rings.polynomial.polynomial_element import is_Polynomial



def recognize_polynomial(pol_n, K=None, N=None, M=None,
                          emb=None, poly_gen=None, solution=None, check=True):
    """
    Given a polynomial pol_n over the real or complex numbers,
    guesses a polynomial over number field K close to it.
    Note, K must be embedded into the real or complex numbers.
    
    INPUT:
    
      - `pol_n` -- a polynomial over a real or complex field `A`
      
      - `K` -- a number field
      
      - `N` (default:None) -- a positive number or None,
        weigh the output minus pol_n by this factor.
        If None, determine one from the precision of pol_n
        (and assume pol_n is accurate to nearly its precision).
        Actually, the value for N=None is pretty arbitrary.
    
      - `M` (default:None) -- a positive number or None,
        weigh the coefficients of the output by this factor.
        If None, determine one from the precision of pol_n
        (and assume pol_n is accurate to nearly its precision).
        Actually, the value for M=None is pretty arbitrary.
        
      - `emb` (default:None) -- an embedding of `K` into a complex
        or real field, use this to embed `K` into `A`. If None,
        then try to coerce.
        
      - `poly_gen` (default:None) -- a polynomial ring generator
        to use for the output, or a string, or None. If a string,
        create a new polynomial ring over `K` with that string
        as variable name. If None, use variable name `y`.
    
      - `solution` (default:None) -- a pair (pol, den), where
        pol is a correct polynomial over `K` such that pol_n
        approximates pol/den . For debugging purposes only.
      
    EXAMPLES::
    
        sage: from recip import *
        sage: Q.<y> = QQ[]
        sage: K.<a> = NumberField(y^2+y-7)
        sage: P.<x> = K[]
        sage: R.<X> = RR[]
        sage: emb = K.embeddings(RR)[0]
        sage: pol_n = emb(a+4)*X^2 + emb(5*a-7/32)*X + emb((7*a-8)/29)
        sage: recognize_polynomial(pol_n, emb=emb)
        (a + 4)*y^2 + (5*a - 7/32)*y + 7/29*a - 8/29
        
        sage: P.<x> = QQ[]
        sage: K.<a> = NumberField(x^4+13*x^2+41)
        sage: emb = K.embeddings(CC)[0]
        sage: P.<X> = PolynomialRing(CC)
        sage: pol_n = emb(a/8+a^2/7+53)*X^5 + emb(a+5)*X^4 - 11/13*X + 83; pol_n
        (51.9117094301786 - 0.345009827503822*I)*X^5 + (5.00000000000000 - 2.76007862003058*I)*X^4 - 0.846153846153846*X + 83.0000000000000
        sage: recognize_polynomial(pol_n, K, emb=emb)
        (1/7*a^2 + 1/8*a + 53)*y^5 + (a + 5)*y^4 - 11/13*y + 83
        
    Larger example::

        sage: P.<y> = PolynomialRing(Rationals())
        sage: kr = NumberField(y^4+26*y^2+5,'a')
        sage: poly = y^2 - 3669057/256*y + 3255076125/64
        sage: s = poly.roots(kr)[1][0]; s
        -56133/1024*a^2 + 6608385/1024
        sage: C = ComplexField(500)
        sage: emb = kr.embeddings(C)[1]
        sage: emb(s)
        6464.12192808629278258077515091885499223178761073356317753629949149850565866041944398811655264025009640104066290246943889048395671706574175947746682694
        sage: r = recognize_polynomial(emb(s), kr, emb=emb)[0]
        sage: r == s
        True


    """
    if K == None:
        if emb == None:
            raise ValueError, "K and emb cannot both be None in " \
                              "recognize_polynomial"
        K = emb.domain()

    if K == QQ:
        bas = [1]
    else:
#        bas = [K(b) for b in K.maximal_order().basis()]
        bas = K.power_basis()
    n = len(bas)

    if not is_Polynomial(pol_n):
        A = pol_n.parent()
        pol_n = PolynomialRing(A, 'x')(pol_n)
    else:
        A = pol_n.base_ring()
    d = pol_n.degree()
        
    if N == None:
        N = 2^ZZ(floor(A.precision()*0.85))
#       I'm using powers of 10 for readability while debugging
        N = 10^(floor(log(N).n()/log(10.n())))
    if M == None:
        M = 2^ZZ(floor(A.precision()*0.05))
        M = 10^(floor(log(M).n()/log(10.n())))


    if emb == None:
        emb = K.hom(A)
    else:
        if not emb.domain() == K:
            raise ValueError, "Domain of embedding emb (=%s) is not " \
                              "%s" % (emb, K)
        if not emb.codomain() == A:
            emb2 = emb.codomain().hom(A)
            #emb = K.hom(emb2(emb.im_gens()[0]), A)
            im_gen = emb(emb.domain().gen())
            emb = K.hom(emb2(im_gen), A, check=check)
        
    # At this point, we have:
    #  -  `K` is embedded into the real or complex field `A` via `emb`. 
    #  -  `pol_n` = sum_i a_i X^i for i in range(d+1)
    #  -  `K` = directsum_j `bas[j]`*QQ for j in range(n)
    # Let pol = sum_{i,j} x_{i,j} * bas[j] * X^i for i, j as above
    # with x_{i,j} integers, and let d be an integer.
    #
    # We embed the rank ((d+1)*n+1) lattice of values for x_{i,j} and d
    # into a real vector space of dimension ((n+2)*(d+1)+1) in such a way
    # that the Euclidean norm approximates
    # N*|pol - d * pol_n|_2 + sum x_{i,j}^2 + d^2.
    #
    # We do this by making a matrix of dimensions as follows: 
    X = Matrix(ZZ, n*(d+1)+1, (n+2)*(d+1)+1)
    # the rows correspond to the basis vectors
    # (one of x_{i,j} or d is 1, the others are zero)
    # the columns to the real values of which to sum the squares.
    #
    # We begin with the terms in sum x_{i,j}^2 + d^2:
    for k in range(n*(d+1)):
        X[k, k] = M
    X[n*(d+1), n*(d+1)] = M
    for i in range(d+1):
        # Next, the terms in - N*d*pol_n
        X[n*(d+1), n*(d+1)+1 + 2*i]     = my_round((-N*pol_n[i]).real())
        X[n*(d+1), n*(d+1)+1 + 2*i + 1] = my_round((-N*pol_n[i]).imag())
        for j in range(n):
            # Finally, the terms in N*pol
            X[n*i+j, n*(d+1)+1 + 2*i]     = my_round((N*emb(bas[j])).real())
            X[n*i+j, n*(d+1)+1 + 2*i + 1] = my_round((N*emb(bas[j])).imag())
    if get_verbose():
        print "recognize_polynomial has constructed a matrix, now doing LLL"
        if get_verbose() > 1:
            print X
    
    if solution != None:
        (pol, den) = solution
        if is_Polynomial(pol):
            list_solution = flatten([Sequence(a) for a in Sequence(pol)])+[den]
        else:
            list_solution = Sequence(pol)+[den]
        print list_solution
        print sum([list_solution[k]*vector(X[k]) for k in range(n*(d+1)+1)])

    X = X.LLL()
    
    if get_verbose():
        print "recognize_polynomial finished LLL"
        if get_verbose() > 1:
            print X
    
    if type(poly_gen) == str:
        y = PolynomialRing(K, poly_gen).gen()
    elif poly_gen == None:
        y = PolynomialRing(K, 'y').gen()
    else:
        y = poly_gen
    
                
    for i in X:
        if i[n*(d+1)] != 0:
            if get_verbose():
                print i
            return sum([sum([i[k*n+l]*y**k*bas[l] for k in range(d+1)]) \
                        for l in range(n)]) / i[n*(d+1)]
    raise RuntimeError, "Failed to recognize polynomial. Bug?"


    
def short_interpolation(a, b):
    """
    Given lists a and b of the same lengths, returns
    the polynomial
        P = sum_i b[i] prod_{j!=i} (x-a[j])
    which satisfies P(a[i]) = b[i] * H'(a[i]),
    where H = prod (x-a[i]).
    
    EXAMPLES::

        sage: from recip import *
        sage: p = short_interpolation([3],[3]); p
        3
        sage: p.parent()
        Univariate Polynomial Ring in x over Integer Ring
        sage: p = short_interpolation([1,3],[2,3.]); p
        5.00000000000000*x - 9.00000000000000
        sage: p.parent()
        Univariate Polynomial Ring in x over Real Field with 53 bits of precision
    """
    universe = Sequence(a+b).universe()
    if not len(a) == len(b):
        raise ValueError, "non-equal lengths in _short_interpolation"
    P = PolynomialRing(universe, 'x')
    x = P.gen()
    return P(sum([b[i] * prod([x-a[j] for j in range(len(a)) if j != i]) \
                for i in range(len(a))]))


def my_round(x):
    """
    Rounds x to the nearest element in ZZ.
    """
    if x in ZZ:
        return ZZ(x)
    return ZZ(x.round())
    



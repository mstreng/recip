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

This file gives the Igusa invariants in terms of theta constants

Main functions:

 * igusa_modular_forms() -- gives the modular forms [h4,h6,h10,h12] in terms
   of theta constants.

 * igusa_modular_forms_to_absolute(i) -- converts a sequence i of values of
   [h4,h6,h10,h12] to the values of the absolute Igusa invariants
   from the author's thesis: h4h6/h10, h12h4^2/h10^2, h4^5/h10^2.
   
 * rosenhain_invariants(g) (for g=2 only) -- gives a triple of Rosenhain
   invariants in terms of theta constants


"""

def igusa_modular_forms_to_absolute(i):
    """
    Converts igusa modular forms to absolute Igusa invariants.
    
    INPUT:
    
      a list [h4,h6,h10,h12]
    
    OUTPUT:
    
      a list [h4*h6/h10, h12*h4^2/h10^2, h4^5/h10^2]
    
    See igusa_modular_forms().
    
    EXAMPLES::
    
        sage: from recip import *
        sage: igusa_modular_forms_to_absolute([1,2,3,4])
        [2/3, 4/9, 1/9]
    """
    [I4, I6prime, I10, h12] = i
    I2 = h12/I10
    i1 = I4*I6prime/I10
    i2 = I2*I4**2/I10
    i3 = I4**5/I10**2
    return [i1, i2, i3]


def igusa_homogeneous_to_absolute(i, prime=True):
    """
    Given a list a list [I2, I4, I6', I10] or [I2, I4, I6, I10] (depending on
    whether prime is True), returns the list
    [h4*h6/h10, h12*h4^2/h10^2, h4^5/h10^2]
    = [I4*I6'/I10, I2*I4^2/I10, I4^5/I10^2]
    
    EXAMPLES::
    
        sage: from recip import *
        sage: k = [1,2,3]
        sage: j = igusa_absolute_to_homogeneous(k); j
        [2, 3, 3, 9]
        sage: igusa_homogeneous_to_absolute(j)
        [1, 2, 3]
        sage: j = igusa_absolute_to_homogeneous(k, False); j
        [2, 3, 0, 9]
        sage: igusa_homogeneous_to_absolute(j, False)
        [1, 2, 3]
        
        sage: l = [0,1,2]
        sage: igusa_homogeneous_to_absolute(igusa_absolute_to_homogeneous(l))
        [0, 1, 2]
        sage: l = [1,0,2]
        sage: igusa_homogeneous_to_absolute(igusa_absolute_to_homogeneous(l))
        [1, 0, 2]
        sage: l = [1,2,0]
        sage: igusa_homogeneous_to_absolute(igusa_absolute_to_homogeneous(l))
        Traceback (most recent call last):
        ...
        ValueError: Invalid input: no homogeneous invariants exist for this triple of absolute invariants
        sage: l = [0,0,1]
        sage: igusa_homogeneous_to_absolute(igusa_absolute_to_homogeneous(l))
        [0, 0, 1]
    """
    if prime:
        [I2, I4, I6prime, I10] = i
    else:
        [I2, I4, I6, I10] = i
        I6prime = (I2*I4-3*I6)/2
    h12 = I2*I10
    return igusa_modular_forms_to_absolute([I4, I6prime, I10, h12])
    

def igusa_absolute_to_homogeneous(i, prime=True):
    """
    Given a list [h4*h6/h10, h12*h4^2/h10^2, h4^5/h10^2]
    = [I4*I6'/I10, I2*I4^2/I10, I4^5/I10^2], returns
    a list [I2, I4, I6', I10] or [I2, I4, I6, I10] (depending on whether
    prime is True, and uniquely determined only up to weighted scaling).
    """
    [i1, i2, i3] = i
    if i3 == 0:
        if not (i1 == 0 and i2 == 0):
            raise ValueError, "Invalid input: no homogeneous invariants exist for this triple of absolute invariants"
        raise ValueError, "Invalid input: homogeneous invariants not well-defined as all we know is I4=0"
    # We scale so that I4^2/I10=1
    I2 = i2        # = (I4^2/I10) * I2
    I4 = i3        # = (I4^2/I10)^2*I4
    I10 = I4**2    # = (I4^2/I10)^-1 * I4^2
    I6prime = i1*I10/I4
    if prime:
        return [I2, I4, I6prime, I10]
    I6 = (I2*I4-2*I6prime)/3
    return [I2, I4, I6, I10]

def igusa_invariants_absolute():
    """
    Returns the absolute Igusa invariants according to
    the choices made in the author's thesis.
    
    This is slow, because a they are completely written out as
    a rational function in theta constants.
    
    When evaluating modular forms, it may be faster to evaluate
    igusa_modular_forms first, and then apply igusa_modular_forms_to_absolute
    to the output.
    
    EXAMPLES::
    
        sage: from recip import *
        sage: invs = igusa_invariants_absolute() # long time
        sage: len(invs) # long time
        3
    """
    return igusa_modular_forms_to_absolute(igusa_modular_forms())

def igusa_modular_forms():
    """
    Returns the homogeneous Igusa invariants h4=I4, h6=I6', h10=I10, h12=I10*I2.

    EXAMPLES:
    
    Check that the Igusa invariants are indeed sums over orbits of a product
    of theta series::
    
        sage: from recip import *
        sage: gens = symplectic_generators(2)
        sage: th = theta_ring(2,2)[0].gens()
        sage: h4 = my_sum(ThetaModForm(th[0]**8).orbit(gens)) # long time
        sage: h6 = my_sum(ThetaModForm((th[0]*th[1]*th[2])**4).orbit(gens)) # long time
        sage: h10 = my_sum(ThetaModForm((th[0]*th[1]*th[2]*th[3]*th[4]*th[6]*th[8]*th[9]*th[12]*th[15])**2).orbit(gens))
        sage: h12 = my_sum(ThetaModForm((th[4]*th[6]*th[8]*th[9]*th[12]*th[15])**4).orbit(gens)) # long time
        sage: igusa_modular_forms() == [h4, h6, h10, h12] # long time
        True
    """
    gens = symplectic_generators(2)
    th = theta_ring(2,2)[0].gens()
    h4 = ThetaModForm(sum([th[i]**8 for i in [0,1,2,3,4,6,8,9,12,15]]))
    th0 = th[0]
    th1 = th[1]
    th2 = th[2]
    th3 = th[3]
    th4 = th[4]
    th6 = th[6]
    th8 = th[8]
    th9 = th[9]
    th12= th[12]
    th15= th[15]
    h6 = ThetaModForm(th0**4*th1**4*th2**4 + th0**4*th1**4*th3**4 + th0**4*th2**4*th3**4 + th1**4*th2**4*th3**4 - th0**4*th2**4*th4**4 + th1**4*th3**4*th4**4 - th0**4*th2**4*th6**4 + th1**4*th3**4*th6**4 - th0**4*th4**4*th6**4 - th1**4*th4**4*th6**4 - th2**4*th4**4*th6**4 - th3**4*th4**4*th6**4 - th0**4*th1**4*th8**4 + th2**4*th3**4*th8**4 + th0**4*th4**4*th8**4 + th3**4*th4**4*th8**4 - th1**4*th6**4*th8**4 - th2**4*th6**4*th8**4 - th0**4*th1**4*th9**4 + th2**4*th3**4*th9**4 - th1**4*th4**4*th9**4 - th2**4*th4**4*th9**4 + th0**4*th6**4*th9**4 + th3**4*th6**4*th9**4 - th0**4*th8**4*th9**4 - th1**4*th8**4*th9**4 - th2**4*th8**4*th9**4 - th3**4*th8**4*th9**4 + th1**4*th2**4*th12**4 - th0**4*th3**4*th12**4 + th0**4*th4**4*th12**4 + th1**4*th4**4*th12**4 - th2**4*th6**4*th12**4 - th3**4*th6**4*th12**4 + th0**4*th8**4*th12**4 + th2**4*th8**4*th12**4 + th4**4*th8**4*th12**4 + th6**4*th8**4*th12**4 - th1**4*th9**4*th12**4 - th3**4*th9**4*th12**4 + th4**4*th9**4*th12**4 + th6**4*th9**4*th12**4 + th1**4*th2**4*th15**4 - th0**4*th3**4*th15**4 - th2**4*th4**4*th15**4 - th3**4*th4**4*th15**4 + th0**4*th6**4*th15**4 + th1**4*th6**4*th15**4 - th1**4*th8**4*th15**4 - th3**4*th8**4*th15**4 + th4**4*th8**4*th15**4 + th6**4*th8**4*th15**4 + th0**4*th9**4*th15**4 + th2**4*th9**4*th15**4 + th4**4*th9**4*th15**4 + th6**4*th9**4*th15**4 - th0**4*th12**4*th15**4 - th1**4*th12**4*th15**4 - th2**4*th12**4*th15**4 - th3**4*th12**4*th15**4)
    h10 = ThetaModForm(prod([th[i]**2 for i in [0,1,2,3,4,6,8,9,12,15]]))
    h12 = ThetaModForm(th1**4*th2**4*th4**4*th6**4*th8**4*th9**4 + th0**4*th3**4*th4**4*th6**4*th8**4*th9**4 + th1**4*th2**4*th3**4*th4**4*th8**4*th12**4 + th0**4*th1**4*th3**4*th6**4*th8**4*th12**4 + th0**4*th2**4*th3**4*th4**4*th9**4*th12**4 + th0**4*th1**4*th2**4*th6**4*th9**4*th12**4 + th0**4*th1**4*th2**4*th4**4*th8**4*th15**4 + th0**4*th2**4*th3**4*th6**4*th8**4*th15**4 + th0**4*th1**4*th3**4*th4**4*th9**4*th15**4 + th1**4*th2**4*th3**4*th6**4*th9**4*th15**4 + th0**4*th1**4*th4**4*th6**4*th12**4*th15**4 + th2**4*th3**4*th4**4*th6**4*th12**4*th15**4 + th0**4*th2**4*th8**4*th9**4*th12**4*th15**4 + th1**4*th3**4*th8**4*th9**4*th12**4*th15**4 + th4**4*th6**4*th8**4*th9**4*th12**4*th15**4)
    return [h4, h6, h10, h12]    
    

def rosenhain_invariants(g):
    """
    Returns a (2g-1)-tuple of Rosenhain invariants e_1,...,e_{2g-1}.
    
    For g<=2, if the e_i(tau) are distinct from each other and from 0,1,inty,
    then CC^g / tau*ZZ^g+ZZ^g is the Jacobian of
    y^2 = x(x-1)*prod_i (x-e_i(tau)).
    
    EXAMPLES::

        sage: from recip import *
        sage: rosenhain_invariants(2)
        [t0^2*t1^2/(t2^2*t3^2), t1^2*t12^2/(t2^2*t15^2), t0^2*t12^2/(t3^2*t15^2)]
        sage: rosenhain_invariants(1)
        Traceback (most recent call last):
        ...
        NotImplementedError: Sorry, Rosenhain invariants currently only implemented for g=2
    """
    if g!=2:
        raise NotImplementedError, "Sorry, Rosenhain invariants currently only implemented for g=2"
    t1 = ThetaModForm([0,0,0,0], den=2, g=2)
    t2 = ThetaModForm([0,0,1/2,1/2], den=2, g=2)
    t3 = ThetaModForm([0,0,1/2,0], den=2, g=2)
    t4 = ThetaModForm([0,0,0,1/2], den=2, g=2)
    e1 = t1**2*t3**2*t2**-2*t4**-2
    t8 = ThetaModForm([1/2,1/2,0,0],g=2,den=2)
    t10 = ThetaModForm([1/2,1/2,1/2,1/2],g=2,den=2)
    e2 = t3**2*t8**2/t4**2/t10**2
    e3 = t1**2*t8**2/t2**2/t10**2
    return [e1,e2,e3]
    
def h6_thesis():
    """
    Reconstructs h6 from the formula in the author's thesis.
    For verification of the thesis and testing of the code only.
    
    EXAMPLE::
    
        sage: from recip import *
        sage: h6_thesis() == igusa_modular_forms()[1]
        True
    """
    C = cartesian_product_iterator([[0,1/2] for i in range(4)])
    T = [vector(c) for c in C if ZZ(4*(c[0]*c[2]+c[1]*c[3])) % 2 == 0]
    Z = [C for C in subsets(T) if len(C) == 3 and 1/2 *(2*sum(C) % 2) in T and not 1/2*(2*sum(C) % 2) in C]
    th = theta_ring(2,2)[0].gens()
    ret = 0
    for [b, c, d] in Z:
        b1 = 2*b[0]
        b2 = 2*b[1]
        b3 = 2*b[2]
        b4 = 2*b[3]
        c1 = 2*c[0]
        c2 = 2*c[1]
        c3 = 2*c[2]
        c4 = 2*c[3]
        d1 = 2*d[0]
        d2 = 2*d[1]
        d3 = 2*d[2]
        # d4 = 2*d[3] doesn't occur in the formula
        e = b1 + b2 + c1 + c2 + d1 + d2 + b1*c1 + b2*c2 + b4*c2 + b1*c3\
            - b2*c4+b1*d1 - b3*d1 + c1*d1 + b2*d2 + c2*d2 + c4*d2 + c1*d3\
            - b2*b3*c1+b2*b4*c2 - b1*b2*c3 - b2*b3*d1 - b3*c1*d1 - b1*c3*d1\
            - b2*c3*d1 - b2*b4*d2-b4*c2*d2 - b1*b2*d3 - b1*c1*d3 - b2*c1*d3
        ret = ret + (-1)**ZZ(e) * th[c_to_num(b, 2)]**4 \
                                * th[c_to_num(c, 2)]**4 \
                                * th[c_to_num(d, 2)]**4
    return ThetaModForm(ret)

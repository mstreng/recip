"""
RECIP -- REpository of Complex multIPlication SageMath code.

this file:
cm_types.sage

This file contains functions and classes for CM-fields and CM-types

See the file README.txt for version information, instructions, and references.

#*****************************************************************************
# Copyright (C) 2010 -- 2020
# Marco Streng <marco.streng@gmail.com>
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

from sage.rings.number_field.number_field_base import is_NumberField
from sage.rings.number_field.number_field_element import is_NumberFieldElement
from sage.rings.number_field.number_field import NumberField_absolute
from sage.rings.polynomial.polynomial_element import is_Polynomial
from sage.rings.all import (CC, ZZ, AA)
from sage.rings.complex_field import is_ComplexField
from sage.rings.complex_interval_field import is_ComplexIntervalField
from sage.structure.element import is_RingElement
from sage.rings.qqbar import QQbar
#from sage.misc.cachefunc import cached_method
from sage.structure.sage_object import SageObject
from sage.rings.rational_field import QQ
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
from sage.rings.number_field.number_field_element import       NumberFieldElement_absolute

from time import perf_counter


def CM_Field(field, name=None, relative_name=None, check=True, embedding=None):
    r"""
    Constructs a CM_Field object, which is an absolute number field that is
    a CM-field, together with useful data such as its complex conjugation
    automorphism, its maximal real subfield, and an optional embedding into a
    complex field.

    INPUT:
        
    - ``polynomial`` -- a polynomial defining a CM-field `K`.
      Alternatively, a totally negative element `r` of a totally
      real number field `K_0`is taken to mean `K = K_0(sqrt{r})`, while a
      triple ``[D, A, B]`` is taken to mean `K = QQ[x]/(x^4+A*x^2+B)`, a
      quartic CM-field with totally real subfield `K_0 = QQ(sqrt{D})`.
    - ``name`` -- the name of the generator of `K` as an absolute
      number field
    - ``relative_name`` -- the name of the generator of the CM-field as a
      relative quadratic extension of a totally real field.
    - ``embedding`` -- (default: None) an embedding of ``K`` into a complex
      field, given by the image of a generator. If a complex field is given,
      then take an embedding such that the generator has maximal imaginary
      part.

    OUTPUT:

    - a CM_Field object corresponding to the input data

    EXAMPLES::

        sage: from recip import *
        
        sage: CM_Field([5, 165, 5445])
        CM Number Field in alpha with defining polynomial x^4 + 165*x^2 + 5445

        sage: CM_Field(-3)
        CM Number Field in alpha with defining polynomial x^2 + 3

        sage: var("x")
        x
        sage: CM_Field(x^2+1)
        CM Number Field in alpha with defining polynomial x^2 + 1

        sage: P.<Y> = QQ[]
        sage: CM_Field(Y^8+1)
        CM Number Field in alpha with defining polynomial Y^8 + 1

        sage: CM_Field(CM_Field(-7), 'b')
        CM Number Field in b with defining polynomial x^2 + 7

        sage: CM_Field(NumberField(x^4+1,'gamma'))
        CM Number Field in gamma with defining polynomial x^4 + 1
    r"""
    if is_NumberFieldElement(field) or field in QQ:
        x = PolynomialRing(field.parent().fraction_field(), 'x').gen()
        poly = x**2-field
    elif isinstance(field, list) and len(field) == 3:
        x = PolynomialRing(QQ, 'x').gen()
        D = QQ(field[0])
        A = QQ(field[1])
        B = QQ(field[2]) 
        if not ((A**2-4*B)/D).is_square():
            raise ValueError( "(A^2-4B)/D must be square in "                                "CM_field([D,A,B]), but [D,A,B] = %s" %                                field)
        poly = x**4 + A*x**2 + B
        if name == None:
            name = 'alpha'
        if relative_name == None:
            relative_name = name

        return CM_Field_quartic([D,A,B], name=name, check=check,                                      relative_name=relative_name)
    

    elif is_NumberField(field):
        if field.is_absolute():
            if name == None:
                name = str(field.gen())
        else:
            if relative_name == None:
                relative_name = str(field.gen())
        poly = field.relative_polynomial()
    elif is_Polynomial(field):
        poly = field
    else:
        try:
            poly = PolynomialRing(QQ, 'x')(field)
        except TypeError:
            raise TypeError( "Incorrect type in CM_Field")
    
    if _is_accepted_complex_field(embedding):
        embedding = highest_root(poly, embedding)
    else:
        try:
            embedding = embedding(embedding.domain().gen())
        except AttributeError:
            pass

    if name == None:
        name = 'alpha'
    if relative_name == None:
        relative_name = name

    base_ring = poly.parent().base_ring()
    
    if base_ring == QQ or base_ring == ZZ:
        if poly.degree() == 4 and poly[1]==0 and poly[3]==0 and                  poly[0] in ZZ and poly[2] in ZZ and poly[4] == 1 and                  not poly[0].is_square():
            A = ZZ(poly[2])
            B = ZZ(poly[0])
            D = ZZ(QuadraticField(A^2-4*B, 'a').discriminant())
            return CM_Field_quartic([D,A,B], name=name, check=check,                                      relative_name=relative_name)
            
        return CM_Field_absolute(poly, name=name, relative_name=relative_name,                                  check=check, embedding=embedding)

    raise NotImplementedError( "CM-fields from relative fields not yet "                                 "implemented")

    
#class CM_FieldElement_absolute(NumberFieldElement_absolute):
#    def __init__(self, parent, f):
#        r"""
#        INPUT:
#        
#        -  ``parent`` - a number field
#        
#        -  ``f`` - defines an element of a number field.
#        r"""
#        NumberFieldElement_absolute.__init__(self, parent, f)


class CM_Field_absolute(NumberField_absolute):
    _real_field = None
    
    # image of generator of _real_field in self, given as sequence:
    _real_gen_in_self = None
    
    # complex conjugate of generator of self, given as sequence:
    _complex_conjugation = None
    
    _relative_name = None
    
    # image of preferred embedding:
    _embedding = None
    
    def __init__(self, absolute_polynomial, relative_field=None,
                 real_field=None, name=None, latex_name=None,
                 relative_name=None, check=True, embedding=None):
        r"""
        This is called by the function CM_Field and initializes a CM_Field_absolute object.

        INPUT:

          - `absolute_polynomial` -- a polynomial over `QQ` defining a
            CM-field.
          - `relative_field` (default: None) -- the same CM-field as a purely
            imaginary quadratic extension of a totally real field. If None,
            then automatically computed.
          - `real_field` (default: None) -- the maximal totally real subfield
            of the CM-field. If None, then automatically computed.
          - `name` -- the name of the generator of ``self``.
          - `latex_name` -- a latex representation of the generator.
          - `relative_name` -- the name of the generator of the relative
            extension mentioned above, if it still needs to be created.
          - `check` (default: True) -- whether to check if the input satisfies
            the requirements.
          - `embedding` (default: None) -- a root of `absolute_polynomial` in
            ambient field.`self` is considered to be embedded in that ambient
            field by identifying its generator with `embedding`. If None, then
            `self` is not embedded.
        
        EXAMPLES::

            sage: from recip import *
            sage: CM_Field(x^4+5*x^2+1) # implicit doctest
            CM Number Field in alpha with defining polynomial x^4 + 5*x^2 + 1
        r"""
        NumberField_absolute.__init__(self, absolute_polynomial, name=name,
                                      latex_name=latex_name,
                                      embedding=embedding)
#        self._element_class = CM_FieldElement_absolute
#        self._zero_element = self(0)
#        self._one_element = self(1)

        if check and not self.is_totally_imaginary():
            raise ValueError( "field defined by %s is not totally "                                "imaginary" % absolute_polynomial)
        self._relative_name = relative_name
        if embedding != None:
            self._embedding = embedding
        if relative_field != None:
            raise NotImplementedError( "relative_field not yet implemented")

    def _init_CM_field_data(self):
        r"""
        Computes the complex conjugation automorphism, the real subfield,
        and the embedding of the real field into self.
        r"""
        complex_conjugation_hom, self._real_field, real_to_self_hom =                                _CM_field_data(self, check=True)
        self._complex_conjugation = complex_conjugation_hom.im_gens()[0].list()
        self._real_gen_in_self = real_to_self_hom.im_gens()[0].list()

    def abs_to_rel(self):
        r"""
        Returns the identity map from self to self.relative_field().

        EXAMPLES::

            sage: from recip import *
            sage: K = CM_Field(x^4+5*x^2+1)
            sage: K.abs_to_rel()
            Ring morphism:
              From: CM Number Field in alpha with defining polynomial x^4 + 5*x^2 + 1
              To:   Number Field in alpha0 with defining polynomial x^2 - 1/2*alpha1 over its base field
              Defn: alpha |--> alpha0
        r"""
        return self.relative_field().structure()[1]

    def a_CM_type(self):
        r"""
        Returns a CM-type of self.
        r"""
        return CM_Type(self.gen())

    def real_field(self):
        r"""
        Returns the maximal totally real subfield of self.

        EXAMPLES::

            sage: from recip import *
            sage: CM_Field([5, 165, 5445]).real_field()
            Number Field in alpha1 with defining polynomial x^2 + 330*x + 21780
            sage: CM_Field(-3).real_field()
            Number Field in alpha0 with defining polynomial x
            
        A genus-three example::
        
            sage: K = CM_Field((x^2+1)*(x^2+3)*(x^2+5)+1)
            sage: K.real_field()
            Number Field in alpha1 with defining polynomial x^3 - 9*x^2 + 23*x - 16

        r"""
        if self._real_field is None:
            self._init_CM_field_data()
        return self._real_field

    def real_generator_in_self(self):
        r"""
        Returns a generator of the real quadratic field expressed
        in terms of the power basis of self.

        EXAMPLES::

            sage: from recip import *
            sage: CM_Field([5, 165, 5445]).real_generator_in_self()
            [0, 0, 2, 0]
        r"""
        if self._real_gen_in_self is None:
            self._init_CM_field_data()
        return self._real_gen_in_self

    def real_to_self(self, x=None):
        r"""
        Given x in self.real_field(), returns x as an element of self.
        If x is None, then returns the inclusion homomorphism from
        self.real_field() to self.

        EXAMPLES::

            sage: from recip import *
            sage: K = CM_Field([5, 165, 5445])
            sage: K0 = K.real_field()
            sage: K.real_to_self(K0.gen())
            2*alpha^2

            sage: K.real_to_self()
            Ring morphism:
              From: Number Field in alpha1 with defining polynomial x^2 + 330*x + 21780
              To:   CM Number Field in alpha with defining polynomial x^4 + 165*x^2 + 5445
              Defn: alpha1 |--> 2*alpha^2
            sage: CM_Field(-3).real_to_self()
            Ring morphism:
              From: Number Field in alpha0 with defining polynomial x
              To:   CM Number Field in alpha with defining polynomial x^2 + 3
              Defn: 0 |--> 0
        r"""
        m = self.real_field().hom(self(self.real_generator_in_self()), self)
        if x is None:
            return m
        return m(x)

    def _complex_conjugation_elt(self):
        r"""
        Returns the complex conjugate of self.gen().
        r"""
        if self._complex_conjugation is None:
            self._init_CM_field_data()
        return self._complex_conjugation

    def complex_conjugation(self):
        r"""
        Returns the complex conjugation automorphism of self.

        EXAMPLES::

            sage: from recip import *
            sage: CM_Field([5, 165, 5445]).complex_conjugation()
            Ring endomorphism of CM Number Field in alpha with defining polynomial x^4 + 165*x^2 + 5445
              Defn: alpha |--> -alpha
        r"""
        return self.hom(self(self._complex_conjugation_elt()), self)

    def embedding(self, codomain=None):
        r"""
        Returns the chosen embedding of self into an ambient field or the given codomain.

        INPUT: codomain (default:None) - If specified, the output is an embedding of self into codomain.
                                         Otherwise, the output is an embedding of self into a number field
                                         or complex field.

        Reflex fields Kr of quartic non-biquadratic CM fields K(a) come with
        a preferred embedding into the Complex Lazy Field. The embedding is the
        type trace T of a is a specific element of i*RR.
        One of T and -T is in i*RR_{>0}, call it S.
        Kr is embedded into CC by sending its generator to S or S/2.
        
        EXAMPLES:
        
        When creating a CM field without specifying an embedding, the
        output is None::

            sage: from recip import *
            sage: CM_Field([5, 165, 5445]).embedding() is None
            True
            sage: CM_Field(x^4+7*x^2+8, embedding = QQbar).embedding()
        
        Reflex fields of CM types with values in CC are naturally
        embedded into CC::
        
            sage: from recip import *
            sage: K = CM_Field([764, 28, 5])
            sage: Phi = K.Phi()
            sage: Kr = Phi.reflex_field()
            sage: Kr.embedding()
            Ring morphism:
              From: CM Number Field in alphar with defining polynomial x^4 + 56*x^2 + 764
              To:   Complex Lazy Field
              Defn: alphar |--> 5.69843276304982?*I
        
            sage: Phiprime = K.Phiprime()
            sage: Krprime = Phiprime.reflex_field()
            sage: Krprime.embedding()
            Ring morphism:
              From: CM Number Field in alpharp with defining polynomial x^4 + 56*x^2 + 764
              To:   Complex Lazy Field
              Defn: alpharp |--> 4.850552962807480?*I
        
            sage: Phibar = Phi.complex_conjugate()
            sage: Krbar = Phibar.reflex_field()
            sage: Krbar is Kr
            True
            sage: Krbar.embedding() # same as for Kr
            Ring morphism:
              From: CM Number Field in alphar with defining polynomial x^4 + 56*x^2 + 764
              To:   Complex Lazy Field
              Defn: alphar |--> 5.69843276304982?*I
            sage: alpha = K.gen()
            sage: Phi.type_norm(alpha+1)
            1/2*alphar^2 + alphar + 15        
            sage: Phibar.type_norm(alpha+1)
            1/2*alphar^2 - alphar + 15
            sage: Kr.embedding()(Phi.type_norm(alpha+1))
            -1.23606797749979? + 5.69843276304982?*I
            sage: Krbar.embedding()(Phibar.type_norm(alpha+1))
            -1.23606797749979? - 5.69843276304982?*I

        The embedding of the reflex field into CC is important because the subfields
        of CC generated by complex multiplication are abelian over the reflex field,
        but not necessarily over conjugates of the reflex field.::

            sage: K = CM_Field([661, 29, 45])
            sage: Phi = K.Phi()
            sage: Zs = list(K.period_matrices_iter(CM_types=[Phi]))
            sage: h = igusa_modular_forms()
            sage: h4 = h[0]
            sage: h6 = h[1]
            sage: h10 = h[2]
            sage: prec = 300
            sage: vals = [h4(Z,prec=prec)*h6(Z,prec=prec)/h10(Z,prec=prec) for Z in Zs]
            sage: C = ComplexField(prec)
            sage: x = polygen(C)
            sage: p = prod([x-v for v in vals])
            sage: p.degree()
            3
            sage: Kr = Phi.reflex_field()
            sage: emb = Kr.embedding(C)
            sage: p = prod([x-v for v in vals])
            sage: rec = recognize_polynomial(p, Kr, emb=emb)
            sage: rec
            y^3 + (-35066459955240/841*alphar^2 - 1487392805457720/841)*y^2 + (597493296539484231514360320/841*alphar^2 + 25343519362621477629659712000/841)*y - 71865977015979239777634263040/841*alphar^2 - 3048296592726424575499072757760/841
            sage: H = NumberField(rec, 'beta')
            sage: H.is_galois_relative()     # C_3 extension
            True
            sage: H.relative_discriminant()  # unramified
            Fractional ideal (1)

            sage: wrong_emb = [phi for phi in Kr.embeddings(C)                                     if abs(abs(phi(Kr.gen()))-abs(emb(Kr.gen()))) > 1][0]
            sage: emb
            Ring morphism:
              From: CM Number Field in alphar with defining polynomial x^4 + 58*x^2 + 661
              To:   Complex Field with 300 bits of precision
              Defn: alphar |--> 6.51278802549251852882256731449156226316290238405637430796238283914922070774072807476615517*I
            sage: wrong_emb
            Ring morphism:
              From: CM Number Field in alphar with defining polynomial x^4 + 58*x^2 + 661
              To:   Complex Field with 300 bits of precision
              Defn: alphar |--> -3.94760587381785896307942602760979688120927892548204030675719951858721000896415640033478443*I
            sage: rec_wrong = recognize_polynomial(p, Kr, emb=wrong_emb)
            sage: H_wrong = NumberField(rec_wrong, 'beta')
            sage: H_wrong.is_galois_relative()    # not Galois
            False
            sage: H_wrong.relative_discriminant() # ramified
            Fractional ideal (661, 1/12*alphar^2 + 719/12)

        In fact, the two polynomials above are real and are conjugate::

            sage: Kr0 = Kr.relative_field().base_field()
            sage: a = [a for a in Kr0.automorphisms() if a(Kr0.gen()) != Kr0.gen()][0]
            sage: [a(Kr.self_to_real(c)) for c in rec] == [Kr.self_to_real(c) for c in rec_wrong]
            True

        r"""
        if self._embedding is None:
            return None
        if codomain is None:
            return self.hom(self._embedding, self._embedding.parent(), check=False)
        return self.hom(codomain(self._embedding), codomain, check=False)

    def embed(self, embedding=QQbar):
        r"""
        Changes the chosen embedding of self into an ambient field.

        INPUT:

          - `embedding` (default: QQbar) -- either a complex field,
          or an embedding of self into a complex field. If a complex field,
          then the embedding is chosen such that the generator of `self`
          has the maximal imaginary part.

        EXAMPLES:

            sage: from recip import *
            sage: K = CM_Field(-7)
            sage: K.embedding()
            sage: K.embed(CC)
            sage: K.embedding()
            Ring morphism:
              From: CM Number Field in alpha with defining polynomial x^2 + 7
              To:   Complex Field with 53 bits of precision
              Defn: alpha |--> 2.64575131106459*I
            sage: K.embed(QQbar(-sqrt(-7)))
            sage: K.embedding()
            Ring morphism:
              From: CM Number Field in alpha with defining polynomial x^2 + 7
              To:   Algebraic Field
              Defn: alpha |--> -2.645751311064591?*I

        r"""
        if _is_accepted_complex_field(embedding):
            self._embedding = highest_root(self.polynomial(), embedding)
        else:
            try:
                embedding = embedding.im_gens()[0]
            except AttributeError:
                pass
            self._embedding = embedding

    def _repr_(self):
        r"""
        Returns a string representation of `self`.

        EXAMPLE::

            sage: from recip import *
            sage: K = CM_Field(x^4+8*x^2+3)
            sage: K # implicit doctest
            CM Number Field in alpha with defining polynomial x^4 + 8*x^2 + 3
            sage: K.__repr__()
            'CM Number Field in alpha with defining polynomial x^4 + 8*x^2 + 3'
        r"""
        # The following lines seem completely natural to me, but don't work:
        # ret = NumberField_absolute._repr_(self)
        # return "CM " + ret
        #return "CM Number Field in %s with defining polynomial %s" %          #       (self._names[0], self.polynomial())
        return "CM Number Field in %s with defining polynomial %s" %                  (self._names + (self.polynomial(),))

    def CM_types(self, L=None, equivalence_classes=False):
        r"""
        Returns the CM-types Phi of self with values in L.

        Not cached.

        INPUT:

        - ``L`` -- (default: QQbar) A number field or CM-field containing a
          normal closure of ``self``, or a complex field.
        - ``equivalence_classes`` --  (default: False) whether to give
          only one element of each equivalence class.

        OUTPUT:

        The CM-types of Phi of self with values in L if L is a CM-field.
        If L is a complex field, then a normal closure
        N of self is constructed and embedded into L, and CM-types
        with values in N are given.

        EXAMPLES::

            sage: from recip import *
            sage: K = CM_Field(x^4 + 8*x^2 + 3)
            sage: K.is_galois()
            False
            sage: types = K.CM_types(QQbar, equivalence_classes=True); types
            [CM-type of CM Number Field in alpha with defining polynomial x^4 + 8*x^2 + 3
              with values in CM Number Field in b with defining polynomial x^8 + 80*x^6 + 2070*x^4 + 18800*x^2 + 32761
              given by alpha |--> [7/37467*b^7 + 445/49956*b^5 + 305/12489*b^3 - 273439/149868*b, 7/74934*b^7 + 445/99912*b^5 + 305/24978*b^3 - 123571/299736*b], CM-type of CM Number Field in alpha with defining polynomial x^4 + 8*x^2 + 3
              with values in CM Number Field in b with defining polynomial x^8 + 80*x^6 + 2070*x^4 + 18800*x^2 + 32761
              given by alpha |--> [7/37467*b^7 + 445/49956*b^5 + 305/12489*b^3 - 273439/149868*b, -7/74934*b^7 - 445/99912*b^5 - 305/24978*b^3 + 123571/299736*b]]
            sage: types[0].domain() is K
            True
            
            sage: K = CM_Field((x^2+2)^2-2, 'a')
            sage: K.is_galois()
            True
            sage: Phi = K.CM_types()[0]; Phi
            CM-type Phi of CM Number Field in a with defining polynomial x^4 + 4*x^2 + 2
            sage: Phi.domain() is K
            True
            sage: phi = [phi for phi in Phi][0]
            sage: phi
            Ring morphism:
              From: CM Number Field in a with defining polynomial x^4 + 4*x^2 + 2
              To:   Complex Lazy Field
              Defn: a |--> 1.847759065022574?*I
            sage: phi.domain() is K
            True

        r"""
        if L is None:
            L = QQbar
        if get_verbose():
            print("Computing CM-types")
            t = clock()
        
        if _is_accepted_complex_field(L):
            N = self.galois_closure('b')
            N._embedding = N.embeddings(L)[0].im_gens()[0]
        else:
            N = L
        if equivalence_classes:
            aut = self.automorphisms()
        else:
            aut = [self.hom(self)]
        emb = self.embeddings(N)
        if not len(emb) == self.degree():
            raise ValueError( "The field L (=%s) does not allow enough "                                "embeddings of self (=%s)" % (L, self))
        emb_pairs = []
        g = self.gen()
        embs_done = []
        c = self.complex_conjugation()
        for phi in emb:
            if not phi(g) in embs_done:
                emb_pairs += [phi]
                embs_done += [phi(g), phi(c(g))]
        Phis = []
        Phis_done = []
        from sage.misc.mrange import cartesian_product_iterator
        for i in cartesian_product_iterator([[0,1] for k in                                               range(len(emb_pairs))]):
            Phi = [self.hom(emb_pairs[k]((c**i[k])(g)),N) for k in                                               range(len(emb_pairs))]
            from sage.sets.set import Set
            if len(aut) == 1:
                Phis += [Phi]
            elif not Set([phi(g) for phi in Phi]) in Phis_done:
                Phis += [Phi]
                Phis_done += [Set([phi(a(g)) for phi in Phi]) for a in aut]
        ret = [CM_Type(Phi) for Phi in Phis]
        
        if get_verbose():
            print("Finished computing CM-types in %s seconds" % (clock()-t))
        return ret

    def is_totally_imaginary_element(self, a):
        r"""
        Return True if and only if a is totally imaginary.
        
        EXAMPLES::

            sage: from recip import *
            sage: K = CM_Field([5,5,5])
            sage: a = K.gen()
            sage: K.is_totally_imaginary_element(a)
            True
            sage: K.is_totally_imaginary_element(a^2)
            False
            sage: K.is_totally_imaginary_element(a+1)
            False

        r"""
        self_to_rel = self.relative_field().structure()[1]
        return self_to_rel(a)[1] != 0 and self_to_rel(a^2)[1] == 0
    
    def is_totally_real_element(self, a):
        r"""
        Return True if and only if a is totally real.

        EXAMPLES::
        
            sage: from recip import *
            sage: K = CM_Field([5,5,5])
            sage: a = K.gen()
            sage: K.is_totally_real_element(a)
            False
            sage: K.is_totally_real_element(a^2)
            True
            sage: K.is_totally_real_element(a^2+a)
            False

        r"""
        self_to_rel = self.relative_field().structure()[1]
        return self_to_rel(a)[1] == 0
    
    def is_totally_positive_real_element(self, a):
        r"""
        Return True if and only if a is totally real and totally positive.
        r"""
        self_to_rel = self.relative_field().structure()[1]
        a = self_to_rel(a)
        if a[1] != 0:
            return False
        return a[0].is_totally_positive()

    def g(self):
        r"""
        Returns half the degree of self.
        
        EXAMPLES::
        
            sage: from recip import *
            sage: CM_Field(-3).g()
            1
            sage: CM_Field(x^4+8*x^2+8).g()
            2
        r"""
        return self.real_field().degree()

    def relative_field(self):
        r"""
        Returns self as a relative number field over self.real_field().

        EXAMPLE::

            sage: from recip import *
            sage: K = CM_Field(x^4+5*x^2+1)
            sage: K.relative_field()
            Number Field in alpha0 with defining polynomial x^2 - 1/2*alpha1 over its base field
        r"""
        return self.relativize(self.real_generator_in_self(), names=self._relative_name)

    def rel_to_self(self, x):
        r"""
        Given x in self.relative_field(), returns x as an element of self.

        EXAMPLE::

            sage: from recip import *
            sage: K = CM_Field(x^4+5*x^2+1)
            sage: K.rel_to_self('alpha0')
            alpha
            sage: K.rel_to_self('alpha1')
            2*alpha^2
        r"""
        return self.relative_field().structure()[0](x)

    def self_to_rel(self, x):
        r"""
        Given x in self, returns x as an element of self.relative_field().

        EXAMPLE::

            sage: from recip import *
            sage: K = CM_Field([5,5,5])
            sage: alpha = K.gen()
            sage: K.self_to_rel(alpha)
            alpha0
        r"""
        return self.relative_field().structure()[1](x)

    def self_to_real(self, x):
        r"""
        Given x in self such that x = xbar, returns x as an element of
        self.real_field().
        r"""
        (a0, a1) = list(self.self_to_rel(x))
        assert a1 == 0
        return a0

    def relative_trace(self, x):
        r"""
        Returns x+xbar as an element of self.real_field().
        r"""
        y = x + self.complex_conjugation()(x)
        return self.self_to_real(y)

    def relative_norm(self, x):
        r"""
        Returns x*xbar as an element of self.real_field().
        r"""
        y = x * self.complex_conjugation()(x)
        return self.self_to_real(y)

    def relative_different(self):
        r"""
        Returns the relative different of self over its real field, as an
        ideal of (the maximal order of) self.
        r"""
        d = self.relative_field().relative_different()
        return self.ideal([self.rel_to_self(g) for g in d.gens()])

    #@cached_method
    def galois_closure(self, names=None, map=False, embed=True):
        r"""
        Return a CM number field K that is the Galois closure of self, i.e., is
        generated by all roots of the defining polynomial of self, and
        possibly an embedding of self into K.
     
        INPUT:
     
        - ``names`` -- variable name for Galois closure
     
        - ``map`` (default: False) -- also return an embedding of self into K

        - ``embed`` (default: True) -- if self is embedded, then give
        K an embedding that extends the embedding of self.
    
        EXAMPLES::

            sage: from recip import *
            sage: K = CM_Field(x^4 + 8*x^2 + 3)
            sage: L, m = K.galois_closure('b', map=True); L
            CM Number Field in b with defining polynomial x^8 + 80*x^6 + 2070*x^4 + 18800*x^2 + 32761
            sage: m
            Ring morphism:
              From: CM Number Field in alpha with defining polynomial x^4 + 8*x^2 + 3
              To:   CM Number Field in b with defining polynomial x^8 + 80*x^6 + 2070*x^4 + 18800*x^2 + 32761
              Defn: alpha |--> -7/74934*b^7 - 445/99912*b^5 - 305/24978*b^3 + 123571/299736*b

            sage: K = CM_Field(x^4 + 8*x^2 + 3, embedding=QQbar) # long time, does not currently work very well: embedding is not there
            sage: L = K.galois_closure(embed=True); L # long time
            CM Number Field in alpha0 with defining polynomial x^8 + 80*x^6 + 2070*x^4 + 18800*x^2 + 32761
            sage: L.embedding() # long time, embedding wasn't there for K, so will not appear for L
            
        Genus-three::
        
            sage: K = CM_Field((x^2+1)*(x^2+3)*(x^2+5)+1)
            sage: K.galois_closure()
            CM Number Field in alpha0 with defining polynomial x^48 + 72*x^46 + 80*x^45 + 2692*x^44 + 4612*x^43 + 64124*x^42 + 123852*x^41 + 1028650*x^40 + 1829576*x^39 + 10447880*x^38 + 12019876*x^37 + 50258300*x^36 - 43026036*x^35 - 136889692*x^34 - 1331961552*x^33 - 2963236937*x^32 - 7718599196*x^31 - 7146191020*x^30 + 8504396980*x^29 + 87525433872*x^28 + 282027249844*x^27 + 605609380956*x^26 + 815516286320*x^25 - 88439203430*x^24 - 4254091272044*x^23 - 15185634342676*x^22 - 35530163136508*x^21 - 61057325330100*x^20 - 69796449648272*x^19 - 6160328327104*x^18 + 226759213095652*x^17 + 760552021937945*x^16 + 1718448648945196*x^15 + 3142781220480156*x^14 + 4919956841221648*x^13 + 6743690042716028*x^12 + 8183680911040904*x^11 + 8840183773407848*x^10 + 8531995881160800*x^9 + 7396811998421180*x^8 + 5807303468842192*x^7 + 4187178898937392*x^6 + 2799780355636992*x^5 + 1732733524962960*x^4 + 957303364563328*x^3 + 439000398925408*x^2 + 150222835195840*x + 29629511650768

        r"""
        if get_verbose():
            t = clock()
            print("Computing Galois closure")
        if names == None:
            names = str(self.gen())+'0'
        K, m = NumberField_absolute.galois_closure(self, names=names, map=True)
        embedding = None
        if embed and self.embedding() != None:
            ambient_ring = self.embedding().codomain()
            for r in K.absolute_polynomial().roots(ambient_ring):
                h = K.hom(r[0], ambient_ring)
                if h(m(self.gen())) == self.embedding()(self.gen()):
                    embedding = r[0]
                    break
            if embedding == None:
                raise ValueError( "Cannot find embedding of galois closure that extends embedding %s."                                    "Possibly due to precision, an incorrect ambient field, or a bug")
        L = CM_Field(K, name = names, embedding = embedding, check=False)
        if get_verbose():
            print("Finished computing Galois closure in %s seconds" % (clock() - t))
        if map:
            #m = self.hom(m(self.gen()), L)
            # m(self.gen()) is in the field K, not L (but they have the same
            # names and defining equation)
            m = self.hom(L(m(self.gen()).polynomial()), L)
            return L, m
        return L

    def class_group_iter(self):
        r"""
        Iterates over the class group, but returns ideals instead
        of ideal classes, and starts with the unit ideal before
        computing the class group. This can save time if the unit
        ideal is all you need.
        r"""
        yield self.ideal(1)
        for c in self.class_group():
            if not c.is_principal():
                yield c.ideal()

    def _principally_polarized_ideal_class_representives_iter(self, n=0, special_g2=True):
        r"""
        Iterates all pairs (A,xi)
        where A is a fractional ideal
        of the maximal order of K=self and xi is a totally imaginary generator
        of (A*Abar*Different(K))^-1, up to the action of K^* given by
        x(A,xi) = (xA, xi/(x*xbar)).

        EXAMPLES::

            sage: from recip import *
            sage: K = CM_Field(-71)
            sage: list(K._principally_polarized_ideal_class_representives_iter())
            [(Fractional ideal (1), -1/71*alpha), (Fractional ideal (1), 1/71*alpha), (Fractional ideal (2, 1/2*alpha - 1/2), -1/142*alpha), (Fractional ideal (2, 1/2*alpha - 1/2), 1/142*alpha), (Fractional ideal (4, 1/2*alpha + 3/2), -1/284*alpha), (Fractional ideal (4, 1/2*alpha + 3/2), 1/284*alpha), (Fractional ideal (3, 1/2*alpha - 1/2), -1/213*alpha), (Fractional ideal (3, 1/2*alpha - 1/2), 1/213*alpha), (Fractional ideal (3, 1/2*alpha + 1/2), -1/213*alpha), (Fractional ideal (3, 1/2*alpha + 1/2), 1/213*alpha), (Fractional ideal (4, 1/2*alpha + 5/2), -1/284*alpha), (Fractional ideal (4, 1/2*alpha + 5/2), 1/284*alpha), (Fractional ideal (2, 1/2*alpha + 1/2), -1/142*alpha), (Fractional ideal (2, 1/2*alpha + 1/2), 1/142*alpha)]

            sage: K = CM_Field(x^4 + 8*x^2 + 3)
            sage: list(K._principally_polarized_ideal_class_representives_iter())
            [(Fractional ideal (1), -1/156*alpha^3 - 17/156*alpha),
             (Fractional ideal (1), 5/156*alpha^3 + 7/156*alpha),
             (Fractional ideal (1), 1/156*alpha^3 + 17/156*alpha),
             (Fractional ideal (1), -5/156*alpha^3 - 7/156*alpha),
             (Fractional ideal (3, 1/2*alpha^3 + 7/2*alpha + 1),
              1/52*alpha^3 + 25/156*alpha),
             (Fractional ideal (3, 1/2*alpha^3 + 7/2*alpha + 1),
              -1/78*alpha^3 - 2/39*alpha),
             (Fractional ideal (3, 1/2*alpha^3 + 7/2*alpha + 1),
              -1/52*alpha^3 - 25/156*alpha),
             (Fractional ideal (3, 1/2*alpha^3 + 7/2*alpha + 1),
              1/78*alpha^3 + 2/39*alpha)]

            
        And in genus three::
        
            sage: K = CM_Field((x^2+1)*(x^2+3)*(x^2+5)+1)
            sage: list(K._principally_polarized_ideal_class_representives_iter())
            [(Fractional ideal (1), 1/916*alpha^5 + 35/916*alpha^3 + 17/916*alpha),
             (Fractional ideal (1), 29/916*alpha^5 + 99/916*alpha^3 + 35/916*alpha),
             (Fractional ideal (1), 27/916*alpha^5 + 29/916*alpha^3 + 1/916*alpha),
             (Fractional ideal (1), -133/916*alpha^5 - 533/916*alpha^3 - 429/916*alpha),
             (Fractional ideal (2, 1/8*alpha^5 - 1/2*alpha^4 + 17/8*alpha^3 - 7/2*alpha^2 + 47/8*alpha - 9/2),
              49/1832*alpha^5 + 341/1832*alpha^3 + 375/1832*alpha),
             (Fractional ideal (2, 1/8*alpha^5 - 1/2*alpha^4 + 17/8*alpha^3 - 7/2*alpha^2 + 47/8*alpha - 9/2),
              47/1832*alpha^5 + 271/1832*alpha^3 + 341/1832*alpha),
             (Fractional ideal (2, 1/8*alpha^5 - 1/2*alpha^4 + 17/8*alpha^3 - 7/2*alpha^2 + 47/8*alpha - 9/2),
              -51/1832*alpha^5 - 411/1832*alpha^3 - 409/1832*alpha),
             (Fractional ideal (2, 1/8*alpha^5 - 1/2*alpha^4 + 17/8*alpha^3 - 7/2*alpha^2 + 47/8*alpha - 9/2),
              -105/1832*alpha^5 - 469/1832*alpha^3 - 411/1832*alpha)]

        r"""
        
        if get_verbose():
            t = clock()
            print("Computing representatives of all principally polarized ideals")
        if self.g() > 2 or not special_g2:
            # The following only works for g >= 2
            D = self.different()
            c = self.complex_conjugation()
            K0 = self.real_field()

            
            # We compute a list unit_gens that is a set of representatives
            # of O_{K}^* / N_{K/K_0}(O_{K}^*)
            unit_group = self.unit_group()
            root_of_unity = unit_group.torsion_generator()
            fu = list(unit_group.fundamental_units())
            unit_gens = list(unit_group.gens())
            n = len(unit_gens)
            assert unit_gens == [root_of_unity] + fu
            reduce_by = [[root_of_unity.order()] + [0 for i in range(n-1)]]
            for u in fu:
                reduce_by.append(unit_group(self(u)*c(self(u))).list())
            M = matrix(reduce_by)
            # The map ZZ^n --> O_{K}^* / N_{K/K_0}(O_{K}^*) given by unit_gens
            # has kernel exactly ZZ^n times M
            M = M.echelon_form()
            # At this stage, M is an upper-triangular matrix with the same
            # property as above.
            d = M.diagonal()
            for pivot in d:
                assert pivot != 0
            d = [range(pivot) for pivot in d]
            # A complete set of representatives of ZZ^n / (ZZ^n times M) is
            # given by the following iterator.
            unit_it = cartesian_product_iterator(d)
            unit_list = []
            for exps in unit_it:
                unit_list.append(self(prod([unit_gens[i]**exps[i] for i in range(n)])))
            
            # We iterate over the ideal classes
            for A in self.class_group_iter():
                ideal_xi = (A*c(A)*D)**-1
                if ideal_xi.is_principal():
                    b = ideal_xi.gens_reduced()[0]
                    for u in unit_list:
                        xi = u*b
                        if self.is_totally_imaginary_element(xi):
                            yield A, xi
            return
        else:
            # The following only works for g <= 2
            #K = self
            D = self.different()
            c = self.complex_conjugation()
            K0 = self.real_field()
            u = self.unit_group()
            torsion_gen = self(u.torsion_generator())
            k = torsion_gen.multiplicative_order()
    #        A_xi = []
            # this is a silly 2^(g-1) = g-1 thing
            m = self.real_to_self()
            assert m.domain() is K0
            fu = [1] + [m(K0(ug)) for ug in K0.unit_group().fundamental_units()]
            fu = fu + [-fuelt for fuelt in fu]
            for A in self.class_group_iter():
    #            A = a.ideal()
                ideal_xi = (A*c(A)*D)**-1
                if ideal_xi.is_principal():
                    xi = ideal_xi.gens_reduced()[0]
                    # Now A*Abar*D = xi**-1, so xi is unique up to units.
                    # We need also xi totally imaginary, so replace xi by u*xi
                    # with u a unit.
                    # For g = 2, my thesis shows units are real * torsion (for g=1,
                    # this is trivial)
                    # Multiplication by the number real does not help at all, hence
                    # can take u = torsion_gen**l.
                    # Also -1 does not help either, so can take l in range(k/2).
                    for l in range(k/2):
                        # The following checks in a silly way whether 
                        # (torsion_gen**l*xi)**2 is totally real.
                        # TODO: improve the algorithm.
                        K0b, mb = self.subfield((torsion_gen**l*xi)**2, 'b')
                        if K0b.degree() == len(K0b.embeddings(AA)) and                                         (-K0b.gen()).is_totally_positive():
                            for u in fu:
                                yield (A,torsion_gen**l*xi*u)
    
    def one_period_matrix(self, CM_types=None, reduced=200):
        r"""
        Returns one isomorphism class of principally polarized
        abelian varieties with one of the given CM-types, given by a
        period matrix that remembers a triple (Phi, A, xi) as in my thesis.
        
        INPUT:
        
         - CM_types (default=None) -- a CM-type, None, or tuple of CM-types.
           If CM_types == None, then take all CM-types with values in QQbar up
           to equivalence.
         - reduced (default=100) -- if False, then don't reduce the period
           matrices. If reduced is a non-zero integer, reduce numerically with
           at least that number of digits precision.

        EXAMPLES::

            sage: from recip import *
            sage: x = var('x')
            sage: K = CM_Field(x^4+68*x^2+578)
            sage: K.one_period_matrix()
            Period Matrix
            [3.8092528978903?*I 1.5778442128152?*I]
            [1.5778442128152?*I 6.9649413235207?*I]
        r"""
        i = self.period_matrices_iter(CM_types, reduced)
        Z = next(i)
        i.close()
        return Z

    def minimal_DAB(self):
        r"""
        Returns the unique triple [D,A,B] for this field
        with D the discriminant of the real quadratic subfield
        and self isomorphic to `\QQ[X]/(X^4+AX^2+B)`.
        r"""
        D = self.real_field().discriminant()
        a = self.gen()
        cc = self.complex_conjugation()
        a = a - cc(a)
        a = a*a.denominator()
        [B,t,A,u,v] = a.minpoly().list()
        assert v == 1 and not u and not t
        return DAB_to_minimal([D,A,B])
    
    def period_matrices_iter(self, CM_types=None, reduced=200):
        r"""
        Iterates over the set of isomorphism classes of principally polarized
        abelian varieties with the given CM-types, each given by a
        period matrix that remembers a triple (Phi, A, xi) as in my thesis.
        
        INPUT:
        
         - CM_types (default=None) -- a CM-type, None, or tuple of CM-types.
           If CM_types == None, then take all CM-types with values in QQbar up
           to equivalence.
         - reduced (default=100) -- if False, then don't reduce the period
           matrices. If reduced is a non-zero integer, reduce numerically with
           at least that number of digits precision.

        EXAMPLES::

            sage: from recip import *
            sage: x = var('x')
            sage: K = CM_Field(x^4+68*x^2+578)
            sage: len(list(K.period_matrices_iter())) # long time, 3 seconds
            8
            sage: next(K.period_matrices_iter())
            Period Matrix
            [3.8092528978903?*I 1.5778442128152?*I]
            [1.5778442128152?*I 6.9649413235207?*I]
        r"""
        if reduced < 100 and reduced != False:
            reduced = 100
        if CM_types is None:
            CM_types = self.CM_types(equivalence_classes=True)
        elif not isinstance(CM_types, list):
            CM_types = [CM_types]

        A_xis = self._principally_polarized_ideal_class_representives_iter()
        for (A,xi) in A_xis:
            for Phi in CM_types:
                if isinstance(Phi, str):
                    if not Phi[:3] == 'Phi':
                        raise ValueError( "Unknown CM-type given by string: "                                            "'%s'" % Phi)
                    Phi = Phi + ' '
                    for Psi in self.CM_types():
                        if Phi in Psi._repr_():
                            Phi = Psi
                            break
                    else:
                        raise ValueError( "Unknown CM-type given by string: "                                            "'%s'" % Phi)
                # TODO: maybe use some variable precision instead of CC
                if Phi.is_totally_positive_imaginary(xi):
                    Z = PeriodMatrix(Phi,A,xi)
                    if reduced > 0:
                        # Somehow, reducing twice works much better
                        Z = Z.reduce(reduced).reduce(reduced)
                    yield Z

    def period_matrices(self, *args, **kwargs):
        r"""
        Returns the list of isomorphism classes of principally polarized
        abelian varieties with the given CM-types, each given by a
        period matrix that remembers a triple (Phi, A, xi) as in my thesis.

        INPUT:

         - CM_types (default=None) -- a CM-type, None, or tuple of CM-types.
           If CM_types == None, then take all CM-types with values in QQbar up
           to equivalence.
         - reduced (default=100) -- if False, then don't reduce the period
           matrices. If reduced is a non-zero integer, reduce numerically with
           at least that number of digits precision.

        EXAMPLES::

            sage: from recip import *
            sage: CM_Field(-3).period_matrices()
            [Period Matrix
             [-1/2*b + 1/2]]

        r"""
        return list(self.period_matrices_iter(*args, **kwargs))


def CM_Type(embeddings, reflex_field=None, reflex_to_codomain=None, check=True):
    r"""
    Create a CM-type using one of the following two possibilities:
    
    INPUT:

        - ``embeddings`` -- A set of embeddings of a CM-field to another
          CM-field

        - ``reflex_field`` -- A Number Field or CM-field object giving the
          reflex field of the CM-type. Automatically constructed if not
          supplied.

        - ``reflex_to_codomain`` -- An embedding of ``reflex_field`` into the
          codomain of the embeddings. Required if ``reflex_field`` is supplied,
          and not used otherwise.

        - ``check`` -- Whether to check if the input is valid.

    OUTPUT:

        The CM-type object given by the set of embeddings.
        
    INPUT:
    
        - ``embeddings`` -- A totally imaginary element xi of a CM-field.
        
        - ``reflex_field`` -- None
        
        - ``reflex_to_codomain`` -- A field, the codomain of the CM-type.
        
        - ``check`` -- Whether to check if the input is valid.

    OUTPUT:

        The CM-type that maps xi to the positive imaginary axis.
    
    EXAMPLES::
    
        sage: from recip import *
        sage: K = CM_Field([5,5,5])
        sage: alpha = K.gen()
        sage: s = 5+2*alpha^2
        sage: l = [alpha, s*alpha, -alpha, -s*alpha]
        sage: Phis = [CM_Type(a) for a in l]; Phis
        [CM-type Phi of CM Number Field in alpha with defining polynomial x^4 + 5*x^2 + 5,
         CM-type Phibar' of CM Number Field in alpha with defining polynomial x^4 + 5*x^2 + 5,
         CM-type Phibar of CM Number Field in alpha with defining polynomial x^4 + 5*x^2 + 5,
         CM-type Phi' of CM Number Field in alpha with defining polynomial x^4 + 5*x^2 + 5]
        sage: Matrix([[1 if Phi.is_totally_positive_imaginary(a) else 0 for Phi in Phis] for a in l])
        [1 0 0 0]
        [0 1 0 0]
        [0 0 1 0]
        [0 0 0 1]

        
    r"""
    if isinstance(embeddings, list):
        return CM_Type_embeddings(embeddings, reflex_field, reflex_to_codomain, check)
    if is_NumberFieldElement(embeddings):
        xi = embeddings
        if not reflex_field is None:
            raise ValueError()
        if xi.parent().g() > 2:
            if not reflex_to_codomain is None:
                raise ValueError()
            return CM_Type_xi(xi)
        codomain = reflex_to_codomain
        # TODO: direct implementation for CM_Type_quartic
        ret = [Phi for Phi in xi.parent().CM_types(L=codomain, equivalence_classes=False)
                   if Phi.is_totally_positive_imaginary(xi, check=check)]
        if len(ret) != 1:
            raise RuntimeError( "Found %s CM_types instead of 1 for xi=%s, probably due to a precision error" % (len(ret), xi))
        return ret[0]

class CM_Type_base(SageObject):
    r"""
    Base class for CM-types.
    r"""

#    def __init__(self, args1*, args2**) # TODO: hoe werkt dit ookalweer?
#        SageObject.__init__(self, args1*, args2**)

    def _type_norm_ideal(self, A, reflex = True):
        r"""
        INPUT:

        - ``A`` -- an ideal of self.domain()
        - ``reflex`` (default:True) -- True or False

        OUTPUT:
        
        - B : the type norm of A, which is an ideal of
          self.reflex_field() (if reflex is True) or
          self.codomain() (if reflex is False).
          
        EXAMPLES::
        
            sage: from recip import *
            sage: x = PolynomialRing(QQ,'x').gen()
            sage: K = CM_Field(x^4+15*x^2+7)
            sage: Phi = K.CM_types()[0]
            sage: Psi = Phi.reflex()
            sage: Kr = Phi.reflex_field()
            sage: A = Kr.class_group().gens()[0].ideal(); A
            Fractional ideal (3, alphar - 1)
            sage: Psi.type_norm(A)
            Fractional ideal (3, alpha^3 + 13*alpha - 1)
            
        Check that the bug with incorrect gcd's is fixed::
        
            sage: x = PolynomialRing(QQ,'x').gen()
            sage: l = CM_Field(x^4+186*x^2+5)
            sage: Psi = l.CM_types()[0]
            sage: B = l.ideal(5).factor()[0][0]
            sage: Psi.type_norm(B^2)
            Fractional ideal (25, 23/2*alphar^2 - 1/2*alphar + 1063/2)
        r"""
        den = A.denominator()
        if den != 1:
            D = den.norm()
            return self._type_norm_ideal(A*D) / self._type_norm_element(D)
        (a, b) = A.gens_two()
        if a == 0:
            a = b
            b = 0
 
        x = self.type_norm(a, reflex)
        if b == 0:
            return x.parent().ideal(x)
        
        # Now b/A and a/A are coprime integral ideals,
        # a is a rational integer, and A is integral.
        # Sanity check:
        if not a in ZZ:
            raise RuntimeError( "First element of gens_two is not a "                                  "rational integer for integral ideal %s" % A)
        z = (b/A).idealcoprime(a.parent().ideal(a))
        # Now z*b/A is an integral (according to the Pari documentation)
        # ideal coprime to a.
        y = self.type_norm(z*b, reflex)
        B = x.parent().ideal(x, y)
        # Let C be the type norm of A and v any valuation of the reflex field.
        #  * If v(a) > 0, then for any pull-back w of v also w(a) > 0,
        #    hence w(z*b) = w(A) and w(a) >= w(A).
        #    For the type norm, we get
        #    v(y) = v(C) and v(x) >= v(A), so v(B) = v(A).
        #  * If v(a) = 0, then for any pull-back w of v also w(z*b) >= w(A) = 0
        #    and w(a) = 0. For the type norm, we get
        #    v(x) = v(C) = 0 and v(y) >= v(C) = 0, so v(B) = v(A) = 0.
        if not B.norm()**A.number_field().degree() ==                 (A.norm()**B.number_field().degree())**self.g():
            raise RuntimeError( "Bug in type norm of ideals. Norm of type "                                  "norm is incorrect: %s instead of %s" %                              (B.norm(), (A.norm()**(B.number_field().degree()/                              A.number_field().degree()))**self.g()))
        return B

    def type_norm(self, x, reflex = True):
        r"""
        INPUT:

        - ``x`` -- an element or ideal of self.domain()
        - ``reflex`` -- True or False

        OUTPUT:
        
        - y : the type norm of x, which is an element or ideal of
          self.reflex_field() (if reflex is True) or
          self.codomain() (if reflex is False).
          
        EXAMPLES::
        
            sage: from recip import *
            sage: K = CM_Field(x^4+15*x^2+7)
            sage: Phi = K.CM_types()[0]
            sage: Phi.type_norm(K.gen())
            1/2*alphar^2 + 15/2
        r"""
        if x in self.domain():
            return self._type_norm_element(x, reflex)
        if is_NumberFieldFractionalIdeal(x):
            if x.number_field() == self.domain():
                return self._type_norm_ideal(x, reflex)
            if x.number_field() == self.reflex_field():
                raise ValueError( "x (%s) is an ideal of the reflex field "                                    "of %s, instead of an ideal of the "                                    "domain" % (x, self))
            raise ValueError( "x (=%s) is an ideal of %s instead of %s" %                                (x, x.number_field(), self.domain()))
        if x in self.codomain():
            raise ValueError( "x (=%s) is an element of the codomain of %s, "                                "instead of its domain" % (x, self))
        if is_NumberFieldElement(x):
            raise ValueError( "x (%s) is an element of the codomain of %s, "                                "instead of an element of the domain" % (x, self))
        raise TypeError( "x (%s) must be an element or ideal of a number "                           "field, but is of type %s" % (x, type(x)))

    def __iter__(self):
        r"""
        Returns an iterator that lists the embeddings of which this CM-type
        consists.
        
        EXAMPLE::

            sage: from recip import *
            sage: k = CM_Field(x^4+15*x^2+7)
            sage: Phi = k.CM_types()[0]
            sage: list(Phi)
            [Ring morphism:
              From: CM Number Field in alpha with defining polynomial x^4 + 15*x^2 + 7
              To:   Complex Lazy Field
              Defn: alpha |--> 3.810227607874509?*I, Ring morphism:
              From: CM Number Field in alpha with defining polynomial x^4 + 15*x^2 + 7
              To:   Complex Lazy Field
              Defn: alpha |--> 0.694381434221063?*I]
        r"""
        return self.embeddings().__iter__()

    def domain(self):
        return self._domain


class CM_Type_xi(CM_Type_base):
    r"""
    CM-type specified by an element xi.
    r"""
    _xi = None
    _domain = None
    
    def __init__(self, xi, reflex_field=None, reflex_to_codomain=None, check=True):
        r"""
        Initializes this CM-type.
        r"""
        self._xi = xi
        self._domain = xi.parent()

    @cached_method
    def embeddings(self, complex_field=None):
        if complex_field is None:
            raise ValueError( "Please specify a complex field for the embeddings of this CM-type")
        if complex_field in ZZ:
            complex_field = ComplexField(complex_field)
        i = complex_field.gen()
        if not i**2 == -1:
            raise TypeError( "Invalid complex field")
        ret = [phi for phi in self.xi().parent().embeddings(complex_field) if phi(self.xi())/i > 0]
        if len(ret) != self.g():
            raise ValueError( "Insufficient precision for embeddings of this CM-type")
        return ret
        
    def g(self):
        return self.domain().g()

    def g_reflex(self):
        raise NotImplementedError()

    def reflex_field(self):
        r"""
        Returns the reflex field of this CM-type.
        
        EXAMPLES::

            sage: from recip import *
            sage: K = CM_Field([197, 15, 7])
            sage: Phi = K.Phi()
            sage: Phi.reflex_field()
            CM Number Field in alphar with defining polynomial x^4 + 30*x^2 + 197
        r"""
        raise NotImplementedError( "Please implement reflex_field for CM types of type " + str(type(self)))

    def reflex_to_codomain(self):
        r"""
        Returns the embedding of self.reflex_field() into self.codomain()

        EXAMPLES::

            sage: from recip import *
            sage: K = CM_Field([197, 15, 7])
            sage: Phi = K.Phi()
            sage: Phi.reflex_to_codomain()
            Ring morphism:
              From: CM Number Field in alphar with defining polynomial x^4 + 30*x^2 + 197
              To:   Complex Lazy Field
              Defn: alphar |--> 4.504609042095571?*I
        r"""
        raise NotImplementedError( "Please implement reflex_to_codomain for CM types of type " + str(type(self)))
    
    def reflex(self):
        r"""
        Returns the reflex type of `self`.
        
        Assumes self is primitive (later versions can be more general).
        
        Do not override this method, override :meth:`_compute_reflex` instead.
        
        EXAMPLES::
        
            sage: from recip import *
            sage: K = CM_Field(x^4+15*x^2+7)
            sage: Phi = K.CM_types()[0]
            sage: Psi = Phi.reflex(); Psi
            CM-type Phi of CM Number Field in alphar with defining polynomial x^4 + 30*x^2 + 197
            sage: Psi.reflex()
            CM-type Phi of CM Number Field in alpha with defining polynomial x^4 + 15*x^2 + 7
            sage: Phi == Psi.reflex()
            True
        r"""
        raise NotImplementedError()
    
    def __eq__(self, other):
        r"""
        EXAMPLES::
        
            sage: from recip import *
            sage: K = CM_Field((x^2+1)*(x^2+3)*(x^2+5)+1)
            sage: CM_Type(K.gen()) == CM_Type(K.gen()^5)
            True
            sage: CM_Type(K.gen()) == CM_Type(K.gen()^7)
            False

        r"""
        
        return self.domain() == other.domain() and                 self.domain().is_totally_positive_real_element(self.xi() / other.xi())

    def xi(self):
        return self._xi

    def codomain(self):
        raise TypeError( "CM-type given by xi has not one specified Sage-codomain. The codomain is the field of complex numbers.")

    def _type_norm_element(self, x, reflex = True):
        r"""
        INPUT:

        - ``x`` -- an element of self.domain()
        - ``reflex`` (default:True) -- True or False

        OUTPUT:
        
        - ``y`` -- the type norm of x, which is an element of
          self.reflex_field() (if reflex is True) or
          self.codomain() (if reflex is False).
          
        EXAMPLES::
        
            sage: from recip import *
            sage: K = CM_Field(x^4+15*x^2+7)
            sage: Phi = K.CM_types()[0]
            sage: Phi.type_norm(K.gen())
            1/2*alphar^2 + 15/2
        r"""
        raise NotImplementedError()

    def _repr_(self):
        r"""
        _repr_ and _str_ return string representations of self.

        EXAMPLE::
        
            sage: from recip import *
            sage: K = CM_Field((x^2+1)*(x^2+3)*(x^2+5)+1)
            sage: K.a_CM_type() # indirect
            CM-type of CM Number Field in alpha with defining polynomial x^6 + 9*x^4 + 23*x^2 + 16 sending alpha to the positive imaginary axis
            sage: CM_Type(K.gen()^7)
            CM-type of CM Number Field in alpha with defining polynomial x^6 + 9*x^4 + 23*x^2 + 16 sending -9*alpha^5 - 23*alpha^3 - 16*alpha to the positive imaginary axis

        r"""
        return "CM-type of " + self.domain().__repr__()                 + " sending " + str(self._xi)                 + " to the positive imaginary axis"
               
    _str_ = _repr_

    def is_totally_positive_imaginary(self, xi, check=True):
        r"""
        Given a totally imaginary element `xi` and a CM-type ``self``, returns
        true if and only if `phi(xi))` is on the positve imaginary axis
        for all `phi` in ``self``.
        r"""
        return self.domain().is_totally_positive_real_element(xi/self._xi)


class CM_Type_embeddings(CM_Type_base):
    r"""
    An object representing a CM-type.
    r"""
    # The following are CM_Field objects:
    _domain = None
    _reflex_field = None
    _codomain = None
    # This is a map from _reflex_field to _codomain
    _reflex_to_codomain = None
    # This is a set of maps from _domain to _codomain
    _embeddings = None
    # weak reference to another CM-type
    _reflex = None

    def __init__(self, embeddings, reflex_field=None, reflex_to_codomain=None,
                 check=True):
        r"""
        Create a CM_Type object corresponding to ``embeddings```.

        INPUT:

        - ``embeddings`` -- A set of embeddings of a CM-field to another
          CM-field

        - ``reflex_field`` -- A Number Field or CM-field object giving the
          reflex field of the CM-type. Automatically constructed if not
          supplied.

        - ``reflex_to_codomain`` -- An embedding of ``reflex_field`` into the
          codomain of the embeddings. Required if ``reflex_field`` is supplied,
          and not used otherwise.

        - ``check`` -- Whether to check if the input is valid.

        OUTPUT:

        A CM-type object given by the set of embeddings.
        
        EXAMPLES::
        
            sage: from recip import *
            sage: K = CM_Field(-3)
            sage: CM_Type([K.hom(K)])
            CM-type of CM Number Field in alpha with defining polynomial x^2 + 3
              with values in CM Number Field in alpha with defining polynomial x^2 + 3
              given by alpha |--> [alpha]
            
            sage: K = CM_Field(CyclotomicField(5))
            sage: CM_Type(K.Hom(K)[0:2])
            CM-type of CM Number Field in zeta5 with defining polynomial x^4 + x^3 + x^2 + x + 1
              with values in CM Number Field in zeta5 with defining polynomial x^4 + x^3 + x^2 + x + 1
              given by zeta5 |--> [zeta5, zeta5^2]
        r"""
        self._embeddings = embeddings
        domain = embeddings[0].domain()
        self._domain = domain
        codomain = embeddings[0].codomain()
        self._codomain = codomain

        if check:
            for e in embeddings:
                if e.domain() != domain or e.codomain() != codomain:
                    raise ValueError( "Different domain or codomain in "                                        "embedding %s and %s" %                                        (embeddings[0], e))
            if 2*len(embeddings) != domain.degree():
                raise ValueError( "Incorrect number of embeddings")
            c = self._domain.complex_conjugation()
            g = domain.gen()
            for i in range(len(embeddings)):
                for j in range(i):
                    if embeddings[i](g) == embeddings[j](g) or                         embeddings[i](c(g)) == embeddings[j](g):
                        raise ValueError( "Embeddings not pairwise distinct "                                            "or not pairwise non-complex-"                                            "conjugate")
        if reflex_field != None:
            if reflex_to_codomain == None:
                raise ValueError( "reflex_to_codomain must be specified if "                                    "reflex_field is specified")
            if isinstance(reflex_field, CM_Field_absolute):
                self._reflex_field = reflex_field
            else:
                self._reflex_field = CM_Field(reflex_field)
            if check and (reflex_to_codomain.domain() != reflex_field or                            reflex_to_codomain.codomain() != codomain):
                raise ValueError( "reflex_to_codomain (=%s) is not a map "                                    "from reflex (=%s) to codomain (=%s)" %                                    (reflex_to_codomain, reflex, codomain))
            self._reflex_to_codomain = reflex_to_codomain
            if check:
                # The next line actually calculates the reflex field, which
                # may be a bit too much work. TODO: faster way to check if the
                # reflex field is correct.
                reflex_field2 = self.reflex_field()
                reflex_to_codomain2 = self.reflex_to_codomain()
                if reflex_field2.degree() != reflex_field.degree():
                    raise ValueError( "Reflex field supplied has incorrect "                                        "degree")
                try:
                    g = reflex_to_codomain2(reflex_field2.gen())
                    if reflex_to_codomain(inverse_field_hom(                                            reflex_to_codomain, g)) != g:
                        raise RuntimeError( "inverse field hom computed "                                              "incorrectly, bug in code")
                except ValueError:
                    raise ValueError( "Incorrect reflex field supplied")
                self._reflex_field = reflex_field
                self._reflex_to_codomain = reflex_to_codomain

    def embeddings(self):
        return self._embeddings
        
    def g(self):
        return self.domain().g()

    def g_reflex(self):
        return self.reflex_field().g()


    def reflex_field(self):
        r"""
        Returns the reflex field of this CM-type.
        
        EXAMPLES::

            sage: from recip import *
            sage: K = CM_Field(x^4+15*x^2+7)
            sage: Phi = K.CM_types()[0]
            sage: Phi.reflex_field()
            CM Number Field in alphar with defining polynomial x^4 + 30*x^2 + 197
        r"""
        if self._reflex_field == None:
            L = self._codomain
            # TODO:
            # I need to write down a proof somewhere that the following set
            # generates the reflex field:
            reflex_gens = [self.type_norm(self.domain().gen()+i,
                   reflex = False) for i in range(2*len(self.embeddings())+1)]
            # TODO: The following is very slow, but can be improved later:
            t = [a for a in L.subfields() if a[0].degree() % 2 == 0] # degree of a CM field is even
            t.sort(key=(lambda x : x[0].degree()))
            for Fm in t:
                F = Fm[0]; m = Fm[1]
                try:
                    for g in reflex_gens:
                        if m(inverse_field_hom(m, g)) != g:
                            raise RuntimeError( "inverse field hom computed incorrectly, bug in code")
                    #reflex_field = F
                    #reflex_to_codomain = m
                    break
                except ValueError:
                    pass
            if F.degree() == L.degree():
                self._reflex_field = self._codomain
                self._reflex_to_codomain = L.hom(L)
            else:
                if self._codomain.embedding() == None:
                    emb = None
                else:
                    emb = self._codomain.embedding()(m(F.gen()))
                F_CM = CM_Field(F, name=str(self.domain().gen())+'r',
                                embedding=emb)
                self._reflex_field = F_CM
                self._reflex_to_codomain = F_CM.hom(m(F.gen()),
                                                    self.codomain())
        return self._reflex_field

    def reflex_to_codomain(self):
        r"""
        Returns the embedding of self.reflex_field() into self.codomain()

        EXAMPLES::

            sage: from recip import *
            sage: K = CM_Field([197, 15, 7])
            sage: Phi = K.Phi()
            sage: Phi.reflex_to_codomain()
            Ring morphism:
              From: CM Number Field in alphar with defining polynomial x^4 + 30*x^2 + 197
              To:   Complex Lazy Field
              Defn: alphar |--> 4.504609042095571?*I
        r"""
        # The following line makes sure that _reflex_to_codomain is computed
        K = self.reflex_field()
        if self._reflex_to_codomain == None:
            raise RuntimeError( "_reflex_to_codomain not yet computed, "                                  "which should have happened; bug in code")
        return self._reflex_to_codomain

    def reflex(self):
        r"""
        Returns the reflex type of `self`.
        
        Assumes self is primitive (later versions can be more general).
        
        Do not override this method, override :meth:`_compute_reflex` instead.
        
        EXAMPLES::
        
            sage: from recip import *
            sage: K = CM_Field(x^4+15*x^2+7)
            sage: Phi = K.CM_types()[0]
            sage: Psi = Phi.reflex(); Psi
            CM-type Phi of CM Number Field in alphar with defining polynomial x^4 + 30*x^2 + 197
            sage: Psi.reflex()
            CM-type Phi of CM Number Field in alpha with defining polynomial x^4 + 15*x^2 + 7
            sage: Phi == Psi.reflex()
            True
        r"""
        ret = self._reflex
        if not ret is None: # if a weak reference was created
            ret = ret()
            if not ret is None: # if the weak reference has not been destroyed
                return ret

        import weakref
        # The following two lines assume self is primitive
        ret = self._compute_reflex()
        ret._reflex = weakref.ref(self)

        self._reflex = weakref.ref(ret)
        return ret
    
    def _compute_reflex(self):
        r"""
        Returns the reflex type of `self`.
        
        Assumes self is primitive (later versions can be more general).
        
        EXAMPLES::
        
            sage: from recip import *
            sage: K = CM_Field(x^4+15*x^2+7)
            sage: Phi = K.CM_types()[0]
            sage: Psi = Phi.reflex(); Psi
            CM-type Phi of CM Number Field in alphar with defining polynomial x^4 + 30*x^2 + 197
            sage: Psi.reflex()
            CM-type Phi of CM Number Field in alpha with defining polynomial x^4 + 15*x^2 + 7
            sage: Phi == Psi.reflex()
            True
        r"""
        
        
        L, m = self.domain().galois_closure(map=True)
        Lpr = self.codomain()
        mpr = self.reflex_to_codomain()
        K = self.domain()
        Kpr = self.reflex_field()
        if not L.degree() == Lpr.degree():
            raise NotImplementedError( "reflex not implemented for this "                                         "CM-type")
        g = K.gen()
        Phigs = [e(g) for e in self.embeddings()]
        S = [h for h in L.Hom(Lpr) if h(m(g)) in Phigs]
        if not 2*len(S) == L.degree():
            raise RuntimeError( "bug in reflex")
        gpr = mpr(Kpr.gen())
        Phigprs_with_duplicates = [inverse_field_hom(h, gpr) for h in S]
        Phigprs = []
        for g in Phigprs_with_duplicates:
            if not g in Phigprs:
                Phigprs = Phigprs + [g]
        Phipr = [Kpr.hom(g, L) for g in Phigprs]
        # The following line assumes self is primitive
        return CM_Type(Phipr, reflex_field=K, reflex_to_codomain=m)
    
    def __eq__(self, other):
        r"""
        EXAMPLES::
        
            sage: from recip import *
            sage: K = CM_Field([5,5,5])
            sage: l = K.CM_types()
            sage: [a == b for a in l for b in l]
            [True, False, False, False, False, True, False, False, False, False, True, False, False, False, False, True]
            
            sage: L = CM_Field(x^4+15*x^2+7)
            sage: any([L.CM_types()[0] == a for a in l])
            False
        r"""
        return self.domain() == other.domain() and                 self._prime == other._prime and                 self._conjugate == other._conjugate
    
    def codomain(self):
        return self._codomain


    def _type_norm_element(self, x, reflex = True):
        r"""
        INPUT:

        - ``x`` -- an element of self.domain()
        - ``reflex`` (default:True) -- True or False

        OUTPUT:
        
        - ``y`` -- the type norm of x, which is an element of
          self.reflex_field() (if reflex is True) or
          self.codomain() (if reflex is False).
          
        EXAMPLES::
        
            sage: from recip import *
            sage: K = CM_Field(x^4+15*x^2+7)
            sage: Phi = K.CM_types()[0]
            sage: Phi.type_norm(K.gen())
            1/2*alphar^2 + 15/2
        r"""
        y = prod([e(x) for e in self._embeddings])
        if not y.norm()**x.parent().degree() ==                 (x.norm()**y.parent().degree())**self.g():
            raise RuntimeError( "norm of type norm is incorrect: %s "                                  "instead of %s" % (y.norm(),               (x.norm()**(y.parent().degree()/x.parent().degree()))**self.g()))
        if reflex:
            return inverse_field_hom(self.reflex_to_codomain(), y)
        return y
    
    def _repr_(self):
        r"""
        _repr_ and _str_ return string representations of self.

        EXAMPLE::
        
            sage: from recip import *
            sage: K = CM_Field(CyclotomicField(5))
            sage: Phi = CM_Type(K.Hom(K)[0:2])
            sage: Phi
            CM-type of CM Number Field in zeta5 with defining polynomial x^4 + x^3 + x^2 + x + 1
              with values in CM Number Field in zeta5 with defining polynomial x^4 + x^3 + x^2 + x + 1
              given by zeta5 |--> [zeta5, zeta5^2]
        r"""
        g = self.domain().gen()
        return "CM-type of " + self.domain().__repr__()                 + "\n  with values in " + self.codomain().__repr__()                 + "\n  given by " + g.__repr__() + " |--> "                 + [h(g) for h in self.embeddings()].__repr__()
               
    _str_ = _repr_

    def is_totally_positive_imaginary(self, xi, check=True):
        r"""
        Given a totally imaginary element `xi`, returns
        true if and only if `phi(xi))` is on the positve imaginary axis
        for all `phi` in the CM-type ``self``.

        INPUT:
        
          - `xi` -- an element of the domain of ``self``
          - `check` -- If `check` is True and `xi` is not totally imaginary,
            then the output is False.
            If `check` is False, then save time by assuming that `xi` is totally
            imaginary.
            

        r"""
        if check and not self.domain().complex_conjugation()(xi) == -xi:
            return False
    
        C = self.codomain()
        if C is CLF or is_ComplexField(C) or is_ComplexIntervalField(C) or C is QQbar:
            Cemb = C
        else:
            Cemb = C.embedding()
            C = Cemb.codomain()
            if not (C is CLF or is_ComplexField(C) or is_ComplexIntervalField(C) or C is QQbar):
                raise RuntimeError()
        return all([(Cemb(phi(xi))/C.gen()) > 0 for phi in self])


class CM_Field_quartic(CM_Field_absolute):
    
    def __init__(self, DAB, relative_field=None,
                 real_field=None, name=None, latex_name=None,
                 relative_name=None, check=True, embedding=None):
        r"""
        This is called by the function CM_Field and initializes a CM_Field_absolute object.

        INPUT:

          - `DAB` -- a triple of positive integers D,A,B such that A^2 > 4B,
            D is a fundamental discriminant and (A^2-4B)/D is a square.
            The CM-field is QQ[X]/(X^4+A*X^2+B).
          - `relative_field` (default: None) -- the same CM-field as a purely
            imaginary quadratic extension of a totally real field. If None,
            then automatically computed.
          - `real_field` (default: None) -- the maximal totally real subfield
            of the CM-field. If None, then automatically computed.
          - `name` -- the name of the generator of ``self``.
          - `latex_name` -- a latex representation of the generator.
          - `relative_name` -- the name of the generator of the relative
            extension mentioned above, if it still needs to be created.
          - `check` (default: True) -- whether to check if the input satisfies
            the requirements.
          - `embedding` (default: None) -- a root of `absolute_polynomial` in
            ambient field.`self` is considered to be embedded in that ambient
            field by identifying its generator with `embedding`. If None, then
            `self` is not embedded.
        
        EXAMPLES::

            sage: from recip import *
            sage: CM_Field([21, 5,1]) # indirect doctest
            CM Number Field in alpha with defining polynomial x^4 + 5*x^2 + 1
        r"""
        x = QQ['x'].gen()
        A = DAB[1]
        B = DAB[2]
        absolute_polynomial = x^4 + A*x^2 + B
        self._DAB = DAB
        CM_Field_absolute.__init__(self, absolute_polynomial,
                 relative_field=relative_field,
                 real_field=real_field, name=name, latex_name=latex_name,
                 relative_name=relative_name, check=check, embedding=embedding)
#        self._element_class = CM_FieldElement_quartic
#        self._zero_element = self(0)
#        self._one_element = self(1)

    def CM_types(self, L=None, equivalence_classes=False):
        r"""
        Returns the CM-types Phi of self with values in L.
        
        Not cached.

        INPUT:

        - ``L`` -- (default: None) A number field or CM-field containing a
          normal closure of ``self``, or a complex field, or None, for using
          a complex lazy field, which tends to work best in practice.
          Actually all other options are buggy (TODO).
        - ``equivalence_classes`` --  (default: False) whether to give
          only one element of each equivalence class.

        OUTPUT:

        The CM-types of Phi of self with values in L if L is a CM-field.
        If L is a complex field, then a normal closure
        N of self is constructed and embedded into L, and CM-types
        with values in N are given.

        EXAMPLES::

            sage: from recip import *
            sage: K = CM_Field(x^4 + 8*x^2 + 3)
            sage: K.CM_types(equivalence_classes=True)
            [CM-type Phi of CM Number Field in alpha with defining polynomial x^4 + 8*x^2 + 3, CM-type Phi' of CM Number Field in alpha with defining polynomial x^4 + 8*x^2 + 3]
        r"""
#        ret = self._CM_types
#        if not ret is None: # if a weak reference was created
#            ret = ret()
#            if not ret is None: # if the weak reference has not been destroyed
#                return ret
#
#        ret = self._compute_CM_types(L, equivalence_classes)
#        self._CM_types = weakref.ref(ret) # Can't create weak reference to "list"!
#        return ret
        if not L is None:
            return CM_Field_absolute.CM_types(self, L, equivalence_classes)
        
        Phi = CM_Type_quartic(self._DAB, scale_down=True, prime=False, domain=self)
        
        ret = [Phi]
        
        if not (equivalence_classes and self.is_galois()):
            ret.append(CM_Type_quartic(self._DAB, scale_down=True, prime=True, domain=self))
        if not equivalence_classes:
            ret = ret + [Phi.complex_conjugate() for Phi in ret]
        return ret
    
    def Phi(self):
        return CM_Type_quartic(self._DAB, scale_down=True, prime=False, domain=self)

    def Phiprime(self):
        return CM_Type_quartic(self._DAB, scale_down=True, prime=True, domain=self)

    def DAB(self):
        return self._DAB
    
    def g(self):
        r"""
        Returns half the degree of self.
        
        EXAMPLES::
        
            sage: from recip import *
            sage: CM_Field(x^4+8*x^2+8).g()
            2
        r"""
        return ZZ(2)

    def is_D4(self):
        r"""
        Returns True if and only if self is quartic non-Galois.
        r"""
        return not self.is_galois()
        
    def is_totally_imaginary_element(self, xi):
        r"""
        Returns True if and only if ``self`` is totally imaginary.
        r"""
        xi = self(xi)
        l = xi.list()
        return l[0] == 0 and l[2] == 0
        
    def is_totally_real(self, xi):
        r"""
        Returns True if and only if ``self`` is totally real.
        r"""
        xi = self(xi)
        l = xi.list()
        return l[1] == 0 and l[3] == 0
    
    def is_totally_real_positive(self, xi, check_real=True):
        r"""
        Returns True if and only if ``self`` is totally real and totally
        positive.
        
        TODO: examples
        r"""
        xi = self(xi)
        l = xi.list()
        if check_real and not self.is_totally_real(xi):
            return False
        a = l[0]
        b = l[2]
        [D,A,B] = self.DAB()
        # 2*xi = 2*a + b * (-A + sqrt(A^2 - 4*B)) = (2*a-b*A) + b*sqrt(A^2-4*B)
        m = 2*a-b*A
        if m <= 0:
            return False
        b = abs(b)
        # now 2*xi = m +/- b*sqrt(A^2-4*B) with m positive
        return m^2 > b^2*(A^2-4*B)


class CM_Type_quartic(CM_Type_embeddings):
    r"""
    An object representing a CM-type.
    r"""
    
    _DAB = None
    _DAB_reflex = None
    
    def __init__(self, DAB, prime=False, conjugate=False, scale_down=False,
                 domain=None, reflex_field=None, check=True):
        r"""
        Create a CM_Type object.
        
        Let phi_1 : K.gen() --> i*sqrt((A+sqrt(A^2-4*B))/2)
        and phi_2 : K.gen() --> i*sqrt((A-sqrt(A^2-4*B))/2).
        (both positive imaginary)
        
        Let Phi = (phi_1, phi_2) and Phi' = (phi_1, \overline{\phi_2}).
        
        INPUT:

        - ``DAB`` -- As in the input for CM-Field
        
        - ``conjugate`` -- If True, return Phibar or Phi'bar, otherwise return
          Phi or Phi'
          
        - ``prime`` -- If True, return Phi' or Phi'bar, otherwise return
          Phi or Phibar
          
        - ``scale_down`` -- TODO: explain
        
        - ``domain`` -- Either None or a CM-Field corresponding to DAB
        
        - ``reflex_field`` -- The reflex field.
        
        - ``check`` -- Whether to check if the input is valid.

        OUTPUT:

        The CM-type object given by the input.
        
        EXAMPLES (TODO)::
        
            sage: from recip import *
            
        r"""
        self._DAB = DAB
        if domain is None:
            domain = CM_Field(DAB)
        elif check:
            if domain.DAB() != DAB:
                raise ValueError( "Incorrect domain %s in construction of "                                    "CM-type with DAB = %s" % (domain, DAB))
        self._domain = domain
        [D,A,B] = DAB
        
        # Reflex field computation:
        # alpha is a root of X^4 + A*X^2 + B
        # Phi consists of the two embeddings sending alpha to the positive imaginary axis
        # Let s = 1 for Phi and s = -1 for Phi'.
        # Let beta = Tr_{this_type}(alpha)
        #          = i(sqrt((A + sqrt(A^2-4B))/2) + s * sqrt((A - sqrt(A^2-4B))/2))
        #          = i(S), where S is positive and
        #      S^2 = ((A + sqrt(A^2-4B))/2) + ((A - sqrt(A^2-4B))/2)
        #            + s * 2 * sqrt((A + sqrt(A^2-4B))/2 * ((A - sqrt(A^2-4B))/2)
        #          = A + s * sqrt((A^2 - (A^2-4B))) = A + 2 * s * sqrt(B)
        # So  beta = i*sqrt(A+2*s*sqrt(B)).
        # It is a root of X^4 + 2*A*X^2 + (A^2 - 4*B)
        # With alphar = beta, this gives the following:
        DAB_reflex = [QuadraticField(B,'a').discriminant(), 2*A, A^2-4*B]
        if scale_down and A%2 == 0 and (A^2-4*B) % 16 == 0:
            # With alphar = beta/2, we get [D, Ar / 4, Br / 16] = [D, A/2, (A^2-4B)/16]
            DAB_reflex[1] = ZZ(A/2)
            DAB_reflex[2] = ZZ(DAB_reflex[2]/16)
            self._scale_down = True
        else:
            self._scale_down = False
        
        self._DAB_reflex = DAB_reflex
        if not reflex_field is None:
            if reflex_field.DAB() != DAB_reflex:
                raise NotImplementedError()
        else:
            name = domain._names[0]
            embedding = 2*RLF(B).sqrt() # embedding of 2sqrt(B)
            if prime:
                embedding = -embedding
            embedding = CLF.gen()*sqrt(embedding + A) # embedding of i*sqrt(A+2*s*sqrt(B))
            if conjugate:
                pass # don't do embedding = -embedding,
                     # the generator of Kr is alphar = Tr_{Phi}(alpha) = Tr_{Phibar}(alpha)bar, not Tr_{Phibar}(alpha)
            if self._scale_down:
                embedding = embedding / 2
            reflex_name = name + 'r'
            if prime:
                reflex_name = reflex_name + 'p'
            reflex_field = CM_Field_quartic([QQ(c) for c in DAB_reflex], name=reflex_name,
                                            relative_name=reflex_name, embedding=embedding)
        
        self._reflex_field = reflex_field
        
        self._prime = prime
        self._conjugate = conjugate


    def embeddings(self):
        if self._embeddings is None:
            [D,A,B] = self._DAB
            s = RLF(A**2-4*B).sqrt()
            alpha = [CLF.gen()*((A+t)/2).sqrt() for t in [s,-s]]
            if self._prime:
                alpha[1] = -alpha[1]
            if self._conjugate:
                alpha = [-a for a in alpha]
            self._embeddings = [self._domain.hom(t, CLF, check=False)                                  for t in alpha]
        return self._embeddings
        
    def g(self):
        return ZZ(2)

    def g_reflex(self):
        return ZZ(2)


    def reflex_field(self):
        r"""
        Returns the reflex field of this CM-type.
        
        EXAMPLES::

            sage: from recip import *
            sage: K = CM_Field(x^4+15*x^2+7)
            sage: Phi = K.CM_types()[0]
            sage: Phi.reflex_field()
            CM Number Field in alphar with defining polynomial x^4 + 30*x^2 + 197
        r"""
        return self._reflex_field


    def reflex_to_codomain(self):
        r"""
        Returns the embedding of self.reflex_field() into self.codomain()
        r"""
        ret = self.reflex_field().embedding()
        if ret is None:
            raise NotImplementedError()
        return ret

    def _compute_reflex(self):
        r"""
        Returns the reflex type of `self`.
        
        EXAMPLES::
        
            sage: from recip import *
            sage: K = CM_Field(x^4+15*x^2+7)
            sage: Phi = K.CM_types()[0]
            sage: Psi = Phi.reflex(); Psi
            CM-type Phi of CM Number Field in alphar with defining polynomial x^4 + 30*x^2 + 197
            sage: Psi.reflex()
            CM-type Phi of CM Number Field in alpha with defining polynomial x^4 + 15*x^2 + 7
            sage: Phi == Psi.reflex()
            True
        r"""
        return CM_Type_quartic(self._DAB_reflex,                                 scale_down=(not self._scale_down),
                               conjugate=self._conjugate,                                 domain=self._reflex_field, reflex_field=self._domain)
    
    def complex_conjugate(self):
        return CM_Type_quartic(self._DAB, scale_down=self._scale_down,
                               conjugate=(not self._conjugate),
                               prime=self._prime, reflex_field=self._reflex_field, domain=self._domain)
                     
    def domain(self):
        return self._domain

    def codomain(self):
        return self.embeddings()[0].codomain()

    def _type_norm_element(self, x, reflex = True):
        r"""
        INPUT:

        - ``x`` -- an element of self.domain()
        - ``reflex`` (default:True) -- True or False

        OUTPUT:
        
        - ``y`` -- the type norm of x, which is an element of
          self.reflex_field() (if reflex is True) or
          self.codomain() (if reflex is False).
          
        EXAMPLES::
        
            sage: from recip import *
            sage: x = QQ['x'].gen()
            sage: K = CM_Field(x^4+5*x^2+5)
            sage: a = K.gen()
            sage: l = [0,1,2,a,a^2,a^3,a+1,a^2+1,a^2+a,a^3+1,a^3+a,a^3+a^2]
            sage: Phi = K.CM_types()[0]
            sage: n = Phi.type_norm
            sage: [(x,y) for x in l for y in l if n(x)*n(y)!=n(x*y)]
            []
            sage: c = Phi.reflex_field().complex_conjugation()
            sage: all([n(x)*c(n(x)) == K(x).norm() for x in l])
            True
            sage: Psi = Phi.reflex()
            sage: Psi.reflex_field() == K
            True
            sage: a = Psi.domain().gen()
            sage: l = [0,1,2,a,a^2,a^3,a+1,a^2+1,a^2+a,a^3+1,a^3+a,a^3+a^2]
            sage: n = Psi.type_norm
            sage: [(x,y) for x in l for y in l if n(x)*n(y)!=n(x*y)]
            []
            sage: c = Psi.reflex_field().complex_conjugation()
            sage: all([n(x)*c(n(x)) == Psi.domain()(x).norm() for x in l])
            True
            
            sage: x = QQ['x'].gen()
            sage: K = CM_Field(x^4+7*x^2+3)
            sage: a = K.gen()
            sage: l = [0,1,2,a,a^2,a^3,a+1,a^2+1,a^2+a,a^3+1,a^3+a,a^3+a^2]
            sage: Phi = K.CM_types()[0]
            sage: n = Phi.type_norm
            sage: [(x,y) for x in l for y in l if n(x)*n(y)!=n(x*y)]
            []
            sage: c = Phi.reflex_field().complex_conjugation()
            sage: all([n(x)*c(n(x)) == K(x).norm() for x in l])
            True
            sage: Psi = Phi.reflex()
            sage: Psi.reflex_field() == K
            True
            sage: a = Psi.domain().gen()
            sage: l = [0,1,2,a,a^2,a^3,a+1,a^2+1,a^2+a,a^3+1,a^3+a,a^3+a^2]
            sage: n = Psi.type_norm
            sage: [(x,y) for x in l for y in l if n(x)*n(y)!=n(x*y)]
            []
            sage: c = Psi.reflex_field().complex_conjugation()
            sage: all([n(x)*c(n(x)) == Psi.domain()(x).norm() for x in l])
            True
            
        r"""
        if not reflex:
            raise NotImplementedError()
        
        [a0,a1,a2,a3] = self._domain.coerce(x).list()
                    
        # Now x = a0+a1*alpha+a2*alpha^2+a3*alpha^3
        # N_self(alpha) = (-1)^prime alpha*phi(alpha)
        #               = -(-1)^prime sqrt(B)
        # beta = tr_self(alpha)
        #      = (-1)^conjugate ( alpha + (-1)^prime phi(alpha))
        # is a root of X^4 + 2*A*X^2 + (A^2-4B)
        # if prime is True, it is in the wrong field for our convention
        # of letting DAB represent the field with the largest (in imaginary
        # part) root
        [D,A,B] = self._DAB
        y = [a0**2+a1**2*A/2+a2**2*B+a3**2*B*A/2-a0*a2*A-a1*a3*A**2/2,
             a0*a1-3*A*a0*a3/2+a1*a2*A/2+a2*a3*B,
             a1**2/2+a3**2*B/2-A*a1*a3/2,
             -a0*a3/2+a1*a2/2]
        
        if self._scale_down:
            y = [y[i]*2**i for i in range(4)]
        if self._conjugate:
            y = [y[0],-y[1],y[2],-y[3]]
        return self._reflex_field(y)

    def type_trace(x):
        r"""
        Returns the type trace of x.
        r"""
        [a0,a1,a2,a3] = self._domain.coerce(x).list()
        # The type trace of 1 is 2
        # The type trace of alpha is beta
        # The type trace of alpha^2 = (-A - sqrt(A^2-4B))/2 is -A
        # The type trace of alpha^3 = (-A*alpha - sqrt(A^2-4B)*alpha)/2 is
        # -A*beta/2 - (1/2)*sqrt(A^2-4B)*(alpha-alpha').
        # This second term has square (1/4)*(A^2-4B)(-A - 2sqrt(B))
        # which is (1/4)(A^2-4B)^2 / (-A + 2 sqrt(B))
        # beta^4 + 2Abeta^2 + (A^2-4B) = 0, so (beta^2-A)^2 = 4B-A^2,
        # As alpha^3 is totally imaginary, so is its type trace.
        if self._scale_down:
            y = [y[i]*2**i for i in range(4)]
        if self._conjugate:
            y = [y[0],-y[1],y[2],-y[3]]
        return self._reflex_field(y)
        
    def _repr_(self):
        r"""
        _repr_ and _str_ return string representations of self.

        EXAMPLE::
        
            sage: from recip import *
            sage: K = CM_Field(CyclotomicField(5))
            sage: Phi = CM_Type(K.Hom(K)[0:2])
            sage: Phi
            CM-type of CM Number Field in zeta5 with defining polynomial x^4 + x^3 + x^2 + x + 1
              with values in CM Number Field in zeta5 with defining polynomial x^4 + x^3 + x^2 + x + 1
              given by zeta5 |--> [zeta5, zeta5^2]
        r"""
        prime = "'" if self._prime else ""
        bar = "bar" if self._conjugate else ""
        return "CM-type Phi%s%s of %s" % (bar, prime, self.domain().__repr__())
               
    _str_ = _repr_
    
    def is_totally_positive_imaginary(self, xi, check=True):
        r"""
        Given a totally imaginary element `xi` and a CM-type ``self``, returns
        true if and only if `psi(xi)` is on the positve imaginary axis
        for all `psi` in ``self``.
        
        EXAMPLES::
        
            sage: from recip import *
            sage: K = CM_Field([5,5,5])
            sage: alpha = K.gen()
            sage: s = 5+2*alpha^2
            sage: l = [alpha, s*alpha, -alpha, -s*alpha]
            sage: Phis = [CM_Type(a) for a in l]; Phis
            [CM-type Phi of CM Number Field in alpha with defining polynomial x^4 + 5*x^2 + 5,
             CM-type Phibar' of CM Number Field in alpha with defining polynomial x^4 + 5*x^2 + 5,
             CM-type Phibar of CM Number Field in alpha with defining polynomial x^4 + 5*x^2 + 5,
             CM-type Phi' of CM Number Field in alpha with defining polynomial x^4 + 5*x^2 + 5]
            sage: Matrix([[1 if Phi.is_totally_positive_imaginary(a) else 0 for Phi in Phis] for a in l])
            [1 0 0 0]
            [0 1 0 0]
            [0 0 1 0]
            [0 0 0 1]

        r"""
        K = self.domain()
        if check and not K.is_totally_imaginary_element(xi):
            raise ValueError( "xi=%s is not totally imaginary" % self)
        
        xi = K(xi)
        prime = self._prime
        bar = self._conjugate
        alpha = K.gen()
        # alpha^2 is (-A+sqrt(A^2-4*B))/2
        [D,A,B] = self._DAB
        s = 2*alpha^2+A
        
        # s is a square root of A^2-4*B, sent to negative, positive
        
        # Let phi_1 : alpha --> i*sqrt((A+sqrt(A^2-4*B))/2), s --> negative
        # and phi_2 : alpha --> i*sqrt((A-sqrt(A^2-4*B))/2), s --> positive.
        
        # Let Phi = (phi_1, phi_2) and Phi' = (phi_1, \overline{\phi_2}).
        
        # Phi maps alpha to the positive imaginary axis
        # Phi' maps -alpha*s to the positive imaginary axis
        # Phibar maps -alpha ...
        # Phi'bar maps alpha*s ....
            # This yields the following check,
            # but TODO: do it in a faster and more direct way

        xi = -xi*alpha # -alpha has the same sign as 1/alpha
        if bar:
            xi = -xi
        if prime:
            xi = -xi*s # s has the same sign as 1/s
        
        assert K.is_totally_real(xi)
        return K.is_totally_real_positive(xi)


def highest_root(poly, field):
    r"""
    Returns a root of poly in field that is the highest in
    the upper half plane, i.e., that has the maximal imaginary part.

    EXAMPLES::

        sage: from recip import *
        sage: P.<x> = QQ[]
        sage: highest_root(x^7+1, CC)
        0.222520933956314 + 0.974927912181824*I
        sage: highest_root(x^6+1, QQbar)
        1*I
    r"""
    highest = None
    for r in poly.roots(field):
        if highest == None or highest.imag() < r[0].imag():
            highest = r[0]
    return highest


def _CM_field_data(K, check=True):
  r"""
  Internal function which returns a triple (c,K0,m), where c is complex
  conjugation on the CM-field K, K0 is the maximal totally real subfield, and
  m is the embedding K0 --> K
  
  EXAMPLES::
  
      sage: from recip import *
      sage: K = CM_Field((x^2+10)^2-5) # indirect doctest
  r"""
  if check:
        assert K.is_totally_imaginary()
  a = K.automorphisms()
  subfields = K.subfields()
  g = K.gen()
  for s in subfields:
    if 2*s[0].degree() == K.degree() and s[0].is_totally_real():
      g0 = s[1](s[0].gen())
      for c in a:
        if c(g) != g and c(g0) == g0:
          return c, s[0], s[1]
  raise ValueError( "Not a CM field")


def _is_accepted_complex_field(K):
    r"""
    Returns true if ``K`` is a complex field accepted by
    the CM-field code.
    r"""
    return K == CC or K == QQbar or is_ComplexField(K) or             K is CLF or is_ComplexIntervalField(K)


def _is_numerical_complex_field(K):
    r"""
    Returns true if ``K`` is a complex field for doing numerical evaluation.
    r"""
    return K == CC or is_ComplexField(K) or                        is_ComplexIntervalField(K)


def inverse_field_hom(m, y):
    r"""
    Given a field homomorphism m and y in the image of m, returns x such that
    m(x) = y

    EXAMPLE::

        sage: from recip import *
        sage: K.<a> = QuadraticField(2)
        sage: L.<b> = CyclotomicField(8)
        sage: m = K.hom(b+b^-1, K)
        sage: inverse_field_hom(m, b^3+b^5 + 7)
        -a + 7
    r"""
    K = m.domain()
    if K.degree() == 1:
        try:
            return K(QQ(y))
        except TypeError:
            raise ValueError( "y (=%s) not in the image QQ of m (=%s)" % (y,m))
    L = m.codomain()
    Lrel = L.relativize(m(K.gen()), 'c')
    (Lrel_to_L, L_to_Lrel) = Lrel.structure()
    y_in_Lrel_base = L_to_Lrel(y)[0]
    if not Lrel_to_L(y_in_Lrel_base) == y:
        raise ValueError( "y (=%s) not in the image of m (=%s)" % (y, m))
    Lrel_base_to_K = Lrel.base_field().hom(K.gen(), K)
    x = Lrel_base_to_K(y_in_Lrel_base)
    if not m(x) == y:
        raise RuntimeError( "Assertion m(x) == y failed in "                    "inverse_field_hom for m = %s, y = %s, x = %s" % (m,y,ret))
    return x


def DAB_to_minimal(DAB, K0=None, m=None, check=False):
    r"""
    Given a triple DAB, finds a triple representing the same field, but such
    that A is minimal (and B is minimal with that A).
    
    K0, if given, should be the quadratic field with discriminant D.

    m, if given, must be a positive integer satisfying m^2 * D = A^2 - 4B
    r"""
    [D,A,B] = DAB
    if K0 is None:
        K0 = QuadraticField(D,'a')
    if m is None:
        m = ZZ(sqrt((A**2-4*B)/D))
    if check:
        assert m**2*D == A**2-4*B
    sqrtD = K0.gen()
    
    # We compute an ideal b such that beta is in b if and only if
    # beta*alpha is integral, where alpha = sqrt(-A+m*sqrtD)
    # Note first that beta*alpha is integral if and only if
    # beta^2*alpha^2 is integral. So let's compute alpha^2:
    alphasq = (-A+m*sqrtD)/2
    a = K0.ideal(alphasq)
    f = a.factor()
    b = K0.ideal(1)
    for (p,e) in f:
        b = b*p**-((e/2).floor())
    # Now beta is in b if and only if beta^2*alpha^2 is integral.
    bas = b.basis()
    
    # Next, we create a matrix such that the square of the Euclidean norm of
    # beta*alpha for beta = x*bas[0]+y*bas[2] is (x,y) M (x,y)^t
    #M = Matrix([[(beta1*beta2*alphasq).trace() for beta1 in bas]      #             for beta2 in bas])
    Q = QuadraticForm(ZZ, 2, [-(bas[0]**2*alphasq).trace(),                               -2*(bas[0]*bas[1]*alphasq).trace(),                               -(bas[1]**2*alphasq).trace()])
    R, T = Q.reduced_binary_form1()
    (a,b,c) = R.coefficients()
    assert a == Q(T.column(0))
    assert c == Q(T.column(1))
    s = [T.column(0)]
    if a == c:
        s.append(T.column(1))
        if b == a:
            s.append(T.column(0)-T.column(1))
        elif b == -a:
            s.append(T.column(0)+T.column(1))
    
    # Now s contains, up to sign, all elements beta with alpha' = beta*alpha of minimal
    # Euclidean norm. The square of this Euclidean norm is a = -trace(alpha'^2) = A'
    
    DABS = []
    for beta in s:
        b = -ZZ((((beta[0]*bas[0]+beta[1]*bas[1])**2*alphasq * 2 + a)**2 - a**2) / 4)
        DABS.append([D,a,b])
    DABS.sort()
    return DABS[0]


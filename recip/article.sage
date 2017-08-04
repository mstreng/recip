"""
RECIP -- REpository of Complex multIPlication SageMath code.

this file:
article.sage

This file contains code used for creating examples of class invariants
using Shimura's reciprocity law in [Streng12].

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

EXAMPLE:

Our first example is the example in Section 6.1 of [Streng12]. This is the same
example presented at Geocrypt 2011 and ECC 2011.

We first load the package, and create the CM-field
`\QQ[X]/(X^2 + 27 X^2 + 52)`, which has real quadratic subfield of discriminant
521. Then we create one period matrix, corresponding to a CM abelian variety.::

    sage: from recip import *
    sage: K = CM_Field([521,27,52])
    sage: Z = K.one_period_matrix('Phi'); Z
    Period Matrix
    [ -0.3803191? + 0.9248553?*I   0.4251995? + 0.3851337?*I]
    [0.42519943? + 0.38513368?*I 0.26063818? + 1.36829357?*I]

At some point, we will need the absolute Igusa invariants (we use the small
choice from the author's thesis), and the real quadratic subfield of the reflex
field, which is `\QQ(\sqrt{52})`, where the 52 is the same 52 as in
[521,27,52]. We simplify things a bit by noting `52 = 2^2 * 13` and let
`a = \sqrt{13}`.::

    sage: I = igusa_invariants_absolute()
    sage: Kr0 = QuadraticField(13,'a')

Next, we take the image of a set of generators of `(O_{K^r}/8O_{K^r})^*`
under the reciprocity map `r` from the article.::

    sage: gammas = reciprocity_map_image(Z, 8)

We also create the matrix U from Section 2.8 of [1].::
    
    sage: U = Z.complex_conjugation_symplectic_matrix(8)

In Section 6.1 of [1], we see first that gamma breaks the 8th powers of the
ten even theta constants up into 4 orbits, and then that U interchanges two
of these orbits, which brings it down to 3 orbits. The following command
gives a visual presentation of these three orbits. The 2 in the command is the
even number `D` with `8 = N = 2D^2`.::

    sage: visualize(gammas + [U], 2) # random
    On the 8-th powers, the action is (note that 16 means 0):
    1: (1,6,16)(2,3,4)(8,15,9)
    2: (1,6,16)(2,3,4)(8,15,9)
    3: ()
    4: ()
    5: ()
    6: ()
    7: (1,4)(2,16)(3,6)(8,15)
    Permutation Group with generators [(), (1,4)(2,16)(3,6)(8,15), (1,6,16)(2,3,4)(8,15,9)] of order 6
    The action on the orbit [0, 1, 2, 3, 4, 6] is as follows
    [            0|            1             2             3             4             5             6             7]
    [-------------+-------------------------------------------------------------------------------------------------]
    [           t0|(-zeta8^2)*t1 (-zeta8^2)*t1            t0            t0            t0            t0            t2]
    [           t1| (zeta8^3)*t6  (zeta8^3)*t6           -t1            t1           -t1            t1            t4]
    [           t2|          -t3           -t3           -t2            t2           -t2            t2            t0]
    [           t3|   (zeta8)*t4    (zeta8)*t4            t3            t3            t3            t3            t6]
    [           t4|(-zeta8^2)*t2 (-zeta8^2)*t2            t4            t4            t4            t4            t1]
    [           t6|           t0            t0           -t6            t6           -t6            t6            t3]
    [      (zeta8)|   (-zeta8^3)    (-zeta8^3)      (-zeta8)       (zeta8)      (-zeta8)       (zeta8)    (-zeta8^3)]
    The action on the orbit [8, 9, 15] is as follows
    [             0|             1              2              3              4              5              6              7]
    [--------------+--------------------------------------------------------------------------------------------------------]
    [            t8|   (zeta8)*t15    (zeta8)*t15            -t8             t8            -t8             t8 (-zeta8^3)*t15]
    [            t9| (-zeta8^2)*t8  (-zeta8^2)*t8             t9             t9             t9             t9     (zeta8)*t9]
    [           t15| (-zeta8^2)*t9  (-zeta8^2)*t9           -t15            t15           -t15            t15  (-zeta8^3)*t8]
    [       (zeta8)|    (-zeta8^3)     (-zeta8^3)       (-zeta8)        (zeta8)       (-zeta8)        (zeta8)     (-zeta8^3)]
    The action on the orbit [12] is as follows
    [             0|             1              2              3              4              5              6              7]
    [--------------+--------------------------------------------------------------------------------------------------------]
    [           t12|   (zeta8)*t12   (-zeta8)*t12            t12           -t12           -t12           -t12 (-zeta8^3)*t12]
    [       (zeta8)|    (-zeta8^3)     (-zeta8^3)       (-zeta8)        (zeta8)       (-zeta8)        (zeta8)     (-zeta8^3)]
    [t12^3/(t8*t9*t15), zeta8*t12^3/(t8*t9*t15), i*t12^3/(t8*t9*t15), zeta8^3*t12^3/(t8*t9*t15)]

Next, we create our modular functions. The 2,2 in the following command stands
for `g=2` followed by `D=2`. In the current implementation, this really is only
for `g` and `D` small, as we currently use a dense polynomial ring in `D^{2g}`
variables.::

    sage: R = theta_ring(2,2)[0]
    sage: zeta8 = R.base_ring().gen()

We create the modular functions
`f = \zeta_8^k (\theta_{12}^3/(\theta_8\theta_9\theta_15))^n`
for `n=1,2` and `k=0,1,2,3` and compute their orbits under U and the image gammas
of the reciprocity map `r`.::

    sage: f0 = R.gens()[12]^3 / prod([R.gens()[i] for i in [8,9,15]])
    sage: T = [ThetaModForm(zeta8^k*f0^n) for k in range(4) for n in [1,2]]
    sage: fs = [t for t in T if t.is_fixed_by(gammas+[U])]; fs
    [(zeta8^2)*t12^6/(t8^2*t9^2*t15^2)]
    sage: f = fs[0]

So `k=n=2` yields the smallest `f` that works. Now we need to compute the
Galois orbit of f(Z) in order to find the minimal polynomial of f(Z).

We need to iterate over the ideals of the reflex field modulo the ideal group
of `H_{\Phi,O}(1)` of Theorem 2.2 of [1]. In this example, this group is just
the group of principal ideals (see Section 6.1 of [1]), so we iterate over the
ideal class group of the reflex field. ::

    sage: shrecip = [Z.Shimura_reciprocity(A.ideal(), n=8, m=1, transformation=False, period_matrix=True) for A in Z.CM_type().reflex_field().class_group()]
    sage: d = len(shrecip); 7
    7

Each of the 7 elements of shrecip is a list consisting of (U, u), where
`g^u ( U ) = g(Z)^a` for one of the 7 ideal classes `a` and for all g in F_N.
First, we create a list of functions `g^u` where `g` ranges over `f` and the
Igusa invariants. ::
    
    sage: funcs = [[f**u for (U,u) in shrecip]] + [[j for k in range(d)] for j in I]

Next, we evaluate these functions. For this, we use a naive implementation of
theta constant evaluation (just evaluate their Fourier expansions). Because
of the sparsity of these expansions, this is only quasi-quadratic. This is
still quite slow and we are looking forward to the implementation by Enge
and Thom\'e of the quasi-linear method of Dupont. In the mean time, we just
wait a few minutes. ::

    sage: rts = [[funcs[j][k](shrecip[k][0], prec=1000) for k in range(d)] for j in range(4)] # long time almost 3 minutes

That's that, we now have all values `g^u(U)` for all `g` under consideration.
Putting these together into minimal polynomials and Hecke representations
('intpols' below) is quick and easy. ::
    
    sage: X = PolynomialRing(rts[0][0].parent(), 'X').gen() # long time (not really, the only thing that takes time is the previous one)
    sage: minpols_num = [prod([X-r for r in rts[k]]) for k in [0,1]] # long time
    sage: minpols = [recognize_polynomial(p, Kr0, N=None, emb=None) for p in minpols_num] # long time
    sage: intpols_num = [[short_interpolation(rts[k], rts[j]) for j in range(k+1,4)] for k in [0,1]] # long time
    sage: intpols = [[recognize_polynomial(p, Kr0, N=None, emb=None) for p in s] for s in intpols_num] # long time
    sage: without_clinv = [minpols[1]] + intpols[1]; without_clinv # long time
    [y^7 + (-155205162116358647755/10201*a + 559600170220938887110/10201)*y^6 + (-152407687697460195175920750535594152550/10201*a + 549513732768094956258970636118192859400/10201)*y^5 + (-2201909580030523730272623848434538048317834513875/20402*a + 7939097894735431844153019089320973153011210882125/20402)*y^4 + (-1047175262927393182849164587480891367594710449395570625/10201*a + 3775644104882200832865729346429752069380200097845736875/10201)*y^3 + (-907392914800494855136752991106041311116404713247380607234375/20402*a + 3271651681305911192688931423723753094763461200379169938284375/20402)*y^2 + (-15014166049656519860045880222971244113390650525905069987454062500/10201*a + 54134345550367190785605984445586939893083531851405365978411062500/10201)*y - 320854170291151322128777010521751890513120770505490537777676328984375/20402*a + 1156856162931200670387093211443242850125709667683265459917987279296875/20402, (155942160719197448511497600/104060401*a - 562257456480820026589520000/104060401)*y^6 + (10915460249997911281051048769982462340888000/104060401*a - 39356251626656444452197645346830542588488000/104060401)*y^5 + (168373146277549827762748743327083206237502305644120000/104060401*a - 607078012314904622487588715272472722561309748288920000/104060401)*y^4 + (238652435808138594034697534364809573290025398381044081800000/104060401*a - 860473594320621938090309645042531340247319597576659017800000/104060401)*y^3 + (104322262281490071026402121264030948196570701298335683621237000000/104060401*a - 376139265828315347221671304362880264396213048032027706818473000000/104060401)*y^2 + (3422978757984824090435381776507567613874776553028721788254983445000000/104060401*a - 12341725426738324424199494569900641042064837165414213925317002945000000/104060401)*y + 25444851830157171979850471655958457967719020294854199175794590509050000000/104060401*a - 91742717970215413696542092180092916555240944945393570493326703047900000000/104060401, (-4012937435827235689317264963498305932800000/104060401*a + 14468851690104080323524823696416410496000000/104060401)*y^6 + (-1506913225655983606143240710922345533620775435750856464000000/104060401*a + 5433252902777486004872984777262721832123795770308636080000000/104060401)*y^5 + (-10885562196559519508593359565510553063925803986258652826001713920000000/104060401*a + 39248452661947760486383130049374702797762681330867052835728546240000000/104060401)*y^4 + (-10353823389206156431412892798015311801020716070867604781938124565025200000000/104060401*a + 37331241126881141754569229980807920080905304920194122898918893817217200000000/104060401)*y^3 + (-4485870855187365031380849955017652628982072100341118267226436318108114676000000000/104060401*a + 16174037383487540399908140658405308702506346444494976170181843406011298948000000000/104060401)*y^2 + (-148450817277784374908896669100295472510110826966451685368788333734587930637950000000000/104060401*a + 535247033579587071899540245431137800567789741426378635706026376130489367192830000000000/104060401)*y - 1586204110477319767185402578315040039982207174583408730603126895573612283837343100000000000/104060401*a + 5719140253677722869121206183461376151381382068131288660643535348602757766642039900000000000/104060401]
    sage: with_clinv = [minpols[0]] + intpols[0]; with_clinv # long time
    [y^7 + (30056912/91809*a - 105080560/91809)*y^6 + (-2510103169664/826281*a + 9050615361664/826281)*y^5 + (-31191346986373120/7436529*a + 112462009885048832/7436529)*y^4 + (-2349120383562514432/66928761*a + 8469874588158623744/66928761)*y^3 + (-26197067707249590272/22309587*a + 94454871140377034752/22309587)*y^2 + (250917334141632512/66928761*a - 904696010264018944/66928761)*y - 499961036800/91809*a + 1800799256576/91809, (155205162116358647755/10201*a - 559600170220938887110/10201)*y^6 + (4293512808155549781824000/10201*a - 15480480581665978617394880/10201)*y^5 + (4796992628983857601067816320/10201*a - 17295802891824100262459263360/10201)*y^4 + (120945583333941400455909591040/30603*a - 436075502251428602426380912640/30603)*y^3 + (1341042908932829083414600581120/10201*a - 4835198970754700340106176593920/10201)*y^2 + (-38533828940561016639837962240/91809*a + 138935696085150872191306301440/91809)*y + 167832146481204715187077120/275427*a - 605127409809396328300544000/275427, (155942160719197448511497600/104060401*a - 562257456480820026589520000/104060401)*y^6 + (-166301726476865527526007193600/312181203*a + 599609402010406725455728179200/312181203)*y^5 + (-2476327951737479042110141713203200/936543609*a + 8928527404854192205498136257740800/936543609)*y^4 + (-552147133863772061563401183356518400/25286677443*a + 1990794802746309324199714339363225600/25286677443)*y^3 + (-56403828236193801907781802283289804800/75860032329*a + 203366894838060336544487139786201497600/75860032329)*y^2 + (20008870365250292462430263758028800/8428892481*a - 72143008066021761377656153466470400/8428892481)*y - 784327759750570090436933294489600/227580096987*a + 2827933954549445333497544074854400/227580096987, (-4012937435827235689317264963498305932800000/104060401*a + 14468851690104080323524823696416410496000000/104060401)*y^6 + (-16711094246362976345862986728316071181619200000/104060401*a + 60252707174372962283239786551247279979110400000/104060401)*y^5 + (-14069073960477216641209061802589478970032128000000/104060401*a + 50726767562795827750926416934161464151375872000000/104060401)*y^4 + (-119107880103422386089716488866303676654026752000000/104060401*a + 429449569024706497123435462920155844683995545600000/104060401)*y^3 + (-11770178875083208903917412174969509904249782272000000/312181203*a + 42437983455495566570747427020049615148124471296000000/312181203)*y^2 + (12526191764620323384312023796102221506858188800000/104060401*a - 45163826693633325731560527040984621868751257600000/104060401)*y - 163671667363475421464773039962382760987852800000/936543609*a + 590126589019696595537103376578535359512576000000/936543609]

As visually very clear, the quadruple of polynomials giving `f(Z)` and
expressing `i_n(Z)` in it takes up less space than the triple of polynomials
giving `i_1(Z)` and expressing `i_n(Z)` in it, so this class invariant really
gives an improvement.

Next, we do another check on the output. It is not proven, but the
smoothness of the denominators give us a lot of confidence. ::

    sage: [(p*101^2*2^28*3^8).denominator() for p in minpols] # long time
    [1, 1]
    sage: [[(p*3^7*101^4).denominator() for p in s] for s in intpols] # long time
    [[1, 1, 1], [1, 1]]

We can actually make the polynomials smaller still, as our smaller class
invariant did introduce quite a large power of 2. ::
 
    sage: X = minpols[0].parent().gen() # long time
    sage: minpols[0](2^3*X)/2^21 # long time
    y^7 + (3757114/91809*a - 13135070/91809)*y^6 + (-39220362026/826281*a + 141415865026/826281)*y^5 + (-60920599582760/7436529*a + 219652363056736/7436529)*y^4 + (-573515718643192/66928761*a + 2067840475624664/66928761)*y^3 + (-799471060401904/22309587*a + 2882533909313264/22309587)*y^2 + (957173668448/66928761*a - 3451141396576/66928761)*y - 238400/91809*a + 858688/91809

EXAMPLE:

Our second example concerns the record field of Enge and Thom{\'e}. It is
mentioned in Section 6.2 of [1] and was presented also at Geocrypt and ECC.
We repeat more or less the same steps, but in less detail. ::

    sage: K = CM_Field([709, 310, 17644])
    sage: O = K.maximal_order()
    sage: Phi = K.CM_types()[3]
    sage: d = K.different().gens_reduced()[0]
    sage: xi = d^-1
    sage: [phi(xi) for phi in Phi]
    [0.00122524459004800?*I, 0.00216656865333580?*I]
    sage: Z = PeriodMatrix(Phi, 1*O, xi)
    sage: gammas = reciprocity_map_image(Z, 8) # long time, less than 7 minutes, note: requires sage-4.7.1 or higher, see trac ticket 11234
    sage: c = Z.complex_conjugation_symplectic_matrix(8) # long time
    sage: visualize(gammas + [c], 2) # long time
    On the 8-th powers, the action is (note that 16 means 0):
    1: ()
    2: ()
    3: ()
    4: ()
    5: ()
    6: ()
    7: ()
    8: ()
    9: ()
    Permutation Group with generators [()] of order 1
    The action on the orbit [1] is as follows
    [           0|           1            2            3            4            5            6            7            8            9]
    [------------+--------------------------------------------------------------------------------------------------------------------]
    [          t1|          t1           t1 (zeta8^2)*t1 (zeta8^2)*t1           t1           t1           t1           t1           t1]
    [     (zeta8)|    (-zeta8)      (zeta8)      (zeta8)      (zeta8)      (zeta8)      (zeta8)      (zeta8)      (zeta8)   (-zeta8^3)]
    The action on the orbit [2] is as follows
    [            0|            1             2             3             4             5             6             7             8             9]
    [-------------+-----------------------------------------------------------------------------------------------------------------------------]
    [           t2|           t2 (-zeta8^2)*t2            t2           -t2            t2            t2            t2            t2            t2]
    [      (zeta8)|     (-zeta8)       (zeta8)       (zeta8)       (zeta8)       (zeta8)       (zeta8)       (zeta8)       (zeta8)    (-zeta8^3)]
    The action on the orbit [3] is as follows
    [            0|            1             2             3             4             5             6             7             8             9]
    [-------------+-----------------------------------------------------------------------------------------------------------------------------]
    [           t3|           t3  (zeta8^2)*t3 (-zeta8^2)*t3  (zeta8^2)*t3            t3            t3            t3            t3            t3]
    [      (zeta8)|     (-zeta8)       (zeta8)       (zeta8)       (zeta8)       (zeta8)       (zeta8)       (zeta8)       (zeta8)    (-zeta8^3)]
    The action on the orbit [4] is as follows
    [         0|         1          2          3          4          5          6          7          8          9]
    [----------+--------------------------------------------------------------------------------------------------]
    [        t4|        t4        -t4         t4         t4         t4         t4         t4         t4         t4]
    [   (zeta8)|  (-zeta8)    (zeta8)    (zeta8)    (zeta8)    (zeta8)    (zeta8)    (zeta8)    (zeta8) (-zeta8^3)]
    The action on the orbit [6] is as follows
    [           0|           1            2            3            4            5            6            7            8            9]
    [------------+--------------------------------------------------------------------------------------------------------------------]
    [          t6|         -t6 (zeta8^2)*t6           t6          -t6           t6           t6           t6           t6           t6]
    [     (zeta8)|    (-zeta8)      (zeta8)      (zeta8)      (zeta8)      (zeta8)      (zeta8)      (zeta8)      (zeta8)   (-zeta8^3)]
    The action on the orbit [8] is as follows
    [         0|         1          2          3          4          5          6          7          8          9]
    [----------+--------------------------------------------------------------------------------------------------]
    [        t8|        t8        -t8        -t8        -t8         t8         t8         t8         t8         t8]
    [   (zeta8)|  (-zeta8)    (zeta8)    (zeta8)    (zeta8)    (zeta8)    (zeta8)    (zeta8)    (zeta8) (-zeta8^3)]
    The action on the orbit [9] is as follows
    [            0|            1             2             3             4             5             6             7             8             9]
    [-------------+-----------------------------------------------------------------------------------------------------------------------------]
    [           t9|          -t9           -t9 (-zeta8^2)*t9 (-zeta8^2)*t9            t9            t9            t9            t9            t9]
    [      (zeta8)|     (-zeta8)       (zeta8)       (zeta8)       (zeta8)       (zeta8)       (zeta8)       (zeta8)       (zeta8)    (-zeta8^3)]
    The action on the orbit [12] is as follows
    [         0|         1          2          3          4          5          6          7          8          9]
    [----------+--------------------------------------------------------------------------------------------------]
    [       t12|       t12        t12       -t12       -t12        t12        t12        t12        t12        t12]
    [   (zeta8)|  (-zeta8)    (zeta8)    (zeta8)    (zeta8)    (zeta8)    (zeta8)    (zeta8)    (zeta8) (-zeta8^3)]
    The action on the orbit [15] is as follows
    [             0|             1              2              3              4              5              6              7              8              9]
    [--------------+--------------------------------------------------------------------------------------------------------------------------------------]
    [           t15|          -t15  (zeta8^2)*t15  (zeta8^2)*t15 (-zeta8^2)*t15            t15            t15            t15            t15            t15]
    [       (zeta8)|      (-zeta8)        (zeta8)        (zeta8)        (zeta8)        (zeta8)        (zeta8)        (zeta8)        (zeta8)     (-zeta8^3)]
    The action on the orbit [0] is as follows
    [         0|         1          2          3          4          5          6          7          8          9]
    [----------+--------------------------------------------------------------------------------------------------]
    [        t0|        t0         t0         t0         t0         t0         t0         t0         t0         t0]
    [   (zeta8)|  (-zeta8)    (zeta8)    (zeta8)    (zeta8)    (zeta8)    (zeta8)    (zeta8)    (zeta8) (-zeta8^3)]

The prime 2 splits sufficiently much that Rosenhain invariants should always work.::

    sage: [e.orbit(gammas+[c]) == [e] for e in rosenhain_invariants(2)] # long time
    [True, True, True]
    
Here's a nice other small invariant::

    sage: R = theta_ring(2,2)[0]
    sage: zeta8 = R.base_ring().gen()    
    sage: ThetaModForm(R.gens()[4]/R.gens()[12]).orbit(gammas) # long time
    [t4/t12, (-t4)/t12]
    
Now the invariants of the talk at Geocrypt and ECC::

    sage: t = ThetaModForm(R.gens()[0]*R.gens()[8]/R.gens()[4]/R.gens()[12])
    sage: U = ThetaModForm(R.gens()[2]/R.gens()[6])^2
    sage: V = ThetaModForm(R.gens()[8]/R.gens()[12])^2
    sage: W = ThetaModForm(R.gens()[0]/R.gens()[4])^2
    sage: [a.orbit(gammas+[c]) for a in [t,U,V,W]] # long time
    [[t0*t8/(t4*t12)], [t2^2/t6^2], [t8^2/t12^2], [t0^2/t4^2]]
    sage: u = U*V
    sage: v = U*W
    sage: t^2/V*U == v
    True
    sage: [t,u,v] # these were presented first
    [t0*t8/(t4*t12), t2^2*t8^2/(t6^2*t12^2), t0^2*t2^2/(t4^2*t6^2)]
    sage: [t,U,V] # maybe these are smaller?
    [t0*t8/(t4*t12), t2^2/t6^2, t8^2/t12^2]

EXAMPLE:

In the first example, the prime 2 splits as a product of 3 ideals in K,
and in the second example, it splits completely. Here's a third example,
one where 2 is inert in the real quadratic subfield, and then ramifies in K.::

    sage: K = CM_Field([93, 11, 7])
    sage: Z = K.one_period_matrix('Phi'); Z
    Period Matrix
    [ 0.19618739074462? + 1.22627074062326?*I -0.50000000000000? + 0.43440480214564?*I]
    [-0.50000000000000? + 0.43440480214564?*I  -0.39237478148923? + 2.4525414812466?*I]
    sage: I = igusa_invariants_absolute()
    sage: Kr0 = QuadraticField(7,'a')
    sage: gammas = reciprocity_map_image(Z, 8)
    sage: c = Z.complex_conjugation_symplectic_matrix(8)
    sage: visualize(gammas + [c], 2)
    On the 8-th powers, the action is (note that 16 means 0):
    1: (1,8)(2,3)(4,12)(9,16)
    2: ()
    3: ()
    4: ()
    5: ()
    Permutation Group with generators [(), (1,8)(2,3)(4,12)(9,16)] of order 2
    The action on the orbit [1, 8] is as follows
    [            0|            1             2             3             4             5]
    [-------------+---------------------------------------------------------------------]
    [           t1|           t8           -t1            t1            t1 (-zeta8^2)*t1]
    [           t8| (zeta8^2)*t1           -t8            t8            t8           -t8]
    [      (zeta8)|     (-zeta8)      (-zeta8)       (zeta8)       (zeta8)    (-zeta8^3)]
    The action on the orbit [2, 3] is as follows
    [            0|            1             2             3             4             5]
    [-------------+---------------------------------------------------------------------]
    [           t2|(-zeta8^3)*t3  (zeta8^2)*t2            t2           -t2            t2]
    [           t3|  (-zeta8)*t2 (-zeta8^2)*t3            t3           -t3 (-zeta8^2)*t3]
    [      (zeta8)|     (-zeta8)      (-zeta8)       (zeta8)       (zeta8)    (-zeta8^3)]
    The action on the orbit [4, 12] is as follows
    [            0|            1             2             3             4             5]
    [-------------+---------------------------------------------------------------------]
    [           t4|(zeta8^3)*t12           -t4           -t4            t4           -t4]
    [          t12|(-zeta8^3)*t4           t12          -t12           t12          -t12]
    [      (zeta8)|     (-zeta8)      (-zeta8)       (zeta8)       (zeta8)    (-zeta8^3)]
    The action on the orbit [6] is as follows
    [           0|           1            2            3            4            5]
    [------------+----------------------------------------------------------------]
    [          t6|         -t6 (zeta8^2)*t6          -t6          -t6          -t6]
    [     (zeta8)|    (-zeta8)     (-zeta8)      (zeta8)      (zeta8)   (-zeta8^3)]
    The action on the orbit [0, 9] is as follows
    [            0|            1             2             3             4             5]
    [-------------+---------------------------------------------------------------------]
    [           t0|(-zeta8^2)*t9            t0            t0            t0            t0]
    [           t9|           t0            t9            t9            t9 (-zeta8^2)*t9]
    [      (zeta8)|     (-zeta8)      (-zeta8)       (zeta8)       (zeta8)    (-zeta8^3)]
    The action on the orbit [15] is as follows
    [             0|             1              2              3              4              5]
    [--------------+--------------------------------------------------------------------------]
    [           t15|(-zeta8^2)*t15  (zeta8^2)*t15           -t15           -t15 (-zeta8^2)*t15]
    [       (zeta8)|      (-zeta8)       (-zeta8)        (zeta8)        (zeta8)     (-zeta8^3)]

So the orbit lengths are one and two. We try the length one orbits.::

    sage: R = theta_ring(2,2)[0]
    sage: zeta8 = R.base_ring().gen()    
    sage: u = ThetaModForm(R.gens()[6]/R.gens()[15]); u.orbit(gammas)
    [t6/t15, (-t6)/((-zeta8^2)*t15), t6/(-t15), (-t6)/((zeta8^2)*t15)]

    sage: T = [ThetaModForm(zeta8^j*(R.gens()[6]/R.gens()[15])^k) for j in range(4) for k in range(1,5)]
    sage: [t for t in T if t.is_fixed_by(gammas+[c])] # long time
    [t6^4/t15^4]

So we get a fourth power only, which does not sound very appealing. We compute
class polynomials anyway.::

    sage: i = u^4
    sage: shrecip = [Z.Shimura_reciprocity(A.ideal(), n=8, m=1, transformation=True, period_matrix=True) for A in Z.CM_type().reflex_field().class_group()]
    sage: d = len(shrecip)
    sage: funcs = [[i**u for (U,M,u) in shrecip]] + [[j for k in range(d)] for j in I]
    sage: rts = [[funcs[j][k](shrecip[k][0], prec=1000) for k in range(d)] for j in range(4)] # long time
    sage: X = PolynomialRing(rts[0][0].parent(), 'X').gen() # long time (not really, the only thing that takes time is the previous line)
    sage: minpols_num = [prod([X-r for r in rts[k]]) for k in [0,1]] # long time
    sage: minpols = [recognize_polynomial(p, Kr0, N=None, emb=None) for p in minpols_num] # long time
    sage: intpols_num = [[short_interpolation(rts[k], rts[j]) for j in range(k+1,4)] for k in [0,1]] # long time
    sage: intpols = [[recognize_polynomial(p, Kr0, N=None, emb=None) for p in s] for s in intpols_num] # long time
    sage: without_clinv = [minpols[1]] + intpols[1]; without_clinv # long time
    [y^4 + (-184495806827624430836280/1369*a + 488130022800109977252840/1369)*y^3 + (-739840169612916553137222630166664400/1874161*a + 1957433098731623057993140146384945000/1874161)*y^2 + (593838754022306022349746857388701866092000/1874161*a - 1571149662015479076302089785396028286292000/1874161)*y + 256630599890114720811619355839188289194840000/1874161*a - 678980746118563319720403133738116657470190000/1874161, (-270221655434999771119223055600/91833889*a + 714939299145193541163344096400/91833889)*y^3 + (-3727191323824921637165372056452606086820000/125720594041*a + 1408745904514050491219230782857892353304000/17960084863)*y^2 + (334216521381289621432992250952777444793660600000/17960084863*a - 884253799623993790114727254619647998017361480000/17960084863)*y + 22673920090914327025893245880450680647891945200000/2565726409*a - 59989553807510358276014064831271118637809400000000/2565726409, (-35321155883012350910345181934457626254613984500000/91833889*a + 93450994485796704365921195510596579651462382700000/91833889)*y^3 + (-141640129436536150758065554738066084110751954721920700997000000/125720594041*a + 53534936879439118823725249011262523955822295914798533958000000/17960084863)*y^2 + (16241229513156186987673214123175291768619277258202566137305570000000/17960084863*a - 42970254277733904075124574567929009538309041327433138746030350000000/17960084863)*y + 1002676338528120243794307573904061269871520624535404681015965100000000/2565726409*a - 2652832237234217402171958143339240543598075536513690929627587600000000/2565726409]
    sage: with_clinv = [minpols[0]] + intpols[0]; with_clinv # long time
    [y^4 + (72*a + 44)*y^3 + (2584*a + 6446)*y^2 + (-31208*a - 82532)*y + 103880*a + 274841, (184495806827624430836280/1369*a - 488130022800109977252840/1369)*y^3 + (-27029176139134935951408360/1369*a + 71512478207111453319866520/1369)*y^2 + (-71829208614459623075013240/1369*a + 190042222864442073538449960/1369)*y + 7388424786263310283746600/1369*a - 19547934564957035488314840/1369, (-270221655434999771119223055600/91833889*a + 714939299145193541163344096400/91833889)*y^3 + (39588264074431100940626955279600/91833889*a - 104740701577697557750895005726800/91833889)*y^2 + (105204600549025167276809030084400/91833889*a - 39763601404659714394158812593200/13119127)*y - 10821451235753605818217228858800/91833889*a + 28630868794703343139798634374800/91833889, (-35321155883012350910345181934457626254613984500000/91833889*a + 93450994485796704365921195510596579651462382700000/91833889)*y^3 + (5174652802249064760968244928094007370259990829300000/91833889*a - 13690844435854520719993743433514496879702946983900000/91833889)*y^2 + (13751481500094966173119928209618461888179517492900000/91833889*a - 5197571458279531383913697384533310446636633505700000/13119127)*y - 1414491245594569190315248249844039393557305072900000/91833889*a + 3742392067521217234569095755708719599036994549900000/91833889]
    sage: [(37^6*p).denominator() for p in minpols] # long time
    [1, 1]
    sage: [[(7^2*37^6*p).denominator() for p in s] for s in intpols] # long time
    [[1, 1, 1], [1, 1]]

So even this class invariant gives an improvement. We could also use the
orbits of length two. There are so many of them, we have a 5-dimensional
lattice P of functions whose 8th powers are class invariants, and some
sublattice C of actual class invariants. This sublattice  C is the fixed
sublattice of a group of order 64::

    sage: len(group_generators_to_list(gammas))
    64

So we can compute this sublattice, and put an L1 norm on it, given by the
height. Then the shortest vector is a good class invariant.

EXAMPLE:

The three fields above were arbitrary fields, but it happened to be so that
the prime two was never inert. And the inert case is the case that is
classically not used (as far as I know) for class invariants. Let's be bold and
try it anyway.::

    sage: from recip import *
    sage: K = CM_Field([5, 17, 61])
    sage: Z = K.one_period_matrix('Phi')
    sage: I = igusa_invariants_absolute()
    sage: Kr0 = QuadraticField(61,'a') # the same 61 as a few lines above
    sage: gammas = reciprocity_map_image(Z, 8)
    sage: c = Z.complex_conjugation_symplectic_matrix(8)
    sage: visualize(gammas + [c], 2) # output is random
    On the 8-th powers, the action is (note that 16 means 0):
    1: (1,8,4,9,3)(2,6,16,15,12)
    2: ()
    3: ()
    4: ()
    5: (1,3)(2,16)(8,9)(12,15)
    Permutation Group with generators [(), (1,3)(2,16)(8,9)(12,15), (1,8,4,9,3)(2,6,16,15,12)] of order 10
    The action on the orbit [1, 3, 4, 8, 9] is as follows
    [            0|            1             2             3             4             5]
    [-------------+---------------------------------------------------------------------]
    [           t1|  (-zeta8)*t8            t1  (zeta8^2)*t1           -t1            t3]
    [           t3|  (-zeta8)*t1            t3  (zeta8^2)*t3           -t3            t1]
    [           t4|           t9            t4  (zeta8^2)*t4           -t4            t4]
    [           t8|   (zeta8)*t4            t8  (zeta8^2)*t8           -t8 (-zeta8^3)*t9]
    [           t9| (zeta8^2)*t3           -t9           -t9            t9 (-zeta8^3)*t8]
    [      (zeta8)|     (-zeta8)      (-zeta8)     (zeta8^3)      (-zeta8)    (-zeta8^3)]
    The action on the orbit [0, 2, 6, 12, 15] is as follows
    [             0|             1              2              3              4              5]
    [--------------+--------------------------------------------------------------------------]
    [            t0| (zeta8^2)*t15             t0             t0             t0             t2]
    [            t2|           -t6             t2            -t2             t2             t0]
    [            t6|            t0            -t6  (-zeta8^2)*t6            -t6             t6]
    [           t12|    (zeta8)*t2            t12           -t12            t12    (zeta8)*t15]
    [           t15| (zeta8^2)*t12           -t15 (-zeta8^2)*t15           -t15    (zeta8)*t12]
    [       (zeta8)|      (-zeta8)       (-zeta8)      (zeta8^3)       (-zeta8)     (-zeta8^3)]

The orbits of length 5 are not surprising, given that the kernel of the norm
map `F_{16}^*\rightarrow F_{4}^*` (in which mu lies) has order 5. That they
remain of length 5 even when complex conjugation is taken into account is
useful, as now we can take the quotient of two orbits.::

    sage: R = theta_ring(2,2)[0]
    sage: zeta8 = R.base_ring().gen()    
    sage: e = [-1, 1, -1, 1, 1, 0, -1, 0, 1, 1, 0, 0, -1, 0, 0, -1]
    sage: f = prod([R.gens()[k]^e[k] for k in range(16)])
    sage: T = [ThetaModForm(zeta8^j*f^k) for j in range(4) for k in range(1,5)]
    sage: [t for t in T if t.is_fixed_by(gammas+[c])] # long time
    [t1^2*t3^2*t4^2*t8^2*t9^2/(t0^2*t2^2*t6^2*t12^2*t15^2), t1^4*t3^4*t4^4*t8^4*t9^4/(t0^4*t2^4*t6^4*t12^4*t15^4)]
    sage: i = ThetaModForm(f^2)

Does not look like the most promising class invariant we have so far, but let's
see. It will not be much work anyway::

    sage: Z.CM_type().reflex_field().class_number()
    1
    sage: funcs = [i] + I
    sage: vals = [funcs[j](Z, prec=300) for j in range(4)] # long time

Now to recognize them as elements of Kr0, i.e., constant polynomials over that
field.::

    sage: P = PolynomialRing(vals[0].parent(), 'X') # long time (not really, the only thing that takes time is the previous line)
    sage: vals = [P(v) for v in vals] # long time
    sage: [recognize_polynomial(p, Kr0, N=None, emb=None) for p in vals] # long time
    [-1/2*a - 5/2, -1232*a + 40208, -65536*a + 323584, 19038208*a + 5531598848]
    
So this is again a very significant improvement in height, though class number
one may not convince anyone. So the next example repeats this trick with a
larger class number.

EXAMPLE:

This example is similar to the previous one, but has a larger class number.::

    sage: l = list(iterate_DAB(30, 30))
    sage: l = [(DAB, CM_Field(DAB)) for DAB in l]
    sage: l = [(DAB, K, K.ideal(2).factor()) for (DAB, K) in l]
    sage: l = [(DAB, K) for (DAB, K, f) in l if len(f) == 1 and f[0][1] == 1]
    sage: l = [(DAB, K.class_number()/K.real_field().class_number()) for (DAB, K) in l]
    sage: [(DAB, h) for (DAB, h) in l if h > 3]
    [([21, 25, 109], 4), ([29, 25, 149], 5)]
    sage: K = CM_Field([29, 25, 149])
    sage: Z = K.one_period_matrix('Phi')
    sage: I = igusa_invariants_absolute()
    sage: Kr0 = QuadraticField(149, 'a') # same 149 as three lines above
    sage: gammas = reciprocity_map_image(Z, 8)
    sage: c = Z.complex_conjugation_symplectic_matrix(8)
    sage: visualize(gammas + [c], 2)
    On the 8-th powers, the action is (note that 16 means 0):
    ...
    orbit [0, 1, 8, 12, 15]
    ...
    orbit [2, 3, 4, 6, 9]
    ...
    sage: R = theta_ring(2,2)[0]
    sage: zeta8 = R.base_ring().gen()    
    sage: e = [1, 1, -1, -1, -1, 0, -1, 0, 1, -1, 0, 0, 1, 0, 0, 1]
    sage: f = prod([R.gens()[k]^e[k] for k in range(16)])
    sage: T = [ThetaModForm(zeta8^j*f^k) for j in range(4) for k in range(1,5)]
    sage: [t for t in T if t.is_fixed_by(gammas+[c])] # long time
    [t0^2*t1^2*t8^2*t12^2*t15^2/(t2^2*t3^2*t4^2*t6^2*t9^2), t0^4*t1^4*t8^4*t12^4*t15^4/(t2^4*t3^4*t4^4*t6^4*t9^4)]
    sage: i = ThetaModForm(f^2)
    sage: shrecip = [Z.Shimura_reciprocity(A.ideal(), n=8, m=1, transformation=True, period_matrix=True) for A in Z.CM_type().reflex_field().class_group()]
    sage: d = len(shrecip)
    sage: funcs = [[i**u for (U,M,u) in shrecip]] + [[j for k in range(d)] for j in I]
    sage: rts = [[funcs[j][k](shrecip[k][0], prec=1000) for k in range(d)] for j in range(4)] # long time
    sage: X = PolynomialRing(rts[0][0].parent(), 'X').gen() # long time (not really, the only thing that takes time is the previous line)
    sage: minpols_num = [prod([X-r for r in rts[k]]) for k in [0,1]] # long time
    sage: minpols = [recognize_polynomial(p, Kr0, N=None, emb=None) for p in minpols_num] # long time
    sage: intpols_num = [[short_interpolation(rts[k], rts[j]) for j in range(k+1,4)] for k in [0,1]] # long time
    sage: intpols = [[recognize_polynomial(p, Kr0, N=None, emb=None) for p in s] for s in intpols_num] # long time
    sage: without_clinv = [minpols[1]] + intpols[1]; without_clinv # long time
    [y^5 + (-2418282952857233369424/17183321251*a + 29519365977317103700176/17183321251)*y^4 + (18657571968827814174248448/17183321251*a - 225676698228500892126561792/17183321251)*y^3 + (-1772614250943575113773715636224/17183321251*a + 7004816258842691941289549414400/17183321251)*y^2 + (6965625693308274607389206861316096/17183321251*a + 91862218076978319776745999932325888/17183321251)*y - 928521990810347943418484426183616233472/17183321251*a - 11333618011834127247599451076308548714496/17183321251, (67843696939489141596895950537494814720/235459998456026925009842449*a - 828207718037464759587741441013340700672/235459998456026925009842449)*y^4 + (-12257312691937195084721183229684104498774016/235459998456026925009842449*a + 148581254610940719416097733527420783734685696/235459998456026925009842449)*y^3 + (-304967001050729084816324463146683300169933389824/235459998456026925009842449*a + 4110452058255955280174707428489140970398718885888/235459998456026925009842449)*y^2 + (-8432023593406349586645449420190347006046711192748032/235459998456026925009842449*a - 102307072679418149121133287754598365849133118759370752/235459998456026925009842449)*y + 110919181216682781503633791589528805897246728327647985664/235459998456026925009842449*a + 1353920838328456722297644139832071467949603951611574484992/235459998456026925009842449, (102663533639339234156654511318396526323989643264/106591217046639622005361*a - 1253168132873027139152980151877076335103743787008/106591217046639622005361)*y^4 + (-788468755675957676359334818168981567390513562124288/106591217046639622005361*a + 9624490169348234463004722185505927893482045282713600/106591217046639622005361)*y^3 + (49806658346906063464978622851939772703084171975318306816/106591217046639622005361*a - 607977589824269924761626249311396725726253950565401755648/106591217046639622005361)*y^2 + (11892116465623212227219255052847525175123656411025718116352/106591217046639622005361*a - 145041760857136573554966147904668265057791140629466427949056/106591217046639622005361)*y + 137170910741158730621452073998807343840543596150445335117824/106591217046639622005361*a - 16890597616368525196940772118962655099105914652834138630914048/106591217046639622005361]
    sage: with_clinv = [minpols[0]] + intpols[0]; with_clinv # long time
    [y^5 + (83075497337/2686562*a - 1014014112993/2686562)*y^4 + (155082250987/1343281*a - 1892895942857/1343281)*y^3 + (254067967826/1343281*a - 3101290832596/1343281)*y^2 + (90778603497/2686562*a - 1108093940637/2686562)*y + 151403184505/2686562*a - 1848111391279/2686562, (2418282952857233369424/17183321251*a - 29519365977317103700176/17183321251)*y^4 + (2813376427149574723392/5352182029*a - 34344388705129573639104/5352182029)*y^3 + (281063334391599905888832/326483103769*a - 3431622775810662832295808/326483103769)*y^2 + (49757608975106394502848/326483103769*a - 607385797316372031731136/326483103769)*y + 1349249041139772592368/5352182029*a - 16469684540226527710608/5352182029, (67843696939489141596895950537494814720/235459998456026925009842449*a - 828207718037464759587741441013340700672/235459998456026925009842449)*y^4 + (255818617576482504504098928281040027648/235459998456026925009842449*a - 3124006514663560468514058645007974825984/235459998456026925009842449)*y^3 + (426630644608252179688506017417084018688/235459998456026925009842449*a - 5214155527658281710830445005988032077824/235459998456026925009842449)*y^2 + (97631372488854352146543063142996082688/235459998456026925009842449*a - 1191877658411774755275584523540831043584/235459998456026925009842449)*y + 136660950448855546231471256489227014144/235459998456026925009842449*a - 1668160417416832978803397540365875564544/235459998456026925009842449, (102663533639339234156654511318396526323989643264/106591217046639622005361*a - 1253168132873027139152980151877076335103743787008/106591217046639622005361)*y^4 + (383293048180275507494037346322869656147361923072/106591217046639622005361*a - 4678687905836451780642635647866873199548037201920/106591217046639622005361)*y^3 + (627959867612408002969363243533912652017668063232/106591217046639622005361*a - 7665227029706998584461581144192425449695548014592/106591217046639622005361)*y^2 + (112185462741310392997036458228281628247338909696/106591217046639622005361*a - 1369398089836923540612806919781494379528780382208/106591217046639622005361)*y + 187105762277224037380424840603213297961640820736/106591217046639622005361*a - 2283916893259042203585521764488246375504219832320/106591217046639622005361]
    sage: [(2*17^2*19^2*29^2*61^2*p).denominator() for p in minpols] # long time
    [1, 1]
    sage: [[(17^4*19^4*29^4*47^2*61^4*p).denominator() for p in s] for s in intpols] # long time
    [[1, 1, 1], [1, 1]]

So the total output size increases a bit (in terms of the number of characters)
with this class invariant, but let's look at the largest coefficient, which
determines the precision for the theta constants::

    sage: h1 = max([max([a.global_height() for a in p]) for p in with_clinv]) # long time
    sage: h1 # long time
    76.4813952507977
    sage: h2 = max([max([a.global_height() for a in p]) for p in without_clinv]) # long time
    sage: h2 # long time
    114.465280863815
    sage: h1/h2 # long time
    0.66816...

So we get an improvement in size of over 30%, even in the inert case!  
"""

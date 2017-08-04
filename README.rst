=========================
REpository of Complex multiplication SageMath code
=========================
.. image:: https://travis-ci.org/mstreng/recip.svg?branch=master
    :target: https://travis-ci.org/mstreng/recip


Once the Travis CI set-up has been completed, the documentation for the package can be found at https://mstreng.github.io/recip/doc/html/

Installation
------------

Local install from source
^^^^^^^^^^^^^^^^^^^^^^^^^

Download the source from the git repository::

    $ git clone https://github.com/mstreng/recip.git

Change to the root directory and run::

    $ sage -pip install --upgrade --no-index -v .

For convenience this package contains a [makefile](makefile) with this
and other often used commands. Should you wish too, you can use the
shorthand::

    $ make install

Usage
-----

Once the package is installed, you can use it in Sage with::

    sage: from recip import *
    sage: CM_Field([5,5,5])
    CM Number Field in alpha with defining polynomial x^4 + 5*x^2 + 5


Source code
-----------

All source code is stored in the folder ``recip`` using the same name as the
package. This is not mandatory but highly recommended for clarity. All source folder
must contain a ``__init__.py`` file with needed includes.

Tests
-----

This package is configured for tests written in the documentation
strings, also known as ``doctests``. For examples, see this
[source file](recip/ultimate_question.py). See also
[SageMath's coding conventions and best practices document](http://doc.sagemath.org/html/en/developer/coding_basics.html#writing-testable-examples).
With additional configuration, it would be possible to include unit
tests as well.

Once the package is installed, one can use the SageMath test system
configured in ``setup.py`` to run the tests::

    $ sage setup.py test

This is just calling ``sage -t`` with appropriate flags.

Shorthand::

    $ make test

Documentation
-------------

The documentation of the package can be generated using Sage's
``Sphinx`` installation::

    $ cd docs
    $ sage -sh -c "make html"

Shorthand::

    $ make doc


#*****************************************************************************
# Copyright (C) 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017 Marco Streng
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

RECIP -- REpository of Complex multIPlication SageMath code.

This started out as code meant for computing with Shimura's RECIProcity law,
but grew into a collection of much of the SageMath code written by me for my
research.

Version 0.2.5 (tested with SageMath 8.0, short tests only).

When using this package in a publication, it is highly likely that it is
approprate certain publications. Please CITE the relevant JOURNAL publications,
as well as giving the URL of this repository.

Here is a list of functionalities of this repository, together with the
publications that should be cited when you use them, and the name of the file
that has examples.

 * Igusa class polynomials (proven correct)
   See both "Igusa class polynomials (not proven correct)" and
   "Denominators of Igusa class polynomials" below.

 * Non-maximal orders of CM-fields and their polarized ideal classes and Igusa
   class polynomials.
   cite [BissonStreng] (code is written for, part of, and based on, this publication)
   see orders.sage for examples

 * (n,n)-isogenies between polarized ideal classes
   cite [BLS]
   see bls.sage for examples

 * Computations related to Shimura's reciprocity law
   cite [Streng12] (code is written for, part of, and based on, this publication)
   see article.sage for examples

 * Igusa class polynomials (not proven correct)
   cite [Streng14], [vWamelen], [Weng] (code is based on these publications)

 * Denominators of Igusa class polynomials
   cite [BouyerStreng] (code is written for, and hence part of, this publication)
   and depending on how the code is used, and on the kind of quartic CM-field,
   also cite one or more of:
   [BouyerStreng], [GL], [LV], [Yang] (large parts of the code are based on these)
   see denominators.sage for examples

Here is a list of SageMath programs written by my students and me that is not part
of this repository.

 * Height reduction of binary forms and hyperelliptic curves.
   (with Florian Bouyer)
   https://bitbucket.org/mstreng/reduce
   cite [BouyerS] (code is written for, part of, and based on, this publication)

 * Solving conics and Mestre's algorithm
   (with Florian Bouyer)
   now part of the standard SageMath functionality

 * Hilbert modular polynomials
   (by Chloe Martindale)
   contact her if you are interested

 * CM class number one for genus 2 and 3
   (by Pınar Kılıçer)
   contact her if you are interested

To use the latest version of this package directly from the web, start SageMath
and type::

    sage: load("https://bitbucket.org/mstreng/recip/raw/master/recip_online.sage")

To use this package offline, download it first and extract it to some
directory, say "somewhere_on_my_drive/recip", then start SageMath and type::

    sage: load_attach_path("somewhere_on_my_drive/recip")
    sage: load("recip.sage")

[ABLPV]  -  Comparing arithmetic intersection formulas for denominators of
            Igusa class polynomials -- Jacqueline Anderson, Jennifer S.
            Balakrishnan, Kristin Lauter, Jennifer Park, and Bianca Viray

[BissonS] - On polarised class groups of orders in quartic CM fields --
            Gaetan Bisson and Marco Streng
            http://arxiv.org/abs/1302.3756

[BLS]    -  Abelian surfaces admitting an (l,l)-endomorphism -- Reinier Broker,
            Kristin Lauter, and Marco Streng
            accepted for publication by Journal of Algebra, 2013
            http://arxiv.org/abs/1106.1884

[BouyerS] - Examples of CM curves of genus 2 defined over the reflex field --
            Florian Bouyer and Marco Streng
            http://arxiv.org/abs/1307.0486

[GJLSVW] -  Igusa class polynomials, embeddings of quartic CM fields, and
            arithmetic intersection theory -- Helen Grundman, Jennifer
            Johnson-Leung, Kristin Lauter, Adriana Salerno, Bianca Viray, and
            Erika Wittenborn
            http://arxiv.org/abs/1006.0208

[GL]     -  Genus 2 curves with complex multiplication -- Eyal Goren and
            Kristin Lauter

[LV]     -  An arithmetic intersection formula for denominators of Igusa class
            polynomials -- Kristin Lauter and Bianca Viray
            arXiv:1210.7841v1

[Yang]   -  Arithmetic interseciton on a Hilbert modular surface and the
            Faltings height -- Tonghai Yang
            http://www.math.wisc.edu/~thyang/general4L.pdf

[recip]  -  recip, SageMath package for explicit complex multiplication -- Marco
            Streng
            https://bitbucket.org/mstreng/recip/

[Streng12]-  An explicit version of Shimura's reciprocity law for Siegel
            modular functions -- Marco Streng
            arXiv:1201.0020

[Streng14]-  Computing Igusa Class Polynomials
            Mathematics of Computation, Vol. 83 (2014), pp 275--309


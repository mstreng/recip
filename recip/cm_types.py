"""
Recip -- Sage package for using Shimura's reciprocity law

See the file README.txt for version information and instructions.

#*****************************************************************************
# Copyright (C) 2010, 2011, 2012, 2013 Marco Streng <marco.streng@gmail.com>
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

from os import path

from __init__ import *

directory = "/".join(path.abspath(__file__).split("/")[:-1]) + "/"

from sage.all import *
from sage.structure.sage_object import load

files = [   "cm_types.sage"]

for f in files:
    load(directory+f)


def add_one(a):
    """
    This is a funtion that does not
    """
    return a+1

from os import path 

from __init__ import *

directory = "/".join(path.abspath(__file__).split("/")[:-1]) + "/"

from sage.all import *

from sage.structure.sage_object import load
files = ["theta_trans.sage"]
for f in files:
    load(directory + f)

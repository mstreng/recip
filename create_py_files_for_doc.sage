





files = [   ("basic", "Basic"),
            ("polynomials", "Polynomials"),
            ("symplectic_matrices", "Symplectic matrices"),
            ("period_matrices", "Period matrices"),
            ("theta_trans", "Theta transformation formulas"),
            ("cm_types", "CM types and CM fields"),
            ("find_class_invariant", "Finding class invariants"),
            ("igusa_invariants", "Igusa invariants"),
            ("list_fields", "List of fields"),
            ("denominators", "Denominators"),
            ("class_polynomials", "Class polynomials"),
            ("orders", "Orders"),
            ("bissonstreng", "Bisson-Streng"),
            ("bls", "Broker-Lauter-Streng")]

print __file__

for (f, t) in files:
    os.chdir("recip")
    a = open(f + ".py", 'w')
    a.write("from os import path \n")
    a.write("\n")
    a.write("from __init__ import *\n")
    a.write("\n")
    a.write('directory = "/".join(path.abspath(__file__).split("/")[:-1]) + "/"\n')
    a.write("\n")
    a.write("from sage.all import *\n\n")
    a.write("from sage.structure.sage_object import load")
    a.write("\n")
    a.write('files = ["' + f + '.sage"]\n')
    a.write("for f in files:\n")
    a.write("    load(directory + f)\n")
    a.close()

    os.chdir("../docs/source")
    b = open(f + ".rst", 'w')
    b.write("..nodoctest\n\n")
    b.write(t + "\n")
    b.write("===============================\n\n")
    b.write(".. automodule:: recip." + f + "\n")
    b.write("   :members:\n")
    b.write("   :undoc-members:\n")
    b.write("   :show-inheritance:\n")
    b.close()
    os.chdir("../..")



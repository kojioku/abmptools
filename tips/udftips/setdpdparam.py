import UDFManager
import sys


def getatomtype(uobj):
    atomtype = uobj.get("Molecular_Attributes.Atom_Type[].Name")
    return atomtype


def getbondtype(uobj):
    bondtype = uobj.get("Molecular_Attributes.Bond_Potential[].Name")
    return bondtype


def getangletype(uobj):
    angletype = uobj.get("Molecular_Attributes.Angle_Potential[].Name")
    return angletype


def putbondtype(uobj):
    '''
    This function sets the bond potential type to Harmonic and sets R0 value
    '''

    print("Setting bond potential type to Harmonic")
    bondname = uobj.get("Molecular_Attributes.Bond_Potential[].Name")
    print(len(bondname))

    for i, bondname in enumerate(bondname):
        print("set bond name:", bondname, "type: Harmonic", "R0:", 0.86)
        uobj.put("Harmonic", "Molecular_Attributes.Bond_Potential[].Potential_Type", [i])
        uobj.put(0.86, "Molecular_Attributes.Bond_Potential[].R0", [i])
    return uobj


def putangletype(uobj):
    '''
    This function sets the angle potential type to Theta and sets theta0 and K values
    '''

    print("Setting angle potential type to Theta")
    anglename = (uobj.get("Molecular_Attributes.Angle_Potential[].Name"))
    print(len(anglename))

    for i, anglename in enumerate(anglename):
        print("set angle name:", anglename, "type: Theta", "theta0:", 0, "K:", 4.0)
        uobj.put("Theta", "Molecular_Attributes.Angle_Potential[].Potential_Type", [i])
        uobj.put(0, "Molecular_Attributes.Angle_Potential[].theta0", [i])
        uobj.put(4.0, "Molecular_Attributes.Angle_Potential[].Theta.K", [i])
    return uobj


if __name__ == '__main__':
    if len(sys.argv) != 4:
        print("Usage: python dpdinfo.py <input.udf> <output.udf> <aij.py>")
        sys.exit(1)
    uobj = UDFManager.UDFManager(sys.argv[1])
    atomtype = getatomtype(uobj)
    print(atomtype)
    bondtype = getbondtype(uobj)
    print(bondtype)
    angletype = getangletype(uobj)
    print(angletype)

    # Set bond and angle potential types
    putbondtype(uobj)

    # Set angle potential types
    putangletype(uobj)


    # Load aij data using exec
    print("Loading aij data")
    aij = []
    with open(sys.argv[3], 'r') as f:
        exec(f.read(), globals())

    # Iterate over Interactions.Pair_Interaction
    pair_names = uobj.get("Interactions.Pair_Interaction[].Name")
    for i, pair_name in enumerate(pair_names):
        for row in aij:
            key1 = f"{row[0]}-{row[1]}"  # A-B
            key2 = f"{row[1]}-{row[0]}"  # B-A
            if pair_name == key1 or pair_name == key2:
                print(f"Found {pair_name} and values are {row[2]}")
                uobj.put(row[2], "Interactions.Pair_Interaction[].DPD.a", [i])
                break

    print("Writing to output file:", sys.argv[2])
    uobj.write(sys.argv[2], define=0)
    print("Done")

    # Generate all combinations of atomtype joined by "-"
#     from itertools import permutations
#     atom_combinations = ["-".join(p) for p in permutations(atomtype, 2)]
# 
#     # Compare combinations with bondtype
#     included = [combo for combo in atom_combinations if combo in bondtype]
#     excluded = [combo for combo in atom_combinations if combo not in bondtype]
# 
#     print("Included in Bond Types:", included)
#     print("Excluded from Bond Types:", excluded)


    # udf.put("Molecular_Attributes.Bond_Potential[].Name")
    # udf.put("Harmonic","Molecular_Attributes.Bond_Potential[].Potential_Type", [i])
    # udf.put(0.86,"Molecular_Attributes.Bond_Potential[].R0,", [i])

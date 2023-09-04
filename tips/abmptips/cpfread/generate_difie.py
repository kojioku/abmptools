import abmptools.cpfmanager

import argparse

def get_args():
    parser = argparse.ArgumentParser(description="generate DIFIE")

    # Input CPF name, required
    parser.add_argument("-i", "--input",
                        required=True,
                        help="Input CPF name (e.g., file-xxx-bbb.cpf)")

    # Start, End, Interval, required
    parser.add_argument("-t", "--time",
                        required=True,
                        nargs=3,
                        type=int,
                        help="Start, End, Interval numbers")

    # Zero padding, default is "1"
    parser.add_argument("-z",
                        "--zero-padding",
                        default=1,
                        type=int,
                        help="Number of digits for zero-padding. Default is '1' (no padding).")

    # Representative structure, default is first structure (Start number)
    parser.add_argument("-s",
                        "--structure",
                        type=int,
                        help="Number of the representative structure for average CPF. Default is the start number.")

    # Target residues, default is removing WAT, HOH at the end of the file
    parser.add_argument("-f",
                        "--fragments",
                        default=None,
                        help="Specify target residues. Default is None")

    # Version, default is "23"
    parser.add_argument("-v",
                        "--version",
                        default=23,
                        help="Specify the output version. Currently, only '23' is available.")

    args = parser.parse_args()

    print("Arguments received:")
    print("Input CPF:", args.input)
    print("Time:", args.time)
    print("Zero-padding:", args.zero_padding)
    print("Representative structure:", args.structure)
    print("Target fragments:", args.fragments)
    print("Version:", args.version)

    return args

if __name__ == "__main__":
    args = get_args()
    cpf = abmptools.CPFManager()
    for i in range(args.time[0], args.time[1]+1, args.time[2]):
        padded = str(i).zfill(args.zero_padding)
        input_cpf = args.input.replace('xxx', padded)
        output_cpf = args.input.replace('xxx', padded + 'out')
        cpf.tgtfrag = cpf.selectfrag(args.fragments)
        cpf = cpf.parse(input_cpf)
        cpf.write('test', output_cpf)

    # cpf = abmptools.CPFManager()
    # cpf.parse('CS4_ligX_md1_12ns_mod_mp2_631gd.cpf')
    # cpf.write('test-from10to23', 'test-i10o23.cpf')


#         cpf.parse('gly5-10.cpf')
#         cpf.write('test-from10to23', 'test-i10o23.cpf')

# cpf10.write('test-from10to23', 'test-i10o23-trim.cpf')

# cpf4 = CPFManager()
# cpf4.parse('gly5-4201.cpf')
# cpf4.write('test-from42to23', 'test-i42o23.cpf')

# print(dir(cpf))
# print(vars(cpf))
# cpf.write('test', 'test.cpf')
# print('cpfver\n', cpf.cpfver)
# print('labels\n', cpf.labels)
# print('atominfo\n', cpf.atominfo.head())
# print('fraginfo\n', cpf.fraginfo)
# print('condition\n', cpf.condition)
# print('static_data\n', cpf.static_data)
# print('mominfo\n', cpf.mominfo)
# print('diminfo\n', cpf.diminfo)
# print('labels\n', cpf.labels)

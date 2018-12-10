#!/usr/env/py36

import argparse
from collections import OrderedDict

# Mutpos: chrom,
#         template,
#         pos,
#         depths,
#         muts,
#         Tcount, Ccount, Gcount, Acount,
#         inscount, delcount, Ncount


def combine(inputs, output):
    combined = OrderedDict()
    for inp in inputs:
        with open(inp, 'r') as fi:
            while True:
                line = fi.readline().split("\t")
                if line == [""]: break

                chrom,ref,pos,depth,mut,Ts,Cs,Gs,As,ins,dels,Ns = line
                if chrom == "Chrom": continue
                new_mp = [chrom,ref,int(depth),int(mut),int(Ts),int(Cs),int(Gs),int(As),int(ins),int(dels),int(Ns)]
                pos = int(pos)
                if pos not in combined.keys():
                    combined[pos] = [new_mp[0],new_mp[1]]+[0]*9
                for i in range(2,len(combined[pos])):
                    combined[pos][i] = combined[pos][i] + new_mp[i]

    combined = OrderedDict(sorted(combined.items(), key=lambda x:x[0]))
    output = output + ".mutpos"
    with open(output, 'w') as fo:
        for key in list(combined.keys()):
            l_out = [str(x) for x in combined[key]]
            l_out.insert(2,str(key))
            fo.write('\t'.join(l_out) + '\n')

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--infile',
                        action ='store',
                        dest = 'inputs',
                        help = 'A comma-delimited list of input files.',
                        required = True)
    parser.add_argument('-o', '--outfile',
                        action = 'store',
                        dest = 'output',
                        help = 'A filename for the output file.',
                        default = 'outfile',
                        required = False)
    args = parser.parse_args()
    inputs = args.inputs.split(',')
    output = args.output

    combine(inputs, output)

if __name__ == "__main__":
    main()

"""
UCM Bash Maker for DCS-3.00

Write a bash script to run the process.
"""

import sys
import os
import re
import argparse
from argparse import ArgumentParser
import time

def main():
    parser = ArgumentParser()
    parser.add_argument("--run_name",
            dest = "run_name",
            help = "An identifier for this particular processing run.",
            required=True
            )
    parser.add_argument("--runIdentifier",
            dest = "runID",
            help = "An identifier for this particular sequencing run.",
            required=True
            )
    parser.add_argument("--sampleIdentifier",
            dest = "samID",
            help = "An identifier for this particular sample.",
            required=True
            )
    parser.add_argument("--inp_path",
            dest = "inp_path",
            help = "Path for the directory of the input FASTQ files.",
            required=True
            )
    parser.add_argument("--endID",
            action = "store",
            dest = "endID",
            help = "Tail ID of the run."
            )
    parser.add_argument("--ref",
            action = "store",
            dest = "ref",
            help = "Path to reference directory plus file for reference genome."
            )
    parser.add_argument("--picard",
            dest = "picard",
            help = "Reference path to picard-tools-1.119 files."
            )
    parser.add_argument("--GATK",
            dest = "GATK",
            help = "Reference path to GATK files."
            )
    parser.add_argument("--DS_path",
            dest = "DS_path",
            help = "Reference path to DCS-3.00 scripts."
            )
    parser.add_argument("--min",
            action = "store",
            dest = "minMem",
            help = "Minimum members for SSCS consensus"
            )
    parser.add_argument("--max",
            action = "store",
            dest = "maxMem",
            help = "Maximum members for SSCS consensus"
            )
    parser.add_argument("--cut",
            action = "store",
            dest = "cutOff",
            help = "Mimimum percent matching for base choice in SSCS consensus "
            )
    parser.add_argument("--Ncut",
            action = "store",
            dest = "Ncut",
            help = "Maximum percent N's allowed"
            )
    parser.add_argument("--rlength",
            action = "store",
            dest = "rlength",
            help = "Length of a single read",
            )
    parser.add_argument("--blength",
            action = "store",
            dest = "blength",
            help = "Length of the barcode sequence on a unprocessed single read. "
            )
    parser.add_argument("--slength",
            action = "store",
            dest = "slength",
            help = "Length of the spacer sequence in a unprocessed single read."
            )
    parser.add_argument("--repFilt",
            action = "store",
            dest = "repFilt",
            help = "Remove tags with homomeric runs of nucleotides of length x."
            )
    parser.add_argument("--softclip_cycles",
            dest = "softclip_cycles",
            help = "Cycles to clip in processing."
            )
    parser.add_argument("--template",
            action = "store",
            dest="template",
            help="Template to use with bash maker.  If not specified, defaults to ucm_bash_template.sh."
            )
    o = parser.parse_args()
    outBash = open(o.run_name + '.script.sh', 'w')

    spath = repr(os.path.realpath(__file__)).replace("'", "").replace("/UCM_BASH_MAKER.py", "")
    if o.template:
        inBash = open(o.template, "r")
    else:
        inBash = open(spath + "/ucm_bash_template.sh", "r")

    for line in inBash:
        outBash.write(line)
        if line.strip() == "###----- NONDEFAULTS -----":
            outBash.write("run_name=" + o.run_name + "\n")
            outBash.write("runID=" + o.runID + "\n")
            outBash.write("samID=" + o.samID + "\n")
            if o.endID:
                outBash.write("endID=" + o.endID + "\n")
                outBash.write("inp_path=" + o.inp_path + "/" + o.samID + "-" + o.endID + "\n")
            else:
                outBash.write("inp_path=" + o.inp_path + "/" + o.samID + "\n")
            if o.ref:
                outBash.write("ref=" + o.ref + "\n")
            if o.picard:
                outBash.write("picard=" + o.picard + "\n")
            if o.GATK:
                outBash.write("GATK=" + o.GATK + "\n")
            if o.DS_path:
                outBash.write("DS_path=" + o.DS_path + "\n")
            if o.minMem:
                outBash.write("minMem=" + o.minMem + "\n")
            if o.maxMem:
                outBash.write("maxMem=" + o.maxMem + "\n")
            if o.cutOff:
                outBash.write("cutOff=" + o.cutOff + "\n")
            if o.Ncut:
                outBash.write("nCutOff=" + o.Ncut + "\n")
            if o.rlength:
                outBash.write("readLength=" + o.rlength + "\n")
            if o.blength:
                outBash.write("barcodeLength=" + o.blength + "\n")
            if o.slength:
                outBash.write("spacerLength=" + o.slength + "\n")
            if o.repFilt:
                outBash.write("repFilt=" + o.repFilt + "\n")
            if o.softclip_cycles:
                outBash.write("softclip_cycles=" + o.softclip_cycles + "\n")

    outBash.write("\n\n#Generated on " + time.ctime(time.time()) + " with the command: \n")
    outBash.write("#python ")
    arguments=sys.argv
    for arg in arguments:
        outBash.write("%s " % arg)
    outBash.write("\n")
    outBash.close()
    inBash.close()

    print("Script made. Run with:\n%s.script.sh\n" % o.run_name)

if __name__ == "__main__":
    main()

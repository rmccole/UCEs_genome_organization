#!/usr/bin/env python
"""
Driver script for enrichment/depletion analysis, for use with files corresponding to
human genome build hg19.

Distributed under the following license:

Copyright 2017 Harvard University, Wu Lab

Licensed under the Apache License, Version 2.0 (the "License"); you may not use this file 
except in compliance with the License. You may obtain a copy of the License at

http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software distributed under the 
License is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, 
either express or implied. See the License for the specific language governing permissions 
and limitations under the License.

"""

import os.path
import sys
import Depletion_enrichment.py as ro3
import argparse


def get_args(strInput=None):
    """

    Collect arguments from command-line, or from strInput if given (only used for debugging)
    """
    parser = argparse.ArgumentParser(description="This program allows you to run the randomoverlaps3.py script against "
                                                 "multiple variant files automatically for the given UCE files. You "
                                                 "must pass at least one UCE file to the script to run. The script will"
                                                 " use the appropriate genome spacing files for each type, which must be"
                                                 " in the same directory")
    parser.add_argument('file', type=argparse.FileType('rU'),
                        help="A file containing a list of paths to the files you want to process, separated by "
                             "newlines")
    parser.add_argument('-c', '--cluster', type=int,
                        help="The cluster size (kb)")
    parser.add_argument('-o', '--output', help="Output file for results [WARNING: Will overwrite any file with the "
                                               "same name in the current directory]")
    parser.add_argument('-a', '--all', type=argparse.FileType('rU'),
                        help="A file containing [a]ll UCEs (exonic + intronic + intergenic)")
    parser.add_argument('-e', '--exonic', type=argparse.FileType('rU'),
                        help="A file containing [e]xonic UCEs")
    parser.add_argument('-i', '--intronic', type=argparse.FileType('rU'),
                        help="A file containing [i]ntronic UCEs")
    parser.add_argument('-t', '--intergenic', type=argparse.FileType('rU'),
                        help="A file containing in[t]ergenic UCEs")
    parser.add_argument('-d', '--debug', action='store_true',
                        help="Set logging level of randomoverlaps3.py to debug")
    if strInput:
        return parser.parse_args(strInput.split())
    else:
        return parser.parse_args()


def get_uces(args):
    """
    
    Check which UCE files have been given and validate that the necessary genome spacing files exist in the current 
    directory
    """
    aUCEFiles = []
    if args.all:
        if os.path.isfile('hg19_nonN_1based_clean.txt'):
            allLen = len(list(args.all))
            aUCEFiles.append(('all', args.all.name, 'hg19.genomic.coordinates.nonN.txt', allLen))
        else:
            print "Cannot find appropriate spacing file for {0}, exiting...".format(args.all.name)
            sys.exit(1)
    if args.intergenic:
        if os.path.isfile('hg19_intergenic_1based.txt'):
            interLen = len(list(args.intergenic))
            aUCEFiles.append(('intergenic', args.intergenic.name, 'hg19_intergenic_1based.txt', interLen))
        else:
            print "Cannot find appropriate spacing file for {0}, exiting...".format(args.intergenic.name)
            sys.exit(1)
    if args.intronic:
        if os.path.isfile('hg19_introns_1based.txt'):
            inLen = len(list(args.intronic))
            aUCEFiles.append(('intronic', args.intronic.name, 'hg19_introns_1based.txt', inLen))
        else:
            print "Cannot find appropriate spacing file for {0}, exiting...".format(args.intronic.name)
            sys.exit(1)
    if args.exonic:
        if os.path.isfile('hg19_exons_1based.txt'):
            exLen = len(list(args.exonic))
            aUCEFiles.append(('exonic', args.exonic.name, 'hg19_exons_1based.txt', exLen))
        else:
            print "Cannot find appropriate spacing file for {0}, exiting...".format(args.exonic.name)
            sys.exit(1)
    if len(aUCEFiles) == 0:
        print "Script must be given at least one valid UCE file"
        sys.exit(1)
    return aUCEFiles


def run(inFile, aUCEs, cluster, debug, output):
    """

    Run randomoverlaps3.py on the given file using the given parameters

    inFile    -- A fileobject containing the test set of intervals
    aUCEs     -- A list of UCE files
    cluster   -- The cluster interval size to be passed to randomoverlaps3.py if given
    debug     -- Boolean for whether randomoverlaps3.py logs to debug or not
    output    -- The name of the output file
    """
    filename = os.path.split(inFile)[1]
    print "Running " + filename
    counter = 0  # Initialize counter so header line is printed only once per run
    if debug:
        log = "debug"
    else:
        log = "warning"
    for tup in aUCEs:
        print "running {}".format(tup[0])
        if cluster:
            aStats = ro3.main(ro3.getArgs("-u {0} -g {1} -i {2} -a {3} -c {4}"
                                          " -d {5}".format(tup[1], tup[2], 1000, inFile, cluster, log), False))
        else:
            aStats = ro3.main(ro3.getArgs("-u {0} -g {1} -i {2} -a {3} "
                                          "-d {4}".format(tup[1], tup[2], 1000, inFile, log), False))
        with open(output, 'a+') as fh:
            if counter == 0:
                fh.write("{0}\t{1}\t{2}\t{3}\n".format(filename, tup[0], tup[3], "\t".join(map(str, aStats))))
                counter += 1  # Print variant file name only once per run
            else:
                fh.write("\t{0}\t{1}\t{2}\n".format(tup[0], tup[3], "\t".join(map(str, aStats))))


def main(args):
    aUCEs = get_uces(args)
    aFiles = [line.strip() for line in args.file]
    # Create output file
    if args.output:
        outFile = args.output
    else:
        outFile = 'results.txt'
        # Write header line once
    header = "CNV Set\tUCE subset\telements\tn\tbp\tmean\ts.d.\tmin\tmax\tKSp-value\tKStestResult\tproportion\tp-value\tObs/Exp\tZtestResult\n"
    with open(outFile, 'w') as fh:  # This also erases any previous output
        fh.write(header)
    for inFile in aFiles:
        if not os.path.isfile(inFile):
            sys.stderr.write("Could not find {0}, skipping...\n".format(inFile))
            continue
        run(inFile, aUCEs, args.cluster, args.debug, outFile)
    print "Wrote results to " + outFile


if __name__ == "__main__":
    args = get_args()
    main(args)
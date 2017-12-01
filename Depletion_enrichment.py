#!/usr/bin/env python

"""
Script to accept a set of uces and pick a set of random regions, matched for length,
in a genome space that is accepted. 

Inputs:
uces: 1-based starts, chr, start end, tab-delimited
genome space where picks are allowed: same as uces
chromosome sizes: chr, size, tab-delimited


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

import sys

if sys.version_info < (2, 7):
    raise Exception("Python 2.7+ is required")

import argparse
import logging
import random
import math
import numpy as np
from scipy import stats

global bVerbose
bVerbose = True

class FoundException(Exception): pass


LOGGING_LEVELS = {'critical': logging.CRITICAL,
                  'error': logging.ERROR,
                  'warning': logging.WARNING,
                  'info': logging.INFO,
                  'debug': logging.DEBUG}


def full_overlap(aIntervalA, aIntervalB):
    """ Returns True if interval A falls completely within interval B otherwise returns False"""
    # Check that both inputs are 3-column intervals
    if not len(aIntervalA) == len(aIntervalB) == 3:
        raise Exception("Regions could not be overlapped")
    if aIntervalA[0] == aIntervalB[0]:
        if aIntervalA[1] >= aIntervalB[1]:
            if aIntervalA[2] <= aIntervalB[2]:
                return True
    else:
        return False


def collapse(aList):
    # Initialize variables
    strChr = iStart = iStop = 0
    aOut = []
    for aInterval in aList:
        # Test if an interval has been stored (always past first loop)
        if strChr:
            # Test if next interval is on a different chr OR if start of
            # next interval is larger than stop of previous interval 
            if strChr != aInterval[0] or aInterval[1] > (iStop + 1):
                # Write interval and clear start/stop coordinates
                aOut.append([strChr, iStart, iStop])
                iStart = iStop = 0
        strChr = aInterval[0]
        # Advance to next interval if iStart is empty
        if not iStart:
            iStart = aInterval[1]
            # If next interval overlaps, adjust stop to larger coordinate
        if aInterval[2] > iStop:
            iStop = aInterval[2]
            # Write last line
    aOut.append([strChr, iStart, iStop])
    return aOut


def partial_overlap(aIntervalA, aIntervalB):
    """
    Returns True if interval A overlaps for at least one bp with interval B,
    then returns the amount of overlap in basepairs

    """
    # Check that both inputs are 3-column intervals
    if not len(aIntervalA) == len(aIntervalB) == 3:
        raise Exception("Regions could not be overlapped")
        # Check that both intervals are on the same chromosome, if not return False
    if aIntervalA[0] == aIntervalB[0]:
        # Unpack coordinates into starts and ends
        iIntervalAStart, iIntervalAEnd = aIntervalA[1:]
        iIntervalBStart, iIntervalBEnd = aIntervalB[1:]
        # Check if start coordinate of interval A lies within interval B
        if point_checker(iIntervalAStart, aIntervalB[1:]):
            if iIntervalAEnd <= iIntervalBEnd:
                iOverlap = iIntervalAEnd - iIntervalAStart + 1
                return True, iOverlap
            else:
                iOverlap = iIntervalBEnd - iIntervalAStart + 1
                return True, iOverlap
        # If not, check if end coordinate of interval A lies within interval B
        elif point_checker(iIntervalAEnd, aIntervalB[1:]):
            if iIntervalAStart >= iIntervalBStart:
                iOverlap = iIntervalAEnd - iIntervalAStart + 1
                return True, iOverlap
            else:
                iOverlap = iIntervalAEnd - iIntervalBStart + 1
                return True, iOverlap
        # If not, check if interval A surrounds interval B
        elif iIntervalAStart <= iIntervalBStart and iIntervalAEnd >= iIntervalBEnd:
            iOverlap = iIntervalBEnd - iIntervalBStart + 1
            return True, iOverlap
        # If not, then intervals do not overlap
        else:
            return False, 0
    else:
        return False, 0


def maxend(aGenomeSpace):
    """
    Return a list of the maximum stop coordinate of for each chromosome given in the genome space file. This assumes
    that the input list is sorted by chr, start, stop

    """
    hEnd = {}
    strChr = stop = 0
    for line in aGenomeSpace:
        if strChr:
            if line[0] == strChr:
                if line[2] > stop:
                    stop = line[2]
            else:
                hEnd[strChr] = stop
                strChr = line[0]
                stop = line[2]
        else:
            strChr = line[0]
            stop = line[2]
    hEnd[strChr] = stop
    return hEnd


def point_checker(iPoint, aInterval):
    """ Returns True if the given point is larger than interval start and smaller than interval end"""
    iIntervalStart, iIntervalEnd = aInterval
    if iIntervalStart <= iPoint <= iIntervalEnd:
        return True
    return False


def random_interval(aInputInterval, aWeightedSpace):
    """

    Generates a random interval matched to the length of the input interval, then checks that interval is in
    defined genome space before returning interval

    """
    iInputLen = aInputInterval[2] - aInputInterval[1]
    iCount = 1
    logging.debug("Looking to match: {}\t{}\t{}".format(*aInputInterval))
    # Get a space from weighted genome space, ensuring interval can fit
    while True:
        aSpace = picker(aWeightedSpace)
        iMinStart = aSpace[1]
        iMaxStart = aSpace[2] - iInputLen
        if iMaxStart >= iMinStart:
            break
        else:
            # Try to find a new space, but don't let loop hang
            iCount += 1
            if iCount > len(aWeightedSpace) * 10:  # Try at least 10x the number of given intervals
                logging.error("Can't place {2}\t{3}\t{4} of length {0} after {1} tries".format(iInputLen, iCount,
                                                                                               *aInputInterval))

                print "Could not find a space large enough to pick a random interval."
                sys.exit(1)
    # Pick random start on space and form random match
    iRandStart = random.randint(iMinStart, iMaxStart)
    iRandStop = iRandStart + iInputLen
    iRandChr = aSpace[0]
    logging.debug("Found match for: {}\t{}\t{}".format(*aInputInterval))
    return [iRandChr, iRandStart, iRandStop]


def overlap(aaIntervals, aaAgainst):
    iOverlapCount = iTotalBPOverlap = 0
    iIntervalCount = 0
    for aTestInterval in aaIntervals:
        iIntervalCount += 1
        logging.debug("Testing interval {}:\t{}".format(iIntervalCount, " ".join(map(str, aTestInterval))))
        for aAgainstInterval in aaAgainst:
            # Error checking to improve efficiency
            if aTestInterval[0] != aAgainstInterval[0]:
                continue
            if aTestInterval[1] > aAgainstInterval[2]:
                continue
            bOverlap, iOverlap = partial_overlap(aTestInterval, aAgainstInterval)
            if bOverlap:
                iOverlapCount += 1
                iTotalBPOverlap += iOverlap
                break
            elif int(aTestInterval[2]) < int(aAgainstInterval[1]):
                break
    return iOverlapCount, iTotalBPOverlap


def norm_distribution(aUCEs, aAgainst, aGenomeSpaceIntervals, iIterations, uceName, againstName):
    # Create list for distribution 
    aOverlapDistribution = []
    bLocPrint = bPrint
    iWrong = 0
    # Loop as many times as specified by iIterations
    for j in xrange(1, (iIterations + 1)):
        logging.debug("Iteration: {}".format(j))
        while True:
            try:
                # Create random match for each UCE
                aRandomMatches = [random_interval(uce, aGenomeSpaceIntervals) for uce in aUCEs]
                logging.debug("Found random match for each UCE in iteration {}".format(j))
                # Check # of matches and # of UCEs are concordant
                if not len(aUCEs) == len(aRandomMatches):
                    raise Exception("Error in creating random matches, could not find a 1 to 1 concordance")
                    # Sort random matches
                aRandomMatches.sort(key=lambda x: (x[0], x[1], x[2]))

                if bVerbose and not bLocPrint:
                    # Print random matches once
                    strRun1RandomFileName = 'run1_randommatches.dist' + str(uceName) + str(againstName) + '.txt'
                    print "Writing file to: " + strRun1RandomFileName
                    sys.stderr.write("Writing matches to " + strRun1RandomFileName + "\n")
                    with open(strRun1RandomFileName, "w") as out:
                        aWriteDistribution = ["\t".join(map(str, line)) for line in aRandomMatches]
                        out.write("\n".join(aWriteDistribution))
                        bLocPrint = True

                # Check that all random matches are non-overlapping
                iRandLen = len(aRandomMatches)
                iCollapsedLen = len(collapse(aRandomMatches))
                if iRandLen == iCollapsedLen:
                    logging.debug("Found {} randoms, {} collapsed randoms in iteration {}".format(iRandLen,
                                                                                                  iCollapsedLen, j))
                    raise FoundException
                else:
                    logging.info("Found {} randoms, {} collapsed randoms in iteration {}, "
                                 "retrying...".format(iRandLen, iCollapsedLen, j))
                    iWrong += 1
                    if iWrong > 100000:  # If randoms keep overlapping
                        print "Cannot find {0} non-overlapping random matches after 1000 tries".format(iRandLen)
                        print "Exiting..."
                        sys.exit(1)
            except NameError:
                print "found it"
            except FoundException:
                # Calculate # of overlaps and bp overlap for all random matches
                iOverlapCount, iTotalBPOverlap = overlap(aRandomMatches, aAgainst)
                logging.debug("Overlaps calculated for iteration {}".format(j))
                aOverlapDistribution.append([iOverlapCount, iTotalBPOverlap])
                break
    logging.info("Found {} instances where randoms overlapped".format(iWrong))
    return aOverlapDistribution


def cluster_distribution(aUCEs, aAgainst, aGenomeSpaceIntervals, iClusterWidth, iIterations, hChrEnds, uceName, againstName):
    try:
        import clustermodule
    except ImportError:
        print "Cannot find clustermodule.py. Ensure file is in working directory, exiting..."
        sys.exit(1)
    bLocPrint = bPrint
    iWrong = 0
    # Cluster UCEs, then associate clusters with UCEs
    aClusteredUCEs = clustermodule.cluster(aUCEs, iClusterWidth, hChrEnds)
    aAssocClusterUCEs = clustermodule.c_trackuces(aClusteredUCEs, aUCEs)
    # Create list for distribution 
    aOverlapDistribution = []
    # Check that there is enough space available to place clusters
    iClusterCoverage = sum([interval_len(line) for line in aClusteredUCEs])
    iSpaceCoverage = sum([interval_len(line[0]) for line in aGenomeSpaceIntervals])
    logging.info("{} cluster coverage, {} space coverage".format(iClusterCoverage, iSpaceCoverage))
    if iSpaceCoverage < iClusterCoverage:
        logging.error("Total coverage of clusters exceeds available space to place clusters")
        sys.exit("Total coverage of clusters exceeds available space to place "
                 "clusters")
        # Loop as many times as specified by iIterations
    for j in xrange(1, (iIterations + 1)):
        logging.debug("Iteration: {}".format(j))
        while True:
            try:
                # Place each cluster in the genome somewhere randomly
                aRandomClusterMatches = []
                for cluster in aAssocClusterUCEs:
                    aRandomCluster = random_interval(cluster[0], aGenomeSpaceIntervals)
                    # Assign associated UCEs to new cluster location
                    for UCE in cluster[1]:
                        strChr = aRandomCluster[0]
                        iStart = aRandomCluster[1] + UCE[0]
                        iStop = iStart + UCE[1]
                        aRandomClusterMatches.append([strChr, iStart, iStop])
                logging.debug("Random clusters created for iteration {}".format(j))
                # Sort clustered random matches
                aRandomClusterMatches.sort(key=lambda x: (x[0], x[1], x[2]))

                if bVerbose and not bLocPrint:
                    # Print random matches once
                    strRun1RandomFileName = 'run1_randommatches.dist' + str(uceName) + str(againstName) + '.txt'
                    print "Writing file to: " + strRun1RandomFileName
                    sys.stderr.write("Writing matches to " + strRun1RandomFileName)
                    with open(strRun1RandomFileName, "w") as out:
                        aWriteDistribution = ["\t".join(map(str, line)) for line in aRandomClusterMatches]
                        out.write("\n".join(aWriteDistribution))
                        bLocPrint = True

                iRandLen = len(aRandomClusterMatches)
                iCollapsedLen = len(collapse(aRandomClusterMatches))
                if iRandLen == iCollapsedLen:
                    logging.debug("Found {} randoms, {} collapsed randoms in iteration {}".format(iRandLen,
                                                                                                  iCollapsedLen, j))
                    raise FoundException
                else:
                    logging.info("Found {} randoms, {} collapsed randoms in iteration {}, "
                                 "retrying...".format(iRandLen, iCollapsedLen, j))
                    iWrong += 1
                    if iWrong > 1000:  # If randoms keep overlapping
                        print "Cannot find {0} non-overlapping random matches after 1000 tries".format(iRandLen)
                        print "Exiting..."
                        sys.exit(1)
            except FoundException:
                # Calculate # of overlaps and bp overlap for clustered random matches
                iOverlapCount, iTotalBPOverlap = overlap(aRandomClusterMatches, aAgainst)
                logging.debug("Overlaps calculated for iteration {}".format(j))
                aOverlapDistribution.append([iOverlapCount, iTotalBPOverlap])
                break

    logging.info("Found {} instances where randoms overlapped".format(iWrong))
    return aOverlapDistribution


def cluster_input(string):
    value = int(string)
    if not value > 0:
        msg = "Cluster width must be greater than 0 kb"
        raise argparse.ArgumentTypeError(msg)
    return value


def cdf(x, mu, sigma):
    y = 0.5 * (1 + math.erf((x - mu) / math.sqrt(2 * sigma ** 2)))
    return y


def Proportion(bp, aOverlapBP):
    "Returns the proportion of random overlaps at or more extreme than the UCE overlaps"
    mean = float(sum(aOverlapBP) / len(aOverlapBP))
    iMean = int(mean)
    #Determine if the uce overlaps value is greater or less than the mean of the randomoverlaps
    if iMean > bp:
        print 'UCE overlaps below random overlaps mean: set may be depleted'
        #Return all values in npArSortedRandomOverlaps that are smaller than or equal to ibpUCEoverlap
        arRandomsLessorEqual = [x for x in aOverlapBP if x <= bp]
        #Calculate the length of this list
        iSmallerOrEqualRandomOverlaps = len(arRandomsLessorEqual)
        proportion = float(float(iSmallerOrEqualRandomOverlaps)/float(len(aOverlapBP)))
    elif iMean <= bp:
        print 'UCE overlaps above random overlaps mean: set may be enriched'
        arRandomsGreaterOrEqual = [x for x in aOverlapBP if x >= bp]
        iRandomsGreaterOrEqualNumber = len(arRandomsGreaterOrEqual)
        proportion = float(iRandomsGreaterOrEqualNumber)/float(len(aOverlapBP))
    return proportion

def KSTest(aOverlapBP):
    "Returns the KS test statistic and p value for rejecting the null hypothesis that aOverlapBP follows a normal distribution with mean and sd equal to those of aOverlapBP"
    mean = float(sum(aOverlapBP) / len(aOverlapBP))
    sd = stdev(aOverlapBP)
    rvNormMatched = stats.norm.rvs(loc=mean, scale=sd, size=1000)
    npArOverlapBP = np.array(aOverlapBP)
    ksStat, KsPval = stats.ks_2samp(npArOverlapBP, rvNormMatched)
    if KsPval <= 0.05:
        strKSresult = "No"
        print 'KS statistic is significant: attention needed'
    else:
        strKSresult = "Yes"
        print 'KS statistic not significant: random overlaps appear normally distributed'
    return ksStat, KsPval, strKSresult

def statistics(aUCEOverlaps, aOverlapDistribution):
    n = aUCEOverlaps[0]
    bp = aUCEOverlaps[1]
    aOverlapBP = zip(*aOverlapDistribution)[1]
    mean = float(sum(aOverlapBP) / len(aOverlapBP))
    sd = stdev(aOverlapBP)
    minimum = min(aOverlapBP)
    maximum = max(aOverlapBP)
    pvalue = cdf(bp, mean, sd)
    obsExp = float(bp) / mean
    proportion = Proportion(bp, aOverlapBP)
    ksStat, ksPval, strKSresult = KSTest(aOverlapBP)
    if pvalue >= 0.975:
        strZtestResult = "Enriched"
    elif pvalue <= 0.025:
        strZtestResult = "Depleted"
    else:
        strZtestResult = "Neither"
    if pvalue > 0.5:
        transformedPvalue = float(1-pvalue)
        return [n, bp, mean, sd, minimum, maximum, ksPval, strKSresult, proportion, transformedPvalue, obsExp, strZtestResult]
    else:
        return [n, bp, mean, sd, minimum, maximum, ksPval, strKSresult, proportion, pvalue, obsExp, strZtestResult]


def stdev(aList):
    dMean = float(sum(aList)) / len(aList)
    dSumOfSquares = sum([((number - dMean) ** 2) for number in aList])
    try:
        dVariance = float(dSumOfSquares) / (len(aList) - 1)
    except ZeroDivisionError:
        print "Cannot calculate statistical variance with only 1 iteration."
        print "Exiting..."
        sys.exit(1)
    return math.sqrt(dVariance)


def writer(aList, uceName, againstName):
    if bVerbose:
        strStatsFileName = 'stats_' + str(uceName) + str(againstName) + '.txt'
        sys.stderr.write("Writing matches to " + strStatsFileName + "\n")
        with open(strStatsFileName, "w") as out:
            out.write("n\tbp\tmean\ts.d.\tmin\tmax\tksPval\tKSresult\tproportion\tp-value\tObs/Exp\tZtestResult\n")
            out.write("\t".join(map(str, aList)))
    print "n\tbp\tmean\ts.d.\tmin\tmax\tproportion\tksPval\tKSresult\tp-value\tObs/Exp\tZtestResult\n"
    print "\t".join(map(str, aList))


def interval_len(aInterval):
    """Given a 3 column interval, return the length of the interval """
    return aInterval[2] - aInterval[1] + 1


def weight(aList):
    aWeightedList = []
    iTotal = sum([interval_len(interval) for interval in aList])  # Get total list coverage
    for interval in aList:
        # Calculate weight (by length of interval) and make tuples of each interval with its weight
        dWeight = float(interval_len(interval)) / iTotal
        aWeightedList.append((interval, dWeight))
    return aWeightedList


def picker(aList):
    x = random.random()
    for elmt, weight in aList:
        if x <= weight:
            return elmt
        x -= weight


def formatInt(aInterval):
    """ Format an 3-column interval correctly """
    return [aInterval[0], int(aInterval[1]), int(aInterval[2])]


def getArgs(strInput=None, verbose=True):
    # Define arguments
    parser = argparse.ArgumentParser(description="This script performs a depletion analysis with the given arguments "
                                                 "and files. The intervals being tested (normally UCEs) are given, as "
                                                 "well as the intervals that are checked for UCE depletion/enrichment. "
                                                 "The script will print the statistical analysis for the run to the "
                                                 "terminal. Unless otherwise specified, interval files should be "
                                                 "provided as arguments. All interval files must be in 1-based and in "
                                                 "3-column format: chr, start, stop")
    parser.add_argument("-u", "--uces", type=argparse.FileType("rU"), required=True,
                        help="The intervals to test (normally UCEs).")
    parser.add_argument("-g", "--genomespace", type=argparse.FileType("rU"), required=True,
                        help="The set of intervals defining the genomic space random sets are to be drawn from")
    parser.add_argument("-a", "--against", type=argparse.FileType("rU"), required=True,
                        help="The set of intervals that are being tested for overlap with UCEs. Total coverage should "
                             "be >= 20 Mb to provide sufficient statistical power.")
    parser.add_argument("-i", "--iterations", type=int, default=1000,
                        help="The number of random sets created to build an expected distribution [default=1000]")
    parser.add_argument("-c", "--cluster", type=cluster_input,
                        help="The maximum size to cluster adjacent intervals (kb)")
    parser.add_argument("-v", "--verbose", action="store_false",
                        help="-v flag prevents the storage of various intermediate files to current directory")
    parser.add_argument("-d", "--debug",
                        help="Debug level [default = None]")
    # Collect arguments
    if strInput:
        if verbose:
            print "Given debug argument string: {0}".format(strInput)
        return parser.parse_args(strInput.split())
    return parser.parse_args()



def main(args):
    # Set debugging level
    if args.debug:
        log_level = LOGGING_LEVELS.get(args.debug.lower(), logging.NOTSET)
        logging.basicConfig(level=log_level, filename=str("debug.log." + str(args.uces.name) + str(args.against.name)), filemode="w",
                            format='%(asctime)s\t%(levelname)s\t%(message)s', datefmt='%Y-%m-%d %H:%M:%S')
    else:
        logging.basicConfig()

    logging.debug(str(args))
    logging.info("Running {} against {} {} times, using {} as the genome space".format(args.uces, args.against,
                                                                                       args.iterations,
                                                                                       args.genomespace))

    # Create interval lists for UCEs, genome space regions and "against" regions
    logging.debug("Reading input files into lists...")
    aUCEs = [formatInt(line.strip().split("\t")) for line in args.uces]
    aAgainst = [formatInt(line.strip().split("\t")) for line in args.against]
    aGenomeSpaceIntervals = [formatInt(line.strip().split("\t")) for line in args.genomespace]
    aGenomeSpaceIntervals.sort(key=lambda x: (x[0], x[1]))
    logging.debug("Lists read and intervals formatted")

    # Weight genome space intervals, only selecting big enough regions if clustered
    if args.cluster:
        # Convert to bp
        iClusterWidth = args.cluster * 1000
        # Get max stop for each chromosome before reducing space
        hEnds = maxend(aGenomeSpaceIntervals)
        aSpace = [line for line in aGenomeSpaceIntervals if interval_len(line) > iClusterWidth]
        logging.info("{} intervals reduced to {} intervals over {} bp".format(len(aGenomeSpaceIntervals), len(aSpace),
                                                                              iClusterWidth))
        if len(aSpace) < 1:
            logging.error("{} intervals over {} bp found, exited".format(len(aSpace), iClusterWidth))
            sys.exit("Cluster size exceeds any one interval in the genome space file, try reducing cluster size")
        aWeightedSpace = weight(aSpace)
    else:
        aWeightedSpace = weight(aGenomeSpaceIntervals)

    # Sort lists
    logging.debug("Sorting lists...")
    aUCEs.sort(key=lambda x: (x[0], x[1], x[2]))
    aWeightedSpace.sort(key=lambda e: e[0])
    aAgainst.sort(key=lambda x: (x[0], x[1], x[2]))
    logging.debug("Lists sorted")

    # Initialize global variables
    global bVerbose
    bVerbose = args.verbose
    global bPrint
    bPrint = False

    # Create distribution of random overlaps, depending on cluster flag
    if args.cluster:
        aOverlapDistribution = cluster_distribution(aUCEs, aAgainst, aWeightedSpace, args.cluster, args.iterations,
                                                    hEnds, args.uces.name, args.against.name)
    else:
        aOverlapDistribution = norm_distribution(aUCEs, aAgainst, aWeightedSpace, args.iterations, args.uces.name, args.against.name)

    logging.debug("Distribution created")
    # Write distribution to file
    if bVerbose:
        strRandomMatchFileName = 'randommatches.dist' + str(args.uces.name) + str(args.against.name) + '.txt'
        print "Writing file to: " + strRandomMatchFileName
        with open(strRandomMatchFileName, "w") as out:
            aWriteDistribution = ["\t".join(map(str, line)) for line in aOverlapDistribution]
            out.write("\n".join(aWriteDistribution))

    # Get UCE overlaps and calculate statistics
    aUCEOverlaps = overlap(aUCEs, aAgainst)
    aStats = statistics(aUCEOverlaps, aOverlapDistribution)
    return aStats


if __name__ == "__main__":
    args = getArgs()
    aStats = main(args)
    writer(aStats, args.uces.name, args.against.name)
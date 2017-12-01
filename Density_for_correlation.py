"""
Script to make windows using pybedtools,
then calculate coverage of each windows by a feature file.

Here, a bed file is provided to give the intervals to be split into
intervals. In this way, nonN regions can be accounted for and
windows containing nonN regions will not have artificially low coverage.


Requires bedtools genome file, in same directory as script.

Requires all bed files to have 0-based starts.

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
import argparse
import pybedtools as pbt
import pandas as pd


def get_args():
    parser = argparse.ArgumentParser(description="Description")
    parser.add_argument("-b", "--bin", type=int)
    parser.add_argument("file", type=argparse.FileType('rU'),
                        help='A file containing a list of paths to the feature files you want to process, '
                        'separated by newlines')
    parser.add_argument("-w", "--windows", type=str, default="hg18.genomic.coordinates.nonN_0based.bed")

    return parser.parse_args()


def getFileNames(args):
    arFilenames = [line.strip() for line in args.file]
    return arFilenames


def makeWindows(args):
    btWindows = pbt.BedTool().window_maker(b=str(args.windows), w = args.bin)
    #Save to a file
    strFilename = str("Intervals_{0}_in_{1}_windows_0based.bed".format(str(args.windows), str(args.bin)))
    btWindows.saveas(strFilename)
    return btWindows


def makeEqualWindows(args):
    btWindows = pbt.BedTool().window_maker(b=str(args.windows), w = args.bin)
    intWindowLength = args.bin
    btEqualWindows = btWindows.filter(lambda x: len(x) == intWindowLength)
    #Turn into a file-based bedtool
    btSavedEqualWindows = btEqualWindows.saveas('EqualWindowsTemp.bed')
    #Get bedtool back
    btEqualWindowsNotStream = getFeatures('EqualWindowsTemp.bed')
    return btEqualWindowsNotStream


def filterBedtool(bedTool, intLength):
    pdFeatures = pd.read_table(bedTool.fn, names=['chr', 'start', 'end'])
    #Calculate length of each feature
    pdFeatures['length'] = pdFeatures['end'] - pdFeatures['start']


def coverage(btWindows, btFeatures):
    btCoverage = btWindows.coverage(btFeatures)
    return btCoverage


def cutToFraction(btObject):
    #Obtains columns 0, 1, 2, and 6 from the output of coverage, i.e. the intervals and the fraction of bp covered
    #in that interval
    btCutToFraction = btObject.cut([0,1,2,6])
    return btCutToFraction


def saveBedTool(btObject, strFilename):
    btObject.saveas(strFilename)


def getFeatures(strFileName):
    btFeatures = pbt.BedTool(strFileName)
    return btFeatures


def getDensity(btWindows, btFeatures):
    btCoverage = coverage(btWindows, btFeatures)
    btCutToFraction = cutToFraction(btCoverage)
    return btCoverage, btCutToFraction


def forEachFeatureFile(strFileName, btWindows):
    #Obtain btFeatures from the input feature file
    btFeatures = getFeatures(strFileName)
    #Obtain full output from coverage (btCoverage) and just the 4 column output with fractions covered (btCutToFraction)
    btCoverage, btCutToFraction = getDensity(btWindows, btFeatures)
    #Save both files with appropriate names
    saveBedTool(btCoverage, str('Full_coverage_{0}.bed'.format(strFileName)))
    saveBedTool(btCutToFraction, str('Density_for_matrix_{0}.bed'.format(strFileName)))
    #return both
    return btCoverage, btCutToFraction


def main():
    print 'Getting arguments'
    args = get_args()
    print 'Making windows'
    btWindows = makeEqualWindows(args)
    #Save btWindows to file
    saveBedTool(btWindows, 'Bins_for_density_0based.bed')

    #Get filenames
    arFilenames = getFileNames(args)
    for filename in arFilenames:
        #Get density for this feature file
        print 'Getting density files for {0}'.format(filename)
        btCoverage, btCutToFraction = forEachFeatureFile(filename, btWindows)


if __name__ == "__main__":
     main()



"""
Script to obtain domains containing UCEs and then measure distance from UCE to nearest TSS, and distance
from domain edge to nearest TSS.

Ruth McCole
January 17th 2018

Copyright 2018 Harvard University, Wu Lab

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

    http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

"""

import argparse
import pandas as pd
import pybedtools as pbt
import seaborn as sns
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from scipy import stats


def get_args():
    parser = argparse.ArgumentParser(description="Description")
    parser.add_argument('-u', "--UCEs", type=str,
                        help='File of UCEs, 0 based, no headers')
    parser.add_argument('-c', "--compareUCEs", type=str,
                        help='File of UCEs to compare to first file, 0 based, no headers, not compatible with'
                             'randoms')
    parser.add_argument('-d', '--domains', type=str,
                        help='File of domains, 0 based, no headers')
    parser.add_argument('-g', '--genes', type=str,
                        help='File of genes, 0 based, with headers')
    parser.add_argument('-f', '--file', type=argparse.FileType('rU'),
                        help='Filename of file containing filenames of all the domains')
    parser.add_argument('-r', '--randoms', type=argparse.FileType('rU'),
                        help='Filename of file containing filenames of all the random regions')
    parser.add_argument('-b', '--bins', type=int, default=10,
                        help='Number of bins for your histogram counts, default is 10')
    parser.add_argument('-x', '--xLim', type=int, nargs=2, default=(0, 900),
                        help='Upper and lower limits of x axis, supply two integers')
    parser.add_argument('-a', '--allDomains', action='store_true', dest='boolA',
                        help='pass -a if you would like to put all the domain types into one big'
                             'set and calculate UCE-TSS and random-TSS distances for this whole set.')

    return parser.parse_args()


def getFeatures(strFileName):
    btFeatures = pbt.BedTool(strFileName)
    return btFeatures


def getDataIntoPandas(strFilename):
    #NB: this assumes your files DO have headers
    print 'Getting data from {0}'.format(strFilename)
    pdData = pd.read_csv(strFilename, sep='\t')
    return pdData


def getDataIntoPandasNoHeader(strFilename, arHeaders):
    #This is for files without headers
    #print 'Getting data from {0}'.format(strFilename)
    pdData = pd.read_csv(strFilename, sep='\t', header=None, names=arHeaders)
    return pdData


def pandaToBedtool(panda):
    arArFeatures = panda.values.tolist()
    btFeatures = getFeatures(arArFeatures)
    return btFeatures


def bedtoolToPanda(btobject):
    pdObject = pd.read_table(btobject.fn, header=None)
    return pdObject


def savePanda(pdData, strFilename):
    pdData.to_csv(strFilename, sep='\t', index=False)


def getFileNames(args):
    arFilenames = [line.strip() for line in args.file]
    return arFilenames


def getTSS(pdGenes):
    #Pull out + and -, select correct fields, put back together
    pdPlusStrand = pdGenes[pdGenes.Gene_strand == '+']
    pdMinusStrand = pdGenes[pdGenes.Gene_strand == '-']

    pdTSSPlus = pdPlusStrand[['Gene_chrom', 'Gene_txStart', 'Gene_name',
                              'GeneSymbol']]

    pdTSSMinus = pdMinusStrand[['Gene_chrom', 'Gene_txEnd', 'Gene_name',
                              'GeneSymbol']]

    #Now create intervals of 1 base so they fit into bedtools format. For plus, add 1 to create end coordinate.
    #For minus, minus 1 to create start coordinate

    pdTSSPlus['TSS_end'] = pdTSSPlus.apply(lambda row: (int(row['Gene_txStart']) + 1), axis=1)
    pdTSSPlus.rename(columns={'Gene_txStart': 'TSS_start'}, inplace=True)


    pdTSSMinus['TSS_start'] = pdTSSMinus.apply(lambda row: (int(row['Gene_txEnd']) - 1), axis=1)
    pdTSSMinus.rename(columns={'Gene_txEnd': 'TSS_end'}, inplace=True)

    #Combine dataframes
    pdConcat = pd.concat([pdTSSPlus, pdTSSMinus], ignore_index=True)

    #Rearrange columns
    pdTSS = pdConcat[['Gene_chrom', 'TSS_start', 'TSS_end', 'Gene_name', 'GeneSymbol']]

    return pdTSS


def createColumnHeaders(pdData, strUniqueWord):
    #Creates column headers containing the string strUniqueWord, and also names the column you want to merge on 'Merge'
    #Get number of columns
    intColumnNumber = pdData.shape[1]
    if intColumnNumber == 3:
        pdData.columns = ['chr', 'start', 'end']
    elif intColumnNumber > 3:
        #Generate list of column names
        arColumnNames = ['chr', 'start', 'end']
        for i in range(3, intColumnNumber):
            strColumnName = strUniqueWord + str(i)
            arColumnNames.append(strColumnName)
        #Assign new column names to dataframe
        pdData.columns = arColumnNames
    else:
        print '{0} columns found: problem'.format(str(intColumnNumber))
    return pdData


def sortFeaturesPandaWithBedtools(pdFeatures):
    #Sort a panda that contains chr, start, and end columns
    #Get headers
    arHeaders = list(pdFeatures)
    #Make bedtool
    btFeatures = pandaToBedtool(pdFeatures)
    #Sort bedtool
    btFeaturesSorted = btFeatures.sort()
    #Make panda
    pdFeaturesSorted = pd.read_table(btFeaturesSorted.fn, header=None)
    #Assign back the headers
    pdFeaturesSorted.columns = arHeaders

    return pdFeaturesSorted


def findIdenticalFeatures(pdFeaturesWithColumnNames):
    pdBoolDupsTrue = pdFeaturesWithColumnNames.duplicated(keep='first')
    pbBoolDupsFalse = ~pdBoolDupsTrue
    pdDups = pdFeaturesWithColumnNames[pdBoolDupsTrue]
    pdNoDups = pdFeaturesWithColumnNames[pbBoolDupsFalse]

    return pdNoDups, pdDups


def removeOldNewDups(pdConcat):
    #Discover duplicates just for breakpoint coordinates, ignorning source
    boolTrueDup = pdConcat.duplicated(subset=['Gene_chrom', 'TSS_start', 'TSS_end'], keep='first')
    #Remove these rows
    pdConcatNoOldNewDups = pdConcat[~boolTrueDup]

    #Also save the rows you removed
    pdDupsFromFullDataset = pdConcat[boolTrueDup]

    return pdConcatNoOldNewDups, pdDupsFromFullDataset


def giveIDs(pdNoDups, strIDWord):
    #First reindex pdNoDups, so that the adding of the ID column with assign works properly, if duplicates were dropped
    pdNoDups.reset_index(inplace=True, drop=True)

    #Make a list of the IDs to add, as a pandas dataframe
    #First find number of rows
    intNoRows = pdNoDups.shape[0]
    arIDColumn = []

    for int in range(1, intNoRows+1):
        strID = strIDWord + str(int)
        arIDColumn.append(strID)
    pdIDs = pd.DataFrame(arIDColumn)

    #Add these IDs as a new column
    pdWithIDs = pdNoDups.assign(ID=pdIDs)

    return pdWithIDs


def giveUniqueIDs(pdDomains):
    # Give columns names
    pdFeaturesWithColumnNames = createColumnHeaders(pdDomains, 'column')

    # Sort
    pdSorted = sortFeaturesPandaWithBedtools(pdFeaturesWithColumnNames)

    # Find and remove identical features
    pdNoDups, pdDups = findIdenticalFeatures(pdSorted)

    # Give unique IDs
    pdWithIDs = giveIDs(pdNoDups, 'Domain')

    return pdWithIDs


def getData():
    #Obtain UCEs, domains, and genes
    args = get_args()
    strUCEFilename = args.UCEs
    strGenesFilename = args.genes
    intBins = args.bins
    tupXlim = args.xLim


    pdUCEs = getDataIntoPandasNoHeader(strUCEFilename,
                                       ['chr_UCEs', 'start_0based_UCEs', 'end_UCEs', 'type', 'UCE_ID'])


    pdGenes = getDataIntoPandas(strGenesFilename)
    #Remove dots from pdGenes headers
    pdGenes.columns = ['Gene_name', 'Gene_chrom', 'Gene_strand', 'Gene_txStart',
                       'Gene_txEnd', 'GeneSymbol']
    # Get TSS
    pdTSS = getTSS(pdGenes)
    # Deduplicate TSS
    pdTSSNoFeatureDups, pdTSSFeatureDups = removeOldNewDups(pdTSS)


    #Dictionary to hold domains
    dictDomains = {}

    if args.domains:
        strDomainsFilename = args.domains

        pdDomains = getDataIntoPandasNoHeader(strDomainsFilename,
                                              ['chr_domains', 'start_0based_domains', 'end_domains'])

        pdDomainsWithIDs = giveUniqueIDs(pdDomains)

        #Add name of file plus domains to dictionary
        dictDomains[strDomainsFilename] = pdDomainsWithIDs


    elif args.file:
        #You are processing multiple domain files

        arFilenames = getFileNames(args)
        #print arFilenames

        #Want to include a name to go with each domain for file saving

        arPdTheseDomains = []


        for strFilename in arFilenames:
            pdTheseDomains = getDataIntoPandasNoHeader(strFilename,
                                              ['chr_domains', 'start_0based_domains', 'end_domains'])
            arPdTheseDomains.append(pdTheseDomains)

            pdTheseDomainsWithIDs = giveUniqueIDs(pdTheseDomains)

            dictDomains[strFilename] = pdTheseDomainsWithIDs

        if args.boolA:
            #Squish all the domains together
            pdAllDomainsTogether = pd.concat(arPdTheseDomains, ignore_index=True)
            #Give them uniqueIDs
            pdAllDomainsWithIDs = giveIDs(pdAllDomainsTogether, 'All_domains_')
            dictDomains['All'] = pdAllDomainsWithIDs

        print 'Domains of the following types will be analyzed'
        for key, value in dictDomains.iteritems():
            print key

    else:
        print 'You must specify either a single domain file with -d or a filename containing a list of the' \
              'names of such files using -f'

    return dictDomains, pdUCEs, pdTSSNoFeatureDups, intBins, tupXlim


def AddKbMbColumns(pdDistances, strDistanceColumn):
    #Obtain the distances scaled to kb
    pdDistances['Distance_Kb'] = pdDistances.apply(lambda row: (float(row[strDistanceColumn])/1000), axis=1)
    #Obtain the distances in Mb
    pdDistances['Distance_Mb'] = pdDistances.apply(lambda row: (float(row[strDistanceColumn])/1000000), axis=1)

    return pdDistances


def getRandoms(args):
    arRandomFilenames = [line.strip() for line in args.randoms]

    arPdRandoms = []


    intRandomIteration = 0

    for strRandomFilename in arRandomFilenames:
        pdTheseRandomRegions = getDataIntoPandasNoHeader(strRandomFilename, ['chr_random', 'start_0based_random',
                                                                             'end_0based_random', 'type'])
        #Give random IDs
        intRandomIteration += 1
        pdTheseRandomRegionsIDs = giveIDs(pdTheseRandomRegions, 'Random_set_{0}_'.format(intRandomIteration))

        # Keep as individual pandas, and also combine into one giant panda
        arPdRandoms.append(pdTheseRandomRegionsIDs)

        pdAllRandoms = pd.concat(arPdRandoms, ignore_index=True)

    return pdAllRandoms, arPdRandoms


def getUCEDomainsAndTSS(dictDomains, pdUCEs, pdTSS, strRandomOrUCE):
    #Convert UCEs and TSS to bedtools
    btUCEs = pandaToBedtool(pdUCEs)
    btTSS = pandaToBedtool(pdTSS)

    #Iterate over the domains dictionary

    dictDistanceResults = {}

    for strDomainName, pdDomain in dictDomains.iteritems():
        pdUCETSSDistances = getDistancesLimitedToDomains(pdDomain, pdUCEs, pdTSS)

        dictDistanceResults[strDomainName] = pdUCETSSDistances

        savePanda(pdUCETSSDistances, '{0}_to_TSS_distance_{1}.txt'.format(strRandomOrUCE, strDomainName))

        #seriesDistances = pdAllDistancesUCEsToTSS[['Distance_bp']]

    return dictDistanceResults


def getDistancesLimitedToDomains(pdDomain, pdRegions, pdTSS):
    btRegions = pandaToBedtool(pdRegions)
    btTSS = pandaToBedtool(pdTSS)
    btDomains = pandaToBedtool(pdDomain)

    btRegionsInDomains = btRegions.intersect(btDomains, u=True)
    btTSSInDomains = btTSS.intersect(btDomains, u=True)

    btRegionsInDomainsSorted = btRegionsInDomains.sort()
    btTSSInDomainsSorted = btTSSInDomains.sort()

    btDistances = btRegionsInDomainsSorted.closest(btTSSInDomainsSorted, d=True)

    pdDistances = bedtoolToPanda(btDistances)
    pdDistances.columns = ['chr_region', 'start_0based_region', 'end_region', 'type', 'Region_ID', 'Gene_chrom', 'TSS_start', 'TSS_end', 'Gene_name', 'GeneSymbol', 'Distance_bp']

    return pdDistances


def getSeriesWithBinsIndex(npHistDensityY, npHistDensityX):
    #Round the numbers
    npHistXRounded = np.around(npHistDensityX, decimals=3)

    #Get into list
    arHistXRounded = npHistXRounded.tolist()

    #Get rid of the last bin - there are always one more bin edges than there are numbers of bins, and you want
    #to make bin edges into lables (or indexes) for the counts, so you just forget the last one.
    #e.g. bin 0 to 0.1 is now just called bin 0, bin 0.49 to 0.5 is called bin 0.49
    arHistXRounded.pop()

    #Make into a series
    seriesHistDensity = pd.Series(npHistDensityY, index=arHistXRounded)

    return seriesHistDensity


def getNumpyHist(seriesDistances, intHistBins, tupXLim):
    #Get into list
    arDistances = seriesDistances.values.tolist()

    #Use np.histogram to split into 50 bins
    npHistCountsY, npHistCountsX = np.histogram(arDistances, bins=intHistBins, range=tupXLim, density=False)

    seriesCountsWithBins = getSeriesWithBinsIndex(npHistCountsY, npHistCountsX)

    return seriesCountsWithBins


def prepareForPlot(pdData, intBins, tupMaxMin):
    # Distances must be got into a histogram

    # Get data into series
    seriesKbDistances = pdData['Distance_Kb']
    # bin for histogram
    seriesCountsWithBins = getNumpyHist(seriesKbDistances, intBins, tupMaxMin)

    return seriesCountsWithBins


def getCountsWithBinsForUCEsAndRandoms(dictFullResults, intBins, tupMaxMin):

    #Dictionary for gathering results
    dictCountsWithBins= {}

    for strKeyDomain, tupleResults in dictFullResults.iteritems():
        print 'Getting counts for {0} bins for UCEs and randoms, for domains {1}'.format(intBins, strKeyDomain)
        arPdRandomDistances = tupleResults[1]
        arSeriesCountsForJoin = []

        for pdRandomDistances in arPdRandomDistances:
            pdRandomDistancesWithKb = AddKbMbColumns(pdRandomDistances, 'Distance_bp')

            #bin for histogram
            seriesRandomCountsWithBins = prepareForPlot(pdRandomDistancesWithKb, intBins, tupMaxMin)

            arSeriesCountsForJoin.append(seriesRandomCountsWithBins)

        #Now join the series together on indexes to make a wide dataframe from which you can get summary stats
        pdRandomsForHistCombined = pd.concat(arSeriesCountsForJoin, axis=1)

        #Get UCE-TSS distances as counts for histogram
        pdUCEDistances = tupleResults[0]
        pdUCEDistancesWithKb = AddKbMbColumns(pdUCEDistances, 'Distance_bp')

        seriesUCEDistancesCountsWithBins = prepareForPlot(pdUCEDistancesWithKb, intBins, tupMaxMin)

        #Insert into counts panda
        pdRandomsForHistCombined.insert(loc=0, column=strKeyDomain, value=seriesUCEDistancesCountsWithBins)

        #Add final counts panda to a dictionary with its corresponding domain name as key
        dictCountsWithBins[strKeyDomain] = pdRandomsForHistCombined
        #Save
        pdRandomsForHistCombined.to_csv('Counts_for_hist_{0}.txt'.format(strKeyDomain), sep='\t')

    return dictCountsWithBins


def getRandomDistancesSeparateAndAllTogether(dictDomains, dictUCEDistanceResults, arPdRandoms, pdTSSNoFeatureDups):
    dictFullResults = {}

    for strDomainKey, pdDomain in dictDomains.iteritems():

        arDistanceBpValuesForThisDomain = []
        arPdRandomDistancesForThisDomain = []

        intRandomIteration = 0

        for pdRandom in arPdRandoms:
            ###New bit
            pdDistances = getDistancesLimitedToDomains(pdDomain, pdRandom, pdTSSNoFeatureDups)
            intRandomIteration += 1
            #Save this
            savePanda(pdDistances, '{0}_to_TSS_distance_{1}.txt'.format('Random_set_{0}'.format(intRandomIteration), strDomainKey))

            # Gather all random distances for each domain type, as list of pandas
            arPdRandomDistancesForThisDomain.append(pdDistances)

            #Obtain series of distancebp - later to be fed to AD test stuff
            arDistanceBpForTheseRandoms = pdDistances['Distance_bp'].values.tolist()
            arDistanceBpValuesForThisDomain.extend(arDistanceBpForTheseRandoms)

        #Combine these with UCE distances for final output
        if strDomainKey in dictUCEDistanceResults:
            print 'Collating full results for {0}'.format(strDomainKey)
            tupleDistanceResults = (dictUCEDistanceResults[strDomainKey], arPdRandomDistancesForThisDomain, arDistanceBpValuesForThisDomain)

            dictFullResults[strDomainKey] = tupleDistanceResults

    return dictFullResults


def processRandoms(args, dictDomains, dictUCEDistanceResults, pdTSSNoFeatureDups, intBins, tupXlim):
    pdAllRandoms, arPdRandoms = getRandoms(args)

    # Process randoms, all from one big panda, and get distance results
    print 'Calculating distances from randoms to TSS'
    #dictRandomDistanceResults = getRandomDistances(dictDomains, pdAllRandoms, pdTSSNoFeatureDups)

    # Process randoms one at a time, get results and also produce summary stats for plotting as histogram

    dictFullResults = getRandomDistancesSeparateAndAllTogether(dictDomains, dictUCEDistanceResults, arPdRandoms, pdTSSNoFeatureDups)
    dictCountsWithBins = getCountsWithBinsForUCEsAndRandoms(dictFullResults, intBins, tupXlim)

    #Save counts for plotting script

    return dictFullResults


def testAD(pdUCETSSDistances, pdRandomTSSDistances):
    #Get two lists of values
    #Feed into AD test
    arUCETSSbpDistances = pdUCETSSDistances['Distance_bp'].values.tolist()
    arRandomTSSbpDistances = pdRandomTSSDistances['Distance_bp'].values.tolist()

    floatStat, critical, approxP = stats.anderson_ksamp([arUCETSSbpDistances, arRandomTSSbpDistances])
    return floatStat, critical, approxP


def testADLists(arUCEDistances, arRandomDistances):

    floatStat, critical, approxP = stats.anderson_ksamp([arUCEDistances, arRandomDistances])
    return floatStat, critical, approxP


def getADResultsFromFullDict(dictFullResults):
    arADResults = []
    for key, tupleResults in dictFullResults.iteritems():
        pdUCEDistances = tupleResults[0]
        arUCEDistances = pdUCEDistances['Distance_bp'].values.tolist()
        arRandomDistances = tupleResults[2]

        floatStat, critical, approxP = testADLists(arUCEDistances, arRandomDistances)
        arADResults.append([key, approxP])

    pdADResults = pd.DataFrame(arADResults, columns=['Dataset', 'AD_pValue'])

    savePanda(pdADResults, 'AD_results_for_all_domains.txt')


def main():
    args = get_args()
    #Get and tidy data
    dictDomains, pdUCEs, pdTSSNoFeatureDups, intBins, tupXlim = getData()

    #Get UCE-TSS shortest distances in bp, only considering UCEs and TSS in the same domain
    dictUCEDistanceResults = getUCEDomainsAndTSS(dictDomains, pdUCEs, pdTSSNoFeatureDups, 'UCE')

    if args.compareUCEs:
        pdCompareUCEs = getDataIntoPandasNoHeader(args.compareUCEs,
                                                  ['chr_UCEs', 'start_0based_UCEs', 'end_UCEs', 'type', 'UCE_ID'])
        dictCompareUCEDistance = getUCEDomainsAndTSS(dictDomains, pdCompareUCEs, pdTSSNoFeatureDups, 'UCEs_compare')

    #Get randoms, if supplied
    args = get_args()
    if args.randoms:
        dictFullResults = processRandoms(args, dictDomains, dictUCEDistanceResults, pdTSSNoFeatureDups,
                                                                              intBins, tupXlim)
        #Each item in dictFullResults is a tuple.
        # Item 0 is UCE TSS distances as panda.
        # Item 1 is list of pandas, one for each random set provided
        #Item 2 is a flat list of all random-TSS distances

        getADResultsFromFullDict(dictFullResults)

    else:
        print 'You have chosen not to supply random regions'

if __name__ == "__main__":
    main()
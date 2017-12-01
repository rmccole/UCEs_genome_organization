"""
Script to take in files made by Density_for_correlation.py and compare them for having matching coordinates. Then create a new
table with one set of coordinates and all the density data from the different input files alongside.

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

import pandas as pd
import argparse

def get_args(strInput=None):
    """

    Collect arguments from command-line, or from strInput if given (only used for debugging)
    """
    parser = argparse.ArgumentParser(description="")
    parser.add_argument('file', type=argparse.FileType('rU'),
                        help="A file containing a list of paths to the files you want to process, separated by "
                             "newlines")

    return parser.parse_args()


def getIntoPandas(arFilenames):

    #Get table from each filename into panda dataframe, make a list of these dataframes

    arPandas = []

    for strFilename in arFilenames:
        #Provide a list of column names
        arColumnNames = ['chr', 'start', 'end', strFilename]

        #Import a table
        pd_density = pd.read_csv(strFilename, sep='\t', names = arColumnNames)

        #Make the three coordinates column into one that can be merged on
        #First add a comma to the end of the chr column
        pd_density['chr'] = pd_density['chr'].astype(str) + ','

        pd_density['start'] = pd_density['start'].astype(str) + ','

        #Combine together the first 3 columns into a new column called coordinate
        pd_density['coordinate'] = pd_density['chr'].map(str) + pd_density['start'].map(str) + pd_density['end'].map(str)

        #remove the columns we don't need any more
        pd_density = pd_density.drop(['chr', 'start', 'end'], axis=1)

        #For attractiveness, put the coordinate column first
        pd_density = pd_density[['coordinate', strFilename]]

        #Append this table to the list of tables
        arPandas.append(pd_density)

    return arPandas



def makeOnePanda(arPandas):
    #from http://stackoverflow.com/questions/23668427/pandas-joining-multiple-dataframes-on-columns
    pandaMerged = reduce(lambda left,right: pd.merge(left,right, on = 'coordinate'), arPandas)
    #print pandaMerged

    return pandaMerged

def printLengths(arPandas, pandaMerged):
    arIntLengths = []
    for panda in arPandas:
        strNameColumn = list(panda)[1]
        #print 'strNameColumn'
        #print strNameColumn
        intLength = panda.shape[0]
        arIntLengths.append(intLength)
        print 'The number of rows in table {0} is {1}'.format(strNameColumn, intLength)
    print 'The number of rows in the merged table is ' + str(len(pandaMerged))


def saveMergedPanda(pandaMerged):

    pandaMerged.to_csv('Merged_density_table.txt', sep='\t', index=False, na_rep='NaN')


def main(args):
    arFilenames = [line.strip() for line in args.file]
    arPandas = getIntoPandas(arFilenames)
    pandaMerged = makeOnePanda(arPandas)
    printLengths(arPandas, pandaMerged)
    saveMergedPanda(pandaMerged)

if __name__ == "__main__":
    args = get_args()
    main(args)
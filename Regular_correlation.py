"""
Script to accept a matrix of values where columns are different genome features (UCEs, exons, etc) and rows are
observations of these features in a particular bin of the genome, and return correlation coefficients and
associated p values. You specify one column to be correlated to any other set of columns.

Designed to accept matrices from Matrix_for_correlation.

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
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from scipy.stats import pearsonr
from scipy.stats import spearmanr
import argparse
from collections import OrderedDict
import seaborn as sns


def get_args():
    parser = argparse.ArgumentParser(description='Correlate a variable with column index -f with variables in'
                                                 'other columns with indexes -c, separate -c indexes with spaces')
    parser.add_argument('-m', '--matrix', type=str,
                        help='The filename of the matrix of data you want to correlate')
    parser.add_argument('-c', '--columns', type=int, nargs='+',
                        help='The indexes of columns you want to correlate that will be along the top of the '
                             'results table (the first column in the file is'
                             'column 0)')
    parser.add_argument('-r', '--rows', type=int, nargs='+',
                        help='The indexes of columns you want to correlate that will be the rows of the results table'
                             ' (the first column in the file is'
                             'column 0)')
    parser.add_argument('-t', '--type', type=str, default='p', choices=['p', 's'],
                        help='Specify the type of correlation, with p for Pearson, and s for Spearman, default Pearson')

    return parser.parse_args()


def importMatrix(strFileName):
    pdMatrix = pd.read_csv(strFileName, sep='\t')
    return pdMatrix


def orderedDictColumns(panda):
    #Get a dictionary of column names whose keys are their indexes, e.g. first column is key 0.
    arColNames = list(panda)
    arIndexes = range(len(arColNames))
    #From http://stackoverflow.com/questions/15372949/preserve-ordering-when-consolidating-two-lists-into-a-dict
    odColumns = OrderedDict(zip(arIndexes, arColNames))
    return odColumns


def getColumnNameList(pdMatrix, arIntColumnIndexes):
    odColumns = orderedDictColumns(pdMatrix)
    arXColumnNames = []
    for int in arIntColumnIndexes:
        strColumn = str(odColumns[int])
        arXColumnNames.append(strColumn)
    return arXColumnNames


def getFixedColumn(pdMatrix, intFixedColumnIndex):
    odColumns = orderedDictColumns(pdMatrix)
    strFixedColumnName = str(odColumns[intFixedColumnIndex])
    return strFixedColumnName


def getSlicedMatrix(pdMatrix, arIntColumns):
    #Get column names
    arStrColumnNames = getColumnNameList(pdMatrix, arIntColumns)
    #Slice matrix
    pdSliced = pdMatrix[arStrColumnNames]
    return pdSliced


def getCorrPval(pdMatrix, strFixedColumnName, arColumnNames, strType):
    print 'getCorrPval started'
    #Put fixed column values into numpy array
    npFixed = pdMatrix[[strFixedColumnName]].values

    #Iterate over the columns, getting Pearson correlation coefficient and pvalue

    #Pearsonr takes 2 lists of numbers and returns two values, the Pearson correlation coefficient and associated
    #p value.

    arStrColumnNames = []
    arFtCorr = []
    arFtPval = []

    for strColumnName in arColumnNames:
        print 'Correlating with {0} '.format(strColumnName)

        #Obtain the column as a numpy array
        npColumn = pdMatrix[[strColumnName]].values

        #Obtain correlation coefficient and p value
        if strType == 'p':
            npCorr, npPval = pearsonr(npFixed, npColumn)
        elif strType == 's':
            npCorr, npPval = spearmanr(npFixed, npColumn)

        ftCorr = float(npCorr)
        ftPval = float(npPval)

        #Append the column name, correlation coefficient, and p value
        arStrColumnNames.append(strColumnName)
        arFtCorr.append(ftCorr)
        arFtPval.append(ftPval)

    #Make pandas dataframes of results
    pdCorr = pd.DataFrame(arFtCorr, index=[arStrColumnNames], columns=[strFixedColumnName])

    #Transpose them so the index (row) is the fixed column and the columns are the other columns
    pdFinalCorr = pdCorr.transpose()

    pdPVal = pd.DataFrame(arFtPval, index=[arStrColumnNames], columns=[strFixedColumnName])
    pdFinalPval = pdPVal.transpose()

    return arStrColumnNames, arFtCorr, arFtPval


def rowIterator(arStrRowNames, arStrColumnNames, pdMatrix, strType):
    print 'rowIterator started'
    #Collect a list of the different rows' correlation coefficients
    arArFtCorr = []
    #Collect a list of the different rows' p-values
    arArFtPVal = []
    #Collect the row names
    arUsedStrRowNames = []

    for strRowName in arStrRowNames:
        #Get row data (it is a column in pdMatrix of course!
        npRowData = pdMatrix[[strRowName]].values
        #Put into correlation function
        arStrUsedColumnNames, arFtCorr, arFtPval = getCorrPval(pdMatrix, strRowName, arStrColumnNames, strType)
        arArFtCorr.append(arFtCorr)
        arArFtPVal.append(arFtPval)
        arUsedStrRowNames.append(strRowName)

        #Check all the column names (top of the results table) have been used in the correct order
        assert arStrUsedColumnNames == arStrColumnNames

    return arStrRowNames, arArFtCorr, arArFtPVal


def BuildResultsPandas(arStrRowNames, arStrColumnNames, arArFtCorr, arArFtPVal):
    #Get the results into nice pandas dataframes
    pdCorrResults = pd.DataFrame(arArFtCorr, index=[arStrRowNames], columns=arStrColumnNames)
    pdPValueResults = pd.DataFrame(arArFtPVal, index=[arStrRowNames], columns=arStrColumnNames)

    return pdCorrResults, pdPValueResults


def SaveResults(pdCorrResults, pdPValueResults, strType):
    if strType == 's':
        strCorrelationName = 'Spearman'
    elif strType == 'p':
        strCorrelationName = 'Pearson'
    pdCorrResults.to_csv(('{0}_coeff.txt'.format(strCorrelationName)), sep='\t')
    pdPValueResults.to_csv(('{0}_pval.txt'.format(strCorrelationName)), sep='\t')


def makeScatterPlots(pdMatrix, arColumnNames, arRowNames):
    for strColumnName in arColumnNames:
        for strRowName in arRowNames:
            drawScatter(pdMatrix, strColumnName, strRowName, (8, 5), 'g', True, 2)


def drawScatter(pdMatrix, strXColumnName, strYColumnName, tupWidthTallness, strSaveFormat, boolRegression, intDotSize):
    #Get x and y axis labels
    strXAxisLabel = strXColumnName
    strYAxisLabel = strYColumnName
    #Make figure
    plt.figure(figsize=tupWidthTallness)
    if not boolRegression:
        sns.regplot(data=pdMatrix, x=strXColumnName, y=strYColumnName, fit_reg=False, line_kws={'color': 'green'}, scatter_kws={'s':intDotSize})
        sns.axlabel(strXAxisLabel, strYAxisLabel)
    if boolRegression:
        sns.regplot(data=pdMatrix, x=strXColumnName, y=strYColumnName, fit_reg=True, line_kws={'color': 'green'}, scatter_kws={'s':intDotSize})
        sns.axlabel(strXAxisLabel, strYAxisLabel)
    #Save figure
    strFileName = strXColumnName + "_" + strYColumnName
    if strSaveFormat == 'p':
            strHistFilename = 'Scatter_{0}.pdf'.format(strFileName)
            plt.savefig(strHistFilename, format='pdf', bbox_inches='tight')
    elif strSaveFormat == 'g':
            strHistFilename = 'Scatter_{0}.png'.format(strFileName)
            plt.savefig(strHistFilename, bbox_inches='tight')
    #Close figure
    plt.close()



def main():
    #Get arguments
    args=get_args()
    #Get matrix
    pdMatrix = importMatrix(args.matrix)

    #Get indexes of columns to be correlated and to form columns along the top of the results table
    arIntColumns = args.columns
    #Get list of column names
    arStrColumnNames = getColumnNameList(pdMatrix, arIntColumns)
    print 'The columns of the results table will be ' + str(arStrColumnNames)

    #Get indexes of columns to be correlated an to form the rows of the results table
    arIntRows = args.rows
    #Get list of column names
    arStrRowNames = getColumnNameList(pdMatrix, arIntRows)
    print 'The rows of the results table will be ' + str(arStrRowNames)

    #Get the type of correlation
    strType = args.type
    if strType == 'p':
        print 'Pearson correlation requested'
    elif strType == 's':
        print 'Spearman correlation requested'

    #Get the results in list form:
    arStrRowNames, arArFtCorr, arArFtPVal = rowIterator(arStrRowNames, arStrColumnNames, pdMatrix, strType)

    #Put into pandas:
    pdCorrResults, pdPValueResults = BuildResultsPandas(arStrRowNames, arStrColumnNames, arArFtCorr, arArFtPVal)

    #Save the pandas:
    SaveResults(pdCorrResults, pdPValueResults, strType)

    #Make scatter plots
    makeScatterPlots(pdMatrix, arStrColumnNames, arStrRowNames)


if __name__ == "__main__":
     main()
"""
Script to take matrix of correlation coefficients and produce a heatmap.

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
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns
from collections import OrderedDict


def get_args():
    parser = argparse.ArgumentParser(description="Description")
    parser.add_argument("-c", "--corrFile", type=str,
                        help='The filename of the file containing the partial correlation coefficient matrix')
    parser.add_argument("-p", "--pValuesFile", type=str,
                        help='The filename of the file containing the p values corresponding to the correlation coefficient matrix')
    parser.add_argument("-t", "--top", type=int, nargs="+",
                        help='The indexes of columns you want to correlate that will be the columns of the heatmap'
                             ' (the first column in the file is'
                             'column 0)')
    parser.add_argument("-s", "--side", type=int, nargs="+",
                        help='The indexes of columns you want to correlate that will be the rows of the heatmap'
                             ' (the first column in the file is'
                             'column 0)')
    parser.add_argument("--annot", dest='boolAnnot', action='store_true',
                        help='Heatmaps will be annotated')
    parser.add_argument("--no-annot", dest='boolAnnot', action='store_false',
                        help='Heatmaps will not be annotated')
    parser.set_defaults(boolAnnot=True)
    parser.add_argument('-w', '--width', type=int, default=8,
                        help='Width of figure, default is 8')
    parser.add_argument('-l', '--length', type=int, default=6,
                        help='Length (aka height) of figure, default is 6')

    return parser.parse_args()


def importMatrix(strFileName):
    #Import the file into a panda
    pdMatrix = pd.read_csv(strFileName, sep='\t', header=0, index_col=0)
    return pdMatrix


def orderedDictColumns(panda):
    #Get a dictionary of column names whose keys are their indexes, e.g. first column is key 0.
    arColNames = list(panda)
    arIndexes = range(len(arColNames))
    #From http://stackoverflow.com/questions/15372949/preserve-ordering-when-consolidating-two-lists-into-a-dict
    odColumns = OrderedDict(zip(arIndexes, arColNames))
    return odColumns


def orderedDictRows(panda):
    #Get a dictionary of row names whose keys are their indexes
    #http://stackoverflow.com/questions/18358938/list-of-index-values-in-pandas-dataframe
    arRowNames = panda.index.values.tolist()
    arRowIndexes = range(len(arRowNames))
    odRows = OrderedDict(zip(arRowIndexes, arRowNames))
    return odRows


def getColumnNameList(pdMatrix, arIntColumnIndexes):
    odColumns = orderedDictColumns(pdMatrix)
    arXColumnNames = []
    for int in arIntColumnIndexes:
        strColumn = str(odColumns[int])
        arXColumnNames.append(strColumn)
    return arXColumnNames


def getRowNameList(pdMatrix, arIntRowIndexes):
    odRows = orderedDictRows(pdMatrix)
    arRowNameList = []
    for intRow in arIntRowIndexes:
        strRow = str(odRows[intRow])
        arRowNameList.append(strRow)
    return arRowNameList


def matrixColumnReorder(pdMatrix, arColOrder):
    pdColReorder = pdMatrix[arColOrder]
    return pdColReorder


def matrixRowReorder(pdMatrix, arRowOrder):
    pdRowReorder = pdMatrix.reindex(arRowOrder)
    return pdRowReorder


def matrixReorder(pdMatrix, arColOrder, arRowOrder):
    #Reorder both the rows and columns in one function
    pdColReordered = matrixColumnReorder(pdMatrix, arColOrder)
    pdColRowReordered = matrixRowReorder(pdColReordered, arRowOrder)
    return pdColRowReordered


def CorrHeatMap(pdCorrReordered, boolAnnot, tupFigSize):
    #Produce a heatmap for correlation coefficients from an array

    ###snsMapDivergecolors is for is you have negative values: positives will be red and negatives blue.
    plt.figure(figsize=tupFigSize)
    if boolAnnot==True:
        sns.heatmap(pdCorrReordered, linewidths=.5, center=0, annot=True, cmap='RdBu_r')
    elif boolAnnot==False:
        sns.heatmap(pdCorrReordered, linewidths=.5, center=0, annot=False, cmap='RdBu_r')
    plt.savefig('Corr_coeff_heatmap.pdf', format='pdf', bbox_inches='tight')


def pValHeatMap(pdPValReordered, boolAnnot, tupFigSize):
    plt.figure(figsize=tupFigSize)
    if boolAnnot:
        sns.heatmap(pdPValReordered, linewidths=.5, annot=True, cmap="Blues_r")
    else:
        sns.heatmap(pdPValReordered, linewidths=.5, annot=False, cmap="Blues_r")
    plt.savefig('Corr_PVals_heatmap.pdf', format='pdf', bbox_inches='tight')


def main():
    args = get_args()
    #Get matrixes
    pdCorr = importMatrix(args.corrFile)
    pdPVal = importMatrix(args.pValuesFile)

    #Get column order list
    arColIndexes = args.top
    
    #Get column names
    arColumnNames = getColumnNameList(pdCorr, arColIndexes)

    #Get row order list
    arRowIndexes = args.side

    #Get row names
    arRowNames = getRowNameList(pdCorr, arRowIndexes)

    #Get decision to annotate or not
    boolAnnot = args.boolAnnot

    #Get fig size
    intWidth = args.width
    intLength = args.length

    tupWidthLength = (intWidth, intLength)

    #Rearrange both matrices
    pdCorrReordered = matrixReorder(pdCorr, arColumnNames, arRowNames)
    pdPValReordered = matrixReorder(pdPVal, arColumnNames, arRowNames)

    print 'Subset of correlation coefficients table to be plotted'
    print pdCorrReordered
    print 'Subset of correlation p-values table to be plotted'
    print pdPValReordered


    CorrHeatMap(pdCorrReordered, boolAnnot, tupWidthLength)
    #Then need to clear the matplotlib figure to stop the second heatmap being plotted on top of first
    #http://stackoverflow.com/questions/36018681/stop-seaborn-plotting-multiple-figures-on-top-of-one-another
    plt.clf()

    pValHeatMap(pdPValReordered, boolAnnot, tupWidthLength)


if __name__ == "__main__":
     main()

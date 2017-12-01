"""
Script to take matrices and calculate partial correlations using matlab functions called using the matlab engine for
python.

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
import matlab.engine
from collections import OrderedDict
import numpy as np


def get_args():
    parser = argparse.ArgumentParser(description="Description")
    parser.add_argument("file", type=str,
                        help='The filename of the matrix you wish to process, columns separated by tabs')
    parser.add_argument("-y", "--YIndex", type=int, help = 'The column index (0-based) of the y variable')
    parser.add_argument("-x", "--XIndexes", type=int, nargs="+",
                        help = 'The column indexes of the x variables, separated by spaces')
    parser.add_argument("-z", "--ZIndexes", type=int, nargs="+",
                        help = 'The column indexes (0based) of the z variables, separated by spaces')

    return parser.parse_args()


def importMatrix(strFileName):
    pdMatrix = pd.read_csv(strFileName, sep='\t')
    return pdMatrix


def pandaColumnToList(panda, strColName):
    arCol = panda[strColName].tolist()
    return arCol


def twoPandaColumns(strcolumn1, strcolumn2, panda):
    #Returns a panda with same indexes and column headings, made of just the two columns named.
    pdTwoColumns = panda[[strcolumn1, strcolumn2]]
    return pdTwoColumns


def pdToNp(panda):
    #Turn a panda into a numpy array
    npArray = panda.values
    return npArray


def npToMatDouble(npArray):
    #Turn a numpy array into a matlab double,
    ### NB does not retain shape
    #From http://stackoverflow.com/questions/10997254/converting-numpy-arrays-to-matlab-and-vice-versa
    mlDouble = matlab.double(npArray.tolist())
    return mlDouble


def orderedDictColumns(panda):
    #Get a dictionary of column names whose keys are their indexes, e.g. first column is key 0.
    arColNames = list(panda)
    arIndexes = range(len(arColNames))
    #From http://stackoverflow.com/questions/15372949/preserve-ordering-when-consolidating-two-lists-into-a-dict
    odColumns = OrderedDict(zip(arIndexes, arColNames))
    return odColumns


def xyzColumnsToMlDoubleXYand_Z(pdMatrix, strXName, strYName, strZName):
    """
    From a panda matrix, specify the column names of the three columns you want to be your x, y, and z variables.
    Returns a n x 2 matrix of x and y values, and an n x 1 matrix of z values.
    If you want more than one column to go into z, this is not the right function.
    """

    #First make a panda of just the x and y columns
    pdXY = twoPandaColumns(strXName, strYName, pdMatrix)
    #Make it into a numpy array
    npXY = pdToNp(pdXY)
    #Make it a matlab double
    mlXY = npToMatDouble(npXY)

    #Now similar with z column:
    pdZ = pdMatrix[strZName]
    npZ = pdToNp(pdZ)
    #Reshape to get each value in its own row
    npZreshape = npZ.reshape(len(npZ),1)
    mlZ = npToMatDouble(npZreshape)

    return mlXY, mlZ


def partialCorr(mlXY, mlZ):
    """
    Run matlab partialcorr function, looking at correlation between pairs of values in mlXy, controlling for variables in z
    """
    eng = matlab.engine.start_matlab()
    partial, pval = eng.partialcorr(mlXY, mlZ, 'rows', 'pairwise', 'type', 'Spearman', nargout=2)
    #Partial and pval are nested lists, need to access the 1th element of the 0th list (2nd element of 1st list for human)
    intPartial = partial[0][1]
    intpval = pval[0][1]

    return intPartial, intpval


def partialCorrResults(intYColumnIndex, arIntXColumnIndexes, arIntZColumnIndexes, pdMatrix):
    """
    Returns two pandas, one with partial correlation coefficients, the other with p values.
    X variables are column headings. Z variables are column rows.
    """
    odColumns = orderedDictColumns(pdMatrix)
    #Get y variable name
    strYVariable = str(odColumns[intYColumnIndex])
    print 'Y variable is {0}'.format(strYVariable)

    #Get list of x variable names, these will be the column names for the results pandas:
    arXColumnNames = getColumnNameList(pdMatrix, arIntXColumnIndexes)
    print 'X variables are {0}'.format(str(arXColumnNames))
    #Get list of y variable names, these will be the row indexes for the results pandas:
    arZColumnNames = getColumnNameList(pdMatrix, arIntZColumnIndexes)
    print 'Z variables are {0}'.format(str(arZColumnNames))

    arArRawPartials = []
    arArRawPVals = []

    for x in arIntXColumnIndexes:
        #For each x variable, make a list of partial coefficients, and a list of pvals, one for each z variable
        arPartials = []
        arPVals = []
        for z in arIntZColumnIndexes:
            print 'Getting partial correlation between {0} and {1}, controlling for {2}'.format(odColumns[x], odColumns[intYColumnIndex], odColumns[z])
            mlXY, mlZ = xyzColumnsToMlDoubleXYand_Z(pdMatrix, odColumns[x], odColumns[intYColumnIndex], odColumns[z])
            intPartial, intPValue = partialCorr(mlXY, mlZ)
            print 'Partial correlation coefficient between {0} and {1}, controlling for {2} is {3}'.format(odColumns[x], odColumns[intYColumnIndex], odColumns[z], intPartial)
            print 'P value for partial correlation between {0} and {1}, controlling for {2} is {3}'.format(odColumns[x], odColumns[intYColumnIndex], odColumns[z], intPValue)
            arPartials.append(intPartial)
            arPVals.append(intPValue)


        #When all z variables are finished, append the newly created lists to a final list of lists.
        arArRawPartials.append(arPartials)
        arArRawPVals.append(arPVals)


    #Create the pretty pandas with columns as x variables, row indexes as z variables
    pdPrettyPartials = pd.DataFrame(data=arArRawPartials, columns=arZColumnNames, index=arXColumnNames)
    pdPrettyPVals = pd.DataFrame(data=arArRawPVals, columns=arZColumnNames, index=arXColumnNames)

    #Transpose pandas so that x variables become columns and yvariables become rows
    pdFinalPartials = pdPrettyPartials.transpose()
    pdFinalPVals = pdPrettyPVals.transpose()

    #Save those pretty pandas
    pdFinalPartials.to_csv(('Partial_coeff_columnsXvariables_rowsZvariables_yvariable{0}.txt'.format(strYVariable)), sep='\t')
    pdFinalPVals.to_csv(('Pvalues_columnsXvariables_rowsZvariables_yvariable_{0}.txt'.format(strYVariable)), sep='\t')

    return pdFinalPartials, pdFinalPVals


def getColumnNameList(pdMatrix, arIntColumnIndexes):
    odColumns = orderedDictColumns(pdMatrix)
    arXColumnNames = []
    for int in arIntColumnIndexes:
        strColumn = str(odColumns[int])
        arXColumnNames.append(strColumn)
    return arXColumnNames


def main():
    print 'Getting filename'
    args = get_args()
    print 'Getting matrix'
    pdMatrix = importMatrix(args.file)
    intYIndex = int(args.YIndex)
    arXIndexes = args.XIndexes

    arZIndexes = args.ZIndexes

    partialCorrResults(intYIndex, arXIndexes, arZIndexes, pdMatrix)


if __name__ == "__main__":
     main()
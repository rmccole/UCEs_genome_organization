"""
Script to print file about number of uces and genes within a domain stats

Wren Saylor
Feb 2018

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

In:
primary - uce file
secondary - file with list of domain filenames
tertiary - genes

Out:
stats file with the intersects; uce x domain, gene x domain, domain size
stats file for all the domains; uce size, all domain size, gene size

"""
import argparse
import pandas as pd
import pybedtools as pbt
from itertools import cycle
import numpy as np
from scipy import stats

# set args
def get_args():
	parser = argparse.ArgumentParser(description="Description")
	parser.add_argument("file",type=str,help='the primary element file') # UCEs
	parser.add_argument("-s","--secondaryfeatures",required=True,type=argparse.FileType('rU'),help="a file with a list of file names with the secondary features to query") # Domains
	parser.add_argument("-t","--tertiaryfeature",type=str,help="the tertiary elements file")# Genes
	return parser.parse_args()

# get bt features
def get_bedtools_features(strFileName):
	return pbt.BedTool(strFileName)

# convert bedtool to panda
def convert_bedtools_to_panda(btfeature):
	return pd.read_table(btfeature.fn,header=None)

# print the number of elements in the feature file
def print_number_of_elements(lengthfeatures,file):
	print 'there are {0} elements in {1}'.format(lengthfeatures,file)

# intersect a file by how many times a feature on b is in the interval
def intersect_bedfiles_c_true(afile,bfile):
	return afile.intersect(bfile,c=True)

def run_print_number_file_features(file):
	btfeatures = get_bedtools_features(file)
	pdfeatures = convert_bedtools_to_panda(btfeatures)
	print_number_of_elements(len(pdfeatures),file)
	return label_coordinate_columns(pdfeatures)

# coordinate labels and size
def label_coordinate_columns(pdfeature):
	pdfeature['Size(Kb)'] = pdfeature.loc[:,2].astype(int)-pdfeature.loc[:,1].astype(int)
	pdfeature.columns.values[0]='chr'
	pdfeature.columns.values[1]='start'
	pdfeature.columns.values[2]='end'
	pdfeature['Size(Kb)'] /= 1000.# convert to Kb
	return pdfeature

# create panda for overlap count datasets
def count_overlap_df(secondary,file,label):
	pdintersect = intersect_bedfiles_c_true(secondary,file)
	pdfeatures = convert_bedtools_to_panda(pdintersect)
	pdcoordinates = label_coordinate_columns(pdfeatures)
	pdcoordinates.columns.values[3]='intersect_{0}'.format(label)
	return pdcoordinates

# remove those without elements
def remove_rows_with_no_overlaps(overlaps,column):
	return overlaps[overlaps[column]!=0]

# get stats for list of columns
def panda_describe_multiple_column(pdfeature):
	intersectcols = [col for col in pdfeature.columns if 'intersect' in col]
	statcols=intersectcols
	return pdfeature[statcols].describe()

# select a list of columns where the column name contains a key word
def subset_column_by_keyword(pdfeature,key):
	return pdfeature.filter(regex=key)

# save panda to file with mode 'a' for appending
def save_panda(pdData,strFilename):
	pdData.to_csv(strFilename,sep='\t',index=True,mode='a')#

# get stats for single column
def panda_describe_single_column(pdfeatures,name):
	return pdfeatures[name].describe()

def main():
	args = get_args()
	pfile = args.file
	secondaryfiles = [line.strip() for line in args.secondaryfeatures]
	tfile = args.tertiaryfeature
	labelprimary = run_print_number_file_features(pfile)
	labeltertiary = run_print_number_file_features(tfile)
	lumpsec,lumpsecstats = [],[]
	for sfile in secondaryfiles:
		run_print_number_file_features(sfile)
		secondary = get_bedtools_features(sfile)
		pdprimary = count_overlap_df(secondary,pfile,'UCEs')
		pdtertiary = count_overlap_df(secondary,tfile,'Genes')
		concattotal = pdprimary.join(pdtertiary,rsuffix='_extra')
		concattotal.drop(columns=['chr_extra','start_extra','end_extra','Size(Kb)_extra'],axis=1,inplace=True)
		withuces = remove_rows_with_no_overlaps(concattotal,'intersect_UCEs')
		lumpsec.append(withuces)
		withdescribe = panda_describe_multiple_column(withuces)
		withdescribe.columns = [str(col) + '_{0}'.format(sfile) for col in withdescribe.columns]
		save_panda(withdescribe,'stats_intersect_domain_{0}.txt'.format(pfile))
		lumpsecstats.append(withdescribe)
	seccat = pd.concat(lumpsec)
	secdescribe = panda_describe_multiple_column(seccat)
	secdescribe.columns = [str(col) + '_All_Domains' for col in secdescribe.columns]
	save_panda(secdescribe,'stats_intersect_domain_{0}.txt'.format(pfile))

if __name__ == "__main__":
	main()
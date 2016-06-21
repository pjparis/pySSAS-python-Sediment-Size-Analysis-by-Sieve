#! /Users/paulp/anaconda/bin/python

# -*- coding: utf-8 -*-
"""
SedSASampleScript_2CSV.py

Script: sediment Grain-size Analysis application testing suite: 

Purpose: script calls SedSASample class method Analyze2CSV to process one or more
		 transect sample sets, writing the full results to a comma delimited text file.

User this script as a template

Created on: Thurs June 09 11:48:33 2016

@author: pjp
"""

# import the SedSAample class:
import sys
sys.path.append('/Users/paulp/Documents/projects/SedSAS/' )
import SedSASampleClass as sssc

# ###### USER INPUTS #################################################################
# 1.) Enter the absolute (full) path to the file(s) containing the data to be analyzed:
fp='/Users/paulp/Documents/Projects/SedSAS/'

# 2.) Enter the name of the file(s) containing the data to be analyzed:
fn='USFWS_survey07142015_transect.csv'

# 3.) Set the field delimiter used to separate data columns (fields) in the input file fn:
delim=','

# 4.) Set up the use columns list. Just modify the existing list as needed:
use_cols=[4,5,6,7,8,9,10,11,12,13,14,15,16]

# 5.) Enter the transect name/id (with the quotes around the name):
trList=['T2'] #['T1','T2','T3','T4','T5','T6','T7','T8','T9','T10','T11','T12','T13','T14','T15','T16']  #['T2']

# 6.) Set up the samples list. Modify the existing list as needed:
samples=['S1','S2','S3','S4']

# 7.) Set up the screens list. Modify the existing list as needed to reflect the sieves 
# used. Values represent phi screen units:
scrns=[-1.0,-0.5,0.0,0.5,1.0,1.25,1.5,1.75,2.0,2.5,3.0,3.5,4.0]

# ###### END USER INPUTS #############################################################



# ####################################################################################


# loop thru the list of transects in the trList and process the samples for each:
for tr in trList:

	inputs=(fp,fn,delim,use_cols,tr,samples,scrns)
	mySand=sssc.SedSASample(inputs)
	mySand.Analyze2CSV(tr+'.csv')
	
print('Fin!')
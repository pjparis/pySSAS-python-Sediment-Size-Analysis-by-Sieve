#! /Users/paulp/anaconda/bin/python

# -*- coding: utf-8 -*-
"""
SedSASApplTesting.py

Sediment Grain-size Analysis application testing suite: Calls SedSASample class 
methods to generate statistics and plots and then write these in a fixed report
format using ReportLab

NOTE: the ReportLab Python libraries are required to use the ReportModule. You can get
these from conda

Created on Tue May 17 11:31:33 2016

@author: paulp
"""
import math, sys
sys.path.append('/Users/paulp/Documents/projects/SedSAS/')
import SedSASampleClass as sssc
import SedSASReportModule as ssrm


# ###### USER INPUTS #################################################################
# 1.) Enter the absolute (full) path to the file(s) containing the data to be analyzed:
fp='/Users/paulp/Documents/projects/SedSAS/'

# 2.) Enter the name of the file(s) containing the data to be analyzed:
fn='USFWS_survey07142015_transect.csv'

# 3.) Set the field delimiter used to separate data columns (fields) in the input file fn:
delim=','

# 4.) Set up the use columns list. Just modify the existing list as needed:
use_cols=[4,5,6,7,8,9,10,11,12,13,14,15,16]

# 5.) Enter the transect name/id (with the quotes around the name):
trList=['T2'] 
# ['T1','T2','T3','T4','T5','T6','T7','T8','T9','T10','T11','T12','T13','T14','T15','T16']     

# 6.) Set up the samples list. Modify the existing list as needed:
samples=['S1','S2','S3','S4']

# 7.) Set up the screens list. Modify the existing list as needed to reflect the sieves 
# used. Values represent phi screen units:
scrns=[-1.0,-0.5,0.0,0.5,1.0,1.25,1.5,1.75,2.0,2.5,3.0,3.5,4.0]

# Input values for the report generator methods:
pname='USFWS - Pea Island National Wildlife Refuge'   # project name
smplDate='2/3/2016'                                   # sampling date

# ###### END USER INPUTS #############################################################



# loop thru the list of transects in the trList and process the samples for each:
for tr in trList:
    # construct argument tuple passed to SASample class instance from user inputs.
    inputs=(fp,fn,delim,use_cols,tr,samples,scrns)
    
    # create an SedSASample class instance for the current transect
    mySand=sssc.SedSASample(inputs)
    
    # if you're calling class methods directly (other than Analyze2CSV) you must 
    # call InterpolateQuantileValues() and return the ordered dictionary D:
    D=mySand.InterpolateQuantileValues()
    
    df=mySand.GetDataFrame()        # return a reference to data frame df
    
     # to generate the PDF archive report for the current transect call the report 
     # constructor:
    c=ssrm.SASReportConstructor(pname,tr)

    
    # to loop thru the samples list for each transect
    for s, Q in D.items():
    	# handle situation if the current transect sample set is missing from sample file
    	if( math.isnan( df[s].sum()) == True):  
        	print('Skipping sample:',tr,s, 'because of missing data...')
    	else:
        	gStats=mySand.ComputeGraphicStats(s,Q)     # compute graphic statistics
        	mStats=mySand.ComputeMomentStats(s)        # compute Method Moment statistics
        	modes=mySand.FindSampleModes(s)            # find mode(s)
        	
        	mySand.PLOTDualSampleWtPercents(s, mode='save')
        	
        	# generate a report page for the curent sample s:
        	ssrm.SASReportPageDispatcher(c, pname, smplDate, tr, s, df, D, gStats, mStats)
        	
    # close out and save the report for the current transect tr:
    ssrm.SASReportDestructor(c)
    
print('Fin!')


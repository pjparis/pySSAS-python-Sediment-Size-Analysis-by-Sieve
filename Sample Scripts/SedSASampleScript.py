#! /Users/paulp/anaconda/bin/python

# -*- coding: utf-8 -*-
"""
SedSASampleScript.py

Script: sediment Grain-size Analysis application testing suite: 

Purpose: script intantiates SedSASample class to process one or more
         transect sample sets, in script, with option to report results via console

User this script as a template

Created on: Thurs June 09 12:07:00 2016

@author: pjp
"""

# import the SedSAample class:
import math
import SedSASampleClass as sssc

# ###### USER INPUTS #################################################################
# 1.) Enter the absolute (full) path to the file(s) containing the data to be analyzed:
fp=' replace this text with the full (absolute) path to data file, less the file name '

# 2.) Enter the name of the file(s) containing the data to be analyzed:
fn=' replace this text with the full name of the source data file '

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

# ###### END USER INPUTS #############################################################


# loop thru the list of transects in the trList and process the samples for each:
for tr in trList:
    # construct argument tuple passed to SASample class instance from user inputs.
    inputs=(fp,fn,delim,use_cols,tr,samples,scrns)

    mySand=sssc.SedSASample(inputs)
    # if you're calling class methods directly (other than Analyze2CSV) you must 
    # call InterpolateQuantileValues() and return the ordered dictionary D:
    D=mySand.InterpolateQuantileValues()
    
    df=mySand.GetDataFrame()        # return a reference to data frame df
    
    # to loop thru the samples list for each transect
    for s, Q in D.items():
    	# handle situation if the current transect sample set is missing from sample file
    	if( math.isnan( df[s].sum()) == True):  
        	print('Skipping sample:',tr,s, 'because of missing data...')
    	else:
        	gStats=mySand.ComputeGraphicStats(s,Q)     # compute graphic statistics
        	mStats=mySand.ComputeMomentStats(s)        # compute Method Moment statistics
        	mySand.PrintGraphicStats(s,gStats)		   # write gStats to the console
        	mySand.PrintMomentStats(mStats)            # write mStats to the console
             
        	modes=mySand.FindSampleModes(s)            # find mode(s)
        	mySand.PrintSampleModes(modes)

        	# uncomment the next line to draw a histogram, PDF, and CDF plots
        	# mode=print prints plot to console window; mode=save saves plot to png file
        	#mySand.PLOTDualSampleWtPercents(s, mode='print') 
             
print('Fin!')
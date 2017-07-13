#! /Users/paulp/anaconda/bin/python

# -*- coding: utf-8 -*-
"""
SedSASampleScript.py

Script: sediment Grain-size Analysis application testing suite: 

Purpose: script intantiates SedSASample class to process one or more
         transect sample sets, in script, with option to report results via console

User this script as a template

Created on: Thurs June 09 12:07:00 2016
Modified to reflect changes/upgrades to the SedSAS class on: July 12, 2017

@author: pjp
"""

# ###### USER INPUTS #################################################################
# 1.) Enter the absolute (full) path to the file(s) containing the data to be analyzed:
fp=' replace this text with the full (absolute) path to data file, less the file name '

# 2.) Enter the name of the file(s) containing the data to be analyzed:
fn=' replace this text with the full name of the source data file '

# 3.) Set the field delimiter used to separate data columns (fields) in the input file fn:
### comma is the default
delim=','

# 4.) Enter the transect name/id (with the quotes around the name):
trList=['T2'] 
# ['T1','T2','T3','T4','T5','T6','T7','T8','T9','T10','T11','T12','T13','T14','T15','T16']     

# 6.) Set up the samples list. Modify the existing list as needed:
samples=['S1','S2','S3','S4']

# 7.) Set up the screens list. Modify the existing list as needed to reflect the sieves 
# used. Values represent phi screen units:
scrns=[-1.0,-0.5,0.0,0.5,1.0,1.25,1.5,1.75,2.0,2.5,3.0,3.5,4.0]

# ###### END USER INPUTS #############################################################


# import prerequisites:
import sys, math
import pandas as pd

# adding the sys.path variable to point to the SedSASClass.py file:
sys.path.append('./')

# import SedSAS class:
import SedSASClass

# read input file, create, populate, and wrangle data frame df:
# Note that the fields (columns) dropped here may not match those for the data you're
# currently working with. Adjust these accordingly, comment the line out if no fields
# require removal.
df = pd.read_csv(fp+fn)
df=df.drop(['Pan Weight','Wet Sample Weight'], axis=1)


# loop thru the list of transects in the trList and process the samples for each:
for tr in trList:
    mySand=SedSASClass.SedSAS(df, tr, samples, scrns)
    # if you're calling class methods directly (other than Analyze2CSV) you must 
    # call InterpolateQuantileValues() and return the ordered dictionary D:
    D=mySand.InterpolateQuantileValues()
    
    dfo=mySand.GetDataFrame()        # return a reference to data frame df
    
    # to loop thru the samples list for each transect
    for s, Q in D.items():
    	# handle situation if the current transect sample set is missing from sample file
        # note that missing transects are handled internally in the class (7/12/2017),but
        # this check does no harm.
    	if( math.isnan( dfo[s].sum()) == True):  
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
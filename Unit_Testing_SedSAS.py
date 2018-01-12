#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug  3 12:03:15 2017

@author: paulp
"""


import sys
import numpy as np
import pandas as pd

# adding the sys.path variable to point to the SedSASClass.py file:
sys.path.append('./')

# import SedSAS class:
import SedSASClass_II

###############################################################################
# user inputs:
fp='./misc/'
fn='FWSGrainSizeAnalysis2016_10.xlsx'
delim=','
samples=['S1','S2','S3','S4','S5']  
scrns=[-1.0,-0.5,0.0,0.5,1.0,1.25,1.5,1.75,2.0,2.5,3.0,3.5,4.0]

###############################################################################

# load the first sheet in the Excel ss file:
df = pd.read_excel(fp+fn)

# since we only need a subset of the total file content, and that content is
# not contiguous, we extract what we want/need in two parts, then combine the
# two together.
df1 = df.iloc[:,2:4]
df2 = df.iloc[:,7:21]     # added capture of reported dry weight
df = pd.concat((df1,df2), axis=1)

df=df.rename(columns = {'transect_id':'Transect', 'sample_number':'Sample' })

trList=['T1']   
for tr in trList:
#    try:
    s=samples[0]
    sc = SedSASClass_II.SedSAS(df, tr, s, scrns)
    dfp=sc.GetDataFrames()
    #gs=sc.ComputeFWLogarithmicGraphicStats(s )
    #sc.PrintComputedStatistics(s, gs, units='phi', method='FWlog' )

    #modes=sc.FindSampleModes(s)
    #sc.PrintSampleModes(s, modes)
    #print(dfp[0])
    
    #print(sc.ComputeArithmeticMethodofMomentsStats(s) )
    #print(sc.ComputeGeometricMethodofMomentsStats(s) )
    #print(sc.ComputeLogarithmicMethodofMomentsStats(s) )
    #fwls = sc.ComputeFWGeometricGraphicStats(s)
    #sc.PrintComputedStatistics(s, fwls, units='mm', method='FWgeo' )
    
    sc.Summary(s)
    
    
    #dfp[0]['midpt_phi']=(dfp[0]['phi'][1:] + dfp[0]['phi'][:-1]) / 2
    
# =============================================================================
#     ## compute and add the class size midpoints in phi units
#     x=np.array( dfp[0]['phi'])
#     mp=(x[1:] + x[:-1]) / 2
#     dfp[0]['midpt_phi']=np.append(dfp[0]['phi'][0:1].values-0.5, mp)    
# 
#     ## using the phi midpoints compute and add class size midpoints in mm
#     dfp[0]['midpt_mm']=2**(dfp[0]['midpt_phi']*-1) 
# =============================================================================
    
    ## compute the mean by airthmetic method of moments:
    #x_bar=( dfp[0][s+'wp']*dfp[0]['midpt_phi'] ).sum() / dfp[0][s+'wp'].sum()

#df=df.loc[ (df['Transect'] == tr) & (df['Sample'] == s) ].copy().T
#df=df.drop(['Transect','Sample'], axis=0)
#df.columns=[s]

#df[s+'wp']=(df/df.sum(axis=0))*100 
#df[s+'cwp']=df[s+'wp'].cumsum()
#df['phi']=scrns
#    except:
        #print('Check df for current transect %s recs. Or, hunt for other errors (maybe in SedSAS???)' %tr )

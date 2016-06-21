# -*- coding: utf-8 -*-
"""

SedSASReportModule.py  - with inputs from the methods in SedSASampleClass
this module creates a printable pdf file for printing, archive, etc.

To use the methods contained first load the module to your script, call the
method SASReportConstructor() to initiate, and then call the 
SASReportPageDispatcher() to do build each report page. Call SASReportDistructor()
to clean things up after each transect. 

METHODS to call:
    SASReportConstructor(pname, tr)
        Arguments:
        - pname=the name of the current project as a Python string
        - tr=the current transect's id or designator, as a Python string
        Returns:
        - c=a Report library canvas object (c)
        >>>Call this method once for each transect
        
    SASReportPageDispatcher(c, pname, smplDate, tr, s, df, gstats, mstats)
        Arguments:
        - c=the ReportLib canvas object returned from SASReportConstructor
        - pname=the name of the current project as a Python string
        - smplDate=date of sample as a string
        - tr=current transect's id as a Python string
        - s=current sample id as a Python string
        - df=pandas data frame returned from SedSASampleClass method: GetDataFrame()        
        - gstats=graphic statistics for current sample s as a Python tuple
        - mstats=method of moment statistics for current sample as a Python tuple
        Returns: nothing
        >>>Call this method once for each sample, where there are n samples per transect
        
    SASReportDestructor(c)
        Arguments:
        - c=the ReportLib canvas object
        Returns: the saved report as a PDF document
        
    
Created on Tue May 17 12:03:44 2016

@author: paulp
"""

# load requisite libraries:
import pandas
import numpy as np

# ## [ReportLab] libraries needed for external report generation
from reportlab.pdfgen import canvas
from reportlab.lib.units import inch, cm
from reportlab.lib.pagesizes import letter
from reportlab.platypus import *


# ## ######################################################################################
# ## SAS Report Generator Functions - Creates a printable pdf file for printing, archive, etc. 
# ## ######################################################################################
# ## SASReportConstructor:
def SASReportConstructor(pname,tr):
    file=pname+'_'+tr+'.pdf'               # set the pdf file name : pData[4] refs the transect id
    return(canvas.Canvas(file, pagesize=letter))    # creates (instantiates) the canvas object c
    
        # some parts to use for the ReportGenerator class, if it should ever come to be...
        #self.project=projectData[0]             # project name, will be printed at top of any generated reoprts
        #self.collectionDate=projectData[1]      # date that the sediment sample(s) was/were collected
        #self.analDate=projectData[2]            # date that either the samples were processed or otherwise analyzed
        
        
# ## ######################################################################################
# ## SASReportDispatcher

def SASReportPageDispatcher(c, pname, smplDate, tr, s, df, D, gstats, mstats):   
    textObject=c.beginText()
    textObject.setTextOrigin(21, 770 )
    textObject.setFont("Helvetica", 14)
    textObject.textLine('Sediment Grain Size Analysis Report')
    textObject.setFont('Helvetica', 11)
    textObject.textLine('    Project: '+pname )
    #textObject.textLine('')
    textObject.setFont("Helvetica", 10)
    textObject.textLine('    Date Sample Collected: '+smplDate )
    textObject.textLine('    Transect: '+tr+'     '+'Sample: '+s )
    textObject.textLine('-'*135)
    #textObject.textLine('Sample Grain-Size Data Table                       Graphic Statistics (Folk, 1980)')
    c.drawText(textObject)
    #c.drawString(21,760, '-'*135)

    # call functions that build the report sheet for each sample location
    SASReportDrawSampleSizeDataTable( c,s,df )
    SASReportDrawGraphicStatisticsTable(c, gstats)
    SASReportDrawMomentStatisticsTable(c, mstats)
    SASReportDrawQuantilesList(c, D, s)
    SASReportDrawPlots(c, s, tr)
    SASReportDrawFootnotes(c)
    
    c.showPage()


# ## ######################################################################################
# ## SASReportDestructor
def SASReportDestructor(c):
    c.drawString(270,100,'End of Report')
    c.save()


# ## ######################################################################################
# ## SASReportDrawSampleSizeDataTable
def SASReportDrawSampleSizeDataTable(c,s,df):
    # ## draw the sediment sample grain size data table
    textObject=c.beginText()
    textObject.setFont('Helvetica', 11)
    c.drawString(21,695,'Sample Grain-Size Data Table')

	# create tmp data frame table1df and populate with contents of df
    table1df=pandas.DataFrame( df['phi'])
    table1df['Raw Wt (gm)']= df[s].astype(np.double).round(2) 
    table1df['Wt.%']=df[s+'wp'].astype(np.double).round(2) 
    table1df['Cum. Wt.%']=df[s+'cwp'].astype(np.double).round(2)
        
    #textObject.setTextOrigin(21, 770)
    lista = [table1df.columns[:,].values.astype(str).tolist()] + table1df.values.tolist()
    t1=Table(lista)

    t1.wrapOn(c, 49,49)
    t1.drawOn(c, 21,435)


# ## ######################################################################################    
# ## SASReportDrawGraphicStatisticsTable
def SASReportDrawGraphicStatisticsTable(c, gstats):
    # ## draw the grain size graphic statistics table
    textObject=c.beginText()
    textObject.setFont('Helvetica', 11)
    c.drawString(300,695,'Graphic Statistics (Folk, 1980):')

    textObject.setFont('Helvetica', 10)
    if np.isnan(gstats[0]):
        c.drawString( 320,670,'Graphic Mean:  Out of Range*')
    else:
        c.drawString(320,670,'Graphic Mean:  '+str(round(gstats[0],2))+' phi' )
    if np.isnan(gstats[1]):
        c.drawString(320,650,'Graphic Median (D50):  Out of Range*')
    else:
        c.drawString(320,650,'Graphic Median (D50)   '+str(round(gstats[1],2))+' phi' )
        
    if np.isnan(gstats[2]):
        c.drawString(320,630,'Graphic Standard Deviation:   Out of Range*')
    else:
        c.drawString(320,630,'Graphic Standard Deviation   '+str(round(gstats[2],2))+' phi' )        

    if np.isnan(gstats[3]):
        c.drawString(320,610,'Graphic Skewness:   Out of Range*')
    else:
        c.drawString(320,610,'Graphic Skewness   '+str(round(gstats[3],2)) )
           
    if np.isnan(gstats[4]):
        c.drawString(320,590,'Graphic Kurtosis:   Out of Range*')
    else:
        c.drawString(320,590,'Graphic Kurtosis   '+str(round(gstats[4],2)) )    
    
    

    
# ## ######################################################################################
# ## SASReportDrawMomentStatisticsTable
def SASReportDrawMomentStatisticsTable(c, mstats):
    # ## draw the grain size moment statistics table
    textObject=c.beginText()
    textObject.setFont('Helvetica', 11)
    c.drawString(300,550,'Moment Statistics (Folk, 1980):')

    textObject.setFont('Helvetica', 10)
    c.drawString(320,525,'Mean:  ' + str(round(mstats[0],2))+' phi' )
    c.drawString(320,505,'Standard Deviation:  '+ str(round(mstats[1],2))+' phi' )
    c.drawString(320,480,'Coarse Fraction, undifferentiated: '+ str(round(mstats[2],2))+ ' gm')
    c.drawString(320,460,'Pan Fraction (fines), undifferentiated: '+ str(round(mstats[3],2))+ ' gm')
        
    textObject.setFont('Helvetica', 10)

    
    
# ## ######################################################################################
# ## SASReportDrawQuantilesList
def SASReportDrawQuantilesList(c, D, s):
    # ## draw the quantiles list:
    c.drawString(21,395,'Sample Quantiles (D5, D10, D16, D25, D50, D75, D84, D90, D95) in phi:')
    t2=Table([D[s]])
    t2.wrapOn(c,500,20)
    t2.drawOn(c, 21,375)

    

# ## ######################################################################################
# ## SASReportDrawPlots
def SASReportDrawPlots(c, s, tr):
    c.drawImage(tr+'_'+s+'_dual.png', 7, 90, 21*cm, 10*cm)
    # ## draw the weight percent histogram and frequency curve to the report sheet:
    #c.drawImage(s+'_wf.png', 21, 90, 9*cm, 9*cm)
    # ## draw the cumulative weight curve to the report sheet:
    #c.drawImage(s+'_cwf.png', 300, 90, 9*cm, 9*cm)   
    # delete the temp files to either make way for the next new set, or to clean things up in the end


   

 # ## ######################################################################################
# ## SASReportDrawFootnotes
def SASReportDrawFootnotes(c):
    textObject=c.beginText()
    textObject.setFont('Helvetica', 9)
    c.drawString(21, 49, 'Folk, R. 1980. Petrology of Sedimentary Rocks. Hemphill Publishing Company, Austin, TX. 184p.')
    c.drawString(21,28,'* out of range values indicate that one or more quantiles fall outside the screen range along the')
    c.drawString(28,15, 'cumulative frequency curve.')
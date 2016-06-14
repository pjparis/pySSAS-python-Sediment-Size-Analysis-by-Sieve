__SedSASample Class Description:__


Conduct a classic grain size analysis on sediment samples partitioned using mechanical sieves. Included class SedSASample computes graphic and method of moment descriptive statistics, along modes, for samples. Results can be written out to the console or written to a comma-separated format text file. Included module SedSASReportModule provides means to write results to a persistent, printable PDF report.

What does SedSASample do?
- computes graphic mean, median, standard deviation, skewness and kurtosis for a sample
- compute method of moments mean and standard deviation for a sample
- locates mode(s)
- prints statistical results to console window
- prints statistical results to a comma-separated format text file


Requisite Libraries:
the SedSASample class will attempt to load: math, numpy, pandas, matplotlib, collections, and pprint at instantiation. 

imported as:
import math
import numpy as np
import pandas
import matplotlib.pyplot as plt
import collections
import pprint


Notes:
1.) SedSASample was written and tested in a Python 3.x environment. It should work under 2.7.x conditions, but do test before proceeding with production work.

2.) SedSASample doesn't [yet] do much error handling and so, if something goes wrong expect to see a jumble of obscure tracebacks and other text that may or may not make much sense. Error handling is being folded into the class, but it may never be implemented sufficiently to cover every possible event. 

2.) If you do not already have Python installed (Windows users), or only have the core interpreter (Linux, Unix, and OSX users), suggest installing Continuum Analytic's Anaconda distribution. Anaconda contains a really nice softare application package manager called conda that makes it easy to install and updated a Python environment, or environments. Get it here: https://www.continuum.io/downloads



Included sample scripts:
- SedSASampleScript.py  (example calling individual class functions)
- SedSASampleScript.2CSV.py  (example calling auto process to text file)
- SedSASampleScript_Reporting.py (example using Report Lab to create PDF report)




Input tuple Arguments:
Input arguments are fed to SedSASample via a simple tuple that contains the following data:

index     value
	0		absolute path to the source data file, less the name of the file itself
	1		name of the source data file, including the extension
	2		the field delimiter used in the source file (usually a comma ',')
	3		a list of columns, by position, that exist in the source file that will be 
				used in the analysis. These are usually the individual sample weights. 
				The position is an integer value that identifies the position of the 
				column (data field) in the source file. The columns are numbered (indexed)
				from left to right beginning with 0. So, the first column in the file is
				column 0, the second is column 1, and so on. For example:
				
					to use columns 4 thru 16 the list will look like this:
					
					[4,5,6,7,8,9,10,11,12,13,14,15,16]
	4		the transect id. The transect identification number or value for the 
			current sample set. 
	5		list of sample site ids. The unique id numbers or values for each of the n
			sediment samples associated with the current transect
	6		a list of the sieve screen sizes used in the analysis. For example, if you
			used sieves with phi sizes:
			
			-1.0,-0.5,0.0,0.5,1.0,1.25,1.5,1.75,2.0,2.5,3.0,3.5,4.0
			
			the list is simply:
			
			[-1.0,-0.5,0.0,0.5,1.0,1.25,1.5,1.75,2.0,2.5,3.0,3.5,4.0]
			

Note that the sample scripts contains pre-prepared inputs and assembles these into the requisite tuple, so one need only update the values in the sample scripts to match
	your analysis. 
	
						
					

Private Methods:
__init__
	- parses the input tuple, gathering the tuple's data into six internal variables. 
	- performs some internal computations that prepare the inputs, and ultimately creates the pandas data frame, that will hold the all the data for analysis. The two most salient operations are the computation of the weight percentages and cumulative weight percentages for the sediment sample. Instantiating the class automatically completes these computation with no use intervention necessary. 

Inputs: the input tuple described above

Returns: none


PYDOC DUMP:

NAME
    SedSASampleClass - # -*- coding: utf-8 -*-

CLASSES
    builtins.object
        SedSASample
    
    class SedSASample(builtins.object)
     |  Class SedSASample - Sediment Size Analysis by Sieve
     |  
     |  DESCRIPTION:
     |  Performs various data extractions, manipulations, statistical analyses, and 
     |  plotting to characterize [sand and gravel sized] sediment particle size 
     |  distributions from mechanical dry sieve fractionation.
     |  
     |  From a sample file of sieved sand (ϕ) fractional weights, SedSASample() can:
     |  
     |   - compute weight percentages and cumulative weight percentages for each ϕ
     |     fraction relative to the total sample
     |  
     |   - compute the mean, median, and sorting (standard deviation) in  ϕ units using  
     |     the graphic method described by Folk (1980) the mean and sorting (standard 
     |     deviation) in ϕ units using the method of moments as described by Folk (1980)
     |  
     |   - create and plot a histogram with frequency curve (PDF) overlay
     |  
     |   - create and plot a cumulative frequency curve (CDF) 
     |  
     |  
     |  NOTES: 
     |      1.) SedSASample requires the following Python libraries (they load upon 
     |          instantiation): 
     |              numpy as np
     |              pandas
     |              matplotlib.pyplot as plt
     |              collections
     |  
     |      2.) SedSASample was written using, and for, Python 3.x. It has not been tested
     |          in a 2.x environment, but would be expected to function with little or no
     |          modifications required.
     |  
     |      3.) Cited reference: Folk, R. 1980. Petrology of Sedimentary Rocks. Hemphill 
     |          Publishing Company, Austin, TX. 184p.
     |  
     |  
     |  HISTORY:
     |      - initial code assembly: April/May, 2016  (pjp)
     |      - method FindSampleModes() added 30 May, 2016   (pjp)
     |      - removed pprint statement from method PrintSampleWeightsDataTable()                        
     |           ## pprint.pprint( self.quantilesDict[s] )  1 June, 2016 (pjp)
     |      - added check on undifferntiated coarse and fine fractions 8 June, 2016 (pjp)
     |      - added method PrintSamplesModes  9 June, 2016 (pjp)
     |  
     |  Methods defined here:
     |  Analyze2CSV(self, fn)
     |      performs internal processing (without need for the user to call the individual 
     |      methods explicitly/directly) to compile the raw data, weight percentages, 
     |      cumulative weight percentages, graphic and method of moment statistics, and
     |      the mode(s). All of these are written directly to an output file (user
     |      specified path and file name) formatted as comma separated values. 
     |      
     |      Input args: 
     |          fn = user-defined output file name
     |      
     |      Returns: none
     |  
     |
     |
     |  Analyze2CSVHeader(self)
     |     generates a columns header string (names of each data column written by method  
     |     Analyze2CSV) for output data file created in method Analyze2CSV
     |	
     |		Input Args: none
     |		
     |		Returns: header as a Python string
     |
     |
     |
     |  ComputeGraphicStats(self, s, Q)
     |      computes "graphic" mean, median, standard deviation, skewness, and kurtosis 
     |      (per Folk, 1980) for a given transect sample.
     |      
     |      Input args:
     |          s  sample identifier as a Python string (EX. 'S1')
     |          Q  list of interpolated quantile values (9 items) for current sample
     |          
     |      Returns:
     |          Python tuple containing MN=mean,MD=median,SD=std. dev.,SK=skewness, 
     |            K=kurtosis in phi units.
     |  
     |
     |
     |  ComputeMomentStats(self, s)
     |      computes the sediment sample mean and standard deviation using the arithmetic 
     |      method of moments method as described by Folk (1980) for current sample.
     |      
     |      Input args:
     |          s = sample identifier as a Python string (EX. 'S1')
     |          
     |       Returns:
     |           Python tuple containing momMean=mean,momSD=standard deviation,
     |           ,undifCF=undifferentiated coarse fraction (the weight of material
     |              captured by the coarsest sieve)
     |            undifFF=undifferentiated fine fraction (the weight of material captured 
     |              by the residual pan at the bottom of the seive stack).
     |  
     |
     |
     |  FindSampleModes(self, s)
     |      locates any modal weight percentage values in current sample. 
     |      
     |      Input args: 
     |          s = sample identifier as a Python string (EX. 'S1)
     |      
     |      Returns: Python ordered dictionary of modes (if any) found in the current 
     |                           sample.
     |  
     |
     |
     |  GetDataFrame(self)
     |      returns the sample set data frame df to the caller.
     |      
     |      Input args: none
     |      
     |      Returns: pandas dataframe for the current sample data
     |  
     |
     |
     |  InterpolateQuantileValues(self)
     |      interpolate quantiles in self.quantilesList for each transect sample.
     |      
     |      Input args: none
     |      
     |      Returns:
     |          an ordered dictionary (key=sample identifier; value=interpolated quantile 
     |          list)
     |  
     |
     |
     |  PLOTDualSampleWtPercents(self, s, mode='print')
     |      plots both the weight percentage and cumulative weight percentage histogram 
     |      and curves (histo+PDF, and CDF) side by side and together, by sieve fraction 
     |      for the current transect sample. Can plot to console or to stored PNG file.
     |      
     |      Input args:
     |          s = sample identifier as a Python string (EX. 'S1)
     |          mode = plot destination: mode='print' plot written to console (default)
     |          mode='save' plot written to png output
     |                                         
     |      Returns: none
     |      
     |      Note: 
     |          This is just a convenient combination of methods PLOTSampleWtPercents and 
     |          PLOTSampleCumWtPercents.
     |  
     |
     |
     |  PLOTSampleCumWtPercents(self, s, mode='print')
     |      plots the individual sample cumulative weight percentages, by sieve 
     |      fraction as a cunulative frequency (CDF) curve for the current transect 
     |      sample. Can plot to console or to stored PNG file.
     |      
     |       Input args:
     |           s = sample identifier as a Python string (EX. 'S1)
     |           mode = plot destination: mode='print' plot written to console (default)
     |           mode='save' plot written to png output
     |                                          
     |       Returns: none
     |  
     |
     |
     |  PLOTSampleWtPercents(self, s, mode='print')
     |      plots the individual sample weight percentages, by sieve fraction 
     |      as a histogram and overprinted frequency (PDF) curve for the current sample. 
     |      Can plot to console or to stored PNG file.
     |      
     |      Input args:
     |          s = sample identifier as a Python string (EX. 'S1)
     |          mode = plot destination: mode='print' plot written to console (default)
     |          mode='save' plot written to png output
     |                                         
     |      Returns: none
     |  
     |
     |
     |  PrintGraphicStats(self, s, gStats)
     |      prints the graphic mean, median, standard deviation, skewness, and kurtosis 
     |      for each transect sample in tabular format to the console. 
     |      
     |      Input args:
     |          s = sample identifier as a Python string (EX. 'S1)
     |          gStats = Python tuple containing the graphic mean, median, std. dev.
     |                          skewness and kurtosis.
     |          
     |      Returns: none
     |      
     |       Note: the gStats tuple is generated by method self.ComputeGraphicStats()
     |  
     |
     |
     |  PrintMomentStats(self, mStats)
     |      prints the method of moments computed sediment size mean and standard 
     |      deviation and undifferentiated coarse and fine fractions for the sample in 
     |      tabular format to the console.
     |      
     |      Input args: 
     |          mStats = Python tuple containing method of moments mean, std. dev., 
     |                   undiff. coarse fraction weight, undiff fine fraction weight.
     |          
     |      Returns: none
     |      
     |      Note: the mStats tuple is generated by method self.ComputeMomentStats()
     | 
     |
     |
     |    PrintSampleModes(self, modes)
     |      Prints the mode values located in method FindSampleModes to the console
     |      
     |      Input args:
     |          modes = Python dictionary containing key (mode in phi units) and value
     |                  (mode in weight percent) pairs
     |                  
     |      Returns: none
     |  
     |
     |
     |  PrintSampleWeightsDataTable(self, s)
     |      prints the raw, weight percentage, and cumulative weight percentage values in 
     |      tabular format for each screen bin to the console.
     |      
     |      Input args: 
     |          s = sample identifier as a Python string (EX. 'S1)
     |          
     |      Returns: none
     |  
     |
     |
     |  ReturnQuantile(self, s, d)
     |      uses simple linear interpolation to simulate the graphical estimation of a
     |      single quantile value. Quantile values are required to compute sediment sample 
     |      statistics, in the style of Folk (1980). See self.QuantilesList for more.
     |      
     |      Input Args:
     |          s  sample identifier  as a Python string (EX. 'S1')
     |          d  quantile value to be interpolated as integer
     |          
     |       Returns:
     |           the interpolated quantile value in phi units
     |           
     |      NOTE this is an internal function called by InterpolateQuantileValues(). 
     |      The user can call this method, but this is discouraged...I don't know why.
     |  
     |  __init__(self, userInputs)
     |      Initialize self.  See help(type(self)) for accurate signature.


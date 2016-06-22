## SedSASample Python Class:


The Python class SedSample contains a series of methods that permit the researcher to carry out a classic grain size analysis on sediment samples partitioned using mechanical sieves. Included class SedSASample computes graphic and method of moment descriptive statistics for sieved samples. Results can be written out to the console or written to a comma-separated format text file. Plotting methods allow for the generation of histogram, PMF, and CDF plots to the console or as PNG formatted graphics files. An external Python module SedSASReportModule (included in this repo) provides means to write results to a persistent, printable PDF report.


**Requisite Libraries:**
the SedSASample class will attempt to load: math, numpy, pandas, matplotlib, and collections at instantiation. 

imported as:<br>
import math<br>
import numpy as np<br>
import pandas<br>
import matplotlib.pyplot as plt<br>
import collections<br>
import pprint<br>


**Notes:**
1. SedSASample was written and tested in a Python 3.x environment. It should work under 2.7.x conditions, but do test before proceeding with production work.

2. SedSASample doesn't [yet] do much error handling and so, if something goes wrong expect to see a jumble of obscure tracebacks and other text that may or may not make much sense. Error handling is being folded into the class, but it may never be implemented sufficiently to cover every possible event. Ask questions of those knowledgable around you, as necessary.

3. If you do not already have Python installed (Windows users), or only have the core interpreter (Linux, Unix, and OSX users), suggest installing Continuum Analytic's Anaconda distribution. Anaconda contains a really nice softare application package manager called conda that makes it easy to install and updated a Python environment, or environments. Get it here: https://www.continuum.io/downloads


**Included sample scripts:**
- SedSASampleScript.py  (example calling individual class functions)
- SedSASampleScript.2CSV.py  (example calling auto process to text file)
- SedSASampleScript_Reporting.py (example using Report Lab to create PDF report)



**Input Arguments (that the user must provide):**
Input arguments are fed to SedSASample via a simple Python tuple that contains the following data:

1. a path (absolute) to the source data files, less the name of the file itself, as a Python string
2. the name of the source data file as a Python string
3. the field delimiter used to separate the data columns in the source file, as a Python string (Ex: delim=','  or delim='|' )
4. a Python list of columns by position in the source file that will be read into the class data frame and used for the analysis. These are usually the individual sieved sample weights. Position is an integer value that identifies the position of the 
	column (data field) in the source file. The columns are numbered (indexed)
	from left to right beginning with 0. So, the first column in the file is
	column 0, the second is column 1, and so on. For example:
				
		to use columns 4 thru 16 the list will look like this:
					
			[4,5,6,7,8,9,10,11,12,13,14,15,16]
			
5. a Python list of one or more transect identifiers, as Python strings. For example:

		transect T1 is entered as:  trlist=['T1']
		transects T1 thru T5 are entered as trlist=['T1','T2','T3','T4','T5']
		
6. a Python list of along-transect sample site ids. Each id is entered in the list as a Python string. For example:

		for a given transect there may be four sand samples collected: S1 thru S4. Enter these as: samples=['S1','S2','S3','S4']
					
7. a Python list of sieve screen (mesh) sizes used in the analysis. Screen sizes are included in the list as Python floating point $phi$ values. For example, if you use an array of sieves sized as: -1.0,-0.5,0.0,0.5,1.0,1.25,1.5,1.75,2.0,2.5,3.0,3.5,4.0
			
			the list would look like this:
			
			scrns=[-1.0,-0.5,0.0,0.5,1.0,1.25,1.5,1.75,2.0,2.5,3.0,3.5,4.0]
			

**Note that the sample scripts contains pre-prepared inputs and assembles these into the requisite tuple, so one need only update the values in the sample scripts to match
	your analysis. **
	
						
					

### The SedSASample class internals:

Methods and their purpose:

*__init__*: 
- parses the user-provided input data (see above), gathering it into six internal variables. 
- performs some internal computations that prepare the inputs, and ultimately creates the pandas data frame, that will hold all the data for analysis. The two most salient operations are the computation of the weight percentages and cumulative weight percentages for the sediment sample. Instantiating the class automatically completes these computations. 

	Inputs: the input tuple described above

	Returns: none

<br>
*__Analyze2CSV(self, fn)__*: 
performs internal processing (without need for the user to call the individual methods explicitly/directly) to compile the raw data, weight percentages, cumulative weight percentages, graphic and method of moment statistics, andthe mode(s). All of these are written directly to an output file (user specified path and file name) formatted as comma separated values. 

	Input args:
		fn = user-defined output file name

	Returns: none


<br>
*__Analyze2CSVHeader(self)__*: 
generates a columns header string (names of each data column written by method  Analyze2CSV) for output data file created in method Analyze2CSV

	Input Args: none

	Returns: header as a Python string


<br>
*__ComputeGraphicStats(self, s, Q)__*: 
computes "graphic" mean, median, standard deviation, skewness, and kurtosis (per Folk, 1980) for a given transect sample.
 
	Input args:
		s  sample identifier as a Python string (EX. 'S1')
	Q  list of interpolated quantile values (9 items) for current sample

	Returns: Python tuple containing MN=mean,MD=median,SD=std. dev.,SK=skewness, K=kurtosis in phi units.


<br>
*__ComputeMomentStats(self, s)__*: 
computes the sediment sample mean and standard deviation using the arithmetic method of moments method as described by Folk (1980) for current sample.

	Input args:
		s = sample identifier as a Python string (EX. 'S1')

	Returns: Python tuple containing momMean=mean,momSD=standard deviation,,undifCF=undifferentiated coarse fraction (the weight of material captured by the coarsest sieve) undifFF=undifferentiated fine fraction (the weight of material captured by the residual pan at the bottom of the seive stack).


<br>
*__FindSampleModes(self, s)__*: 
locates any modal weight percentage values in current sample. 

	Input args: 
		s = sample identifier as a Python string (EX. 'S1)

	Returns: Python ordered dictionary of modes (if any) found in the current sample.


<br>
*__GetDataFrame(self)__*: 
returns the sample set data frame df to the caller.

	Input args: none

 	Returns: pandas dataframe for the current sample data


<br>
*__InterpolateQuantileValues(self)__*: 
interpolate quantiles in self.quantilesList for each transect sample.

	Input args: none

	Returns: an ordered dictionary (key=sample identifier; value=interpolated quantile list)


<br>
*__PLOTDualSampleWtPercents(self, s, mode='print')__*: 
plots both the weight percentage and cumulative weight percentage histogram and curves (histo+PDF, and CDF) side by side and together, by sieve fraction for the current transect sample. Can plot to console or to stored PNG file.

	Input args:
		s = sample identifier as a Python string (EX. 'S1)
		mode = plot destination: mode='print' plot written to console (default)
		mode='save' plot written to png output

	Returns: none

	Note: 
	This is just a convenient combination of methods PLOTSampleWtPercents and PLOTSampleCumWtPercents.


<br>
*__PLOTSampleCumWtPercents(self, s, mode='print')__*: 
plots the individual sample cumulative weight percentages, by sieve fraction as a cunulative frequency (CDF) curve for the current transect sample. Can plot to console or to stored PNG file.

	Input args:
		s = sample identifier as a Python string (EX. 'S1)
		mode = plot destination: mode='print' plot written to console (default)
		mode='save' plot written to png output

	Returns: none


<br>
*__PLOTSampleWtPercents(self, s, mode='print')__*: 
plots the individual sample weight percentages, by sieve fraction as a histogram and overprinted frequency (PDF) curve for the current sample. Can plot to console or to stored PNG file.

	Input args:
		s = sample identifier as a Python string (EX. 'S1)
		mode = plot destination: mode='print' plot written to console (default)
		mode='save' plot written to png output

	Returns: none


<br>
*__PrintGraphicStats(self, s, gStats)__*: 
prints the graphic mean, median, standard deviation, skewness, and kurtosis for each transect sample in tabular format to the console. 

	Input args:
		s = sample identifier as a Python string (EX. 'S1)
		gStats = Python tuple containing the graphic mean, median, std. dev.skewness and kurtosis.

	Returns: none

	Note: the gStats tuple is generated by method self.ComputeGraphicStats()


<br>
*__PrintMomentStats(self, mStats)__*: 
prints the method of moments computed sediment size mean and standard deviation and undifferentiated coarse and fine fractions for the sample in tabular format to the console.

	Input args: 
		mStats = Python tuple containing method of moments mean, std. dev., undiff. coarse fraction weight, undiff fine fraction weight.

	Returns: none

	Note: the mStats tuple is generated by method self.ComputeMomentStats()


<br>
*__PrintSampleModes(self, modes)__*: 
Prints the mode values located in method FindSampleModes to the console

	Input args:
		modes = Python dictionary containing key (mode in phi units) and value (mode in weight percent) pairs

	Returns: none


<br>
*__PrintSampleWeightsDataTable(self, s)__*: 
prints the raw, weight percentage, and cumulative weight percentage values in tabular format for each screen bin to the console.

	Input args: 
		s = sample identifier as a Python string (EX. 'S1)

	Returns: none


<br>
*__ReturnQuantile(self, s, d)__*: 
uses simple linear interpolation to simulate the graphical estimation of a single quantile value. Quantile values are required to compute sediment sample statistics, in the style of Folk (1980). See self.QuantilesList for more.

	Input Args:
		s  sample identifier  as a Python string (EX. 'S1')
		d  quantile value to be interpolated as integer

	Returns: the interpolated quantile value in phi units

	NOTE this is an internal function called by InterpolateQuantileValues(). The user can call this method, but this is discouraged...I don't know why.



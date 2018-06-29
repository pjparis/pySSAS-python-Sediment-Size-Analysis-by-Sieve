#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug  3 12:01:40 2017

@author: paulp
"""

# -*- coding: utf-8 -*-

# ## import (if you can) prerequisite Python dependencies:
try:
	import sys
	import os
	import numpy as np
	import pandas
	import matplotlib.pyplot as plt
except(ImportError):
	print('Oops! There was a problem trying to import one of the requisite Python modules:',
	    sys.exc_info()[0], sys.exc_info()[1] )


# STATIC ATTRIBUTES:
# list of quantiles to be estimated (interpolated) as required by Folk and Ward
#  graphical statistics...
#quantilesList=[5,10,16,25,50,75,84,90,95]


class SedSAS(object):
	'''
	Expects input as a Pandas dataframe with two data columns. The first is the aperture 
	for each sieve used in the analysis. The second is the weight of sand/sediment captured
	in each sieve. In addition, a sample ID number is requested. This used to identify the
	particular sample when (as needed) SedSAS must post a message back to the user via the
	system console, or to write one or more records to a new results file.
	''' 
	def __init__(self, df, id='1'):
	
		self.id=id
		self.df=df
		#self.trim_end_vals=trim_end_vals
		
		print('\n\n Processing sample: ', self.id)
		try:
			# Check for all zeros in data--empty dataframe. If yes, then exit gracefully...
			rowsum=self.df['Weight'].sum()
			if( rowsum == 0.0 ):
				print('Oops! The dataframe for sample', str(self.id), 'appears to be empty.')
				print('New SedSASample class instance NOT created.')
			else:
				# 1.) Check the partially-differentiated/partially-determined size
				#     fractions at the coarse and fine ends of the weight classes (bins)
				pdsf = self.ComputePartiallyDeterminedSampleFractions()
				
				# compute the fraction-weight-percent column in dataframe df:
				self.df['Weight_Percent'] = self.df['Weight'] / \
				   self.df['Weight'].sum(axis=0)*100

				# compute the cumulative fractional weight percent column in dataframe df:
				self.df['Cum_Weight_Percent'] = self.df['Weight_Percent'].cumsum()

				## compute and add the class size midpoints in phi units
				self.ComputeWeightClassBinMidpoints()

				# interpolate or extrapolate quantiles:
				#self.Q=self.EstimateQuantileValues()
		except:
			print("Oops! Unexpected error [SedSAS init]:", sys.exc_info()[0], sys.exc_info()[1])

# ####### END OF __init__



### CLASS METHODS CALLED BY THE USER:
	### GENERAL METHODS:
	def ReturnWeightsData(self):
		''' Returns to the caller the working dataframe built by the class constructor
		as a Pandas dataframe object. Data returned are: raw weights, weight percentages,
		cumulative weight percentages, and mid-point class labels (in phi and mm units).
		The function takes no user-supplied arguments 
		'''
		try:
			return(self.df)
		except:
			print("Oops! Unexpected Error: [Your call to ReturnWeightsData]", sys.exc_info()[0])
			print('Maybe the weights table doesn\'t exist? Could be due to a problem that')
			print('occurred when the SedSASClass instance was created.')


	def ReturnQuantiles(self):
		'''returns the sample set quantiles associated with the Folk and Ward
        statistics (5th, 10th, 16th, 25th, 50th, 75th, 84th, 90th, and 95th)

		Input args: none

		Returns: a  Python tuple holding a pair of Python lists. The first
		list are quantiles in phi units. The second, quantiles in millimeters
		'''
		try:
			Quantiles=[3,10,16,25,50,75,84,90,95]   # quantiles in phi
			return( self.EstimateQuantileValues(Quantiles) )
		except:
			print("\n Oops! Unexpected Error: [Your call to ReturnQuantileList]", sys.exc_info()[0])
			print('Maybe the quantile list doesn\'t exist? Could be due to a problem that')
			print('occurred when the SedSASClass instance was created. \n')




	### STATISTICAL/COMPUTATIONAL METHODS ###
	
	## STATISTICS COMPUTATION METHOD 0: COMBINATION OF METHODS 1, 2, 4, and 5 PLUS
	## MODES AND TEXTURES
	def ComputeGSStats(self):
		'''n This is more or less a convenience method that combines four (of the
		most popular?) computational methods together into a single, relatively 
		easy to use spot. The four methods include the inclusive graphic statistics
		of Folk and Ward (logarithmic), along with geometric and
		logarithmic method of moments. Modes are also included. Qualitative 
		descriptors (textures) are included, too.
		
		The method takes no input arguments and returns a Python dictionary housing
		all of the results from the analyses
		'''
		## generate the statistics:
		fwlog=self.ComputeFWLogInclusiveGraphicsStats(return_description=True)
		McLog=self.ComputeMcCammonHiEffGraphicStats(return_description=True)
		momlog=self.ComputeLogarithmicMethodofMomentsStats(return_description=True)
		modes=self.FindSampleModes()
		
		tdict={'FWLogMean':fwlog[0][0],
		       'FWLogSort':fwlog[0][1],
		       'FWLogSkew':fwlog[0][2],
		       'FWLogKurt':fwlog[0][3],
		       'FWLogSizeClass':fwlog[1][0],
		       'FWLogSortCLass':fwlog[1][1],
		       'FWLogSkewClass':fwlog[1][2],
		       'FWLogKurtClass':fwlog[1][3],
		       'McLogMean':McLog[0][0],
		       'McLogSort':McLog[0][1],
		       'McLogSkew':McLog[0][2],
		       'McLogKurt':McLog[0][3],
		       'McLogSizeClass':McLog[1][0],
		       'McLogSortCLass':McLog[1][1],
		       'McLogSkewClass':McLog[1][2],
		       'McLogKurtClass':McLog[1][3],
		       'MoMLogMean':momlog[0][0],
		       'MoMLogSort':momlog[0][1],
		       'MoMLogSkew':momlog[0][2],
		       'MoMLogKurt':momlog[0][3],
		       'MoMLogSizeClass':momlog[1][0],
		       'MoMLogSortCLass':momlog[1][1],
		       'MoMLogSkewClass':momlog[1][2],
		       'MoMLogKurtClass':momlog[1][3]
		       }

		## add mode(s) to the dictionary:
		for n, mode in enumerate(modes):
			if n==0:
				tdict['PrimaryMode_phi'] = mode[0]
				tdict['PrimaryMode_mm'] = mode[1]
			if n==1:
				tdict['SecondaryMode_phi'] = mode[0]
				tdict['SecondaryMode_mm'] = mode[1]
			if n==2:
				tdict['TertiaryMode_phi'] = mode[0]
				tdict['TertiaryMode_mm'] = mode[1]
				
		return(tdict)
		
	
	
	## STATISTICS COMPUTATION METHOD 1: FOLK & WARD 'CLASSIC' LOG INCLUSIVE GRAPHICS:
	
	###def ComputeFWLogInclusiveGraphicsStats()
	def ComputeFWLogInclusiveGraphicsStats(self, return_description=False):
		'''computes "inclusive graphic" mean, standard deviation, skewness, and 
		kurtosis (per Folk and Ward, 1957 and Folk, 1980) for a given transect 
		sample. NOTE this is the 'classic' Folk and Ward particle size analysis
		method (traditionally using hand-drawn probability ordinate cumulative 
        frequency plots).

		Input args: return_description: whether or not to return qualitative
        textural descriptions, as per Folk, 1980 for the 4 statistical
        results. Default is False.

		Returns:
			Python list containing MN=mean, SD=std.deviation, 
			SK=skewness, K=kurtosis in phi units, or unknowns if cumulative
			weight percent[0] is > 50%.
		'''
		try:
			if( self.df['Cum_Weight_Percent'][0] < 50.0 ):
				Quantiles=[5,10,16,25,50,75,84,90,95]
				Q=self.EstimateQuantileValues(Quantiles)  # quantiles in phi

				MN=round( (Q[2]+Q[4]+Q[6])/3, 3)           # for the mean:
				SD=round( ((Q[6]-Q[2])/4)+((Q[8]-Q[0])/6.6), 3)  # sorting  

				# the graphic skewness (SK) and kurtosis (K):
				SK=round( ((Q[6]+Q[2]-2*Q[4])/(2*(Q[6]-Q[2]))) + \
					((Q[8]+Q[0]-2*Q[4])/(2*(Q[8]-Q[0]))), 3)
				K=round( (Q[8]-Q[0])/(2.44*(Q[5]-Q[3])), 3 )

				if(return_description == True):
					descriptions=[]
					descriptions.append( self.ClassifyByGrainSizePhi(MN) )
					descriptions.append( self.ClassifyLogarithmicSorting(SD) )
					descriptions.append( self.ClassifyFWLogSkewness(SK) )
					descriptions.append( self.ClassifyFWKurtosis(K) )
				else:
					descriptions=['-','-','-','-']
					
			else:
				MN='< '+str(self.df['Aperture'][0] )
				SD,SK,K='not computed','not computed','not computed'
				print('Five or more quantiles require extrapolation which is insufficient for reliable value estimation--no FW-Logarithmic statistics were generated')

		except (NameError, AttributeError, KeyError):
			print("\n Oops! Unexpected Error: [Your call to ComputeFWLogarithmicGraphicStats]", 
			   sys.exc_info()[0])
			print('Maybe the problem is due to something that occurred when the SedSASClass')
			print('instance was created. \n')
			MN,SD,SK,K=0,0,0,0
			descriptions=['-','-','-','-']

		return( [MN,SD,SK,K], descriptions )


	## STATISTICS COMPUTATION METHOD 2: McCammon GRAPHIC:
    
	### def ComputeMcCammonHiEffGraphicStats()
	def ComputeMcCammonHiEffGraphicStats(self, return_description=False):
		'''computes high-efficiency "graphic" mean and standard deviation only
		 (per McCammon, 1962) for a given sediment sample. This is an experimental
         high-efficiency method (97% vs 88% for FW) for computing particle-size
         statistical estimators. Might prove an interesting comparison betw these
         and Log MoM results...

		Input args: return_description: whether or not to return qualitative
        textural descriptions, as per Folk, 1980 for the mean and sorting
        results. Default is False.

		Returns:
			Python tuple containing MN=mean, SD=std. deviation in phi.
		'''
		try:
			if( self.df['Cum_Weight_Percent'][0] < 50.0 ):

				# for the mean
				Quantiles=[5,15,25,35,45,55,65,75,85,95]
				Q=self.EstimateQuantileValues(Quantiles)  # quantiles in phi
				
				MN=round( sum(Q)/10,3)  # mean

				# for sorting (standard deviation):
				Quantiles=[3,10,20,30,70,80,90,97]
				Q=self.EstimateQuantileValues(Quantiles)  # quantiles in phi
				SD=round( (Q[4]+Q[5]+Q[6]+Q[7]-Q[0]-Q[1]-Q[2]-Q[3])/9.1, 3 )

				SK,K = '-','-'    # empty skew and kurt
				
				if(return_description == True):
					descriptions=[]
					descriptions.append(self.ClassifyByGrainSizePhi(MN) )
					descriptions.append( self.ClassifyLogarithmicSorting(SD) )
					descriptions.append('-');  descriptions.append('-')
				else:
					descriptions=['-','-','-','-']
				
			else:
				MN='< '+str(self.df['Aperture'][0] )
				SD='not computed'
				SK,K = '-','-'
				print('Five or more quantiles require extrapolation which is insufficient for reliable value estimation--no FW-geometric statistics were generated')

		except (NameError, AttributeError, KeyError):
			print("\n Oops! Unexpected Error: [Your call to ComputeMcCammonHiEffGraphicStats]", 
			   sys.exc_info()[0])
			print('Maybe the problem is due to something that occurred when the SedSASClass')
			print('instance was created. This is often the case.\n')
			MN,SD, SK, K = 0,0,'-','-'
			descriptions=['-','-','-','-']
			
		return( [MN,SD, SK, K], descriptions )


	## STATISTICS COMPUTATION METHOD 3: METHOD OF MOMENTS - ARITHMETIC:
	def ComputeArithmeticMethodofMomentsStats(self):
		'''computes the mean, sorting (standard deviation), skewness, and
		kurtosis (the first four statistical moments) using the method of
		moments technique (based on weight distributions--weight fractions
		captured in each seive). 

		Inputs: None

		Returns: Python list containing the four computed moments

		Units are mm, where appropriate
		Note that this method is not in common use in the geosciences.
		'''
		try:
			## compute the mean (M1):
			M1=round( ( self.df['Weight_Percent']*
               self.df['MidPt_mm'] ).sum() / self.df['Weight_Percent'].sum(), 3)

			## compute the 2nd moment (M2 std. dev.):
			M2=round( ( ((self.df['Weight_Percent']*
               (self.df['MidPt_mm']-M1)**2 ).sum())/ 100)**0.5, 3 )

			## compute the 3rd moment (M3 skewness):
			M3=round( ((self.df['Weight_Percent']*
               (self.df['MidPt_mm']-M1)**3 ).sum())/(100*(M2**3)), 3 )

			## compute the 4th moment (M4 kurtosis):
			M4=round( ((self.df['Weight_Percent']*
               (self.df['MidPt_mm']-M1)**4 ).sum())/(100*(M2**4)), 3 )

		except (NameError, AttributeError, KeyError):
			print("\n Oops! Unexpected Error: [Your call to ComputeArithmeticMethodofMomentsStats]", 
			   sys.exc_info()[0])
			print('Maybe the problem is due to something that occurred when the SedSASClass')
			print('instance was created. This is frequently the case.\n')
			M1,M2,M3,M4=0,0,0,0
			descriptions=[]
			
		return( [M1,M2,M3,M4] )


	## STATISTICS COMPUTATION METHOD 4: METHOD OF MOMENTS - LOGARITHMIC:
	def ComputeLogarithmicMethodofMomentsStats(self, return_description=False):
		'''computes the mean, sorting (standard deviation), skewness, and
		kurtosis (the first four statistical moments) using the logarithmic 
		method of moments technique (based on weight distributions--weight 
		fractions captured in each seive). These are taken from Table 1 in 
        Friedman, 1962.

		Inputs: None

		Returns: Python list containing the four computed moments 

		Size units are in phi, where appropriate
		'''
		try:
			## compute the mean (M1):
			M1=round( ( self.df['Weight_Percent']*
                       self.df['MidPt_phi'] ).sum() / 100.0, 3 )

			## compute the 2nd moment (M2 std. dev.):
			M2=round( ( (self.df['Weight_Percent']*
              (self.df['MidPt_phi']-M1)**2).sum()/100 )**0.5, 3 )

			## compute the 3rd moment (M3 skewness):
			M3=round( (self.df['Weight_Percent']* 
            (self.df['MidPt_phi']-M1)**3 ).sum()/(100*(M2**3)), 3 )

			## compute the 4th moment (M4 kurtosis):
			M4=round( (self.df['Weight_Percent']* 
             (self.df['MidPt_phi']-M1)**4 ).sum()/(100*(M2**4)), 3 )

			if(return_description == True):
				descriptions=[]
				descriptions.append( self.ClassifyByGrainSizePhi(M1) )
				descriptions.append( self.ClassifyLogarithmicSorting(M2) )
				descriptions.append( self.ClassifyMoMLogSkewness(M3) )
				descriptions.append( self.ClassifyMoMKurtosis(M4) )
			else:
				descriptions=['-','-','-','-']

		except (NameError, AttributeError, KeyError):
			print("\n Oops! Unexpected Error: [Your call to ComputeLogarithmicMethodofMomentsStats]", 
			   sys.exc_info()[0])
			print('Maybe the problem is due to something that occurred when the SedSASClass')
			print('instance was created. This is typically the case.\n')
			M1,M2,M3,M4=0,0,0,0
			descriptions=['-','-','-','-']
			
		return( [M1,M2,M3,M4], descriptions )


	## STATISTICS COMPUTATION METHOD X: METHOD OF MOMENTS - GEOMETRIC:
	## TAKEN FROM BLOTT AND PYE, 2001
	## THIS METHOD IS UNDOCUMENTED, AND VERY MUCH EXPERIMENTAL...
	def ComputeGeometricMethodofMomentsStats(self, return_description=False):
		'''computes the mean, sorting (standard deviation), skewness, and
		kurtosis (the first four statistical moments) using the geometric 
		method of moments technique (based on natural log of weight 
		distributions--weight fractions captured in each seive). 

		Inputs: None

		Returns: Python list containing the four computed moments

		Units are mm, where appropriate
		'''
		try:
			## compute the mean (M1):
			M1=round( np.exp( ( self.df['Weight_Percent']*
                    np.log(self.df['MidPt_mm']) ).sum() /
                                     self.df['Weight_Percent'].sum() ), 3 )

			## compute the 2nd moment (M2 std. dev.):
			M2=round( np.exp( ( ((self.df['Weight_Percent']*
                              (np.log(self.df['MidPt_mm'])-
                               np.log(M1))**2 ).sum())/ 100)**0.5 ), 3 )

			## compute the 3rd moment (M3 skewness):
			M3=round( (self.df['Weight_Percent']*((np.log(self.df['MidPt_mm'])-
                   np.log(M1))**3)).sum() / (100*np.log(M2)**3), 3 )

			## compute the 4th moment (M4 kurtosis):
			M4=round( (self.df['Weight_Percent']*(np.log(self.df['MidPt_mm'])-
                  np.log(M1))**4 ).sum()/(100*np.log(M2)**4), 3 )

			if(return_description == True):
				descriptions=[]
				descriptions.append( self.ClassifyByGrainSizeMM(M1) )
				descriptions.append( self.ClassifyGeometricSorting(M2) )
				descriptions.append( self.ClassifyMoMGeoSkewness(M3) )
				descriptions.append( self.ClassifyMoMKurtosis(M4) )
			else:
				descriptions=['-','-','-','-']

		except (NameError, AttributeError, KeyError):
			print("\n Oops! Unexpected Error: [Your call to ComputeGeometricMethodofMomentsStats]", 
			   sys.exc_info()[0])
			print('Maybe the problem is due to something that occurred when the SedSASClass')
			print('instance was created. This is frequently the case.\n')
			M1,M2,M3,M4=0,0,0,0
			descriptions=['-','-','-','-']
			
		return( [M1,M2,M3,M4], descriptions )


	## STATISTICS COMPUTATION METHOD: FINDING THE SAMPLE MODES
	def FindSampleModes(self):
		'''locates any modal weight percentage values in current sample. 

		Input args: None

		Returns: Python tuple of modes (if any) found in the 
			current sample. The 0th item in each tuple contains the mode 
			in phi units. The 1st item is the mode in mm.
           
		# HOW THIS THING WORKS:
		# first time thru it reads the 0th item in the sorted weight list B and 
		# posts this to the modes dictionary (this is the largest or primary 
		# mode). Then, it finds any additional modes by taking each weight 
		# percentage value in the sorted list B, locating that weight value and 
		# its index in the weight percentage column (df[s+'wp'] of the df, then 
		# comparing that weight percent value to those 
		# immediately before and immediately after the current. If both adjacent 
		# weights are less than the current weight we have another mode, which we 
		# insert into the modes dictionary.
		'''
		# initialize an empty modes list:
		modes=[]

		try:
			# sort the weight percentages column in df in descending 
			# order, put into list B:
			B=sorted(list( self.df['Weight_Percent'] ), reverse=True)

			for i in range(len(B)):
				if i==0:                   # for primary mode
					phi=self.df['Aperture'][ list( self.df['Weight_Percent'] ).index(B[0])]
					mm=2**(phi*-1)
					modes.append( (phi,mm) )
				else:
					Amax=B[i]
					Amax_index=list( self.df['Weight_Percent'] ).index(Amax)
				# as long as Amax_index is not allowed to exceed the length of 
				# the w weight percentage column in df, less one. If it equals 
				# the length of the column then, when you offset by +1 in 
				# searching for a mode you generate an out of range index...bomb! 
				# If the index is okay, we can proceed to check for a mode by 
				# looking at wt values surrounding the Amax value.
					if( (Amax_index < len(self.df['Weight_Percent'])-1) and \
						(self.df['Weight_Percent'].iloc[Amax_index-1] < Amax) and \
						(self.df['Weight_Percent'].iloc[Amax_index+1] < Amax) ):
						phi=self.df['Aperture'].iloc[Amax_index]
						mm=2**(phi*-1)
						modes.append( (phi,mm) )
					
			if(len(modes) > 1):
				print('WARNING! Number of modes in sample is', len(modes), 'Textural shape descriptions of a')
				print('multimodal sample distribution are unreliable or possibly even nonsensical.')
				
			return(modes)
		except (NameError, AttributeError, KeyError, ValueError):
			print("\n Oops! Unexpected Error: [Your call to FindSampleModes]", 
			   sys.exc_info()[0])
			print('Maybe the problem is due to something that occurred when the SedSASClass')
			print('instance was created. Perhaps this is the case?\n')
			return( () )



	### PLOTTING METHODS: ###
	def PLOTSampleWtPercents(self, printTo='console'):
		'''plots the individual sample weight percentages, by sieve fraction
		as a histogram and overprinted frequency (PDF) curve for the current
		sample. Can plot to console or to stored PNG file.

		Input args:
			printTo = plot destination: printTo='console' plot written to monitor
				(default)
			printTo='file' plot written to png output file in local directory

		Returns: none
		'''
		plt.style.use('fivethirtyeight')
		fig1 = plt.figure(figsize=(7,5))
		axh1 = fig1.add_subplot(1,1,1)
		fig1.tight_layout()  #pad=2.1, w_pad=0.5, h_pad=2.0)

		bins = np.arange(len(self.df['Aperture']))
		w=0.75
		axh1.bar(bins,self.df['Weight_Percent'], align='center', width=w)
		axh1.plot(bins,self.df['Weight_Percent'],'g--',linewidth=3.5)
		axh1.set_xlim(-0.25,len(self.df['Aperture']))
		#axh1.set_ylim(0,105)
		xTickLabels=list(map(str, self.df['Aperture']))
		xTickLabels[-1]='fines'
		axh1.set_xticks(bins )   #, xTickLabels, rotation=90
		axh1.set_xticklabels(xTickLabels, rotation=90)
		axh1.set_xlabel('Sediment Size (phi)')
		axh1.set_ylabel('Weight Percent')
		axh1.set_title('Histogram and Frequency Curve - Sample: '+str(self.id) )
		
		print('-'*80)
		if(printTo=='console'):
			plt.show()
		if(printTo=='file'):
			print('Writing plot file to local directory...')
			fig_name='HistoPlot_'+str(self.id)+'_wf.png'
			plt.savefig(fig_name, bbox_inches='tight')
			plt.close()


	def PLOTSampleCumWtPercents(self, printTo='console'):
		'''plots the individual sample cumulative weight percentages, by sieve
		fraction as a cunulative frequency (CDF) curve for the current
		transect sample. Can plot to console or to stored PNG file.

		Input args:
			mode = plot destination: mode='print' plot written to console
			(default)
			mode='save' plot written to png output

		Returns: none
		'''
		plt.style.use('fivethirtyeight')
		fig2 = plt.figure(figsize=(7,5))
		axh2 = fig2.add_subplot(1,1,1)
		fig2.tight_layout()  #pad=2.1, w_pad=0.5, h_pad=2.0)

		bins = np.arange(len(self.df['Aperture']))
		#w=0.75
		#plt.bar(bins,self.df[s+'cwp'], align='center', width=w)
		axh2.plot(bins,self.df['Cum_Weight_Percent'],'g--',linewidth=3.5)
		axh2.set_xlim(-0.25,len(self.df['Aperture']))
		#plt.ylim(0,105)
		xTickLabels=list(map(str, self.df['Aperture']))
		xTickLabels[-1]='fines'
		axh2.set_xticks(bins)
		axh2.set_xticklabels(xTickLabels, rotation=90 )
		axh2.set_xlabel('Sediment Size (phi)')
		axh2.set_ylabel('Weight Percent')
		axh2.set_title('Cumulative Frequency Curve-Linear Ordinate Sample: '+str(self.id) )
		fig_name='CumCurve_'+str(self.id)+'_cwf.png'
		print('-'*80)
		if(printTo=='console'):
			plt.show()
		if(printTo=='file'):
			print('Writing plot file to local directory...')
			plt.savefig(fig_name,bbox_inches='tight')
			plt.close()


	def PLOTDualSampleWeightPercents(self, printTo='console'):
		'''plots both the weight percentage and cumulative weight percentage
		histogram and curves (histo+PDF, and CDF) side by side and together,
		by sieve fraction for the current transect sample. Can plot to
		console or to stored PNG file.

		Input args:
			mode = plot destination: mode='print' plot written to console
				(default)
			mode='save' plot written to png output

		Returns: none

		Note: 
			This is just a convenient combination of methods
			- PLOTSampleWtPercents and 
			- PLOTSampleCumWtPercents.
		'''
		#import seaborn as sns
		#sns.set()
		plt.style.use('fivethirtyeight')
		fig3 = plt.figure(figsize=(21,7))
		#axh3 = fig3.add_subplot(1,1,1)
		#fig3.tight_layout()  #pad=2.1, w_pad=0.5, h_pad=2.0)

		bins = np.arange(len(self.df['Aperture']))
		w=0.75
		xTickLabels=list(map(str, self.df['Aperture']))	  # convert list of numeric screen sizes 
													# to list of strings for labels
		xTickLabels[-1]='fines'

		# write header to console:
		#print('-'*90)
		#print('  Transect: '+self.transect+'   Sample: '+s)
		print('')

		#fig1=plt.figure(figsize=(21,7))

		# weight percentages subplot
		axh3=fig3.add_subplot(1,2,1)
		axh3.bar(bins,self.df['Weight_Percent'], align='center', width=w)
		axh3.plot(bins,self.df['Weight_Percent'],'g--',linewidth=3.5)
		axh3.set_xlim(-0.25,len(self.df['Aperture']))
		#axh3.set_ylim(0,105)
		axh3.set_xticklabels(xTickLabels, rotation=90)
		axh3.set_xticks(bins)
		axh3.set_xlabel('Sediment Size (phi)')
		axh3.set_ylabel('Weight Percent')
		axh3.set_title('Histogram and Frequency Curve')

		# cumulative weight percentages subplot:
		axh4=fig3.add_subplot(1,2,2)
		#axh4.bar(bins,self.df[s+'cwp'], align='center', width=w)   # uncomment to draw bars
		axh4.plot(bins,self.df['Cum_Weight_Percent'],'g-', linewidth=2.0)
		axh4.set_xlim(-0.25,len(self.df['Aperture']))
		axh4.set_ylim(0,105)
		axh4.set_xticklabels(xTickLabels, rotation=90)
		axh4.set_xticks(bins)
		axh4.set_yticks([10,20,30,40,50,60,70,80,90,100])
		#axh4.minorticks_on()
		axh4.set_xlabel('Sediment Size (phi)')
		axh4.set_ylabel('Cumulative Weight Percent')
		axh4.set_title('Cumulative Frequency Curve-Linear Ordinate') #: Transect '+transect+'  Sample '+s)
		#axh4.grid(which='both')

		fig_name='HistoCumPlot_'+str(self.id)+'_dual.png'
		if(printTo=='console'):
			plt.show()
		if(printTo=='file'):
			print('Writing plot to local directory...')
			plt.savefig(fig_name, bbox_inches='tight')
			plt.close()

		print('-'*90)

### END OF CLASSES CALLED BY THE USER ###


### CLASS METHODS CALLED INTERNALLY (THESE ARE NOT TYPICALLY CALLED BY THE USER...):
	def ComputePartiallyDeterminedSampleFractions(self):
		'''Finds the amount of the total sediment fraction in the sample that: 1.) failed to 
		pass through the largest aperature sieve in the stack, and/or 2.) the fraction of the
		sample that fell into the pan at the bottom of the stack. (if either of these values
		exceeds 5% of the total sample the user is warned about potential impacts on results)
	
		Requires: a Pandas dataframe containing the sieved sample weights from SEDSASClass
	
		Returns: a Python list containing the computed percent coarse fraction (pcf) and percent
		fine fraction (pff), and optionally, warns the user of excess sample in either 
		fraction on the console.
		'''
		pcf=(self.df['Weight'][0] / self.df['Weight'].sum()) * 100
		if( pcf >= 5.0 ):
			print('WARNING: percent of coarsest fraction in sample', str(self.id),'is:')
			print(str(round(pcf,2)),'percent. This exceeds 5% of total by weight.')
			print('Values in excess of 5% can introduce significant error in some analyses.\n')

		pff=(self.df['Weight'][len(self.df)-1] / self.df['Weight'].sum() ) * 100
		if( pff >= 5.0 ):
			print('WARNING: percent of fine fraction in sample', str(self.id),'is:')
			print(str(round(pff,2)),'percent. This exceeds 5% of total by weight.') 
			print('If this percentage represents the pan fraction the user is warned')
			print('that associated errors may be large enough to compromise results. If the')
			print('pan fraction was not included in sample, this warning may be safely ignored.\n')

		return([pcf, pff])
	
	
	def ComputeWeightClassBinMidpoints(self):
		'''Computes the size midpoint value for each sieve weight class. The weight class
		is defined by the collection of sieve apertures used in the analysis, as recorded
		in the 0th column in df.
	
		Requires: a Pandas dataframe containing the sieve apertures used in the analysis
		in phi units
	
		Returns: -
		'''
		x=np.array( self.df['Aperture'] )
		mp=(x[1:] + x[:-1]) / 2
		self.df['MidPt_phi']=np.append(self.df['Aperture'][0:1].values-0.5, mp)
		self.df['MidPt_mm']=2**(self.df['MidPt_phi']*-1)
		
	
	
	def ExtrapolateQuantileLinearFit(self, q, qnt):
		''' Linearly extrapolates quantile value. With extrapolation there is no lower
		cumulative (the lowest cumulative weight, the weight of that amount of material
		contained in the largest aperture sieve used in the analysis is > q) weight value
		so we must build one into the data artificially. Here we use the simple technique
		GROW_STACK_by_1, which simply adds an additional sieve to the stack that is 1
		size larger than the largest sieve aperture used in the analysis. The lower y is 
		then assumed to be 0.
		'''
		print('WARNING: '+str(qnt)+'% quantile < min cum wt%, extrapolating solution')

		Xk_1=float(q['Aperture'].iloc[0]) - 1.0
		Yk_1=0.0
		Xk=float(q['Aperture'].iloc[0])
		Yk=float(q['Cum_Weight_Percent'].iloc[0])

		return( Xk_1+((Xk-Xk_1)*(qnt-Yk_1)/(Yk-Yk_1)) )


	def InterpolateQuantileLinearFit(self, p, q, qnt):
		''' Uses simple linear interpolation to estimate quantile values between known
		cumulative weight values.
		'''
		Yk_1=float(p['Cum_Weight_Percent'].iloc[[-1]])
		Yk=float(q['Cum_Weight_Percent'].iloc[[0]])
		Xk_1=float(p['Aperture'].iloc[[-1]])
		Xk=float(q['Aperture'].iloc[[0]])
		return( Xk_1+((Xk-Xk_1)*((qnt-Yk_1)/(Yk-Yk_1))) )


	def EstimateQuantileValues(self, quantilesList):
		'''interpolate quantiles in quantilesList for each sample. 
		Calls function ReturnQuantile() for each quantile value to be computed 
		from the seived weights

		Requires: -

		Returns: a Python tuple holding a pair of Python list of the
		interpolated quantiles. The first list is in Phi units, the second in 
		millimeters.
		'''
		qntLPhi=[];   qntLmm=[]

		for quantile in quantilesList:
			p=self.df.query('Cum_Weight_Percent'+'<'+str(quantile) )
			q=self.df.query('Cum_Weight_Percent'+'>'+str(quantile) )

			if(p.empty) & (q.empty):     # if both sub dfs are empty, abort
				return(np.nan)
			if p.empty:            # we must extrapolate to estimate quantile:
				if( float(q['Cum_Weight_Percent'].iloc[0]) > quantilesList[0] ):
					qval=self.ExtrapolateQuantileLinearFit(q, quantile)
			else:                  # otherwise, (linear) interpolate quantile:
				qval=self.InterpolateQuantileLinearFit(p, q, quantile)

			qntLPhi.append( round(qval, 3) )
#			qntLmm.append( round(2**((-1)*qval),3 ) )

		return( (qntLPhi) )


	### TEXTURAL DESCRIPTION METHODS ####
	## WENTWORTH, 1922
	def ClassifyByGrainSizePhi(self, MN):
		'''
		'''
		try:
			if( MN <= -8.0):
				return('Boulders')
			if( MN <= -6.0 and MN > -8.0):
				return('Cobbles')
			if( MN <= -2.0 and MN > -6.0):
				return('Pebbles')
			if( MN <= -1.0 and MN > -2.0):
				return('Gravel Granules')
			if( MN <= 0.0 and MN > -1.0):
				return('Very Coarse Sand')
			if( MN <= 1.0 and MN > 0.0):
				return('Coarse Sand')
			if( MN <= 2.0 and MN > 1.0):
				return('Medium Sand')
			if( MN <= 3.0 and MN > 2.0):
				return('Fine Sand')
			if( MN <= 4.0 and MN > 3.0):
				return('Very Fine Sand')
			if( MN <= 5.0 and MN > 4.0):
				return('Coarse Silt')
			if( MN <= 6.0 and MN > 5.0):
				return('Medium Silt')
			if( MN <= 7.0 and MN > 6.0):
				return('Fine Silt')
			if( MN <= 9.0 and MN > 7.0):
				return('Very Fine Silt')
			if( MN > 9.0):
				return('Clay')
		except (NameError, AttributeError, KeyError) as e:
			print('Oops! An error occured.', sys.exc_info()[0], e)


	def ClassifyLogarithmicSorting(self, SD):
		''' Classifies the resultant sorting (standard deviation) value as computed using 
		LOGARITHMIC methods from Method of Moments and Folk and Ward approaches
		
		Called by STATISTICAL METHOD 1 (FOLK AND WARD LOGARITHMIC) and METHOD 4 (METHOD
		OF MOMENT LOGARITHMIC)
		FROM FOLK AND WARD, 1957
		'''
		if( SD < 0.35):
			return('Very well sorted')
		if( SD >= 0.35 and SD < 0.50):
			return('Well sorted')
		if( SD >= 0.50 and SD < 0.71):
			return('Moderately well sorted')
		if( SD >= 0.71 and SD < 1.00):
			return('Moderately sorted')
		if( SD >= 1.00 and SD < 2.00):
			return('Poorly sorted')
		if( SD >= 2.00 and SD < 4.00):
			return('Very poorly sorted')
		if( SD > 4.00):
			return('Extremely poorly sorted')


	def ClassifyFWLogSkewness(self, SK):
		''' Classifies a resultant skewness value as computed using the
		Folk and Ward (1957) graphic Logarithmic method of computation. 
		
		Called by STATISTICAL METHOD 1 (FOLK AND WARD LOGARITHMIC)
		FROM FOLK AND WARD, 1957
		'''
		if( SK > 0.3 and SK <= 1.00):
			return('Strongly fine-skewed')
		if( SK > 0.1 and SK <= 0.3):
			return('Fine skewed')
		if( SK < 0.1 and SK >= -0.1):
			return('Near symmetrical')
		if( SK < -0.1 and SK >= -0.3):
			return('Coarse skewed')
		if( SK < -0.3 and SK >= -1.00):
			return('Strongly coarse-skewed')


	def ClassifyMoMLogSkewness(self, SK):
		''' Classifies a resultant skewness value as computed using the
		Method of Moments Logarithmic computation. 
		
		Called by STATISTICAL METHOD 4 (METHOD OF MOMENTS LOGARITHMIC)
		FROM BLOTT AND PYE, 2001
		'''
		if( SK > 1.30):
			return('Very fine skewed')
		if( SK > 0.43 and SK <= 1.30):
			return('Fine skewed')
		if( SK > -0.43 and SK <= 0.43):
			return('Symmetrical')
		if( SK > -1.30 and SK <= -0.43):
			return('Coarse skewed')			
		if( SK < -1.30):
			return('Very coarse skewed')


	def ClassifyFWKurtosis(self, K):
		''' Classifies a resultant kurtosis value as computed using the
		Folk and Ward (1957) graphic Logarithmic and Geometric methods of computation. 
		The classification is based on Folk and Ward descriptions for kurtosis
		
		Called by STATISTICAL METHOD 1 (FOLK AND WARD LOGARITHMIC)
		FROM FOLK AND WARD, 1957
		'''
		if( K < 0.67):
			return('Very platykurtic')
		if( K >= 0.67 and K < 0.90):
			return('Platykurtic')
		if( K >= 0.90 and K < 1.11):
			return('Mesokurtic')
		if( K >= 1.11 and K < 1.50):
			return('Leptokurtic')
		if( K >= 1.50 and K < 3.00):
			return('Very leptokurtic')			
		if( K > 3.00):
			return('Extremely leptokurtic')		


	def ClassifyMoMKurtosis(self, K):
		''' Classifies a resultant kurtosis value as computed using the
	Logarithmic and Geometric methods of moments computations. 
		
		Called by STATISTICAL METHOD 4 (METHOD OF MOMENTS LOGARITHMIC):
		FROM BLOTT AND PYE, 2001
		'''
		if( K < 1.70):
			return('Very platykurtic')
		if( K >= 1.70 and K < 2.55):
			return('Platykurtic')
		if( K >= 2.55 and K < 3.70):
			return('Mesokurtic')
		if( K >= 3.70 and K < 7.40):
			return('Leptokurtic')
		if( K > 7.40 ):
			return('Very leptokurtic')
			
### END OF CLASS METHODS CALLED INTERNALLY



# ## ####################################################################################
# ## End of SedSASample class listing 
# ## ####################################################################################






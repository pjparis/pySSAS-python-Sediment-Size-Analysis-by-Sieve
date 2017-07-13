# -*- coding: utf-8 -*-

# ## import prerequisite Python module dependencies:
import math
import numpy as np
import pandas
import matplotlib.pyplot as plt
import collections


class SedSAS(object):
	'''
		Class SedSAS - Sediment Size Analysis by Sieve

		DESCRIPTION:
		Performs various data extractions, manipulations, statistical analyses, and plotting
		to characterize [sand and gravel sized] sediment particle size distributions from 
		mechanical dry sieve fractionation.

		From a sample file of sieved sand (ϕ) fractional weights, SedSASample() can:

		 - compute weight percentages and cumulative weight percentages for each ϕ
		   fraction relative to the total sample

		 - compute the mean, median, and sorting (standard deviation) in  ϕ units using the 
		   graphic method described by Folk (1980) the mean and sorting (standard 
		   deviation) in ϕ units using the method of moments as described by Folk (1980)

		 - create and plot a histogram with frequency curve (PDF) overlay
 
		 - create and plot a cumulative frequency curve (CDF) 


		NOTES: 
			1.) SedSASample requires the following Python libraries (they load upon instantiation): 
					numpy as np
					pandas
					matplotlib.pyplot as plt
					collections

			2.) SedSASample was written using, and for, Python 3.x. It has not been tested
				in a 2.x environment, but would be expected to function with little or no
				modifications required.

			3.) Cited reference: Folk, R. 1980. Petrology of Sedimentary Rocks. Hemphill 
				Publishing Company, Austin, TX. 184p.


		HISTORY:
			- initial code assembly: April/May, 2016  (pjp)
			- method FindSampleModes() added 30 May, 2016   (pjp)
			- removed pprint statement from method PrintSampleWeightsDataTable() \
						## pprint.pprint( self.quantilesDict[s] )  1 June, 2016 (pjp)
			- added check on undifferntiated coarse and fine fractions 8 June, 2016 (pjp)
			- added method PrintSamplesModes  9 June, 2016 (pjp)
			- added method Analyze2Stdout 30 June, 2017 (pjp)
			- modified class input stream to simplify class instantiation July, 2017 (pjp)
			- added assertion handling to constructor (__init__) July, 2017 (pjp)

	'''

	def __init__(self,df,transect,samples,screens):
		 # load the prepared and ready-to-go pandas data frame:
		self.df0=df       
		 # current transect id (as a Python string)
		self.transect=transect	 # was input 4		  
		# samplesList=list of transect sample site ids (usually S1-S4 as Python strings)
		self.samplesList=samples   # was input 5			
		# screens=list of seive screen sizes used in the analysis--transposed (floats)
		self.screens=np.array(screens).T   # was input 6	
		
		# quantilesList=list of quantiles that will be interpolated to compute statistics.
		# These are fixed, based on requirements from Folk, (1980).
		self.quantilesList=[5,10,16,25,50,75,84,90,95] 
		
		# quantilesDict=ordered dict. of interpolated quantiles (keys=sample id.; 
		# values=list of interpolated quantile values (9 items as per quantilesList)) 
		self.quantilesDict=collections.OrderedDict()
		
		# alert user as to the current transect sample(s) being processed:
		print('')
		print('Processing transect:', self.transect )
		print('')
		
		# ####### some internal temporary data frame ops (6 of 'em) that build the final 
		# df data frame:
		#
		# 1. load the source data file into an initial pandas data frame df0, then extract 
		#	only the records (rows) of the current transect tr. A transpose of the 
		#	subset, with column headers added, and then the individual screen weights are 
		#	parsed out for final inclusion in df1. df1 is next, and finally, checked for 
		#	data (no NaNs) and if none are found then the data is passed on to step 2. 
		#	Otherwise, an note is posted to the user (console) informing that there is 
		#	missing data in the current sample and so will be skipped over. Processing is
		#	then passed to the next sample

		# extract recs for current transect and computes transpose
		try:
			df1=self.df0.loc[ self.df0['Transect'] == self.transect ].copy().T
			df1=df1.drop(['Transect','Sample'], axis=0)		  # added 6/29/2017
			df1.columns=self.samplesList					  # add column headers
			# keep only the [sieve] screen weight columns

			# 2. checks the percentage of sand/sediment that was stopped by the largest screen 
			# as well as the percentage captured in the pan at the base of the sieve stack. If 
			# either of the percentages is > 5%, warn the user. If the latter, proceed further
			# to replace existing pan fraction with 0.0 gms. Why do this? See print statements
			# below.
			for s in self.samplesList:
				pcf=df1[s][0]/df1[s].sum()*100	# percent coarse fraction
				pff=df1[s][-1]/df1[s].sum()*100   # percent fine fraction

				if(pcf>5.0):
					print('')
					print('WARNING: percent of unidfferentiated coarse fraction in sample',s,'is:' ,str(round(pcf,2)), '--Values in excess of 5% can introduce significant error in resulting statistical computations.')
				if( pff > 5.0 ):
					print('')
					print('WARNING: percent of unidfferentiated fine fraction in sample',s, 'is:' ,str(round(pff,2)), '--As values in excess of 5% can introduce significant error in resulting statistical computations, this fraction will not be considered in the analysis.')
					# replace the pan fraction (screen [-1]) with 0.0
					df1[s][-1]=0.0	
			
			
			# 3. compute the sample fraction weight percent sums (df2) and cumulative weight
			# percentages (df3) for the transect samples...
			df2=(df1/df1.sum(axis=0))*100		 # compute weight percentages
			df3=df2.cumsum()					 # compute cumulative weight percentages from 
											 # weight percentages
		
			# 4. build unique column names for weight percentage and cumulative weight 
			# percentage values...
			wDfColNamesList=[];   cDfColNamesList=[]
			for s in self.samplesList:
				wDfColNamesList.append(s+'wp')
				cDfColNamesList.append(s+'cwp')
			df2.columns=wDfColNamesList		 # set weight percent column names
			df3.columns=cDfColNamesList		 # set cumulative weight percent column names
		
			# 5. append both the weight percentage and cumulative weight percentage values to 
			# df
			frames=[df1,df2,df3]
			self.df=pandas.concat(frames, axis=1)
		
			# 6. add column of seive screen sizes (phi) to df
			self.df['phi']=self.screens			
			self.foundTransect = True
		except:
			print('Missing data for transect:', self.transect )
			self.foundTransect = False
	# ####### End of internal df calcs. Now, let's get to work!	
	# ####### End of Constructor (__init__)	



	# ####### Class Methods:

	def ReturnQuantile(self,s,d):
		'''uses simple linear interpolation to simulate the graphical estimation of a
		   single quantile value. Quantile values are required to compute sediment sample 
		   statistics, in the style of Folk (1980). See self.QuantilesList for required...
		   
		   Input Args:
			   s  sample identifier  as a Python string (EX. 'S1')
			   d  quantile value to be interpolated as integer
			   
			Returns:
				the interpolated quantile value in phi units
				
		   NOTE this is an internal function called by InterpolateQuantileValues(). 
		   The user can call this method, but this is discouraged...I don't know why.
		'''
		p=self.df.query(s+'cwp'+'<'+str(d) ) # query cumulative weight % df for lower bnd weight % recs
		q=self.df.query(s+'cwp'+'>'+str(d) ) # query cumulative weight % df for lower bnd weight % recs
		if p.empty:						  # if d < largest screen size (phi) in column, return -9
			return(np.nan)
		else:								# otherwise, interpolate...
			A=float(p[s+'cwp'].iloc[[-1]])
			B=float(q[s+'cwp'].iloc[[0]])
			C=float(p['phi'].iloc[[-1]])
			D=float(q['phi'].iloc[[0]])
			return( (((d - A)/(B - A))*( D-C ))+C )



	def InterpolateQuantileValues(self):
		'''interpolate quantiles in self.quantilesList for each transect sample. Calls function 
                   ReturnQuantile() for each quantile value to be computed from the seived weights

		   Input args: none

		   Returns:
			   an ordered dictionary (key=sample identifier; value=interpolated quantile 
			   list)
		'''
		for s in self.samplesList:
			tmpList=[]
			for d in self.quantilesList:
				tmpList.append( round(self.ReturnQuantile(s,d),2) )
			self.quantilesDict[s]=tmpList

		return(self.quantilesDict)



	def GetDataFrame(self):
		'''returns the sample set data frame df to the caller.

		   Input args: none

		   Returns: pandas dataframe for the current sample data
		'''
		return(self.df)



	def ComputeGraphicStats(self,s,Q):
		'''computes "graphic" mean, median, standard deviation, skewness, and kurtosis 
		   (per Folk, 1980) for a given transect sample.

		   Input args:
			   s  sample identifier as a Python string (EX. 'S1')
			   Q  list of interpolated quantile values (9 items) for current sample

		   Returns:
			   Python tuple containing MN=mean,MD=median,SD=std. dev.,SK=skewness,K=kurtosis in phi units.
		'''
		MN=(Q[2]+Q[4]+Q[6])/3	   # for the graphic mean:
		MD=Q[4]					 # for the graphic median:
		SD=((Q[6]-Q[2])/4)+((Q[8]-Q[0])/6.6) # inclusive graphic standard deviation:

		# the graphic skewness (SK) and kurtosis (K):
		SK=((Q[6]+Q[2]-2*Q[4])/(2*(Q[6]-Q[2]))) + ((Q[8]+Q[0]-2*Q[4])/(2*(Q[8]-Q[0])))
		K=(Q[8]-Q[0])/(2.44*(Q[5]-Q[3]))											   

		return( (MN,MD,SD,SK,K) )



	def ComputeMomentStats(self,s):
		'''computes the sediment sample mean and standard deviation using the arithmetic 
		   method of moments method as described by Folk (1980) for current sample.

		   Input args:
			   s = sample identifier as a Python string (EX. 'S1')

			Returns:
				Python tuple containing momMean=mean,momSD=standard deviation,
						,undifCF=undifferentiated coarse fraction (the weight of material
						 captured by the coarsest sieve), undifFF=undifferentiated fine 
						 fraction (the weight of material captured by the residual pan at 
						 the bottom of the seive stack).
		   '''

		# make unique copies of the screens and sample fraction weight lists, then remove 
		# the pan fractions (the undifferentiated fine stuff that ends up in the 
		# "catch-all"  pan at the bottom of the seive stack). Also, remove the coarsest 
		# fraction from the weights:
		S=list(self.df['phi'].copy(deep=True))   # make a copy of the screens column vect.
		S.remove(S[-1])			  # strip off the pan fraction (fines), the last item 
									 # in the screens list
		W=list(self.df[s].copy(deep=True))  # repeat for the fractional weights as you did 
											# with the screens

		# before we strip off the coarse and fine fractions in order to compute the MoM 
		# mean and std dev we want to capture the undifferentiated segment of the sample. 
		# The undifferentiated fraction consists of the sand remaining in the coarsest 
		# seive and that that falls thru into the pan at the bottom of the stack. As these
		# fractions are unbounded on one end we can't interpolate quantile values from
		# them. If, later we add the ability to extrapolate, these fractions will be more 
		# useful.  For now, however, we just report them back to the user.
		undifCF=W[0]
		undifFF=W[-1]

		# now we strip the undiff off of the raw weights list
		W.remove(W[0])			  # strip off the coarsest fraction--it's 
									# undifferentiated (unbounded)
		W.remove(W[-1])			 # strip off the finest (pan) fraction--it's 
									# undifferentiated, too

		# convert the weights list W into an array:
		W=np.array(W)

		# "loop" thru the screens copy list S and compute the screen midpoints. Turn the resulting list into an array
		midPt=np.array( [(float(S[i])+float(S[i+1]))/2 for i in range(len(S)-1)] )

		# compute the mean sediment size using Folk (1980)'s Method of Moments (page 46) 
		momMean=(midPt * W).sum() / W.sum()
		# compute the sediment sample size standard deviation using Folk(1980)'s Method of Moments (page 46)
		momSD=( (W*(momMean-midPt)**2).sum() / W.sum())**0.5

		# we want to provide back to the caller not only the mean and sd but also the 
		# undifferentiated fractions, too
		return( (momMean,momSD,undifCF,undifFF)) 



	def FindSampleModes(self, s):
		'''locates any modal weight percentage values in current sample. 

		   Input args: 
			   s = sample identifier as a Python string (EX. 'S1)

		   Returns: Python ordered dictionary of modes (if any) found in the current 
		   			sample. The dictionary key contains the mode in phi size units. 
		   			The dictionary values are the modes in weight percentages.
		   '''
		# mode dictionary: key=phi size; value=weight percentage
		modes=collections.OrderedDict()

		# sort the weight percentages in descending order
		B=sorted(list( self.df[s+'wp'] ), reverse=True)

		# HOW THIS THING WORKS:
		# first time thru it reads the 0th item in the sorted weight list B and posts this
		# to the modes dictionary (this is the largest or primary mode). Then, it finds 
		# any additional modes by taking each weight percentage value in the sorted list
		# B, locating that weight value and its index in the weight percentage column 
		# (df[s+'wp'] of the df, then comparing that weight percent value to those 
		# immediately before and immediately after the current. If both adjacent weights 
		# are less than the current weight we have another mode, which we insert into the 
		# modes dictionary.
	  
		for i in range(len(self.df[s+'wp'])):
			# first time thru we collect the primary mode. This is the largest
			# weight percentage in the sample, which, in the sorted list B, is
			# the first item in that list...
			if i == 0:
				modes[ self.screens[ list( self.df[s+'wp'] ).index(B[i]) ]] = (B[i])   
			# once the primary mode is plucked from the data search the weights
			# for others...
			else:	 
				Amax=B[i]
				Amax_index=list( self.df[s+'wp'] ).index(Amax)
				# as long as Amax_index is not allowed to exceed the length of the w
				# weight percentage column in df, less one. If it equals the length of
				# the column then, when you offset by +1 in searching for a mode you
				# generate an out of range index...bomb! If the index is okay, we can 
				# proceed to check for a mode by looking at wt values surrounding the
				# Amax value.
				if( (Amax_index < len(self.df[s+'wp'])-1) and \
				(self.df[s+'wp'][Amax_index-1] < Amax) and \
				(self.df[s+'wp'][Amax_index+1] < Amax) ):
						modes[self.screens[Amax_index]]=Amax
					
		return(modes)



	def Analyze2CSVHeader(self):
		''' generates a columns header string (names of each data column written by method Analyze2CSV)
                    for output data file created in method Analyze2CSV
		
			Input Args: none
			
			Returns: header as a Python string
		'''
		ws='Sample,'
		# craft a columns header and write to write string	
		for i in range(3):
			if(i==0):	
				w=['\''+str(h)+'phi_rwts\',' if(h != 4.0) else '\'fines_rwts\'' for h in self.screens ]
				ws=ws+''.join(w)
			if(i==1):
				ws=ws+','
				w=['\''+str(h)+'phi_pwts\',' if(h != 4.0) else '\'fines_pwts\'' for h in self.screens ]
				ws=ws+''.join(w)
			if(i==2):
				ws=ws+','
				w=['\''+str(h)+'phi_cpwts\',' if(h != 4.0) else '\'fines_cpwts\'' for h in self.screens ]
				ws=ws+''.join(w)
		ws=ws+','+'\'g_mean\''+','+'\'g_median\''+','+'\'g_std\''+','+'\'g_skew\''
		ws=ws+','+'\'g_kurt\''+','+'\'m_mean\''+','+'\'m_std\''+','+'\'undif_cor\''
		ws=ws+','+'\'undif_fin\''+','+'\'mode1_phi\''+','+'\'mode1_wt\''
		ws=ws+','+'\'mode2_phi\''+','+'\'mode2_wt\''
		ws=ws+','+'\'mode3_phi\''+','+'\'mode3_wt\''

		return(ws)



	def Analyze2Stdout(self):
		''' performs internal processing (without need for the user to call the individual 
		methods explicitly/directly) to compile the raw data, weight percentages, 
		cumulative weight percentages, graphic and method-of-moment statistics, and
		the mode(s)  All of the these are written directly, formatted to standard out.
		
		Input args:
		
		Returns: none
		'''
		# compute the quantile values from the raw weights:
		D=self.InterpolateQuantileValues()
					
		for s, Q in D.items():
			if( math.isnan( self.df[s].sum()) == True):  #pandas.isnull(self.df[s]) )
				print('Skipping sample:',self.transect,s, 'because of missing data...')
			else:
				print('')
				print('Sample:', str(s))
				self.PrintGraphicStats(s, self.ComputeGraphicStats(s, Q))
				self.PrintMomentStats(self.ComputeMomentStats(s))
				self.PrintSampleModes(self.FindSampleModes(s))

				print('')
		return()
		
		
	def Analyze2CSV(self, fn):
		''' performs internal processing (without need for the user to call the individual 
			methods explicitly/directly) to compile the raw data, weight percentages, 
			cumulative weight percentages, graphic and method of moment statistics, and
			the mode(s). All of these are written directly to an output file (user
			specified path and file name) formatted as comma separated values. 

			Input args: 
				fn = user-defined output file path and name

			Returns: none
		'''
		# call Analyze2CSVHeader to generate a header string for the output file
		hdr=self.Analyze2CSVHeader()

		D=self.InterpolateQuantileValues()
		# to generate the requisite content that will be written to the output file, call
		# and run the following methods automatically:

		# open the target text file and write the header string:
		f=open(fn, 'a' )
		f.write(hdr + '\n')

		# start loop here:
		for s, Q in D.items():
			if( math.isnan( self.df[s].sum()) == True):  #pandas.isnull(self.df[s]) )
				print('Skipping sample:',self.transect,s, 'because of missing data...')
			else:
				ws=s+','
				gStats=self.ComputeGraphicStats(s, Q)
				mStats=self.ComputeMomentStats(s)
				modes=self.FindSampleModes(s)

				# comprehensions extract data frame columns (series), round values to
				# places, convert to strings, and then add (concatenate) to write
				# string ws...
				w=[str(round(h,3))+',' for h in self.df[s]]		# raw wt series to str
				ws=ws+''.join(w)
				w=[str(round(h,3))+',' for h in self.df[s+'wp']]   # wt% series to str
				ws=ws+''.join(w)
				w=[str(round(h,3))+',' for h in self.df[s+'cwp']]  # cumwt% series to str
				ws=ws+''.join(w)

				w=[str(round(h,3))+',' for h in gStats]
				ws=ws+''.join(w)
				w=[str(round(h,3))+',' for h in mStats]
				ws=ws+''.join(w)
				ws=ws.rstrip(',')   # remove that last residual comma

				# get the modes, up to the first three, to add to out file:
				keys=list(modes.keys())
				vals=list(modes.values())

				w=''
				for i in range(3):
					try:
						key=keys[i]; val=vals[i]
						w=str(key)+','+str(round(val,3))
						ws=ws+','+w
					except:
						ws=ws+','

				print('Writing data for sample:', s, 'to file:', fn)
				f.write(ws + '\n')
				#print('')
				#print( ws )

		return()



	def PrintSampleWeightsDataTable(self, s):
		'''prints the raw, weight percentage, and cumulative weight percentage values in 
			tabular format for each screen bin to the console.

			Input args: 
				s = sample identifier as a Python string (EX. 'S1)

			Returns: none
		'''
		# ## print the report header
		print('-'*62)
		print('-'*62)
		print('Sample Weights Table','	 Transect: ',self.transect,'	 Sample: ',s )
		print('-'*62)
		# create temporary data frame table1df and populate with screen values (phi)
		table1df=pandas.DataFrame( self.df['phi'])		 
		# add the raw fraction weights to the df
		table1df['Raw Wt (gm)']=round(self.df[s].astype(float),2)   
		# add the fraction weight percents and cumulative weight percents to df
		table1df['Wt.%']=round(self.df[s+'wp'].astype(float),2)		
		table1df['Cum. Wt.%']=round(self.df[s+'cwp'].astype(float),2) 
		print(table1df)
		print('-'*62)
		print('Sample Quantiles:')
		print('( D5,   D10,  D16,  D25,  D50,  D75,  D84,  D90,  D95)')
		print( self.quantilesDict.get(s) )
		print('-'*62)



	def PrintGraphicStats(self, s, gStats):  
		'''prints the graphic mean, median, standard deviation, skewness, and kurtosis 
			for each transect sample in tabular format to the console. 

			Input args:
				s = sample identifier as a Python string (EX. 'S1)
				gStats = Python tuple containing the graphic mean, median, std. dev.
						skewness and kurtosis.

			Returns: none

			 Note: the gStats tuple is generated by method self.ComputeGraphicStats()
		'''
		print('-'*80)
		print('Graphic-Derived Sample Grain Size Statistics (Folk, 1980), Transect:',self.transect,'  Sample:',s)
		if np.isnan(gStats[0]):
			print('Est. Graphic Mean:','Out of Range*')
		else:
			print('Est. Graphic Mean', round(gStats[0],2),'Phi' )

		if np.isnan(gStats[1]):
			print('Est. Graphic Median:','Out of Range*')
		else:
			print('Est. Graphic Median', round(gStats[1],2),'Phi' )

		if np.isnan(gStats[2]):
			print('Est. Graphic Standard Deviation:','Out of Range*')
		else:
			print('Est. Graphic Standard Deviation', round(gStats[2],2),'Phi' )		

		if np.isnan(gStats[3]):
			print('Est. Graphic Skewness:','Out of Range*')
		else:
			print('Est. Graphic Skewness', round(gStats[3],2) )

		if np.isnan(gStats[4]):
			print('Est. Graphic Kurtosis:','Out of Range*')
		else:
			print('Est. Graphic Kurtosis', round(gStats[4],2) )

		print('-'*80)
		print('* out of range values indicate that one or more quantiles fall outside the screen')
		print('  range along the cumulative frequency curve.'+'\n')
		print('-'*80)



	def PrintMomentStats(self, mStats):
		'''prints the method of moments computed sediment size mean and standard deviation and undifferentiated coarse and fine fractions for the sample in tabular format to the console.

			Input args: 
				mStats = Python tuple containing method of moments mean, std. dev., 
						undiff. coarse fraction weight, undiff fine fraction weight.

			Returns: none

			Note: the mStats tuple is generated by method self.ComputeMomentStats()
		'''
		# ## print the report header
		print('-'*80)
		print('Method of Moments-derived Sample Grain Size Statistics (Folk, 1980)')
		print('')
		print('Mean:', round(mStats[0],2),'phi', ' Standard Deviation:', round(mStats[1],2),'phi')
		print('')
		print('Coarse Fraction, undifferentiated: ', round(mStats[2],2),'gm')
		print('Pan Fraction (fines), undifferentiated: ', round(mStats[3],2),'gm')
		print('-'*80)



	def PrintSampleModes(self, modes):
		'''Prints the mode values located in method FindSampleModes to the console

		   Input args:
			   modes = Python dictionary containing key (mode in phi units) and value
					   (mode in weight percent) pairs

		   Returns: none
		'''
		print('Mode(s):')
		print('Phi', ' Weight %')
		print('-------------------')
		for k,v in modes.items():
			print(k,str(round(v,3)))
		print('-------------------')

		return()



	def PLOTSampleWtPercents(self, s, mode='print'):
		'''plots the individual sample weight percentages, by sieve fraction 
			as a histogram and overprinted frequency (PDF) curve for the current sample. 
			Can plot to console or to stored PNG file.

			Input args:
				s = sample identifier as a Python string (EX. 'S1)
				mode = plot destination: mode='print' plot written to console (default)
											   mode='save' plot written to png output

			Returns: none
		'''
		plt.style.use('ggplot')
		bins = np.arange(len(self.screens))
		w=0.75
		plt.bar(bins,self.df[s+'wp'], align='center', width=w)
		plt.plot(bins,self.df[s+'wp'],'g--',linewidth=3.5)
		plt.xlim(-0.25,len(self.screens))
		#plt.ylim(0,105)
		xTickLabels=list(map(str, self.screens))
		xTickLabels[-1]='fines'
		plt.xticks(bins, xTickLabels )
		plt.xlabel('Sediment Size (phi)')
		plt.ylabel('Weight Percent')
		plt.title('Histogram and Frequency Curve'+'  Transect: '+self.transect+'   Sample: '+s)
		fig_name=self.transect+'_'+s+'_wf.png'
		print('-'*80)
		if(mode=='print'):
			plt.show()
		if(mode=='save'):
			plt.savefig(fig_name)
			plt.close()



	def PLOTSampleCumWtPercents(self, s, mode='print'):
		'''plots the individual sample cumulative weight percentages, by sieve 
		   fraction as a cunulative frequency (CDF) curve for the current transect sample.
		   Can plot to console or to stored PNG file.

			Input args:
				s = sample identifier as a Python string (EX. 'S1)
				mode = plot destination: mode='print' plot written to console (default)
											   mode='save' plot written to png output
											   
			Returns: none
		'''
		plt.style.use('ggplot')
		bins = np.arange(len(self.screens))
		#w=0.75
		#plt.bar(bins,self.df[s+'cwp'], align='center', width=w)
		plt.plot(bins,self.df[s+'cwp'],'g--',linewidth=3.5)
		plt.xlim(-0.25,len(self.screens))
		#plt.ylim(0,105)
		xTickLabels=list(map(str, self.screens))
		xTickLabels[-1]='fines'
		plt.xticks(bins, xTickLabels )
		plt.xlabel('Sediment Size (phi)')
		plt.ylabel('Weight Percent')
		plt.title('Cumulative Frequency Curve-Linear Ordinate'+'  Transect: '+self.transect+'   Sample: '+s)
		fig_name=self.transect+'_'+s+'_cwf.png'
		print('-'*80)
		if(mode=='print'):
			plt.show()
		if(mode=='save'):
			plt.savefig(fig_name)
			plt.close()



	def PLOTDualSampleWtPercents(self, s, mode='print'):
		'''plots both the weight percentage and cumulative weight percentage histogram 
			and curves (histo+PDF, and CDF) side by side and together, by sieve fraction 
			for the current transect sample. Can plot to console or to stored PNG file.

			Input args:
				s = sample identifier as a Python string (EX. 'S1)
				mode = plot destination: mode='print' plot written to console (default)
											   mode='save' plot written to png output
											   
			Returns: none

			Note: 
				This is just a convenient combination of methods PLOTSampleWtPercents and 
				PLOTSampleCumWtPercents.
		'''
		plt.style.use('ggplot')
		bins = np.arange(len(self.screens))
		w=0.75
		xTickLabels=list(map(str, self.screens))   # convert list of numeric screen sizes to list of strings for labels
		xTickLabels[-1]='fines'

		# write header to console:
		#print('-'*90)
		#print('  Transect: '+self.transect+'   Sample: '+s)
		print('')

		fig1=plt.figure(figsize=(12,4))

		# weight percentages subplot
		ax1=fig1.add_subplot(1,2,1)
		ax1.bar(bins,self.df[s+'wp'], align='center', width=w)
		ax1.plot(bins,self.df[s+'wp'],'g--',linewidth=3.5)
		ax1.set_xlim(-0.25,len(self.screens))
		#ax1.set_ylim(0,105)
		ax1.set_xticklabels(xTickLabels)
		ax1.set_xticks(bins)
		ax1.set_xlabel('Sediment Size (phi)')
		ax1.set_ylabel('Weight Percent')
		ax1.set_title('Histogram and Frequency Curve')

		# cumulative weight percentages subplot:
		ax2=fig1.add_subplot(1,2,2)
		#ax2.bar(bins,self.df[s+'cwp'], align='center', width=w)   # uncomment to draw bars
		ax2.plot(bins,self.df[s+'cwp'],'g-', linewidth=2.0)
		ax2.set_xlim(-0.25,len(self.screens))
		ax2.set_ylim(0,105)
		ax2.set_xticklabels(xTickLabels)
		ax2.set_xticks(bins)
		ax2.set_yticks([10,20,30,40,50,60,70,80,90,100])
		#ax2.minorticks_on()
		ax2.set_xlabel('Sediment Size (phi)')
		ax2.set_ylabel('Cumulative Weight Percent')
		ax2.set_title('Cumulative Frequency Curve-Linear Ordinate') #: Transect '+transect+'  Sample '+s)
		#ax2.grid(which='both')

		fig_name=self.transect+'_'+s+'_dual.png'
		if(mode=='print'):
			plt.show()
		if(mode=='save'):
			plt.savefig(fig_name)
			plt.close()

		print('-'*90)

# ## ####################################################################################
# ## End of SedSASample class listing 
# ## ####################################################################################





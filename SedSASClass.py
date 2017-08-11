#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug  3 12:01:40 2017

@author: paulp
"""

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

    def __init__(self,df,transect,sample,screens):
        
         # copy class user inputs to class variables:
         self.df=df
         self.transect=transect
         self.sample=sample
         self.screens=np.array(screens).T
         
         # list of quantiles to be estimated (interpolated) as required by Folk
         # and Ward graphical statistics...
         self.quantilesList=[5,10,16,25,50,75,84,90,95]
         
         ''' Build Main Dataframe df:
         reworks the user-supplied dataframe (as an input to the class) 
         putting it into a format that can be used in the various methods
         in the class. Here we also add weight percent and cumulative weight
         percent columns along with an integer size (phi) reference column.
        '''
         try:
             ## select record for current transect and sample from df:
            self.df=self.df.loc[ (self.df['Transect'] == self.transect) & \
                    (self.df['Sample'] == self.sample) ].copy()

            ## capture the dry weight field attribute:
            self.dry_wt=float( self.df['dry_weight'] )
                                
            self.df=self.df.T
            
            ## drop the Transect, Sample, and dry_weight fields from df:
            self.df=self.df.drop(['Transect','Sample', 'dry_weight'], axis=0)
            
            ## add column labels (names) to df:
            self.df.columns=[self.sample]
            
            ## compute and add the fraction weight percentages
            self.df[self.sample+'wp']=(self.df/self.df.sum(axis=0))*100 
            
            ## compute and add the cumulative fractional weight percentages:
            self.df[self.sample+'cwp']=self.df[self.sample+'wp'].cumsum()
            
            ## add integer field of sieve mesh (screen) sizes represented in df:
            self.df['phi']=self.screens
         except AttributeError as e:
            print('Oops! We had a problem creating/converting the dataframe df. Check the input dataframe formatting...:', e)
            
         try:
            ## compute and add the class size midpoints in phi units
            x=np.array( self.df['phi'])
            mp=(x[1:] + x[:-1]) / 2
            self.df['midpt_phi']=np.append(self.df['phi'][0:1].values-0.5, mp)    

            ## using the phi midpoints compute and add class size midpoints in
            ## mm
            self.df['midpt_mm']=2**(self.df['midpt_phi']*-1) 
         except:
            print('Oops!. We had a problem creating the midpoint columns')
            
         ''' Build the Fractional Dataframes (dfg, dfs, and dff):
         creates subset dataframes, each holding a particle size class
        (gravel, sand, slit/clay). The subset dataframes are used in lieu of
        the main df for computations when there is an excess of material 
        collected in either the coarsest (largest screen) or finest (pan) 
        fractions.
        gravel fraction (dfg), a sand fraction (dfs), and a fines fraction (dff)
        '''
         try:
            ### for the coarse fraction:
            self.dfg=self.df.loc[ self.df['phi'] <= -1 ].copy()
            self.dfg[self.sample+'wp']=(self.dfg/self.dfg.sum(axis=0))*100 
            self.dfg[self.sample+'cwp']=self.dfg[self.sample+'wp'].cumsum()
            
            ### for the sand fraction:
            self.dfs=self.df.loc[ (self.df['phi'] > -1)  &  \
                (self.df['phi'] < 4) ].copy()
            self.dfs[self.sample+'wp']=(self.dfs/self.dfs.sum(axis=0))*100 
            self.dfs[self.sample+'cwp']=self.dfs[self.sample+'wp'].cumsum()
            
            ### for the fine fraction:
            self.dff=self.df.loc[ self.df['phi'] >= 4 ].copy()
            self.dff[self.sample+'wp']=(self.dff/self.dff.sum(axis=0))*100 
            self.dff[self.sample+'cwp']=self.dff[self.sample+'wp'].cumsum()
         except:
            print('ERROR! There was a problem in class method _init_ trying to subset df based on particle size fractions (see Step 7.)')
        
         ## call InterpolateQuantileVales to generate the quantile list that
         ## is passed to the Folk and Ward computation methods, along with a
         ## few other methods in the class:
         try:
             self.Q = self.InterpolateQuantileValues(self.sample)
         except:
             print('Oops! Call to InterpolateQuantileValues failed.')
             print('Quantile list Q not created.')
         
         # alert user as to the current transect sample being processed:
         print('-'*80)
         print('Particle-Size Composition Analysis, Transect:', self.transect, 'Sample:', self.sample )
         print('-'*80)
        
    # ####### End of Constructor (__init__)    



    # ####### Class Methods:


    def GetDataFrames(self):
        '''returns the sample set data frames df,dfg,dfs,and dff to the caller.

           Input args: none

           Returns: a Python tuple containing 3 pandas dataframes: the first [0]
                    is the df for the entire sample dataset, the second [1] is a
                    subset containing only the gravel fraction, the third [2] 
                    holds the sand particle size fraction, and the fourth [3]
                    and final one is for the fine (silt/clay) fraction(s).
        '''
        try:
            return( (self.df,self.dfg,self.dfs,self.dff) )
        except AttributeError as e:
            print('Oops! [GetDataFrames]', e)
            
            
            
    
    def Summary(self, s):
        ''' Summary is a convenience method that calls the class methods:
            - PrintSampleWeightsDataTable
            - PrintUndifferentiatedSizeFractionWeights
            and goes on to compute the difference between the sample dry 
            sediment weight and the summed individual raw weights (in grams)
            reported for each seive mesh (screen) by the investigator in df: 
        '''
        self.PrintSampleWeightsDataTable(s)

        ## compute and report the difference between the dry sample weight
        ## reported in the data and the sum of individual binned (seived) 
        ## weights also in df:
        print('')
        print('Pre and post-process dry sample weights:')
        print('Reported pre-seive total dry sample weight:', \
              round(self.dry_wt, 3),'gms' )
        print('Computed sum of inidividual (Raw Wt (gm)) weights:', \
              round(self.df[s].sum(),3),'gms' )
        print('')
        
        self.PrintUndifferentiatedSizeFractionWeights(s)
            
        return()
        
        
        
        
    def CheckPercentofCoarseFineFractions(self, s):
        ''' checks the percentage of sand/sediment that was stopped by the 
        largest screen as well as the percentage captured in the pan at 
        the base of the sieve stack. If either of the percentages is > 5%,
        warn the user. If the latter, proceed further to replace existing 
        pan fraction with 0.0 gms. Why do this? See print statements below.
    
        Inputs: s - ID referencing the current sample 
    
        Returns: none
    
        NOTE this is an internal function called by ReturnQuantile(). The user 
        can call this method, but this is discouraged...I don't know why.
        '''
        ### check for excess material in the gravel and/or fines fractions:
        pcf=self.df[s][0]/self.df[s].sum()*100    # percent coarse fraction
        pff=self.df[s][-1]/self.df[s].sum()*100   # percent fine fraction
                
        if(pcf > 5.0):
            print('WARNING: percent of unidfferentiated coarse fraction in sample',s,'is:' ,str(round(pcf,2)),'which exceeds 5% of total by weight--As values in excess of 5% here can introduce significant error in resulting  computations, this fraction will not be considered in the analysis. NOTE THAT NO EXTRAPOLATED VALUES ARE USED IN THE FOLK AND WARD (FW) GRAPHiCAL METHOD COMPUTATIONS. \n')
        
        if(pff > 5.0):
            print('WARNING: percent of unidfferentiated fine fraction in sample',s, 'is:' ,str(round(pff,2)),'which exceeds 5% of total by weight--As values in excess of 5% here can introduce significant error in resulting  computations, this fraction will not be considered in the analysis. NOTE THAT NO EXTRAPOLATED VALUES ARE USED IN THE FOLK AND WARD (FW) GRAPHICAL METHOD COMPUTATIONS. \n')
        return()




    def ReturnQuantile(self,s,d):
        '''uses simple linear interpolation to simulate the graphical estimation 
        of a single quantile value. Quantile values are required to compute 
        sediment sample statistics, in the style of Folk (1980). See 
        self.QuantilesList for required...
           
        Input Args:
            s  sample identifier  as a Python string (EX. 'S1')
            d  quantile value to be interpolated as integer
               
        Returns:
            the interpolated quantile value in phi units
                
        NOTE this is an internal function called by 
           InterpolateQuantileValues(). The user can call this method, but this 
           is discouraged...I don't know why.
        '''
        p=self.df.query(s+'cwp'+'<'+str(d) )
        q=self.df.query(s+'cwp'+'>'+str(d) )
        if p.empty:            # if d < largest screen size (phi) in column, return -9
            return(np.nan)
        else:                                # otherwise, interpolate...
            A=float(p[s+'cwp'].iloc[[-1]])
            B=float(q[s+'cwp'].iloc[[0]])
            C=float(p['phi'].iloc[[-1]])
            D=float(q['phi'].iloc[[0]])
            return( (((d - A)/(B - A))*( D-C ))+C )




    def InterpolateQuantileValues(self, s):
        '''interpolate quantiles in self.quantilesList for each transect sample. 
        Calls function ReturnQuantile() for each quantile value to be computed 
        from the seived weights

        Input args: s - sample identifier  as a Python string (EX. 'S1')

        Returns: a Python list of interpolated quantiles 
        '''
        qList=[]
        self.CheckPercentofCoarseFineFractions(s)
        for d in self.quantilesList:
            qList.append( round(self.ReturnQuantile(s,d),2) )

        return(qList)



    ## computational method 1: Folk and Ward 'classi' graphic:
    def ComputeFWLogarithmicGraphicStats(self,s):
        '''computes "graphic" mean, median, standard deviation, skewness, and 
        kurtosis (per Folk and Ward, 1957 and Folk, 1980) for a given transect 
        sample. NOTE this is the 'classic' Folk and Ward (1957) and Folk (1980)
        particle size analysis method.

        Input args:
            s  sample identifier as a Python string (EX. 'S1')
            Q  list of interpolated quantile values (9 items) for current 
               sample in phi units

        Returns:
            Python tuple containing MN=mean, SD=std.deviation, 
            SK=skewness, K=kurtosis in phi units.
        '''
        Q=self.Q
        MN=(Q[2]+Q[4]+Q[6])/3                      # for the graphic mean:
        SD=((Q[6]-Q[2])/4)+((Q[8]-Q[0])/6.6)       # inclusive graphic 
                                                   # standard deviation:

        # the graphic skewness (SK) and kurtosis (K):
        SK=((Q[6]+Q[2]-2*Q[4])/(2*(Q[6]-Q[2]))) + \
           ((Q[8]+Q[0]-2*Q[4])/(2*(Q[8]-Q[0])))
        K=(Q[8]-Q[0])/(2.44*(Q[5]-Q[3]))

        return( (MN,SD,SK,K) )




    def ComputeFWGeometricGraphicStats(self,s):
        '''computes geometric "graphic" mean, median, standard deviation, 
        skewness, and kurtosis (per Folk and Ward, 1957) for a given transect 
        sample.

        Input args:
            s  sample identifier as a Python string (EX. 'S1')
            Q  list of interpolated quantile values (9 items) for current 
               sample in phi units

         Returns:
             Python tuple containing MN=mean, SD=std. deviation,
             SK=skewness, K=kurtosis in mm.
        '''
        ## convert quantiles in phi units to ln(mm) ( np.log(d(mm)=2^-phi) )
        ## where np.log is numpy's natural log (ln) and 2^-phi is the 
        ## conversion from phi to mm:
        
        Qm=[ np.log(2**((-1)*q)) for q in self.Q]
        
        # note that np.log is the natural logarithm
        MN=(Qm[2] + Qm[4] + Qm[6])/3                  # the graphic mean:
        SD=((Qm[6] - Qm[2])/4)+((Qm[8] - Qm[0])/6.6)  # inclusive graphic
                                                      # standard deviation:

        # the graphic skewness (SK):
        SK=(Qm[6] + Qm[2] - 2*Qm[4]) / (2*(Qm[6] - Qm[2])) \
        + (Qm[8] + Qm[0] - 2*(Qm[4])) /(2*(Qm[8] - Qm[0]))
        # and kurtosis (K)
        K=(Qm[8] - Qm[0]) / (2.44*(Qm[5] - Qm[3]))                                               

        return( (MN,SD,SK,K) )




    def ComputeArithmeticMethodofMomentsStats(self, s):
        '''computes the mean, sorting (standard deviation), skewness, and
        kurtosis (the first four statistical moments) using the method of
        moments technique (based on weight distributions--weight fractions
        captured in each seive). 
                           
        Inputs: s sample identifier as a Python string (EX. 'S1')
        
        Returns: Python tuple containing the four computed moments
        
        Units are mm, where appropriate
        '''
        ## compute the mean (M1) by airthmetic MoM:
        M1=( self.df[s+'wp']*self.df['midpt_mm'] ).sum() / self.df[s+'wp'].sum()
        
        ## compute the 2nd moment (M2 std. dev.) by airthmetic MoM:
        M2=( ((self.df[s+'wp']* (self.df['midpt_mm']-M1)**2 ).sum())/ 100)**0.5
        
        ## compute the 3rd moment (M3 skewness) by arithmetic MoM:
        M3=((self.df[s+'wp']* (self.df['midpt_mm']-M1)**3 ).sum())/(100*(M2**3))
        
        ## compute the 4th moment (M4 kurtosis) by arithmetic MoM:
        M4=((self.df[s+'wp']* (self.df['midpt_mm']-M1)**4 ).sum())/(100*(M2**4))
        
        return( (M1,M2,M3,M4) )
        
        
        
        
    def ComputeGeometricMethodofMomentsStats(self, s):
        '''computes the mean, sorting (standard deviation), skewness, and
        kurtosis (the first four statistical moments) using the geometric 
        method of moments technique (based on natural log of weight 
        distributions--weight fractions captured in each seive). 
                           
        Inputs: s sample identifier as a Python string (EX. 'S1')
        
        Returns: Python tuple containing the four computed moments
        
        Units are mm, where appropriate
        '''
        ## compute the mean (M1) by airthmetic MoM:
        M1=np.exp( ( self.df[s+'wp']*np.log(self.df['midpt_mm']) ).sum() / \
                  self.df[s+'wp'].sum() )
        
        ## compute the 2nd moment (M2 std. dev.) by airthmetic MoM:
        M2=np.exp( ( ((self.df[s+'wp']* (np.log(self.df['midpt_mm'])- \
                  np.log(M1))**2 ).sum())/ 100)**0.5 )
        
        ## compute the 3rd moment (M3 skewness) by arithmetic MoM:
        M3=((self.df[s+'wp']* (np.log(self.df['midpt_mm']) - \
                  np.log(M1))**3).sum())/(100*np.log(M2**3))
        
        ## compute the 4th moment (M4 kurtosis) by arithmetic MoM:
        M4=((self.df[s+'wp']* (np.log(self.df['midpt_mm']) - \
                  np.log(M1))**4 ).sum())/(100*np.log(M2**4))
        
        return( (M1,M2,M3,M4) )
        
        
        
        
    def ComputeLogarithmicMethodofMomentsStats(self, s):
        '''computes the mean, sorting (standard deviation), skewness, and
        kurtosis (the first four statistical moments) using the logarithmic 
        method of moments technique (based on weight distributions--weight 
        fractions captured in each seive). 
                           
        Inputs: s sample identifier as a Python string (EX. 'S1')
        
        Returns: Python tuple containing the four computed moments 
        
        Size units are in phi, where appropriate
        '''
        ## compute the mean (M1) by airthmetic MoM:
        M1=( self.df[s+'wp']*self.df['midpt_phi'] ).sum() / self.df[s+'wp'].sum()
        
        ## compute the 2nd moment (M2 std. dev.) by airthmetic MoM:
        M2=( ((self.df[s+'wp']* (self.df['midpt_phi']-M1)**2 ).sum())/ 100)**0.5
        
        ## compute the 3rd moment (M3 skewness) by arithmetic MoM:
        M3=((self.df[s+'wp']* (self.df['midpt_phi']-M1)**3 ).sum())/(100*(M2**3))
        
        ## compute the 4th moment (M4 kurtosis) by arithmetic MoM:
        M4=((self.df[s+'wp']* (self.df['midpt_phi']-M1)**4 ).sum())/(100*(M2**4))
        
        return( (M1,M2,M3,M4) )
        
        


    def FindSampleModes(self, s):
        '''locates any modal weight percentage values in current sample. 

           Input args: 
               s = sample identifier as a Python string (EX. 'S1)
               prt2stdout True/False(default) gives user the option to print 
               results to
               the console (std out) in addition to results in the returned dict.

           Returns: Python ordered dictionary of modes (if any) found in the 
               current sample. The dictionary key contains the mode in phi size 
               units. The dictionary values are the modes in weight percentages.
           '''
        # mode dictionary: key=phi size; value=weight percentage
        modes=collections.OrderedDict()

        # sort the weight percentages in descending order
        B=sorted(list( self.df[s+'wp'] ), reverse=True)

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
      
        for i in range(len(self.df[s+'wp'])):
            # first time thru we collect the primary mode. This is the largest
            # weight percentage in the sample, which, in the sorted list B, is
            # the first item in that list...
            if i == 0:
                modes[ self.screens[ list( self.df[s+'wp'] ).index(B[i]) ]] \
                   = (B[i])   
            # once the primary mode is plucked from the data search the weights
            # for others...
            else:     
                Amax=B[i]
                Amax_index=list( self.df[s+'wp'] ).index(Amax)
                # as long as Amax_index is not allowed to exceed the length of 
                # the w weight percentage column in df, less one. If it equals 
                # the length of the column then, when you offset by +1 in 
                # searching for a mode you generate an out of range index...bomb! 
                # If the index is okay, we can proceed to check for a mode by 
                # looking at wt values surrounding the Amax value.
                if( (Amax_index < len(self.df[s+'wp'])-1) and \
                (self.df[s+'wp'][Amax_index-1] < Amax) and \
                (self.df[s+'wp'][Amax_index+1] < Amax) ):
                        modes[self.screens[Amax_index]]=Amax

        return(modes)



    def PrintSampleWeightsDataTable(self, s):
        '''prints the raw, weight percentage, and cumulative weight percentage 
        values in tabular format for each screen bin to the console.

        Input args: 
            s = sample identifier as a Python string (EX. 'S1)

        Returns: none
        '''
        # ## print the report header
        print('-'*62)
        print('-'*62)
        print('Sample Weights Table','     Transect: ',self.transect,'     Sample: ',s )
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
        print('Sample Quantiles:  Transect:', self.transect, 'Sample:', s)
        print('( D5,  D10, D16, D25, D50, D75, D84, D90, D95)')
        print( self.Q )
        print('-'*62)




    def PrintUndifferentiatedSizeFractionWeights(self, s):
        '''
        prints the undifferentiated particle size fraction weights for the
        current sample. Reported are the undifferentiated coarse (gravel), 
        sand, and slit/clay (fines) particle fractions, along with a total
        weight.
        
        Inputs:  none
        Returns: none
        '''
        print('')
        print('-'*80)
        print('Particle-Size Composition, Transect:', self.transect, \
              'Sample', s,':')
        print('Undifferentiated gravel fraction (< -1 phi):', \
              str(self.dfg[s].sum()), 'gm')
        print('Undifferentiated sand fraction (-1 to 4 phi):', \
              str(self.dfs[s].sum()), 'gm')
        print('Undifferentiated fine fraction (>= 4 phi):', \
              str(self.dff[s].sum()), 'gm')
        print('Total Sample Weight:', self.df[s].sum(), 'gm')
        print('Median particle size:', self.Q[4], 'phi',\
              '  ',round(2**((-1)*self.Q[4]),3), 'mm' )
        print('')




    def PrintComputedStatistics(self, s, stats, units='phi', method='FWlog' ):
        '''prints the computed mean, median, standard deviation, skewness, and 
           kurtosis for each transect sample in tabular format to the console. 

            Input args:
                s = sample identifier as a Python string (EX. 'S1)
                stats = Python tuple containing the graphic mean,  std. 
                deviation, skewness, and kurtosis.
                units - phi (default) or mm. The particle size units--based on 
                        the statistics generated 
                method - the method used to compute particle size statistics.
                         Choose from one of the five class functions:
                             'FWlog' for ComputeFWLogarithmicGraphicStats()
                             'FWgeo' for ComputeFWGeometricGraphicStats()
                             'MoMar' for ComputeArithmeticMethodofMomentsStats
                             'MoMgeo' for ComputeGeometricMethodofMomentsStats
                             'Momlog' for ComputeLogarithmicMethodofMomentsStats
                        default is 'FWlog'

            Returns: none
        '''
        if(method=='FWlog'): method=('Method of Folk and Ward (Logarithmic)')
        if(method=='FWgeo'): method=('Method of Folk and Ward (Geometric)')
        if(method=='MoMar'): method=('Method of Moments (Arithmetic)')
        if(method=='MoMgeo'): method=('Method of Moments (Geometric)')
        if(method=='MoMlog'): method=('Method of Moments (Logarithmic)')
        
        print('-'*80)
        print('Sample Particle Size Statistics, Transect:',self.transect,'  Sample:',s)
        
        print('Computation:', method)
        print('')
        
        if np.isnan(stats[0]):
            print('Estimated Mean:','Out of Range*')
        else:
            print('Estmated Mean', round(stats[0],2), units )

        if np.isnan(stats[1]):
            print('Estmated Standard Deviation:','Out of Range*')
        else:
            print('Estmated Standard Deviation', round(stats[2],2), units )

        if np.isnan(stats[2]):
            print('Estmated Skewness:','Out of Range*')
        else:
            print('Estmated Skewness', round(stats[3],2) )

        if np.isnan(stats[3]):
            print('Estmated Kurtosis:','Out of Range*')
        else:
            print('Estmated Kurtosis', round(stats[4],2) )

        print('-'*80)
        print('* out of range values indicate that one or more quantiles fall outside the range along the cumulative frequency curve.'+'\n')
        print('-'*80)



    def PrintSampleModes(self, s, modes):
        '''Prints the sample modes located by method FindSampleModes to the
        console.
        
        Inputs: s - the ID referencing the current sediment sample
                modes - the mode(s) (up to the first 3) located by method
                FindSampleModes() as a Python dictionary.
        '''
        print('Mode(s) for:', s)
        print('Phi', ' Weight %')
        print('-------------------')
        for k,v in modes.items():
            print(k,str(round(v,3)))
        print('-------------------')
        
        return()



    def Analyze2CSVHeader(self):
        ''' generates a columns header string (names of each data column written 
        by method Analyze2CSV) for output data file created in method Analyze2CSV
        
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
        ''' performs internal processing (without need for the user to call the 
        individual methods explicitly/directly) to compile the raw data, weight 
        percentages, cumulative weight percentages, graphic and method-of-moment 
        statistics, and the mode(s)  All of the these are written directly, 
        formatted to standard out.
        
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
                self.ComputeLogGraphicStats(s, Q, True)
                self.PrintMomentStats(self.ComputeMomentStats(s))
                self.FindSampleModes(s, True)

                print('')
        return()
        
        
        
    def Analyze2CSV(self, fn):
        ''' performs internal processing (without need for the user to call the 
        individual methods explicitly/directly) to compile the raw data, weight 
        percentages, cumulative weight percentages, graphic and method of moment 
        statistics, and the mode(s). All of these are written directly to an 
        output file (user specified path and file name) formatted as comma 
        separated values. 

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
                w=[str(round(h,3))+',' for h in self.df[s]]        # raw wt series to str
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





    def PrintSampleModes(self, s, modes):
        '''Prints the mode values located in method FindSampleModes to the
           console

           Input args:
               modes = Python dictionary containing key (mode in phi units) and
               value (mode in weight percent) pairs

           Returns: none
        '''
        print('Mode(s): Transect', self.transect, 'Sample:', s )
        print('Phi', ' Weight %')
        print('-------------------')
        for k,v in modes.items():
            print(k,str(round(v,3)))
        print('-------------------')

        return()



    def PLOTSampleWtPercents(self, s, mode='print'):
        '''plots the individual sample weight percentages, by sieve fraction
           as a histogram and overprinted frequency (PDF) curve for the current
           sample. Can plot to console or to stored PNG file.

            Input args:
                s = sample identifier as a Python string (EX. 'S1)
                mode = plot destination: mode='print' plot written to console
                (default)
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
           fraction as a cunulative frequency (CDF) curve for the current
           transect sample. Can plot to console or to stored PNG file.

            Input args:
                s = sample identifier as a Python string (EX. 'S1)
                mode = plot destination: mode='print' plot written to console
                (default)
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
        '''plots both the weight percentage and cumulative weight percentage
           histogram and curves (histo+PDF, and CDF) side by side and together,
           by sieve fraction for the current transect sample. Can plot to
           console or to stored PNG file.

            Input args:
                s = sample identifier as a Python string (EX. 'S1)
                mode = plot destination: mode='print' plot written to console
                (default)
                mode='save' plot written to png output
                                               
            Returns: none

            Note: 
                This is just a convenient combination of methods
                - PLOTSampleWtPercents and 
                - PLOTSampleCumWtPercents.
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





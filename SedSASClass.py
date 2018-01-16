#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug  3 12:01:40 2017

@author: paulp
"""

# -*- coding: utf-8 -*-

# ## import prerequisite Python module dependencies:
import numpy as np
import pandas
import matplotlib.pyplot as plt
#from scipy.optimize import curve_fit
#import collections


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
            
            must supply a Transect column containing the transect identifier
            must supply a Sample column containing the sample identifier
            must supply a dry_weight column containing the total computed
               dry weight of the sample (value computed by summing all 
               individual raw weights)

    '''

    def __init__(self, sample_id, df, screens, extrap_method='Linear'):
        
         # copy class user inputs to class variables:
         self.df=df
         self.sample_id=sample_id
         self.screens=np.array(screens).T
         
         # list of quantiles to be estimated (interpolated) as required by Folk
         # and Ward graphical statistics...
         self.quantilesList=[5,10,16,25,50,75,84,90,95]
         
         # construct the master dataframe. This will be central to much of the
         # processing
         self.BuildWorkingDataFrame()
        
         ## call InterpolateQuantileValues to generate the quantile list that
         ## is passed to the Folk and Ward computation methods, along with a
         ## few other methods in the class:
         self.Q = self.InterpolateQuantileValues(extrap_method)
         #try:
         #    
         #except:
         #    print('Oops! Call to InterpolateQuantileValues failed.')
         #    print('Quantile list Q not created.')
         
         # alert user as to the current transect sample being processed:
         print('-'*70)
         print('Particle-Size Composition Analysis. Processing Sample ID:',
               self.sample_id )
         print('-'*70)
         
         ## set a flag to control the creation and writing of a text file 
         ## header line that is used by the method: Analysis2CSV. Initialize
         ## to True
         self.csvHeaderFlag=True
        
    # ####### End of Constructor (__init__)    



    # ####### Class Methods:

    def BuildWorkingDataFrame(self):
        ''' from the user-provide dataframe (as an argument to this function
        and WHICH CONTAINS ONLY THE RAW phi WEIGHTS--use indexing to send 
        only this portion of the original dataframe) and, using the raw 
        weights, computes the individual and cumulative weight percentages, 
        and the phi-class midpoint, adding each of these three to a new
        sample dataframe. This dataframe iis returned to the caller for use
        as input to other class methods.
        
        inputs: dataframe containing only the raw sediment sample weights
        
        returns: new dataframe with, along with a reproduction of the original 
        sieve fraction weights, the individual and cumulative weight percent-
        ages. The phi midpoint for each size class is also included.
        '''
        
        # first make a copy of the transposed user-supplied weight dataframe 
        self.sdf=self.df.T
        self.sdf.columns=['weight']

        ## add integer field of sieve mesh (screen) sizes represented in df:
        self.sdf['phi']=self.screens
        
        ## compute and add the class size midpoints in phi units
        x=np.array( self.sdf['phi'] )
        mp=(x[1:] + x[:-1]) / 2
        self.sdf['midpt_phi']=np.append(self.sdf['phi'][0:1].values-0.5, mp)    

        ## using the phi midpoints compute and add class size midpoints in mm
        self.sdf['midpt_mm']=2**(self.sdf['midpt_phi']*-1)
        
        ## compute and add the fraction weight percentages
        self.sdf['weight_percent']=(self.sdf['weight']/self.sdf['weight'].sum(axis=0))*100 
            
        ## compute and add the cumulative fractional weight percentages:
        self.sdf['cum_weight_percent']=self.sdf['weight_percent'].cumsum()
            
        return()
        
        

    def GetDataFrame(self):
        '''returns the sample set data frame sdf as generated by class method:
            BuildSampleDataFrame() under the class constructor __init__ to
            the caller

           Input args: none

           Returns: a  pandas dataframe.
        '''
        try:
            return( self.sdf )
        except AttributeError as e:
            print('Oops! [GetDataFrame]', e)
            
            
            
            
    def GetQuantileList(self):
        '''returns the sample set qunatiles list as generated by class method:
            InterpolateQuantileValues() under the class constructor __init__ to
            the caller

           Input args: none

           Returns: a  Python list.
        '''
        try:
            return( self.Q )
        except AttributeError as e:
            print('Oops! [GetQuantilList]', e)
            
            
            
            
    def Summary(self, s):
        ''' Summary is a convenience method that calls the class methods:
            - PrintSampleWeightsDataTable
            - PrintUndifferentiatedSizeFractionWeights
            and goes on to report the sample dry sediment weight and the summed 
            individual raw weights (in grams) reported for each seive mesh 
            (screen) by the investigator in df: 
        '''
        self.PrintSampleWeightsDataTable(s)

        ## compute and report the difference between the dry sample weight
        ## reported in the data and the sum of individual binned (seived) 
        ## weights also in df:
        print('')
        print('Pre and post-process dry sample weights:')
        print('Reported pre-seive total dry sample weight:', \
              round(self.dry_wt, 3),'gms' )
        print('Computed sum of individual (Raw Wt (gm)) weights:', \
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
        can call this method, but this is discouraged...
        '''
        ### check for excess material in the gravel and/or fines fractions:
        pcf=self.df[s][0]/self.df[s].sum()*100    # percent coarse fraction
        pff=self.df[s][-1]/self.df[s].sum()*100   # percent fine fraction
                
        if(pcf > 5.0):
            print('WARNING: percent of unidfferentiated coarse fraction in sample',s,'is:' ,str(round(pcf,2)),'which exceeds 5% of total by weight--values in excess of 5% here can introduce significant error in the analysis. \n')
            print('NOTE THAT SOME QUANTILE VALUES WILL BE DETERMINED BY EXTRAPOLATION. \n')
        
        if(pff > 5.0):
            print('WARNING: percent of unidfferentiated fine fraction in sample',s, 'is:' ,str(round(pff,2)),'which exceeds 5% of total by weight--values in excess of 5% here can introduce significant error in the analysis. \n')
            print('NOTE THAT SOME QUANTILE VALUES WILL BE DETERMINED BY EXTRAPOLATION. \n')
        return()



    def ExtrapolateQuantileLinearFit(self, q, qnt):
        '''
        '''
        ### check magnitude of cum wt pct for phi[0]. if it's 
        ### > 30% then we place the x-intercept at pho[0]-2
        ### if it's <= 30% then we place the x-intercept at 
        ### phi[0] - 1. This limits the likelihood of large 
        ### sample fractions in the topmost sieve yielding 
        ### exaggerated quantile estimates:
        print('Warning: Extrapolating', str(qnt), 'percent quantile')
        cwp1=float(q['cum_weight_percent'].iloc[0])
        if( cwp1 > 30.0 ):
            Xk_1=float(q['phi'].iloc[0]) - 2.0
            Yk_1=0.0
            Xk=float(q['phi'].iloc[0])
            Yk=float(q['cum_weight_percent'].iloc[0])
            
        elif( cwp1 > 10.0 and cwp1 <= 30.0 ):
            Xk_1=float(q['phi'].iloc[0]) - 1.0
            Yk_1=0.0
            Xk=float(q['phi'].iloc[0])
            Yk=float(q['cum_weight_percent'].iloc[0])
            
        else:
            Xk_1=float(q['phi'].iloc[0])
            Yk_1=float(q['cum_weight_percent'].iloc[0])
            Xk=float(q['phi'].iloc[1])
            Yk=float(q['cum_weight_percent'].iloc[1])
        #Xe=Xk_1+((Xk-Xk_1)*((qnt-Yk_1)/(Yk-Yk_1)))
        return( Xk_1+((Xk-Xk_1)*((qnt-Yk_1)/(Yk-Yk_1))) )
        

    def poly2f(self,x,a,b,c):
        return( a*x**2+b*x+c )
    

    def poly3f(self, x, a, b, c, d):     # 3rd order fit is really close-overfit?
        return( a*x**3+b*x**2+c*x+d ) 
        
        
        
    def ExtrapolateQuantilePolyFit(self, q, qnt, extrap_method='Poly2'):
        '''Provides a non-linear (polynomial) option to extrapolate quantile 
        values used in the graphic statistical computations. Fits a 2nd or 3rd
        degree polynomial (based on value of degree argument (domain: 2,3)) to
        data in dataframe q and uses the resulting polynomial function to
        extimate quantile values where the minimum measured cumulative weight 
        percentage is greater than the smallest quaantile(s). This typically
        occurs when samples are skewed toward the coarser fractions.
        
        This approach has the potential to provide more accurate results than 
        does a simple linear fit (see self.ExtrapolateQuantileLinearFit) but 
        the user should examine and compare results provided here with the
        linear returns.
        
        Inputs:
            q : dataframe containing all cumulative weights associated with
                seive 'bins' greated than the qnt value
            qnt : the current quantile value to be computed/determined
            degree : the degree or order of the poly function to be fit to the
                data. Acceptable values are 2 or 3.
                
        Returns: the extrapolated quantile value in phi units
        '''
        step=-0.001
        w=0
        X=np.array( q['phi'], dtype='float64').ravel()
        Y=np.array( q['cum_weight_percent'], dtype='float64').ravel()
        
        if(extrap_method=='Poly2'):
            params, pcov = curve_fit(self.poly2f, X, Y)
            A=params[0]; B=params[1]; C=params[2]

            # iterate to converge on closest guess to qnt +/- step:
            for guess in np.arange(X.min(), X.min()-2, step):
                r=A*guess**2+B*guess+C
                if(r <= qnt+0.1) & (r >= qnt-0.1):
                    w=guess
                    
        if(extrap_method=='Poly3'):
            params, pcov = curve_fit(self.poly3f, X, Y)
            A=params[0]; B=params[1]; C=params[2]; D=params[3]

            # iterate to converge on closest guess to qnt +/- step:
            for guess in np.arange(X.min(), X.min()-2, step):
                r=A*guess**3+B*guess**2+C*guess+D
                if(r <= qnt+0.1) & (r >= qnt-0.1):
                    w=guess
                    
        plt.figure()
        plt.plot(X, Y, 'ko', label="Original Data")
        plt.plot(X, self.poly3f(X, *params), 'r-', label="Fitted Curve")
        plt.plot((-3,3),(5,5))     # 5th quantile
        plt.plot((-3,3),(10,10))   # 10th quantile
        plt.ylim(0,100)
        plt.xlim(-3,4)
        plt.legend()
        plt.show()
        
        return( w )
        
        
        
    def InterpolateQuantileLinearFit(self, p, q, qnt):
        '''
        '''
        Yk_1=float(p['cum_weight_percent'].iloc[[-1]])
        Yk=float(q['cum_weight_percent'].iloc[[0]])
        Xk_1=float(p['phi'].iloc[[-1]])
        Xk=float(q['phi'].iloc[[0]])
        return(Xk_1+((Xk-Xk_1)*((qnt-Yk_1)/(Yk-Yk_1))) )
        
        
        
        
    def ReturnQuantile(self, qnt, extrap_method):
        '''Computes current (user-requested) quantile value: qnt. 
        Process simulates the graphical estimation of a single quantile value
        where the quantile, phi units, is read visually from a graph of 
        cumulate weight percentage for the current sample. Quantile values are 
        required to compute sediment sample statistics, in the style of 
        Folk (1980).  
           
        Input Args:
            qnt  quantile value to be interpolated/extrapolated : integer
               
        Returns:
            the interpolated or extrapolated quantile value in phi units
                
        NOTE this is an internal function called by 
           InterpolateQuantileValues(). The user can call this method, but this 
           is discouraged...
        '''
        p=self.sdf.query('cum_weight_percent'+'<'+str(qnt) )
        q=self.sdf.query('cum_weight_percent'+'>'+str(qnt) )
        if(p.empty) & (q.empty):     # if both sub dfs are empty, abort
            return(np.nan)
        if p.empty:            # we must extrapolate to estimate quantile
            if(extrap_method=='Linear'):
                qval=self.ExtrapolateQuantileLinearFit(q, qnt)
            if(extrap_method=='Poly2') or (extrap_method=='Poly3'):
                qval=self.ExtrapolateQuantilePolyFit(q, qnt, extrap_method)

        else:                  # otherwise, (linear) interpolate quantile:
            qval=self.InterpolateQuantileLinearFit(p, q, qnt)
 
        return( qval )


   

    def InterpolateQuantileValues(self, extrap_method):
        '''interpolate quantiles in self.quantilesList for each transect sample. 
        Calls function ReturnQuantile() for each quantile value to be computed 
        from the seived weights

        Input args: extrap_method - the type of extrapolation to use, if 
                    required. Options are: 'Linear' (default), 'Poly2' second
                    degree (quadratic) polynomial fit to cumulative weights, 
                    Poly3 a third order polynomial fit to the cumulative
                    weights in dataframe q (q is extracted from df in method
                    ReturnQuantile())

        Returns: a Python list of interpolated quantiles 
        '''
        qntList=[]
        ###self.CheckPercentofCoarseFineFractions(s)
        for qnt in self.quantilesList:
            qntList.append( round(self.ReturnQuantile(qnt, extrap_method),3))
           
        return(qntList)



    ## computational method 1: Folk and Ward 'classic' graphic:
    def ComputeFWLogarithmicGraphicStats(self):
        '''computes "graphic" mean, standard deviation, skewness, and 
        kurtosis (per Folk and Ward, 1957 and Folk, 1980) for a given transect 
        sample. NOTE this is the 'classic' Folk and Ward (1957) and Folk (1980)
        particle size analysis method (traditionally using hand-drawn
        probability ordinate cumulative frequency plots).

        Input args:
            None

        Returns:
            Python list containing MN=mean, SD=std.deviation, 
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

        return( [MN,SD,SK,K] )




    def ComputeFWGeometricGraphicStats(self):
        '''computes geometric "graphic" mean, standard deviation, 
        skewness, and kurtosis (per Folk and Ward, 1957) for a given transect 
        sample.

        Input args:
            Q  list of interpolated quantile values (9 items) for current 
               sample in phi units

         Returns:
             Python tuple containing MN=mean, SD=std. deviation,
             SK=skewness, K=kurtosis in mm.
        '''
        ## convert quantiles in phi units to ln(mm) ( np.log(d(mm)=2^-phi) )
        ## where np.log is numpy's natural log (ln) and 2^-phi is the 
        ## conversion from phi to mm:
        Qp=self.Q
        Qm=[2**(-1*i) for i in Qp]
                
        # note that np.log is the natural logarithm
        MN=np.exp( (np.log(Qm[2]) + np.log(Qm[4]) + np.log(Qm[6]))/3 ) # the graphic mean:
        SD=np.exp( ((np.log(Qm[2]) - np.log(Qm[6]))/4) + ((np.log(Qm[0]) - \
         np.log(Qm[8]))/6.6) )  # inclusive graphic standard deviation:

        # the graphic skewness (SK):
        A=(np.log(Qm[2]) + np.log(Qm[6]) - 2*np.log(Qm[4])) / (2*(np.log(Qm[6]) - np.log(Qm[2])))
        B=(np.log(Qm[0]) + np.log(Qm[8]) - 2*np.log(Qm[4])) / (2*(np.log(Qm[8]) - np.log(Qm[0])))
        SK=A+B
        
        # and kurtosis (K)
        K=(np.log(Qm[0]) - np.log(Qm[8])) / (2.44*( np.log(Qm[3]) - np.log(Qm[5])))

        return( [MN,SD,SK,K] )




    def ComputeArithmeticMethodofMomentsStats(self):
        '''computes the mean, sorting (standard deviation), skewness, and
        kurtosis (the first four statistical moments) using the method of
        moments technique (based on weight distributions--weight fractions
        captured in each seive). 
                           
        Inputs: None
        
        Returns: Python tuple containing the four computed moments
        
        Units are mm, where appropriate
        '''
        ## compute the mean (M1) by airthmetic MoM:
        M1=( self.sdf['weight_percent']*self.sdf['midpt_mm'] ).sum() / self.sdf['weight_percent'].sum()
        
        ## compute the 2nd moment (M2 std. dev.) by airthmetic MoM:
        M2=( ((self.sdf['weight_percent']*(self.sdf['midpt_mm']-M1)**2 ).sum())/ 100)**0.5
        
        ## compute the 3rd moment (M3 skewness) by arithmetic MoM:
        M3=((self.sdf['weight_percent']*(self.sdf['midpt_mm']-M1)**3 ).sum())/(100*(M2**3))
        
        ## compute the 4th moment (M4 kurtosis) by arithmetic MoM:
        M4=((self.sdf['weight_percent']*(self.sdf['midpt_mm']-M1)**4 ).sum())/(100*(M2**4))
        
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
        M3=(self.df[s+'wp']*((np.log(self.df['midpt_mm']) - \
                  np.log(M1))**3)).sum() / (100*np.log(M2)**3)
        
        ## compute the 4th moment (M4 kurtosis) by arithmetic MoM:
        M4=(self.df[s+'wp']*(np.log(self.df['midpt_mm']) - \
                  np.log(M1))**4 ).sum()/(100*np.log(M2)**4)
        
        return( (M1,M2,M3,M4) )
        
        
        
        
    def ComputeLogarithmicMethodofMomentsStats(self):
        '''computes the mean, sorting (standard deviation), skewness, and
        kurtosis (the first four statistical moments) using the logarithmic 
        method of moments technique (based on weight distributions--weight 
        fractions captured in each seive). 
                           
        Inputs: None
        
        Returns: Python tuple containing the four computed moments 
        
        Size units are in phi, where appropriate
        '''
        ## compute the mean (M1) by airthmetic MoM:
        M1=( self.df[s+'wp']*self.df['midpt_phi'] ).sum() / self.df[s+'wp'].sum()
        
        ## compute the 2nd moment (M2 std. dev.) by airthmetic MoM:
        M2=( (self.df[s+'wp']*((self.df['midpt_phi']-M1)**2)).sum()/ 100)**0.5
        
        ## compute the 3rd moment (M3 skewness) by arithmetic MoM:
        M3=((self.df[s+'wp']* (self.df['midpt_phi']-M1)**3 ).sum())/(100*(M2**3))
        
        ## compute the 4th moment (M4 kurtosis) by arithmetic MoM:
        M4=((self.df[s+'wp']* (self.df['midpt_phi']-M1)**4 ).sum())/(100*(M2**4))
        
        return( (M1,M2,M3,M4) )
        
        


    def FindSampleModesDepr(self, s):
        '''locates modal weight percentage values in current sample. 

           Input args: 
               s = sample identifier as a Python string (EX. 'S1)

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



    def FindSampleModes(self):
        '''locates modal weight percentage values in current sample. 

           Input args: None

           Returns: Python tuple of modes (if any) found in the 
               current sample. The 0th item contains the mode in phi size 
               units. The 1st item is the mode in mm.
           Note that this method is a revised/updates, and hopefully improved,
           version of the existing, and original, FindSampleModes method, which
           is still in service. If this rev'ed version proves to be a better 
           device than the original, the latter will be deprecated and removed.
           Otherwise...
        '''
        #self.df.sort_values( [s+'wp'], axis=1)
        # initialize an empty modes list:
        modes=[]
        
        # sort the weight percentages column in df in descending 
        # order, put into list B:
        B=sorted(list( self.df[s+'wp'] ), reverse=True)
        
        for i in range(len(B)):
            if i==0:                   # for primary mode
                phi=self.screens[ list( self.df[s+'wp'] ).index(B[0])]
                mm=2**(phi*-1)
                modes.append( (phi,mm) )
            else:
                Amax=B[i]
                Amax_index=list( self.df[s+'wp'] ).index(Amax)
                if( (Amax_index < len(self.df[s+'wp'])-1) and \
                   (self.df[s+'wp'][Amax_index-1] < Amax) and \
                   (self.df[s+'wp'][Amax_index+1] < Amax) ):
                    phi=self.screens[Amax_index]
                    mm=2**(phi*-1)
                    modes.append( (phi,mm) )  
        return(modes)
           


    def PrintSampleWeightsDataTable(self):
        '''prints the raw, weight percentage, and cumulative weight percentage 
        values in tabular format for each screen bin to the console.

        Input args: none

        Returns: none
        '''
        # ## print the report header
        print('-'*62)
        print('-'*62)
        print('Sample Weights Table','     Sample: ', self.sample_id )
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
        print('Undifferentiated gravel fraction (<= -1 phi):', \
              str(self.dfg[s].sum()), 'gm')
        print('Undifferentiated sand fraction (-1 to 4 phi):', \
              str(self.dfs[s].sum()), 'gm')
        print('Undifferentiated fine fraction (>= 4 phi):', \
              str(self.dff[s].sum()), 'gm')
        print('Total Sample Weight:', self.df[s].sum(), 'gm')
        print('Median particle size:', self.Q[4], 'phi',\
              '  ',round(2**((-1)*self.Q[4]),3), 'mm' )
        print('')




    def PrintComputedStatistics(self, stats, units='phi', method='FWlog' ):
        '''prints the computed mean, median, standard deviation, skewness, and 
           kurtosis for each transect sample in tabular format to the console. 

            Input args:
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
                             'MoMlog' for ComputeLogarithmicMethodofMomentsStats
                        default is 'FWlog'

            Returns: none
        '''
        if(method=='FWlog'): method=('Method of Folk and Ward (Logarithmic)')
        if(method=='FWgeo'): method=('Method of Folk and Ward (Geometric)')
        if(method=='MoMar'): method=('Method of Moments (Arithmetic)')
        if(method=='MoMgeo'): method=('Method of Moments (Geometric)')
        if(method=='MoMlog'): method=('Method of Moments (Logarithmic)')
        
        print('-'*80)
        print('Sample Particle Size Statistics, Sample:', self.sample_id)
        
        print('Computation:', method)
        print('')
        
        if np.isnan(stats[0]):
            print('Estimated Mean:','Out of Range*')
        else:
            print('Estmated Mean', round(stats[0],2), units )

        if np.isnan(stats[1]):
            print('Estmated Standard Deviation:','Out of Range*')
        else:
            print('Estmated Standard Deviation', round(stats[1],2), units )

        if np.isnan(stats[2]):
            print('Estmated Skewness:','Out of Range*')
        else:
            print('Estmated Skewness', round(stats[2],2) )

        if np.isnan(stats[3]):
            print('Estmated Kurtosis:','Out of Range*')
        else:
            print('Estmated Kurtosis', round(stats[3],2) )

        print('-'*80)
        print('* out of range values indicate that one or more quantiles needed for the computation fall outside  data range along the cumulative frequency curve.'+'\n')
        print('-'*80)




    def PrintSampleModes(self, s, modes):
        '''Prints the sample modes located by method FindSampleModes to the
        console.
        
        Inputs: s = current sample identifier as a Python string (EX. 'S1)
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

### CONSTRUCTION CONTINUES BELOW...




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




    def PLOTDualSampleWeightPercents(self, mode='print'):
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
        plt.style.use('ggplot')
        bins = np.arange(len(self.screens))
        w=0.75
        xTickLabels=list(map(str, self.screens))   # convert list of numeric screen sizes 
                                                   # to list of strings for labels
        xTickLabels[-1]='fines'

        # write header to console:
        #print('-'*90)
        #print('  Transect: '+self.transect+'   Sample: '+s)
        print('')

        fig1=plt.figure(figsize=(21,7))

        # weight percentages subplot
        ax1=fig1.add_subplot(1,2,1)
        ax1.bar(bins,self.sdf['weight_percent'], align='center', width=w)
        ax1.plot(bins,self.sdf['weight_percent'],'g--',linewidth=3.5)
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
        ax2.plot(bins,self.sdf['cum_weight_percent'],'g-', linewidth=2.0)
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

        fig_name=str(self.sample_id)+'_'+'_dual.png'
        if(mode=='print'):
            plt.show()
        if(mode=='save'):
            plt.savefig(fig_name)
            plt.close()

        print('-'*90)

# ## ####################################################################################
# ## End of SedSASample class listing 
# ## ####################################################################################





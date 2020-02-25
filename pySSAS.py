
import numpy as np
import pandas as pd


class Data_Loaders:
    '''
    Data_Loader: class to handle the loading of source data for analysis
    '''
    def __init__(self):
        self.df = 'No dataframe defined'
        
    def dataframe_from_csv(self, src_file):
        '''Expects a csv text file from the caller with first line a column
            header; the remaining lines contain data in specified order (consult
            the README for details.
        '''
        try:
            self.df = pd.read_csv(src_file)
            # call the builder class function to 
            return(self.df)
        except AttributeError as e:
            print('Oops! Something seems to have gone wrong in loading the source data:', e)
            
            
            
    def Return_Dataframe(self):
        return(self.df)
    
    
class Data_Constructors:
    '''
    This could be an internal class that's called by the data_loaders to take
    the raw data as provided by the user and add to the raw weights: the weight
    percentages, the cumulative weight percentages, and 
    '''
    def __init__(self, df, id, screens, extrap_method='linear'):
        self.df = df
        self.id = id
        self.screens = np.array(screens).T
        
        # list of quantiles to be estimated (interpolated) as required by Folk
        # and Ward graphical statistics...
        self.quantilesList=[5,10,16,25,50,75,84,90,95]

# *******
    def BuildWorkingDataFrame(self):
        ''' from the user-provide dataframe (as an argument to this function
        and WHICH CONTAINS ONLY THE RAW phi WEIGHTS--use indexing to send 
        only this portion of the original dataframe) and, using the raw 
        weights, computes the individual and cumulative weight percentages, 
        and the phi-class midpoint, adding each of these three to a new
        sample dataframe. This dataframe is returned to the caller for use
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
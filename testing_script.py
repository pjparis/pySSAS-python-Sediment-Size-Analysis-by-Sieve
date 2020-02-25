from pySSAS import Data_Loaders
from pySSAS import Data_Constructors

data = Data_Loaders()
df = data.dataframe_from_csv('test_data.csv')
#### Interpolating the quantile values in the SedSASample class

To replicate the old-fashioned graphic methods for estimating sediment grain size quantiles for analysis, as per Folk and Ward (1957) we must invoke a double-dose of 2-D, simple linear interpolation. Even though the absissa is in sediment particle size phi units which are based 2 logarithms we treat the x-axis scaling as a linear one. The ordinate on the other hand, representing the sieved fraction cumulative weight percent, is linear throughout. 

-the data is stored in a Python Pandas dataframe

to estimate the $k^{th}$ quantile:

1. search the sample column and extract all records (rows) whose cumulative weights are less than that of the quantile value (D) you are searching: *p=cum_phi.query('S2 <'+str(d) )*

2. search the sample column and extract all records (rows) whose cumulative weights are greater than that of the quantile value (D) you are searching: *q=cum_phi.query('S2>'+str(d) )*

3. for the group of query results (records) that are less than D extract and store the last record (slice position: -1). this would be the record whose cumulative weight percentage is closest to, but still less than, D
4. for the group of query results (records) that are greater than D extract and store the first record (slice position 1). This would be the record whose cumulative weight percentage is closest to, but still greater than, D
- The computation then proceeds as a proportionality equation of the form:


$$( \frac{D_k-CumWt_{lower bnd}}{CumWt_{upper bnd}-CumWt_{lower bnd}} )$$
<br><br>
the resulting quotient is then multiplied by the difference between the bounding Phis:

$$ \Phi_{upper bnd} - \Phi_{lower bnd} $$

finally, the resulting product is added to the lower bounding Phi. 
<br><br>
The entire equation looks like this:


$$ D\Phi=(( \frac{D_k-CumWt_{lower bnd}}{CumWt_{upper bnd}-CumWt_{lower bnd}} )*(\Phi_{upper bnd} - \Phi_{lower bnd}))+\Phi_{lower bnd} $$
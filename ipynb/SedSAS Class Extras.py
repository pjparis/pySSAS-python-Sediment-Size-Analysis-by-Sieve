	def ClassifyGeometricSorting(self, SD):
		''' Classifies the resultant sorting (standard deviation) value as computed using 
		either GEOMETRIC Method from Method of Moments or Folk and Ward approaches
		
		Called by STATISTICAL METHOD 2 (FOLK AND WARD GEOMETRIC) and METHOD 4 (METHOD
		OF MOMENT GEOMETRIC)
		'''
		try:
			if( SD < 1.27):
				return('Very well sorted')
			if( SD >= 1.27 and SD < 1.41):
				return('Well sorted')
			if( SD >= 1.41 and SD < 1.64):
				return('Moderately well sorted')
			if( SD >= 1.64 and SD < 2.00):
				return('Moderately sorted')
			if( SD >= 2.00 and SD < 4.00):
				return('Poorly sorted')
			if( SD >= 4.00 and SD < 16.00):
				return('Very poorly sorted')
			if( SD > 16.00):
				return('Extremely poorly sorted')
		except (NameError, AttributeError, KeyError) as e:
			print('Oops! An error occured.', sys.exc_info()[0], e)
			
			
		def ClassifyByGrainSizeMM(self, MN):
		'''
		'''
		try:
			if( MN >= 256.0):
				return('Boulders')
			if( MN >= 64.0 and MN < 256.0):
				return('Cobbles')
			if( MN >= 4.0 and MN < 64.0):
				return('Pebbles')
			if( MN >= 2.0 and MN < 4.0):
				return('Gravel Granules')
			if( MN >= 1.0 and MN < 2.0):
				return('Very Coarse Sand')
			if( MN >= 0.5 and MN < 1.0):
				return('Coarse Sand')
			if( MN >= 0.25 and MN < 0.5):
				return('Medium Sand')
			if( MN >= 0.125 and MN < 0.25):
				return('Fine Sand')
			if( MN >= 0.0625 and MN < 0.125):
				return('Very Fine Sand')
			if( MN >= 0.03125 and MN < 0.0625):
				return('Coarse Silt')
			if( MN >= 0.015625 and MN < 0.03125):
				return('Medium Silt')
			if( MN >= 0.0078125 and MN < 0.015625):
				return('Fine Silt')
			if( MN >= 0.00390625 and MN < 0.0078125):
				return('Very Fine Silt')
			if( MN < 0.00390625):
				return('Clay')
		except (NameError, AttributeError, KeyError) as e:
			print('Oops! An error occured.', sys.exc_info()[0], e)



	def ClassifyFWGeoSkewness(self, SK):
		''' Classifies a resultant skewness value as computed using the
		Folk and Ward (1957) graphic Geometric method of computation. 
		
		Called by STATISTICAL METHOD 2 (FOLK AND WARD GEOMETRIC)
		'''
		if( SK < -0.3 and SK >= -1.00):
			return('Very fine skewed')
		if( SK < -0.1 and SK >= -0.3):
			return('Fine skewed')
		if( SK >= -0.1 and SK <= 0.1):
			return('Symmetrical')
		if( SK > 0.1 and SK <= 0.3):
			return('Coarse skewed')
		if( SK > 0.3 and SK <= 1.00):
			return('Very coarse skewed')




	def ClassifyMoMGeoSkewness(self, SK):
		''' Classifies a resultant skewness value as computed using the
		Method of Moments Geometric computation. 
		     
		Called by STATISTICAL METHOD 4 (METHOD OF MOMENTS GEOMETRIC)
		'''
		if( SK < -1.30):
			return('Very fine skewed')
		if( SK > -1.30 and SK <= -0.43):
			return('Fine skewed')
		if( SK > -0.43 and SK <= 0.43):
			return('Symmetrical')
		if( SK > 0.43 and SK <= 1.30):
			return('Coarse skewed')			
		if( SK > 1.30):
			return('Very coarse skewed')



		momgeo=self.ComputeGeometricMethodofMomentsStats(return_description=True)


		       'MoMGeoMean':momgeo[0][0],
		       'MoMGeoSort':momgeo[0][1],
		       'MoMGeoSkew':momgeo[0][2],
		       'MoMGeoKurt':momgeo[0][3],
		       'MoMGeoSizeClass':momgeo[1][0],
		       'MoMGeoSortCLass':momgeo[1][1],
		       'MoMGeoSkewClass':momgeo[1][2],
		       'MoMGeoKurtClass':momgeo[1][3],
		       
		       
		       

# 	def Stats2CSV(self, S, fn):
# 		'''Writes statistics results to an ASCII comma separated text file. Multiple calls 
# 		will append subsequent results to existing file. If file does not exist it is
# 		created anew.
# 		S=computed statistics in a Python list or tuple
# 		fn=file path and name into which the 
# 	'''
# 		Qp=self.Q[0]
# 		Sp=S[0]
# 		
# 		# create a file header string for new file:
# 		hdr='id'+','+'Mean'+','+'Median'+','+'Sorting'+','+'Skewness'+','+'Kurtosis'
# 		hdr=hdr+','+'D25'+','+'D75'+'\n'
# 	
# 		if( os.path.exists(fn)):
# 			fp=open(fn, 'a')
# 		else:
# 			fp=open(fn, 'w')
# 			print(hdr)
# 			fp.write(hdr)
# 
# 		write_string=str(self.id)+','+str(Sp[0])+','+str(Qp[4])+','+str(Sp[1])
# 		write_string=write_string+','+str(Sp[2])+','+str(Sp[3])+','+str(Qp[3])
# 		write_string=write_string+','+str(Qp[5])+'\n'
# 		
# 		fp.write(write_string)
# 		
# 		fp.close()

		return()

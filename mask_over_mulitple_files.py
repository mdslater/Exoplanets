  import pyfits
  import numpy as np
  import os
  import sys
  import scipy.ndimage as snd
  import csv
  import astropy.table
  import astropy.units as u
  import math
  import scipy.integrate as integrate
  import scipy.special as special
  from scipy.integrate import quad, dblquad
  from scipy import integrate
  from astropy.time import Time

  #define function that finds radius of star
  def dist(x, y, xx, yy):
	x = x * 1.
	y = y * 1.
	xx = xx * 1.
	yy = yy * 1.

	return np.sqrt((yy-y)**2+(xx-x)**2) 

  with open('cor_Data.txt') as f:
	fitfiles=f.read().split()
		
		
  #set a variable equal to your file#

  data_path=[]
  for fitfile in fitfiles:
	with open(fitfile) as f:
		data_path=fitfile
		hdulist=pyfits.open(data_path) * 1
		dataset=hdulist[0].data * 1
		header=hdulist[0].header
		#create mask#
		mask=dataset * 1
		
		#make mask field all 0s#
		mask=np.multiply(mask,0.0)
		
		thrs=(10 * np.std(dataset))
		
		#define x and y values#
		y_size=len(dataset)
		x_size=len(dataset[0])
		
		#for loop to define each pixel in the mask as either 1 or 0, depending on whether its 			luminosity#
		#is higher or lower than the threshold#
  
		for y in range(y_size):
		 	for x in range(x_size):
		 		if dataset[y,x] > thrs:
		 			mask[y,x]=1
		 
		
		 
		 
		#labels where each pixel above the threshold#
		labels, num=snd.label(mask, np.ones((3,3)))
		
		#defines the center of mass for each star#
		centers = snd.center_of_mass(mask,labels,range(1,num+1))
		
		#defines the x and y values of each center of mass#
		x = np.array(centers)[:,0]
		y = np.array(centers)[:,1]
		
		#convert x and y lists into integer lists
		x=x.astype(int)
		y=y.astype(int)	
		
		#turn x and y into lists from arrays
		x=x.tolist()
		y=y.tolist()
	
		#drop the dead pixels: to-do -> x,y array
		 
		def drop_dead(matrix,x1,y1):			
 			for i in range(len(x1)):		
  				if matrix[x1[i],y1[i]-1]==0 and matrix[x1[i],y1[i]+1]==0 and matrix[x1[i]+1,y1[i]]==0 and matrix[x1[i]-1,y1[i]]==0:	
  					matrix[x1[i],y1[i]] = 0
  			
		drop_dead(mask,x,y)
  
  		
  
		#define radius of stars
		radius=[]
		for i in range(len(x)):
 			x_s=x[i]
 			y_s=y[i]
		 	list_x=[x_s]
		 	list_y=[y_s]
		 	c_list_x=[]
		 	c_list_y=[]
		 	n_list_x=[]
		 	n_list_y=[]
		 	dist_array=[]
		 	for j in range(len(list_x)):
		 		if list_x[j] not in c_list_x:
		 			n_list_x.append(list_x[j])
		 			n_list_y.append(list_y[j])
		 	condition=True
		 	while condition:
		 		for k in range(len(n_list_x)):
		 			dist_array.append(dist(n_list_x[k],n_list_y[k],x_s,y_s))	
		 			if mask[n_list_x[k]+0,n_list_y[k]+1]==1:
		 				list_x.append(n_list_x[k]+0)
		 				list_y.append(n_list_y[k]+1)
		 			if mask[n_list_x[k]+0,n_list_y[k]-1]==1:
		 				list_x.append(n_list_x[k]+0)
		 				list_y.append(n_list_y[k]-1)
		 			if mask[n_list_x[k]+1,n_list_y[k]+0]==1:
		 				list_x.append(n_list_x[k]+1)
		 				list_y.append(n_list_y[k]+0)
		 			if mask[n_list_x[k]-1,n_list_y[k]+0]==1:
		 				list_x.append(n_list_x[k]-1)
		 				list_y.append(n_list_y[k]+0)
		 						
		 			c_list_x.append(n_list_x[k])
		 			c_list_y.append(n_list_y[k])
		 			
		 			
		 		n_list_x=[]
		 		n_list_y=[]
		 	
		 		for i in range(len(list_x)):
		 			if list_x[i] not in c_list_x:
		 				n_list_x.append(list_x[i])
 						n_list_y.append(list_y[i])
 				if len(n_list_x)==0:
 					
 					radius.append(max(dist_array))
 					condition = False
 		
		#delete x,y,r for r=0
		new_x=[]
		new_y=[]
		new_radius=[]
		for i in range(len(x)):
			if radius[i] != 0:
				new_x.append(x[i])
		for i in range(len(y)):
			if radius[i] != 0:
				new_y.append(y[i])

		for i in range(len(radius)):
			if radius[i] != 0:
				new_radius.append(radius[i])
		x=new_x
		y=new_y
		radius=new_radius
		coordinates_and_radii=zip(x,y,radius)	  
		
		
		#calculate for fluxes of each star
		total_flux=[]
		total=0.
		numbers1=[0,-1,-2,-3,-4,-5,-6,-7,1,2,3,4,5,6,7]
		numbers2=[0,-1,-2,-3,-4,-5,-6,-7,1,2,3,4,5,6,7]
		
		for i in range(len(radius)):		
		  	if radius[i] > 0.:	
		  		for n1 in numbers1:
					for n2 in numbers2:
						if dist(x[i],y[i],x[i]+n1,y[i]+n2) <= radius[i]:
							total+=dataset[x[i]+n1,y[i]+n2]
		
			total_flux.append(total)
			total=0.
		
		total_flux=[j for j in total_flux if j != 0]

		#combine all data
		coord_rad_flux=zip(x,y,radius,total_flux)
		#get time and date for each file
		time=[]
		date=[]
		time.append(header["TIME-OBS"])
		date.append(header["DATE-OBS"])
		for i in range(len(time)):		
			time.insert(0,date[0])

				#delete x,y,r for r=0
		new_x=[]
		new_y=[]
		new_radius=[]
		for i in range(len(x)):
			if radius[i] != 0:
				new_x.append(x[i])
		for i in range(len(y)):
			if radius[i] != 0:
				new_y.append(y[i])

		for i in range(len(radius)):
			if radius[i] != 0:
				new_radius.append(radius[i])
		x=new_x
		y=new_y
		radius=new_radius
		coordinates_and_radii=zip(x,y,radius)	  
		
		
		#calculate for fluxes of each star
		total_flux=[]
		total=0.
		numbers1=[0,-1,-2,-3,-4,-5,-6,-7,1,2,3,4,5,6,7]
		numbers2=[0,-1,-2,-3,-4,-5,-6,-7,1,2,3,4,5,6,7]
		
		for i in range(len(radius)):		
		  	if radius[i] > 0.:	
		  		for n1 in numbers1:
					for n2 in numbers2:
						if dist(x[i],y[i],x[i]+n1,y[i]+n2) <= radius[i]:
							total+=dataset[x[i]+n1,y[i]+n2]
		
			total_flux.append(total)
			total=0.
		
		total_flux=[j for j in total_flux if j != 0]

		#combine all data
		coord_rad_flux=zip(x,y,radius,total_flux)
		#get time and date for each file
		time=[]
		date=[]
		time.append(header["TIME-OBS"])
		date.append(header["DATE-OBS"])
		for i in range(len(time)):		
			time.insert(0,date[0])
		

		#take out extra stars
		x_1star=[]
		y_1star=[]
		radius_1star=[]
		fluxes_1star=[]
		for i in range(len(x)):
			if x[i] == 3858 or x[i] == 3858 + 1 or x[i] == 3858 + 2 or x[i] == 3858 + 3 or x[i] == 3858 + 4 or x[i] == 3858 + 5 or x[i] == 3858 - 1 or x[i] == 3858 - 2 or x[i] == 3858 - 3 or x[i] == 3858 - 4 or x[i] == 3858 - 5:
				if y[i] == 1120 or y[i] == 1120 + 1 or y[i] == 1120 + 2 or y[i] == 1120 + 3 or y[i] == 1120 + 4 or y[i] == 1120 + 5 or y[i] == 1120 - 1 or y[i] == 1120 - 2 or y[i] == 1120 - 3 or y[i] == 1120 - 4 or y[i] == 1120 - 5: 
					x_1star.insert(0,x[i])
					fluxes_1star.insert(0,total_flux[i])
					y_1star.insert(0,y[i])
					radius_1star.insert(0,radius[i])
			if x[i] == 3774 or x[i] == 3774 + 1 or x[i] == 3774 + 2 or x[i] == 3774 + 3 or x[i] == 3774 + 4 or x[i] == 3774 + 5 or x[i] == 3774 - 1 or x[i] == 3774 - 2 or x[i] == 3774 - 3 or x[i] == 3774 - 4 or x[i] == 3774 - 5:
				if y[i] == 387 or y[i] == 387 + 1 or y[i] == 387 + 2 or y[i] == 387 + 3 or y[i] == 387 + 4 or y[i] == 387 + 5 or y[i] == 387 - 1 or y[i] == 387 - 2 or y[i] == 387 - 3 or y[i] == 387 - 4 or y[i] == 387 - 5: 
					x_1star.insert(0,x[i])
					fluxes_1star.insert(0,total_flux[i])
					y_1star.insert(0,y[i])
					radius_1star.insert(0,radius[i])	
			if x[i] == 3676 or x[i] == 3676 + 1 or x[i] == 3676 + 2 or x[i] == 3676 + 3 or x[i] == 3676 + 4 or x[i] == 3676 + 5 or x[i] == 3676 - 1 or x[i] == 3676 - 2 or x[i] == 3676 - 3 or x[i] == 3676 - 4 or x[i] == 3676 - 5:
				if y[i] == 3080 or y[i] == 3080 + 1 or y[i] == 3080 + 2 or y[i] == 3080 + 3 or y[i] == 3080 + 4 or y[i] == 3080 + 5 or y[i] == 3080 - 1 or y[i] == 3080 - 2 or y[i] == 3080 - 3 or y[i] == 3080 - 4 or y[i] == 3080 - 5: 
					x_1star.insert(0,x[i])
					fluxes_1star.insert(0,total_flux[i])
					y_1star.insert(0,y[i])
					radius_1star.insert(0,radius[i])
			if x[i] == 3487 or x[i] == 3487 + 1 or x[i] == 3487 + 2 or x[i] == 3487 + 3 or x[i] == 3487 + 4 or x[i] == 3487 + 5 or x[i] == 3487 - 1 or x[i] == 3487 - 2 or x[i] == 3487 - 3 or x[i] == 3487 - 4 or x[i] == 3487 - 5:
				if y[i] == 2118 or y[i] == 2118 + 1 or y[i] == 2118 + 2 or y[i] == 2118 + 3 or y[i] == 2118 + 4 or y[i] == 2118 + 5 or y[i] == 2118 - 1 or y[i] == 2118 - 2 or y[i] == 2118 - 3 or y[i] == 2118 - 4 or y[i] == 2118 - 5: 
					x_1star.insert(0,x[i])
					fluxes_1star.insert(0,total_flux[i])
					y_1star.insert(0,y[i])
					radius_1star.insert(0,radius[i])
			if x[i] == 3482 or x[i] == 3482 + 1 or x[i] == 3482 + 2 or x[i] == 3482 + 3 or x[i] == 3482 + 4 or x[i] == 3482 + 5 or x[i] == 3482 - 1 or x[i] == 3482 - 2 or x[i] == 3482 - 3 or x[i] == 3482 - 4 or x[i] == 3482 - 5:
				if y[i] == 1894 or y[i] == 1894 + 1 or y[i] == 1894 + 2 or y[i] == 1894 + 3 or y[i] == 1894 + 4 or y[i] == 1894 + 5 or y[i] == 1894 - 1 or y[i] == 1894 - 2 or y[i] == 1894 - 3 or y[i] == 1894 - 4 or y[i] == 1894 - 5: 
					x_1star.insert(0,x[i])
					fluxes_1star.insert(0,total_flux[i])
					y_1star.insert(0,y[i])
					radius_1star.insert(0,radius[i])
			if x[i] == 3127 or x[i] == 3127 + 1 or x[i] == 3127 + 2 or x[i] == 3127 + 3 or x[i] == 3127 + 4 or x[i] == 3127 + 5 or x[i] == 3127 - 1 or x[i] == 3127 - 2 or x[i] == 3127 - 3 or x[i] == 3127 - 4 or x[i] == 3127 - 5:
				if y[i] == 1874 or y[i] == 1874 + 1 or y[i] == 1874 + 2 or y[i] == 1874 + 3 or y[i] == 1874 + 4 or y[i] == 1874 + 5 or y[i] == 1874 - 1 or y[i] == 1874 - 2 or y[i] == 1874 - 3 or y[i] == 1874 - 4 or y[i] == 1874 - 5: 
					x_1star.insert(0,x[i])
					fluxes_1star.insert(0,total_flux[i])
					y_1star.insert(0,y[i])
					radius_1star.insert(0,radius[i])
			if x[i] == 3102 or x[i] == 3102 + 1 or x[i] == 3102 + 2 or x[i] == 3102 + 3 or x[i] == 3102 + 4 or x[i] == 3102 + 5 or x[i] == 3102 - 1 or x[i] == 3102 - 2 or x[i] == 3102 - 3 or x[i] == 3102 - 4 or x[i] == 3102 - 5:
				if y[i] == 367 or y[i] == 367 + 1 or y[i] == 367 + 2 or y[i] == 367 + 3 or y[i] == 367 + 4 or y[i] == 367 + 5 or y[i] == 367 - 1 or y[i] == 367 - 2 or y[i] == 367 - 3 or y[i] == 367 - 4 or y[i] == 367 - 5: 
					x_1star.insert(0,x[i])
					fluxes_1star.insert(0,total_flux[i])
					y_1star.insert(0,y[i])
					radius_1star.insert(0,radius[i])
			if x[i] == 2993 or x[i] == 2993 + 1 or x[i] == 2993 + 2 or x[i] == 2993 + 3 or x[i] == 2993 + 4 or x[i] == 2993 + 5 or x[i] == 2993 - 1 or x[i] == 2993 - 2 or x[i] == 2993 - 3 or x[i] == 2993 - 4 or x[i] == 2993 - 5:
				if y[i] == 591 or y[i] == 591 + 1 or y[i] == 591 + 2 or y[i] == 591 + 3 or y[i] == 591 + 4 or y[i] == 591 + 5 or y[i] == 591 - 1 or y[i] == 591 - 2 or y[i] == 591 - 3 or y[i] == 591 - 4 or y[i] == 591 - 5: 
					x_1star.insert(0,x[i])
					fluxes_1star.insert(0,total_flux[i])
					y_1star.insert(0,y[i])
					radius_1star.insert(0,radius[i])
			if x[i] == 2798 or x[i] == 2798 + 1 or x[i] == 2798 + 2 or x[i] == 2798 + 3 or x[i] == 2798 + 4 or x[i] == 2798 + 5 or x[i] == 2798 - 1 or x[i] == 2798 - 2 or x[i] == 2798 - 3 or x[i] == 2798 - 4 or x[i] == 2798 - 5:
				if y[i] == 972 or y[i] == 972 + 1 or y[i] == 972 + 2 or y[i] == 972 + 3 or y[i] == 972 + 4 or y[i] == 972 + 5 or y[i] == 972 - 1 or y[i] == 972 - 2 or y[i] == 972 - 3 or y[i] == 972 - 4 or y[i] == 972 - 5: 
					x_1star.insert(0,x[i])
					fluxes_1star.insert(0,total_flux[i])
					y_1star.insert(0,y[i])
					radius_1star.insert(0,radius[i])
			if x[i] == 2237 or x[i] == 2237 + 1 or x[i] == 2237 + 2 or x[i] == 2237 + 3 or x[i] == 2237 + 4 or x[i] == 2237 + 5 or x[i] == 2237 - 1 or x[i] == 2237 - 2 or x[i] == 2237 - 3 or x[i] == 2237 - 4 or x[i] == 2237 - 5:
				if y[i] == 658 or y[i] == 658 + 1 or y[i] == 658 + 2 or y[i] == 658 + 3 or y[i] == 658 + 4 or y[i] == 658 + 5 or y[i] == 658 - 1 or y[i] == 658 - 2 or y[i] == 658 - 3 or y[i] == 658 - 4 or y[i] == 658 - 5: 
					x_1star.insert(0,x[i])
					fluxes_1star.insert(0,total_flux[i])
					y_1star.insert(0,y[i])
					radius_1star.insert(0,radius[i])
			if x[i] == 2225 or x[i] == 2225 + 1 or x[i] == 2225 + 2 or x[i] == 2225 + 3 or x[i] == 2225 + 4 or x[i] == 2225 + 5 or x[i] == 2225 - 1 or x[i] == 2225 - 2 or x[i] == 2225 - 3 or x[i] == 2225 - 4 or x[i] == 2225 - 5:
				if y[i] == 1300 or y[i] == 1300 + 1 or y[i] == 1300 + 2 or y[i] == 1300 + 3 or y[i] == 1300 + 4 or y[i] == 1300 + 5 or y[i] == 1300 - 1 or y[i] == 1300 - 2 or y[i] == 1300 - 3 or y[i] == 1300 - 4 or y[i] == 1300 - 5: 
					x_1star.insert(0,x[i])
					fluxes_1star.insert(0,total_flux[i])
					y_1star.insert(0,y[i])
					radius_1star.insert(0,radius[i])
			if x[i] == 2190 or x[i] == 2190 + 1 or x[i] == 2190 + 2 or x[i] == 2190 + 3 or x[i] == 2190 + 4 or x[i] == 2190 + 5 or x[i] == 2190 - 1 or x[i] == 2190 - 2 or x[i] == 2190 - 3 or x[i] == 2190 - 4 or x[i] == 2190 - 5:
				if y[i] == 1847 or y[i] == 1847 + 1 or y[i] == 1847 + 2 or y[i] == 1847 + 3 or y[i] == 1847 + 4 or y[i] == 1847 + 5 or y[i] == 1847 - 1 or y[i] == 1847 - 2 or y[i] == 1847 - 3 or y[i] == 1847 - 4 or y[i] == 1847 - 5: 
					x_1star.insert(0,x[i])
					fluxes_1star.insert(0,total_flux[i])
					y_1star.insert(0,y[i])
					radius_1star.insert(0,radius[i])
			if x[i] == 2061 or x[i] == 2061 + 1 or x[i] == 2061 + 2 or x[i] == 2061 + 3 or x[i] == 2061 + 4 or x[i] == 2061 + 5 or x[i] == 2061 - 1 or x[i] == 2061 - 2 or x[i] == 2061 - 3 or x[i] == 2061 - 4 or x[i] == 2061 - 5:
				if y[i] == 1388 or y[i] == 1388 + 1 or y[i] == 1388 + 2 or y[i] == 1388 + 3 or y[i] == 1388 + 4 or y[i] == 1388 + 5 or y[i] == 1388 - 1 or y[i] == 1388 - 2 or y[i] == 1388 - 3 or y[i] == 1388 - 4 or y[i] == 1388 - 5: 
					x_1star.insert(0,x[i])
					fluxes_1star.insert(0,total_flux[i])
					y_1star.insert(0,y[i])
					radius_1star.insert(0,radius[i])
			if x[i] == 1950 or x[i] == 1950 + 1 or x[i] == 1950 + 2 or x[i] == 1950 + 3 or x[i] == 1950 + 4 or x[i] == 1950 + 5 or x[i] == 1950 - 1 or x[i] == 1950 - 2 or x[i] == 1950 - 3 or x[i] == 1950 - 4 or x[i] == 1950 - 5:
				if y[i] == 519 or y[i] == 519 + 1 or y[i] == 519 + 2 or y[i] == 519 + 3 or y[i] == 519 + 4 or y[i] == 519 + 5 or y[i] == 519 - 1 or y[i] == 519 - 2 or y[i] == 519 - 3 or y[i] == 519 - 4 or y[i] == 519 - 5: 
					x_1star.insert(0,x[i])
					fluxes_1star.insert(0,total_flux[i])
					y_1star.insert(0,y[i])
					radius_1star.insert(0,radius[i])
			if x[i] == 1920 or x[i] == 1920 + 1 or x[i] == 1920 + 2 or x[i] == 1920 + 3 or x[i] == 1920 + 4 or x[i] == 1920 + 5 or x[i] == 1920 - 1 or x[i] == 1920 - 2 or x[i] == 1920 - 3 or x[i] == 1920 - 4 or x[i] == 1920 - 5:
				if y[i] == 3668 or y[i] == 3668 + 1 or y[i] == 3668 + 2 or y[i] == 3668 + 3 or y[i] == 3668 + 4 or y[i] == 3668 + 5 or y[i] == 3668 - 1 or y[i] == 3668 - 2 or y[i] == 3668 - 3 or y[i] == 3668 - 4 or y[i] == 3668 - 5: 
					x_1star.insert(0,x[i])
					fluxes_1star.insert(0,total_flux[i])
					y_1star.insert(0,y[i])
					radius_1star.insert(0,radius[i])
			if x[i] == 1834 or x[i] == 1834 + 1 or x[i] == 1834 + 2 or x[i] == 1834 + 3 or x[i] == 1834 + 4 or x[i] == 1834 + 5 or x[i] == 1834 - 1 or x[i] == 1834 - 2 or x[i] == 1834 - 3 or x[i] == 1834 - 4 or x[i] == 1834 - 5:
				if y[i] == 1514 or y[i] == 1514 + 1 or y[i] == 1514 + 2 or y[i] == 1514 + 3 or y[i] == 1514 + 4 or y[i] == 1514 + 5 or y[i] == 1514 - 1 or y[i] == 1514 - 2 or y[i] == 1514 - 3 or y[i] == 1514 - 4 or y[i] == 1514 - 5: 
					x_1star.insert(0,x[i])
					fluxes_1star.insert(0,total_flux[i])
					y_1star.insert(0,y[i])
					radius_1star.insert(0,radius[i])
			if x[i] == 1780 or x[i] == 1780 + 1 or x[i] == 1780 + 2 or x[i] == 1780 + 3 or x[i] == 1780 + 4 or x[i] == 1780 + 5 or x[i] == 1780 - 1 or x[i] == 1780 - 2 or x[i] == 1780 - 3 or x[i] == 1780 - 4 or x[i] == 1780 - 5:
				if y[i] == 3167 or y[i] == 3167 + 1 or y[i] == 3167 + 2 or y[i] == 3167 + 3 or y[i] == 3167 + 4 or y[i] == 3167 + 5 or y[i] == 3167 - 1 or y[i] == 3167 - 2 or y[i] == 3167 - 3 or y[i] == 3167 - 4 or y[i] == 3167 - 5: 
					x_1star.insert(0,x[i])
					fluxes_1star.insert(0,total_flux[i])
					y_1star.insert(0,y[i])
					radius_1star.insert(0,radius[i])
			if x[i] == 1739 or x[i] == 1739 + 1 or x[i] == 1739 + 2 or x[i] == 1739 + 3 or x[i] == 1739 + 4 or x[i] == 1739 + 5 or x[i] == 1739 - 1 or x[i] == 1739 - 2 or x[i] == 1739 - 3 or x[i] == 1739 - 4 or x[i] == 1739 - 5:
				if y[i] == 3328 or y[i] == 3328 + 1 or y[i] == 3328 + 2 or y[i] == 3328 + 3 or y[i] == 3328 + 4 or y[i] == 3328 + 5 or y[i] == 3328 - 1 or y[i] == 3328 - 2 or y[i] == 3328 - 3 or y[i] == 3328 - 4 or y[i] == 3328 - 5: 
					x_1star.insert(0,x[i])
					fluxes_1star.insert(0,total_flux[i])
					y_1star.insert(0,y[i])
					radius_1star.insert(0,radius[i])
			if x[i] == 1715 or x[i] == 1715 + 1 or x[i] == 1715 + 2 or x[i] == 1715 + 3 or x[i] == 1715 + 4 or x[i] == 1715 + 5 or x[i] == 1715 - 1 or x[i] == 1715 - 2 or x[i] == 1715 - 3 or x[i] == 1715 - 4 or x[i] == 1715 - 5:
				if y[i] == 2618 or y[i] == 2618 + 1 or y[i] == 2618 + 2 or y[i] == 2618 + 3 or y[i] == 2618 + 4 or y[i] == 2618 + 5 or y[i] == 2618 - 1 or y[i] == 2618 - 2 or y[i] == 2618 - 3 or y[i] == 2618 - 4 or y[i] == 2618 - 5: 
					x_1star.insert(0,x[i])
					fluxes_1star.insert(0,total_flux[i])
					y_1star.insert(0,y[i])
					radius_1star.insert(0,radius[i])
			if x[i] == 1643 or x[i] == 1643 + 1 or x[i] == 1643 + 2 or x[i] == 1643 + 3 or x[i] == 1643 + 4 or x[i] == 1643 + 5 or x[i] == 1643 - 1 or x[i] == 1643 - 2 or x[i] == 1643 - 3 or x[i] == 1643 - 4 or x[i] == 1643 - 5:
				if y[i] == 747 or y[i] == 747 + 1 or y[i] == 747 + 2 or y[i] == 747 + 3 or y[i] == 747 + 4 or y[i] == 747 + 5 or y[i] == 747 - 1 or y[i] == 747 - 2 or y[i] == 747 - 3 or y[i] == 747 - 4 or y[i] == 747 - 5: 
					x_1star.insert(0,x[i])
					fluxes_1star.insert(0,total_flux[i])
					y_1star.insert(0,y[i])
					radius_1star.insert(0,radius[i])
			if x[i] == 1445 or x[i] == 1445 + 1 or x[i] == 1445 + 2 or x[i] == 1445 + 3 or x[i] == 1445 + 4 or x[i] == 1445 + 5 or x[i] == 1445 - 1 or x[i] == 1445 - 2 or x[i] == 1445 - 3 or x[i] == 1445 - 4 or x[i] == 1445 - 5:
				if y[i] == 2625 or y[i] == 2625 + 1 or y[i] == 2625 + 2 or y[i] == 2625 + 3 or y[i] == 2625 + 4 or y[i] == 2625 + 5 or y[i] == 2625 - 1 or y[i] == 2625 - 2 or y[i] == 2625 - 3 or y[i] == 2625 - 4 or y[i] == 2625 - 5: 
					x_1star.insert(0,x[i])
					fluxes_1star.insert(0,total_flux[i])
					y_1star.insert(0,y[i])
					radius_1star.insert(0,radius[i])
			if x[i] == 1309 or x[i] == 1309 + 1 or x[i] == 1309 + 2 or x[i] == 1309 + 3 or x[i] == 1309 + 4 or x[i] == 1309 + 5 or x[i] == 1309 - 1 or x[i] == 1309 - 2 or x[i] == 1309 - 3 or x[i] == 1309 - 4 or x[i] == 1309 - 5:
				if y[i] == 643 or y[i] == 643 + 1 or y[i] == 643 + 2 or y[i] == 643 + 3 or y[i] == 643 + 4 or y[i] == 643 + 5 or y[i] == 643 - 1 or y[i] == 643 - 2 or y[i] == 643 - 3 or y[i] == 643 - 4 or y[i] == 643 - 5: 
					x_1star.insert(0,x[i])
					fluxes_1star.insert(0,total_flux[i])
					y_1star.insert(0,y[i])
					radius_1star.insert(0,radius[i])
			if x[i] == 1307 or x[i] == 1307 + 1 or x[i] == 1307 + 2 or x[i] == 1307 + 3 or x[i] == 1307 + 4 or x[i] == 1307 + 5 or x[i] == 1307 - 1 or x[i] == 1307 - 2 or x[i] == 1307 - 3 or x[i] == 1307 - 4 or x[i] == 1307 - 5:
				if y[i] == 3239 or y[i] == 3239 + 1 or y[i] == 3239 + 2 or y[i] == 3239 + 3 or y[i] == 3239 + 4 or y[i] == 3239 + 5 or y[i] == 3239 - 1 or y[i] == 3239 - 2 or y[i] == 3239 - 3 or y[i] == 3239 - 4 or y[i] == 3239 - 5: 
					x_1star.insert(0,x[i])
					fluxes_1star.insert(0,total_flux[i])
					y_1star.insert(0,y[i])
					radius_1star.insert(0,radius[i])
			if x[i] == 1004 or x[i] == 1004 + 1 or x[i] == 1004 + 2 or x[i] == 1004 + 3 or x[i] == 1004 + 4 or x[i] == 1004 + 5 or x[i] == 1004 - 1 or x[i] == 1004 - 2 or x[i] == 1004 - 3 or x[i] == 1004 - 4 or x[i] == 1004 - 5:
				if y[i] == 617 or y[i] == 617 + 1 or y[i] == 617 + 2 or y[i] == 617 + 3 or y[i] == 617 + 4 or y[i] == 617 + 5 or y[i] == 617 - 1 or y[i] == 617 - 2 or y[i] == 617 - 3 or y[i] == 617 - 4 or y[i] == 617 - 5: 
					x_1star.insert(0,x[i])
					fluxes_1star.insert(0,total_flux[i])
					y_1star.insert(0,y[i])
					radius_1star.insert(0,radius[i])
			if x[i] == 717 or x[i] == 717 + 1 or x[i] == 717 + 2 or x[i] == 717 + 3 or x[i] == 717 + 4 or x[i] == 717 + 5 or x[i] == 717 - 1 or x[i] == 717 - 2 or x[i] == 717 - 3 or x[i] == 717 - 4 or x[i] == 717 - 5:
				if y[i] == 163 or y[i] == 163 + 1 or y[i] == 163 + 2 or y[i] == 163 + 3 or y[i] == 163 + 4 or y[i] == 163 + 5 or y[i] == 163 - 1 or y[i] == 163 - 2 or y[i] == 163 - 3 or y[i] == 163 - 4 or y[i] == 163 - 5: 
					x_1star.insert(0,x[i])
					fluxes_1star.insert(0,total_flux[i])
					y_1star.insert(0,y[i])
					radius_1star.insert(0,radius[i])
			if x[i] == 611 or x[i] == 611 + 1 or x[i] == 611 + 2 or x[i] == 611 + 3 or x[i] == 611 + 4 or x[i] == 611 + 5 or x[i] == 611 - 1 or x[i] == 611 - 2 or x[i] == 611 - 3 or x[i] == 611 - 4 or x[i] == 611 - 5:
				if y[i] == 2677 or y[i] == 2677 + 1 or y[i] == 2677 + 2 or y[i] == 2677 + 3 or y[i] == 2677 + 4 or y[i] == 2677 + 5 or y[i] == 2677 - 1 or y[i] == 2677 - 2 or y[i] == 2677 - 3 or y[i] == 2677 - 4 or y[i] == 2677 - 5: 
					x_1star.insert(0,x[i])
					fluxes_1star.insert(0,total_flux[i])
					y_1star.insert(0,y[i])
					radius_1star.insert(0,radius[i])
			if x[i] == 542 or x[i] == 542 + 1 or x[i] == 542 + 2 or x[i] == 542 + 3 or x[i] == 542 + 4 or x[i] == 542 + 5 or x[i] == 542 - 1 or x[i] == 542 - 2 or x[i] == 542 - 3 or x[i] == 542 - 4 or x[i] == 542 - 5:
				if y[i] == 1165 or y[i] == 1165 + 1 or y[i] == 1165 + 2 or y[i] == 1165 + 3 or y[i] == 1165 + 4 or y[i] == 1165 + 5 or y[i] == 1165 - 1 or y[i] == 1165 - 2 or y[i] == 1165 - 3 or y[i] == 1165 - 4 or y[i] == 1165 - 5: 
					x_1star.insert(0,x[i])
					fluxes_1star.insert(0,total_flux[i])
					y_1star.insert(0,y[i])
					radius_1star.insert(0,radius[i])
			if x[i] == 475 or x[i] == 475 + 1 or x[i] == 475 + 2 or x[i] == 475 + 3 or x[i] == 475 + 4 or x[i] == 475 + 5 or x[i] == 475 - 1 or x[i] == 475 - 2 or x[i] == 475 - 3 or x[i] == 475 - 4 or x[i] == 475 - 5:
				if y[i] == 991 or y[i] == 991 + 1 or y[i] == 991 + 2 or y[i] == 991 + 3 or y[i] == 991 + 4 or y[i] == 991 + 5 or y[i] == 991 - 1 or y[i] == 991 - 2 or y[i] == 991 - 3 or y[i] == 991 - 4 or y[i] == 991 - 5: 
					x_1star.insert(0,x[i])
					fluxes_1star.insert(0,total_flux[i])
					y_1star.insert(0,y[i])
					radius_1star.insert(0,radius[i])
			if x[i] == 164 or x[i] == 164 + 1 or x[i] == 164 + 2 or x[i] == 164 + 3 or x[i] == 164 + 4 or x[i] == 164 + 5 or x[i] == 164 - 1 or x[i] == 164 - 2 or x[i] == 164 - 3 or x[i] == 164 - 4 or x[i] == 164 - 5:
				if y[i] == 2696 or y[i] == 2696 + 1 or y[i] == 2696 + 2 or y[i] == 2696 + 3 or y[i] == 2696 + 4 or y[i] == 2696 + 5 or y[i] == 2696 - 1 or y[i] == 2696 - 2 or y[i] == 2696 - 3 or y[i] == 2696 - 4 or y[i] == 2696 - 5: 
					x_1star.insert(0,x[i])
					fluxes_1star.insert(0,total_flux[i])
					y_1star.insert(0,y[i])
					radius_1star.insert(0,radius[i])

		#write the fluxes for each fits file onto a csv file
		with open('fluxes.csv','a')as fd:
			writer=csv.writer(fd)
			writer.writerow(time)		
		with open('fluxes.csv','a')as fd:
			writer=csv.writer(fd)
			for val in coord_rad_flux:
				writer.writerow(val)

		
		
	mask=0.
	thrs=0.
	y_size=0
	x_size=0
	total_flux=[]
	time=[]
	date=[]
	x=[]
	y=[]
	centers=[]
	radius=[]
	labels=[]
	
		
		
		
		
		
	











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

	

#set a variable equal to your file#
with open('cor_Data.txt') as f:
	fitfiles=f.read().split()
data_path=[]
for file_nr,fitfile in enumerate(fitfiles):
	with open(fitfile) as f:
		
		#Identify the stars found in the previous image
		x_stars=[]
		y_stars=[]
		if file_nr != 0:		
			x_stars= x
			y_stars= y

		

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
	
		for i in range(len(radius)):		
		  	if radius[i] > 0.:
				numbers=range(-int(radius[i])-2,int(radius[i])+3)
		  		for n1 in numbers:
					for n2 in numbers:
						if dist(x[i],y[i],x[i]+n1,y[i]+n2) <= radius[i]:
							total+=dataset[x[i]+n1,y[i]+n2]
		
			total_flux.append(total)
			total=0.
		
		total_flux=[j for j in total_flux if j != 0]


		#get time and date for each file
		time=[]
		date=[]
		time.append(header["TIME-OBS"])
		date.append(header["DATE-OBS"])
		for i in range(len(time)):		
			time.insert(0,date[0])

		#create list to separate data in fluxes document
		data_separater=[]

		#create information table for each file		
		ids=range(1,len(x)+1)

		if file_nr==0:
			#create master-table
			master_x=x
			master_y=y
			master_radius=radius
			master_flux=total_flux
			sum_elements=[0] * len(master_x)

		
		iter_x=[]
		iter_y= []
		iter_run=range(-5,6)
		for i in range(len(master_x)):	
			iter_x.append(master_x[i]-5)
			iter_x.append(master_x[i]-4)
			iter_x.append(master_x[i]-3)
			iter_x.append(master_x[i]-2)
			iter_x.append(master_x[i]-1)
			iter_x.append(master_x[i])
			iter_x.append(master_x[i]+1)
			iter_x.append(master_x[i]+2)
			iter_x.append(master_x[i]+3)
			iter_x.append(master_x[i]+4)
			iter_x.append(master_x[i]+5)
		for i in range(len(master_y)):	
			iter_x.append(master_y[i]-5)
			iter_x.append(master_y[i]-4)
			iter_x.append(master_y[i]-3)
			iter_x.append(master_y[i]-2)
			iter_x.append(master_y[i]-1)
			iter_x.append(master_y[i])
			iter_x.append(master_y[i]+1)
			iter_x.append(master_y[i]+2)
			iter_x.append(master_y[i]+3)
			iter_x.append(master_y[i]+4)
			iter_x.append(master_y[i]+5)
		#counting the number of times each star is seen in all the fits files			
		for i in range(len(master_x)):
			for j in range(len(x)):		
				for k in iter_run:
					if master_x[i] in range(x[j] - iter_run[k],x[j]+1+iter_run[k]):
						if master_y[i] in range(y[j] - iter_run[k], y[j] + 1 + iter_run[k]):
							sum_elements[i] += 1
							break
	

		#create lists for x, y, radius, and flux for stars that show up in the previous file as 		well as the current file in order to take into account the movement of the sky
		n_x=[]
		n_y=[]
		n_radius=[]
		n_flux=[]

			
		#only taking data from stars that appear in the previous image in order to remove data
		#from background noise
		if file_nr !=0:
			for i in range(len(x)):
				for j in range(len(x_stars)):
					for k in iter_run:
						if x[i] in range(x_stars[j] - iter_run[k],x_stars[j]+1+iter_run[k]):
							if y[i] in range(y_stars[j] - iter_run[k], y_stars[j] + 1 + iter_run[k]):
								n_x.insert(0,x[i])
								n_y.insert(0,y[i])
								n_radius.insert(0,radius[i])
								n_flux.insert(0,total_flux[i])
								break

									
			#add the extra stars into master list, later we will decide they are background 			#noise	
			for i in range(len(x)):
					if x[i] not in iter_x and y[i] not in iter_y:
						master_x.append(x[i])
						master_y.append(y[i])
						master_radius.append(radius[i])
						master_flux.append(total_flux[i])
						sum_elements.append(1)


		
		if file_nr == 0:
			for i in range(len(x)):
				id_flux=[total_flux[i]]
		for i in range(len(x)):			
			print id_flux

			
		#create a master table with info on each star that is found, including the number of 			times it is found out of all the files
		master_table=zip(master_x,master_y,sum_elements)
		





		#combine all data
		coord_rad_flux=zip(n_x,n_y,n_radius,n_flux)
		
		#write the x, y, radius, and flux values for each fits file onto a csv file	
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
	centers=[]
	radius=[]
	labels=[]

  #write the master table of the x values, y values, radii, and fluxes into a csv file	
  with open('master_table.csv','a')as fd:
	writer=csv.writer(fd)
	for val in master_table:
		writer.writerow(val)
		
		
		
	











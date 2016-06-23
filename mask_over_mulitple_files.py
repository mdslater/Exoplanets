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


#define function to write master table
def get_next_pos(master_table_var,shift_y,shift_x,it_max,y1,x1,radius,flux,master_ident,master_y1,master_x1,master_r,master_flux,master_count):

	y1=int(y1+y_shift)
	x1=int(x1+x_shift)
	y_new = 1000.
	x_new = 1000.
	check1= False
	check2 = False
	for i, test_id in enumerate(master_ident):	
		for it_no in range(it_max+1):
			
			if y1 + it_no == master_y1[i] and y_new == 1000.:
				check1 = True
				y_new = (y1 + it_no) * 1.	
				#print 'test1'
				#print y_new
			elif y1 - it_no == master_y1[i] and y_new == 1000.:
				check1 = True
				y_new = (y1 - it_no) * 1.
				#print 'test2'
			if x1 + it_no == master_x1[i] and x_new == 1000.:
				check2 = True
				x_new = (x1 + it_no) * 1.
				#print 'test3'
			elif x1 - it_no == master_x1[i] and x_new == 1000.:
				check2 = True
				x_new = (x1 - it_no) * 1.
				#print 'test4'
			if check1 == True and check2 == True:
				star_id = master_ident[i]
				master_count[i] += 1
		#print 'test5'
				break
		if check1 == True and check2 == True:
			break	
	
		y_new = 1000.
		x_new = 1000.
		check1= False
		check2 = False
		
		
	if y_new == 1000. or x_new == 1000.:	
		print 'test6'
		y_new = int((y1 - y_shift) * 1.)
		x_new = int((x1 - x_shift) * 1.)
		star_id=len(master_ident)+1
		master_ident.append(star_id)	
		master_y1.append(y_new)
		master_x1.append(x_new)
		master_r.append(radius)
		master_flux.append(flux)
		master_count.append(1)
	#master_table_var=zip(master_ident,master_y1,master_x1,master_r,master_flux,master_count)
	master_table_var=[master_ident,master_y1,master_x1,master_r,master_flux,master_count]

	return star_id, y_new, x_new, radius, flux, master_table_var




#define function that finds radius of star
def dist(x, y, xx, yy):
	x = x * 1.
	y = y * 1.
	xx = xx * 1.
	yy = yy * 1.

	return np.sqrt((yy-y)**2+(xx-x)**2) 

id_flux=[]	

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





		if file_nr==0:
			#create master-table
			master_id=range(len(x))
			master_y=y * 1
			master_x=x * 1
			master_radius=radius * 1 
			master_flux=total_flux * 1 
			sum_elements=[1] * len(master_x)
			master_table=zip(master_id,master_y,master_x,master_radius,master_flux,sum_elements)
			max_flux=0.
			for ii,fluxx in enumerate(master_flux):
				if fluxx > max_flux:
					max_flux=fluxx
					y_zerop = y[ii] * 1.
					x_zerop = x[ii] * 1.
					
			y_shift = 0
			x_shift = 0
			
		#search where the pixel with the maximum intensity is located in each file
		#we will use this location to calculate how much each pixel moved relative to the first 		#file
		if file_nr!= 0:
			max_flux=0.
			for i,fluxx in enumerate(total_flux):
				if fluxx > max_flux:
					max_flux=fluxx
					y_zerop_n = int(y[i])
					x_zerop_n = int(x[i])
			
			#calculate the x and y shift relative to the first file	
			y_shift = y_zerop_n-y_zerop
			x_shift = x_zerop_n-x_zerop
		print y_shift, x_shift
		
		




		if file_nr == 0:
			id_flux=[[] for _ in range(len(master_x))]
			id_x=[[] for _ in range(len(master_x))]
			id_y=[[] for _ in range(len(master_x))]
			for i in range(len(master_x)):
				id_flux[i].append(total_flux[i])
				id_x[i].append(range(x[i]-5,x[i]+6))
				id_y[i].append(range(y[i]-5,y[i]+6))







		for i in range(len(y)):
			master_table=get_next_pos(master_table,y_shift,x_shift,5,y[i],x[i],radius[i],total_flux[i],master_id,master_y,master_x,master_radius,master_flux,sum_elements)[5]
			print master_table[5]


		#write the x, y, radius, and flux values for each fits file onto a csv file	
		with open('fluxes.csv','a')as fd:
			writer=csv.writer(fd)
			writer.writerow(time)
		#print master_table	
		

		



		
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



		
		
	











#====================================================
#This script is used to find the total luminosity of each star in a group of images of the same star #field. It then outputs a csv file with the luminosity vs time tables for each star detected. From there, you can graph them and use their graphs to search for any exoplanets. 
#Prepare a list ('cor_Data.txt') of your images that have already been corrected for darks and #flats.
#Create a folder titled "results" in the folder where this python is located. Your resulting CSV tables will be placed in there.
#Run this script in the folder of target images!
#Once the run is complete, you can read the fluxes of each star in the CSV file titled "lists_fluxes_per_star.csv" and the time for each flux in the file titled " flux_times.csv". Each star is represent
#====================================================

#====================================================
#          Loading necessary packages
#====================================================
import pyfits
import numpy as np
import scipy.ndimage as snd
from scipy import stats
import csv
import astropy.table
import astropy.units as u
import math
from astropy.time import Time

#====================================================
#Defining functions that will be uses later on in the script
#====================================================
total_fluxes_list=[]
#define function to write master table
def get_next_pos(master_table_var,shift_y,shift_x,it_max,y1,x1,radius,flux,master_ident,master_y1,master_x1,master_r,master_flux,master_count):
	
	y_corr=int(y1-shift_y)
	x_corr=int(x1-shift_x) 
	
	y_new = 1000.
	x_new = 1000.
	check1= False
	check2 = False
	for i, test_id in enumerate(master_ident):	
		for it_no in range(it_max+1):
			
			if y_corr + it_no == master_y1[i] and y_new == 1000.:
				check1 = True
				y_new = int(y_corr + it_no)	
				
			elif y_corr - it_no == master_y1[i] and y_new == 1000.:
				check1 = True
				y_new = int(y_corr - it_no)
				
			if x_corr + it_no == master_x1[i] and x_new == 1000.:
				check2 = True
				x_new = int(x_corr + it_no)
				
			elif x_corr - it_no == master_x1[i] and x_new == 1000.:
				check2 = True
				x_new = int(x_corr - it_no)
				
			if check1 == True and check2 == True:
				star_id = master_ident[i]
				master_count[i] += 1
		
				break
		if check1 == True and check2 == True:
			break	
	
		y_new = 1000.
		x_new = 1000.
		check1= False
		check2 = False
		
		
	if y_new == 1000. or x_new == 1000.:	
		
		star_id=len(master_ident)+1
		master_ident.append(star_id)	
		master_y1.append(y_corr)
		master_x1.append(x_corr)
		master_r.append(radius)
		master_flux.append(flux)
		master_count.append(1)
	
	master_table_var=[master_ident,master_y1,master_x1,master_r,master_flux,master_count]

	return star_id, y_new, x_new, master_table_var, master_ident, master_y1, master_x1, master_r, master_flux, master_count


#define function that creates a list of fluxes over time per star

def flux_tables(shift_y,shift_x,it_max,y1,x1,flux,master_ident,master_y1,master_x1,time,radius):
	y_corr=int(y1-shift_y)
	x_corr=int(x1-shift_x) 
	y_new = 1000.
	x_new = 1000.
	check1= False
	check2 = False
	for i in range(len(master_ident)):
		for it_no in range(it_max+1):
			if y_corr + it_no == master_y1[i] and y_new == 1000.:
				check1 = True
				y_new = int(y_corr + it_no)	
				
			elif y_corr - it_no == master_y1[i] and y_new == 1000.:
				check1 = True
				y_new = int(y_corr - it_no)
				
			if x_corr + it_no == master_x1[i] and x_new == 1000.:
				check2 = True
				x_new = int(x_corr + it_no)	
				
			elif x_corr - it_no == master_x1[i] and x_new == 1000.:
				check2 = True
				x_new = int(x_corr - it_no)
				break

		if check1 == True and check2 == True:
			star_fluxes[i].append(flux)
			flux_time[i].append(time)
			star_radii[i].append(radius)
			
			break




#define function that finds radius of star
def dist(x, y, xx, yy):
	x = x * 1.
	y = y * 1.
	xx = xx * 1.
	yy = yy * 1.

	return np.sqrt((yy-y)**2+(xx-x)**2) 

id_flux=[]
star_fluxes= [[] for m in range(10000)]	
flux_time= [[] for m in range(10000)]
star_radii= [[] for m in range(10000)]

#set a variable equal to your file#
with open('cor_Data.txt') as f:
	fitfiles=f.read().split()
data_path=[]
for file_nr,fitfile in enumerate(fitfiles):
	with open(fitfile) as f:
		print 'Current image:',
		print fitfile
		
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
		thrs = 165
		

		
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
		#Radius checked with mask.fits file, they match up perfectly.

		for i in range(len(radius)):
			radius[i] = radius[i] * 2.3
			radius[i] = int(round(radius[i]))



		#resultant radii will have decimal points, round them to nearest integer.
		for i in range(len(radius)):
			radius[i] = int(round(radius[i]))
 		

			
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
				numbers=range(-int(radius[i])-10,int(radius[i])+10)
		  		for n1 in numbers:
					for n2 in numbers:
						if dist(x[i],y[i],x[i]+n1,y[i]+n2) <= radius[i]:
							total+=dataset[x[i]+n1,y[i]+n2]
		
			total_flux.append(total)
			total=0.
	
		total_flux=[j for j in total_flux if j != 0]	
  		#fluxes checked using ds9, they match up perfectly.
	
		if file_nr==0:
			master_y=[]
			master_x=[]
			master_radius=[] 
			master_flux=[]
			#create master-table
			master_id=range(len(x))
			for i in range(len(x)):
				master_y.append(y[i] * 1)
				master_x.append(x[i] * 1)
				master_radius.append(radius[i] * 1) 
				master_flux.append(total_flux[i] * 1)
			sum_elements=[0] * len(master_x)
			master_table=zip(master_id,master_y,master_x,master_radius,master_flux,sum_elements)			
			max_flux = 0.
			for ii,fluxx in enumerate(master_flux):
				if fluxx > max_flux:
					max_flux=fluxx
					y_zerop = y[ii] * 1.
					x_zerop = x[ii] * 1.
					
			y_shift = 0.
			x_shift = 0.
			
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
		
		#change the radius list for file_nr != 0 so that the same stars detected before have 			#the same radius. We do this so that fluxuations in the radius of the same star does 			#not cause any fluctiations in the fluxes.


		check11 = False
		check22 = False
		for j in range(len(x)):
			x_corrected = x[j] - x_shift
			y_corrected = y[j] - y_shift
			for i in range(len(master_x)):
				for it_no in range(11):
					if x_corrected + it_no == master_x[i] and check11 == False:
						check11 = True
					elif x_corrected - it_no == master_x[i] and check11 == False:
						check11 = True
					if y_corrected + it_no == master_y[i] and check22 == False:
						check22 = True
					elif y_corrected - it_no == master_y[i] and check22 == False:
						check22 = True
					if check11 == True and check22 == True: 
						radius[j] = master_radius[i]				
	
				check11 = False
				check22 = False

			


		
		#calculate for fluxes of each star now using the corrected radii
		total_flux=[]
		total=0.
	
		for i in range(len(radius)):		
		  	if radius[i] > 0.:
				numbers=range(-int(radius[i])-10,int(radius[i])+10)
		  		for n1 in numbers:
					for n2 in numbers:
						if dist(x[i],y[i],x[i]+n1,y[i]+n2) <= radius[i]:
							total+=dataset[x[i]+n1,y[i]+n2]
		
			total_flux.append(total)
			total=0.
		
		total_flux=[j for j in total_flux if j != 0]	

		






		#get time and date for each file
		mjd=[]
		time=[]
		date=[]
		time.append(header["TIME-OBS"])
		date.append(header["DATE-OBS"])
		mjd=zip(time,date)










		for i in range(len(y)):
			master_table=get_next_pos(master_table,y_shift,x_shift,10,y[i],x[i],radius[i],total_flux[i],master_id,master_y,master_x,master_radius,master_flux,sum_elements)[5]


		if file_nr == 0:
			caly = y
			calx = x
			calrad = radius

		#calculate for fluxes of each star that shows up in the first image. We then use the 			#sum of these fluxes for each image to calibrate the fluxes for each image
		cal_total_flux=[]
		total=0.
		
		for i in range(len(calrad)):		
		  	if calrad[i] > 0.:
				numbers=range(-int(calrad[i])-10,int(calrad[i])+10)
		  		for n1 in numbers:
					for n2 in numbers:
						if dist(calx[i],caly[i],calx[i]+n1,caly[i]+n2) <= calrad[i]:
							total+=dataset[calx[i]+x_shift+n1,caly[i]+y_shift+n2]
		
			cal_total_flux.append(total)
			total=0.
		sum_cal_total_flux = np.sum(cal_total_flux)

		





		 
		#calibrate fluxes for each image by dividing the fluxes by the sum from above
		sum_total_flux = np.median(total_flux)
		
		for i in range(len(total_flux)):
			total_flux[i]=total_flux[i]/sum_cal_total_flux

		
			
		

		for i in range(len(y)):		
			flux_tables(y_shift,x_shift,10,y[i],x[i],total_flux[i],master_id,master_y,master_x,time,radius[i])


		with open('./results/y,x,radii,flux.csv','a')as fd:
			writer=csv.writer(fd)
			for val in mjd:
				writer.writerow(val)
		image_info = zip(y,x,radius,total_flux)
		with open('./results/y,x,radii,flux.csv','a')as fd:
			writer=csv.writer(fd)
			for val in image_info:
				writer.writerow(val)
		
	


		






corr_star_fluxes = []
for i in range(len(star_fluxes)):
	if len(star_fluxes[i]) > 3:
		corr_star_fluxes.append(star_fluxes[i])
corr_flux_times = []
for i in range(len(flux_time)):
	if len(flux_time[i]) > 3:
		corr_flux_times.append(flux_time[i])
corr_star_radii = []
for i in range(len(star_radii)):
	if len(star_radii[i]) > 3:
		corr_star_radii.append(star_radii[i])

master_table = zip(master_id,master_y,master_x,master_radius,master_flux,sum_elements)

#write the master table of the x values, y values, radii, and fluxes into a csv file	
with open('./results/master_table.csv','a')as fd:
	writer=csv.writer(fd)
	for val in master_table:
		writer.writerow(val)
with open('./results/lists_fluxes_per_star.csv','a')as fd:
	writer=csv.writer(fd)
	for val in corr_star_fluxes:
		writer.writerow(val)
with open('./results/flux_times.csv','a')as fd:
	writer=csv.writer(fd)
	for val in corr_flux_times:
		writer.writerow(val)
with open('./results/star_radii.csv','a')as fd:
	writer=csv.writer(fd)
	for val in corr_star_radii:
		writer.writerow(val)


	



		
		
	











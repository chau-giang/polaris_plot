
# import astropy.constants as ac
import numpy as np
import matplotlib.pyplot as plt
import os
import astropy.units as au
import astropy.constants as ac
from astropy.io import fits
import math 
import scipy
from scipy.interpolate import interp2d
from matplotlib.pyplot import cm
import matplotlib
import scipy
from scipy.interpolate import interp1d
#plt.close('all')
 

class plot_from_file:
	def __init__(self, filename, direction):
		self.filename = filename
		self.direction = direction
				
	# Plot physical properties from the "physical_parameter file"
	def plot_Td_d(self):
		data = np.genfromtxt(self.filename)
		d = data[:, 0]/ac.au.to('m').value
		Td = data[:, 2]
		plt.plot(d, Td)
		plt.xscale('log')
		plt.yscale('log')
		if self.direction  == 0:  # along x axis
			plt.xlabel('x axis [au]')
		else:
			plt.xlabel('z axis [au]')
		plt.ylabel('Dust temperature (K)')
		plt.show()
    
    
	def plot_Tg_d(self):
		data = np.genfromtxt(self.filename)
		d = data[:, 0]/ac.au.to('m').value
		Tg = data[:, 3]
		plt.plot(d, Tg)
		plt.xscale('log')
		plt.yscale('log')
		if self.direction  == 0:  # along x axis
	   		plt.xlabel('x axis [au]')
		else:
			plt.xlabel('z axis [au]')
		plt.ylabel('Gas temperature (K)')
		plt.show()
	
	
	#===============================================================
	# Plot grain size from the "grain_size_file"			
	def plot_aalign_d(self):
		data = np.genfromtxt(self.filename)
		d = data[:, 0]/ac.au.to('m').value
		a_align = data[:, 1]*1e6
		plt.plot(d, a_align)
		plt.xscale('log')
		plt.yscale('log')
		if self.direction == 0:
			plt.xlabel('x axis [au]')
		else:
			plt.ylabel('z axis [au]')
		plt.ylabel(r'$\sf{\rm a_{\rm align} (\mu m)}$')
		#plt.show()

        
	def plot_adisr_d(self, sav=True):
		data = np.genfromtxt(self.filename)
		d = data[:, 0]/ac.au.to('m').value
		a_disr = data[:, 2]*1e6
		a_disr_max = data[:, 3]*1e6

		# remove the part outside the disruption region, at which
		# a_disr = a_disr_max
		index = np.where(a_disr == a_disr_max)[0]
		index = index[1:]
		a_disr[index] = np.nan
		a_disr_max[index] = np.nan

		plt.plot(d, a_disr)
		plt.plot(d, a_disr_max)
		plt.xscale('log')
		plt.yscale('log')
		if self.direction == 0:
			plt.xlabel('x axis [au]')
		else:
			plt.ylabel('z axis [au]')
		plt.ylabel(r'$\sf{\rm a_{\rm disr} - a_{\rm disr,max} (\mu m)}$')
		#plt.show()

	def plot_amax_aJ_lowJ_d(self, sav=True):
		data = np.genfromtxt(self.filename)
		d = data[:, 0]/ac.au.to('m').value
		a_min_aJ_lowJ = data[:, 4]*1e6
		a_max_aJ_lowJ = data[:, 5]*1e6

		# remove the part outside the disruption region, at which
		# a_min_aJ_lowJ = a_max_aJ_lowJ
		index = np.where(a_min_aJ_lowJ == a_max_aJ_lowJ)[0]
		index = index[1:]
		a_min_aJ_lowJ[index] = np.nan
		a_max_aJ_lowJ[index] = np.nan

		plt.plot(d, a_min_aJ_lowJ)
		plt.plot(d, a_max_aJ_lowJ)
		plt.xscale('log')
		plt.yscale('log')
		if self.direction == 0:
			plt.xlabel('x axis [au]')
		else:
			plt.ylabel('z axis [au]')
		plt.ylabel(r'$\sf{\rm a_{\rm min,aJ}^{\rm low-J} - a_{\rm max,aJ}^{\rm low-J} (\mu m)}$')
		#plt.show()
        
	def plot_amax_aJ_highJ_d(self, sav=True):
		data = np.genfromtxt(self.filename)
		d = data[:, 0]/ac.au.to('m').value
		a_min_aJ_highJ = data[:, 6]*1e6
		a_max_aJ_highJ = data[:, 7]*1e6

		# remove the part outside the disruption region, at which
		# a_min_aJ_highJ = a_max_aJ_highJ
		index = np.where(a_min_aJ_highJ == a_max_aJ_highJ)[0]
		index = index[1:]
		a_min_aJ_highJ[index] = np.nan
		a_max_aJ_highJ[index] = np.nan

		plt.plot(d, a_min_aJ_highJ)
		plt.plot(d, a_max_aJ_highJ)
		plt.xscale('log')
		plt.yscale('log')
		if self.direction == 0:
			plt.xlabel('x axis [au]')
		else:
			plt.ylabel('z axis [au]')
		plt.ylabel(r'$\sf{\rm a_{\rm min,aJ}^{\rm high-J} - a_{\rm max,aJ}^{\rm high-J} (\mu m)}$')
		#plt.show()


	def plot_adg_50_d(self, sav=True):
		data = np.genfromtxt(self.filename)
		d = data[:, 0]/ac.au.to('m').value
		a_min_dg_50 = data[:, 8]*1e6
		a_max_dg_50 = data[:, 9]*1e6

		# remove the part outside the disruption region, at which
		# a_min_aJ_highJ = a_max_aJ_highJ
		index = np.where(a_min_dg_50 == a_max_dg_50)[0]
		index = index[1:]
		a_min_dg_50[index] = np.nan
		a_max_dg_50[index] = np.nan

		plt.plot(d, a_min_dg_50)
		plt.plot(d, a_max_dg_50)
		plt.xscale('log')
		plt.yscale('log')
		if self.direction == 0:
			plt.xlabel('x axis [au]')
		else:
			plt.ylabel('z axis [au]')
		plt.ylabel(r'$\sf{\rm a_{\rm min,JB}^{\rm DG,50} - a_{\rm max,JB}^{\rm DG,50} (\mu m)}$')
		#plt.show()


	def plot_adg_100_d(self, sav=True):
		data = np.genfromtxt(self.filename)
		d = data[:, 0]/ac.au.to('m').value
		a_min_dg_100 = data[:, 8]*1e6
		a_max_dg_100 = data[:, 9]*1e6

		# remove the part outside the disruption region, at which
		# a_min_aJ_highJ = a_max_aJ_highJ
		index = np.where(a_min_dg_100 == a_max_dg_100)[0]
		index = index[1:]
		a_min_dg_100[index] = np.nan
		a_max_dg_100[index] = np.nan

		plt.plot(d, a_min_dg_100)
		plt.plot(d, a_max_dg_100)
		plt.xscale('log')
		plt.yscale('log')
		if self.direction == 0:
			plt.xlabel('x axis [au]')
		else:
			plt.ylabel('z axis [au]')
		plt.ylabel(r'$\sf{\rm a_{\rm min,JB}^{\rm DG,100} - a_{\rm max,JB}^{\rm DG,100} (\mu m)}$')
		#plt.show()

 
 
class plot_line_from_fits_file:
	def __init__(self, filename, index, direction, plane, colors, labels, rmax):
		self.filename = filename       # filename
		self.index = index			   # index of output parameter: (0): Td, (1): Tg, (2): aalig, (3): U
									   # (4): cos(psy), (5): gamma_rad
		self.direction = direction	   # along horizontal (0) or vertical (1) direction of the midplane
		self.plane = plane			   # xy (0) or xz (1) plane
		self.colors = colors		   # color of line
		self.labels = labels		   # label of line
		self.rmax = rmax 			   # outer boundary of grid

	def line_from_fits(self):
		# To take data along the horizontal or vertical direction of the midplane as the function of distances
		hdulist = fits.open(self.filename)
		hdu = hdulist[0]
		nr_pixel = int(hdu.header[3] / 2)
		
		data = hdu.data[self.index, self.plane, :, :]
		if self.direction == 0: # take data from the horizontal direction
			data = data[nr_pixel, 0:nr_pixel]
		else: # take data from the vertical direction
			data = data[0:nr_pixel, nr_pixel-4]


		d = np.linspace(0, self.rmax, nr_pixel)[::-1]

		if self.index == 2: #align
			data = data*1e6
		if self.index == 4: #psi
			data = np.arccos(data)*180/math.pi
			

		plt.plot(d, data, color = self.colors, label = self.labels)
		
		y_label = np.array(['Dust temperature (K)', 'Gas temperature (K)', r'$\sf{\rm a_{\rm align} (\mu m)}$', 
							'Radiation field strength U', r'$\sf{\rm \cos(\theta)}$', r'$\sf{\gamma_{\rm rad}}$'])
		plt.ylabel(y_label[self.index])
		plt.xlabel('Radius from the center region')
		if (self.index == 2) or (self.index == 3):
			plt.yscale('log')
		plt.legend()
		#plt.show()
    

 
class plot_fits_file:
	def __init__(self, filename, index, plane, title_local):
		self.filename = filename
		self.index = index
		self.plane = plane
		self.title_local = title_local

	def map(self):
		hdulist = fits.open(self.filename)
		hdu = hdulist[0]
		nr_pixel = int(hdu.header[3] / 2)
				
		if self.index == 2: # grain alignment size
		    plt.imshow(hdu.data[self.index, self.plane, :, :]*1e6, origin = 'lower', norm = matplotlib.colors.LogNorm())
		elif self.index == 3: # radiation field strength
		    plt.imshow(hdu.data[self.index, self.plane, :, :], origin = 'lower', norm = matplotlib.colors.LogNorm())
		elif self.index == 4: # angle between radiation and magnetic field
		    plt.imshow(np.arccos(hdu.data[self.index, self.plane, :, :])*180/math.pi, origin = 'lower')
		else:
		    plt.imshow(hdu.data[self.index, self.plane, :, :], origin = 'lower')
		    
		xaxis = np.array(['-15000', '-10000', '-5000', '0', '5000', '10000', '15000'])
		plt.xticks(np.linspace(0, nr_pixel*2, 7), xaxis)
		plt.yticks(np.linspace(0, nr_pixel*2, 7), xaxis)
		if self.plane == 0: # on x-y plane
		    plt.xlabel('x [au]')
		    plt.ylabel('y [au]')
		else: # on x-z plane
		    plt.xlabel('x [au]')
		    plt.ylabel('z [au]')
		    
		title_prop = np.array([r'$\sf{\rm T_{\rm d} (K)}$', r'$\sf{\rm T_{\rm g} (K)}$', r'$\sf{\rm a_{\rm align} (\mu m)}$', 
							r'$\sf{\rm u_{\rm rad}/u_{\rm ISRF}}$', r'$\sf{\rm \cos(\theta)}$', r'$\sf{\gamma_{\rm rad}}$'])					
		plt.title(str(title_prop[self.index]) + '-' + str(self.title_local))
		plt.colorbar()
		plt.show()
		

class polarized_thermal:
    def __init__(self, filename, color_name, label_name):   
        self.filename = filename
        self.color_name = color_name
        self.label_name = label_name

    def index_radius(self, distance, nr_pixel): # find the grid in pencil beam (1 grid cell) at chosen distance
            x = np.linspace(0, nr_pixel*2, nr_pixel*2)   
            X, Y = np.meshgrid(x, x)
            R = np.sqrt((X-nr_pixel)**2 + (Y-nr_pixel)**2)
            
            index = np.where(R < distance)		# remove grid inside the distance
            coor_x = index[0]
            coor_y = index[1]
            R[coor_x, coor_y] = np.nan
            
            index = np.where(R > distance+1)	#remove grid outside the distance
            coor_x = index[0]
            coor_y = index[1]
            R[coor_x, coor_y] = np.nan	
            
            index = np.where(np.isnan(R) == False)
            coor_x = index[0]
            coor_y = index[1]
            return coor_x, coor_y
        


    def P_I_matrix(self):           	# get polarization degree, optical depth, column density and
        hdulist = fits.open(self.filename)      # fits file of Stoke parameter and optical depth
        hdu = hdulist[0]
        nr_pixel = np.int(hdu.data.shape[2]/2)

        max_grid = 0
        for i in range(0, nr_pixel):
                coor_x, coor_y = self.index_radius(i, nr_pixel)
                if (len(coor_x) > max_grid):
                        max_grid = len(coor_x)
        
        I_matrix = np.zeros([max_grid, nr_pixel])  # intensity maxtrix
        P_matrix = I_matrix.copy()                 # polarization matrix: row:  distance from the center to the outer
                                                   #				column: data in all grid in this distance
        tau_matrix = I_matrix.copy()
        NH_matrix = I_matrix.copy()
        nr_point = np.zeros(nr_pixel)
        
        for i in range(0, nr_pixel):	
            coor_x, coor_y = self.index_radius(i, nr_pixel)
            
            I = hdu.data[0, 0, coor_x, coor_y]
            Q = hdu.data[1, 0, coor_x, coor_y]
            U = hdu.data[2, 0, coor_x, coor_y]
            P = np.sqrt(Q**2 + U**2)/I
            tau = hdu.data[4, 0, coor_x, coor_y]
            NH = hdu.data[5, 0, coor_x, coor_y]*1e-4 #cm-2
            cout = 0
            for j in range(0, len(P)):
                I_matrix[j,i] = I[j]
                P_matrix[j,i] = P[j]				
                tau_matrix[j,i] = tau[j]
                NH_matrix[j,i] = NH[j]
                cout += 1
            nr_point[i] = cout
                    
        return I_matrix, P_matrix, tau_matrix, NH_matrix, nr_point 


    # Function to plot the polarization degtee of the fix wavelength with normlalized intensity
    def mean_value(self, I_matrix, P_matrix, tau_matrix, NH_matrix, nr_point):
        I_mean = np.zeros(len(I_matrix[0, :]))
        P_mean = I_mean.copy()
        tau_mean = I_mean.copy()
        NH_mean = I_mean.copy()
        for i in range(0, len(I_mean)):
            I_mean[i] = sum(I_matrix[:, i])/nr_point[i]
            P_mean[i] = sum(P_matrix[:, i])/nr_point[i]
            tau_mean[i] = sum(tau_matrix[:, i])/nr_point[i]
            NH_mean[i] = sum(NH_matrix[:, i])/nr_point[i]
        return I_mean, P_mean, tau_mean, NH_mean		
        


    def plot_P_I(self):
        I_matrix, P_matrix, tau_matrix, NH_matrix, nr_point = self.P_I_matrix()
        I_mean, P_mean, tau_mean, NH_mean = self.mean_value(I_matrix, P_matrix, tau_matrix, NH_matrix, nr_point)

        fig,ax = plt.subplots(figsize = (6, 4.5))
        ax.scatter(I_matrix/np.max(I_matrix), P_matrix*100, alpha = 0.01, s = 0.2, color = self.color_name)
        ax.plot(I_mean/max(I_mean), P_mean*100, color = self.color_name, linewidth = 2, label = self.label_name)


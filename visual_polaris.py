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
import struct
 
class remove_index:
	def index_distance(d, x1, x2=[]):
        # only take data from the center to the outer 
		index = np.where(d >= 0)[0]
		d = d[index]
		x1 = x1[index]
		if len(x2)==0:
			return d, x1
		else:
			x2 = x2[index]
		return d, x1, x2

	def index_no_data(x1, x2):
		# remove data outside the active region
		index = np.where(x1 == x2)[0]
		#index = index[1:]
		x1[index] = np.nan
		x2[index] = np.nan
		return x1, x2


def find_index_nearest(array, value):
    """
    Args:
        array: list of value
        value: value to find

    Returns:
        index: index where the value in array is closest to the input value
    """

    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return idx


class plot_from_file:
	def __init__(self, filename, direction):	
		self.filename = filename   # link to file data
		self.direction = direction
				
	def plot_data(self, data_x, data_y, label):
		plt.plot(data_x, data_y)
		plt.xscale('log')
		plt.yscale('log')
		if 'x' in self.direction:
			plt.xlabel('x axis [au]')
		else:
			plt.ylabel('z axis [au]')
		plt.ylabel(label)
		plt.show()

	def plot_size_min_max(self, data_x, data_y_min, data_y_max, label):
		plt.plot(data_x, data_y_min)
		plt.plot(data_x, data_y_max)
		plt.xscale('log')
		plt.yscale('log')
		if 'x' in self.direction:
			plt.xlabel('x axis [au]')
		else:
			plt.ylabel('z axis [au]')
		plt.ylabel(label)
		plt.show()

	# Plot physical properties from the "physical_parameter file"
	def plot_nH_d(self,   draw=True):
		data = np.genfromtxt(self.filename)
		d = data[:, 0]/ac.au.to('m').value
		nH = data[:, 1]/ac.m_p.to('kg').value
		d, nH = remove_index.index_distance(d, nH)

		if draw: 
			self.plot_data(d, nH, r'$\sf{\rm n_{\rm H} (cm^{-3})$')
		else:   
			return d, nH


	def plot_Td_d(self,   draw=True):
		data = np.genfromtxt(self.filename)
		d = data[:, 0]/ac.au.to('m').value
		Td = data[:, 2]
		d, Td = remove_index.index_distance(d, Td)

		if draw: 
			self.plot_data(d, Td, r'$\sf{\rm n_{\rm H} (cm^{-3})$')
		else:    
			return d, Td


	def plot_Tg_d(self,   draw=True):
		data = np.genfromtxt(self.filename)
		d = data[:, 0]/ac.au.to('m').value
		Tg = data[:, 3]
		d, Tg = remove_index.index_distance(d, Tg)

		if draw: 
			self.plot_data(d, Tg, 'Gas temperature (K)')
		else:    
			return d, Tg


	#===============================================================
	# Plot grain size from the "grain_size_file"			
	def plot_aalign_d(self,   draw=True):
		data = np.genfromtxt(self.filename)
		d = data[:, 0]/ac.au.to('m').value
		a_align = data[:, 1]*1e6
		d, a_align = remove_index.index_distance(d, a_align)

		if draw:
			self.plot_data(d, a_align, r'$\sf{\rm a_{\rm align} (\mu m)}$')
		else:   
			return d, a_align

	def plot_adisr_d(self,   draw=True):
		data = np.genfromtxt(self.filename)
		d = data[:, 0]/ac.au.to('m').value
		a_disr = data[:, 2]*1e6
		a_disr_max = data[:, 3]*1e6
		d, a_disr, a_disr_max = remove_index.index_distance(d, a_disr, a_disr_max)
		a_disr, a_disr_max = remove_index.index_no_data(a_disr, a_disr_max)

		if draw: 
			self.plot_size_min_max(d, a_disr, a_disr_max, r'$\sf{\rm a_{\rm disr} - a_{\rm disr,max} (\mu m)}$')
		else:
			return d, a_disr, a_disr_max
		
		
	def plot_aaJ_lowJ_d(self,   draw=True):
		data = np.genfromtxt(self.filename)
		d = data[:, 0]/ac.au.to('m').value
		a_min_aJ_lowJ = data[:, 4]*1e6
		a_max_aJ_lowJ = data[:, 5]*1e6
		d, a_min_aJ_lowJ, a_max_aJ_lowJ = remove_index.index_distance(d, a_min_aJ_lowJ, a_max_aJ_lowJ)
		a_min_aJ_lowJ, a_max_aJ_lowJ = remove_index.index_no_data(a_min_aJ_lowJ, a_max_aJ_lowJ)
		
		if draw: 
			self.plot_size_min_max(d, a_min_aJ_lowJ, a_max_aJ_lowJ, r'$\sf{\rm a_{\rm min,aJ}^{\rm low-J} - a_{\rm max,aJ}^{\rm low-J} (\mu m)}$')
		else:    
			return d, a_min_aJ_lowJ, a_max_aJ_lowJ
		
	def plot_aaJ_highJ_d(self,   draw=True):
		data = np.genfromtxt(self.filename)
		d = data[:, 0]/ac.au.to('m').value
		a_min_aJ_highJ = data[:, 6]*1e6
		a_max_aJ_highJ = data[:, 7]*1e6
		d, a_min_aJ_highJ, a_max_aJ_highJ = remove_index.index_distance(d, a_min_aJ_highJ, a_max_aJ_highJ)
		a_min_aJ_highJ, a_max_aJ_highJ = remove_index.index_no_data(a_min_aJ_highJ, a_max_aJ_highJ)
		
		if draw:
			self.plot_size_min_max(d, a_min_aJ_highJ, a_max_aJ_highJ, r'$\sf{\rm a_{\rm min,aJ}^{\rm high-J} - a_{\rm max,aJ}^{\rm high-J} (\mu m)}$')
		else:   
			return d, a_min_aJ_highJ, a_max_aJ_highJ



	def plot_adg_50_d(self,   draw=True):
		data = np.genfromtxt(self.filename)
		d = data[:, 0]/ac.au.to('m').value
		a_min_dg_50 = data[:, 8]*1e6
		a_max_dg_50 = data[:, 9]*1e6
		d, a_min_dg_50, a_max_dg_50 = remove_index.index_distance(d, a_min_dg_50, a_max_dg_50)
		a_min_dg_50, a_max_dg_50 = remove_index.index_no_data(a_min_dg_50, a_max_dg_50)

		if draw:
			self.plot_size_min_max(d, a_min_dg_50, a_max_dg_50, r'$\sf{\rm a_{\rm min,JB}^{\rm DG,50} - a_{\rm max,JB}^{\rm DG,50} (\mu m)}$')
		else:   
			return d, a_min_dg_50, a_max_dg_50
		
	def plot_adg_100_d(self,   draw=True):
		data = np.genfromtxt(self.filename)
		d = data[:, 0]/ac.au.to('m').value
		a_min_dg_100 = data[:, 10]*1e6
		a_max_dg_100 = data[:, 11]*1e6
		d, a_min_dg_100, a_max_dg_100 = remove_index.index_distance(d, a_min_dg_100, a_max_dg_100)
		a_min_dg_100, a_max_dg_100 = remove_index.index_no_data(a_min_dg_100, a_max_dg_100)
		
		if draw: 
			self.plot_size_min_max(d, a_min_dg_100, a_max_dg_100, r'$\sf{\rm a_{\rm min,JB}^{\rm DG,100} - a_{\rm max,JB}^{\rm DG,100} (\mu m)}$')
		else:    
			return d, a_min_dg_100, a_max_dg_100



	def amax_JB_pm(self, f_p, factor, B, draw=True, s = 0.5, n_atomic = 1.332e+28):
        # f_p: iron fraction inside grains
		d, nH = self.plot_nH_d(False) 
		d, Td = self.plot_Td_d(False)
		d, Tg = self.plot_Tg_d(False)
		dens = nH * Td * np.sqrt(Tg)
		amax_JB_pm = 4.52e-10 / factor * s * n_atomic * f_p * B / dens #m
		amax_JB_pm = amax_JB_pm * 1e6

		if draw:
			self.plot_data(d, amax_JB_pm, r'$\sf{\rm a_{\rm max,JB} - pm }$')
		else:
			return d, amax_JB_pm

	def amax_JB_spm(self, Ncl, phi_sp, factor, B, draw=True, s = 0.5):
		# f_p: iron fraction inside grains
		d, nH = self.plot_nH_d(False) 
		d, Td = self.plot_Td_d(False)
		d, Tg = self.plot_Tg_d(False)
		dens = nH * Td * np.sqrt(Tg)
		amax_JB_spm = 3.775e19 / factor * s * Ncl * phi_sp * B / dens #m 
		amax_JB_spm = amax_JB_spm * 1e6 #um
		
		if draw:
			self.plot_data(d, amax_JB_spm, r'$\sf{\rm a_{\rm max,JB} - spm }$')
		else:
			return d, amax_JB_spm

class binary_file:
	def __init__(self, filename):
		self.filename = filename
   
	def exp_list(self, start, stop, total_number, base):
		"""Calculates exponential distribution between two values.

		Args:
			start (float): starting value
			stop (float): last value
			total_number (int): total amount of distributed values.
			base (float): distribution factor


		Returns:
			number_list (list): Distributed numbers
		"""

		number_list = np.zeros(total_number + 1)
		number_list[0] = start
		number_list[total_number] = stop

		if base > 1:
			dx = (stop - start) * (base - 1.0) / (pow(base, total_number) - 1)

			for i_x in range(0, total_number):
				number_list[i_x] = start + dx * \
					(pow(base, i_x) - 1) / (base - 1.0)
		else:
			raise ValueError('only positive exp bases are allowed!')
		return number_list


	def radial_list(self, dictionary):
		"""
		Args:
			Dictionary of the grid

		Returns:
			Position of the center of the cell in the radial direction
		"""
		radial_list = self.exp_list(dictionary['Rmin'], dictionary['Rmax'], dictionary['Nr'], dictionary['fr'])/au.au.to('m')
		length_cell = radial_list[1:] - radial_list[0:len(radial_list)-1]
		radial_list = radial_list[0:len(radial_list)-1] + length_cell/2
		return radial_list


	def read_binary_grid_file(self, save_radiation=True, full_temperature=True, one_cell_azimuthal=True):
		"""
		Args: binary file
			  is binary file has only one cell in the azimuthal direction

		Returns: dictionary and matrix of data of size: nr x nth x nph x data_len of 
		"""

		with open(self.filename, 'rb') as fp:
	 
		# Basic information of grid
			grid_id, = (struct.unpack('<H', fp.read(2)))
			parameter_size, = struct.unpack('<H', fp.read(2))
			ids = struct.unpack('<' + parameter_size*'H', fp.read(2*parameter_size))
			Rmin, = struct.unpack('<d', fp.read(8))
			Rmax, = struct.unpack('<d', fp.read(8))
			Nr, = struct.unpack('<H', fp.read(2))
			Nph, = struct.unpack('<H', fp.read(2))
			Nth, = struct.unpack('<H', fp.read(2))
			fr, = struct.unpack('<d', fp.read(8))
			fph, = struct.unpack('<d', fp.read(8))
			fth, = struct.unpack('<d', fp.read(8))

			# Take data in all cells in grid and reshape it into 4D matrix
			data = struct.unpack('<' + parameter_size * Nr * Nth * Nph * 'd', fp.read(8 *  parameter_size * Nr * Nth * Nph))  
			center_point = struct.unpack('<' + parameter_size * 'd', fp.read(8 * parameter_size))
			matrix = np.reshape(data, (Nr, Nph, Nth, parameter_size))

			# Dichtionary
			parameter = dict()
			parameter['parameter_size_cell'] = parameter_size
			parameter['Rmin'] = Rmin
			parameter['Rmax'] = Rmax
			parameter['Nr'] = Nr
			parameter['Nph'] = Nph
			parameter['Nth'] = Nth
			parameter['fr'] = fr
			parameter['fph'] = fph
			parameter['fth'] = fth

			parameter['m_H'] = 0
			parameter['T_d'] = 1
			parameter['B_x'] = 2
			parameter['B_y'] = 3
			parameter['B_z'] = 4
			parameter['T_g'] = 5
			
			if full_temperature:
				parameter['nr_grain_size'] = 160 
				array = np.linspace(6, 325, parameter['nr_grain_size']*2)
				array = array.reshape(parameter['nr_grain_size'], 2)
				parameter['T_d_a'] = array[:, 0]
				parameter['abs_ini'] = array[:, 1] # take odd order

				if one_cell_azimuthal:
					parameter['a_align'] =  326
					parameter['a_min_aJ_lowJ'] =  327
					parameter['a_max_aJ_lowJ'] =  328
					parameter['a_min_aJ_highJ'] =  329
					parameter['a_max_aJ_highJ'] =  330
					parameter['a_min_JB_DG_50'] =  331
					parameter['a_max_JB_DG_50'] =  332
					parameter['a_min_JB_DG_100'] =  333
					parameter['a_max_JB_DG_100'] =  334
					parameter['anisotropic_degree'] =  335
					parameter['cos(psi)'] =  336

					if save_radiation:
						parameter['nr_wave'] = 100
						array = np.linspace(337, parameter['parameter_size_cell']-1, parameter['nr_wave']*4)
						array = array.reshape(parameter['nr_wave'], 4)
						parameter['urad'] = array[:, 0]
						parameter['ux'] = array[:, 1]
						parameter['uy'] = array[:, 2]
						parameter['uz'] = array[:, 3]
	                        
			else:
				if one_cell_azimuthal:
					parameter['a_align'] =  6
					parameter['a_min_aJ_lowJ'] =  7
					parameter['a_max_aJ_lowJ'] =  8
					parameter['a_min_aJ_highJ'] =  9
					parameter['a_max_aJ_highJ'] =  10
					parameter['a_min_JB_DG_50'] =  11
					parameter['a_max_JB_DG_50'] =  12
					parameter['a_min_JB_DG_100'] =  13
					parameter['a_max_JB_DG_100'] =  14
					parameter['anisotropic_degree'] =  15
					parameter['cos(psi)'] =  16
    
					if save_radiation:
						parameter['nr_wave'] = 100
						array = np.linspace(17, parameter['parameter_size_cell']-1, parameter['nr_wave']*4)
						array = array.reshape(parameter['nr_wave'], 4)
						parameter['urad'] = array[:, 0]
						parameter['ux'] = array[:, 1]
						parameter['uy'] = array[:, 2]
						parameter['uz'] = array[:, 3]

		return parameter, matrix

	def wavelength_list(self):
		"""
		Returns: list of wavelengths
		"""
		
		location = os.getcwd()
		wave_list = np.genfromtxt(os.path.join(location, 'wave.dat'))
		return wave_list

class fits_file_MCMC:
    def __init__(self, filename):
        self.filename = filename

    def read_fits(self):
        hdulist = fits.open(self.filename)
        hdu = hdulist[0]
       
        dat = dict()
        dat_ = hdu.data
        dat['dust_temperature'] = dat_[0,:,:,:]
        dat['gas_temperature'] = dat_[1,:,:,:]
        dat['a_align'] = dat_[2,:,:,:]
        dat['urad'] = dat_[3,:,:,:]
        dat['cos(spi)'] = dat_[4,:,:,:]
        dat['anisotropic_degree'] = dat_[5,:,:,:]
        dat['amax_aJ_lowJ'] = dat_[7,:,:,:]
        dat['amax_aJ_highJ'] = dat_[9,:,:,:]
        dat['initial_abs_rate'] = dat_[14,:,:,:]
        dat['plane_xy'] = 0
        dat['plane_xz'] = 1

        # Construct cell center position (in AU)
        h = hdu.header
        dat['domain'] = dict()
        dat['domain']['Nx'] = h['NAXIS1']
        dat['domain']['Ny'] = h['NAXIS2']
        dat['domain']['dx'] = h['CDELT1B']
        dat['domain']['dy'] = h['CDELT2B']
        dat['x'] = h['CRVAL1B'] + np.arange(dat['domain']['Nx'])*dat['domain']['dx']
        dat['y'] = h['CRVAL2B'] + np.arange(dat['domain']['Ny'])*dat['domain']['dy']
        dat['domain']['xmin'] = h['CRVAL1B'] - 0.5*h['CDELT1B']
        dat['domain']['ymin'] = h['CRVAL2B'] - 0.5*h['CDELT2B']
        dat['domain']['xmax'] = dat['x'].max() + 0.5*h['CDELT1B']
        dat['domain']['ymax'] = dat['y'].max() + 0.5*h['CDELT2B']

        return dat

    def plot_line_from_fits(self, parameter, plane, direction, colors, labels):
        # Take data along the horizontal or vertical direction of the midplane as the function of distances
        dat = self.read_fits()
        data = dat[parameter]
        data_plane = data[dat[plane]]

        if 'horizontal' in direction: # take data from the horizontal direction
            data_direction = data_plane[int(dat['Ny']/2),  0:dat['Nx']]
            d = dat['x']
        else: # take data from the vertical direction
            data_direction = data[0: dat['Ny'], int(dat['Nx']/2)]
            d = dat['y']


        grain_size = ['a_align', 'amax_aJ_lowJ', 'amax_aJ_highJ']
        if any(x in parameter for x in grain_size): 
            data_direction = data_direction*1e6
        if ('cos(psi)' in parameter): 
            data_direction = np.arccos(data_direction)*180/math.pi
                

        plt.plot(d, data_direction, color = colors, label = labels)
        plt.xlim(dat['domain']['xmin'], dat['domain']['xmax'])
        plt.xlabel('Radius from the center region [au]')
        plt.ylabel(parameter)
        if any(x in parameter for x in grain_size): 
            plt.yscale('log')
        plt.legend()
        plt.show()
    
 
    def plot_map_from_fits(self, parameter, plane):
        dat = self.read_fits()
        data = dat[parameter]
        data_plane = data[dat[plane]]
 
        grain_size = ['a_align', 'amax_aJ_lowJ', 'amax_aJ_highJ']
        if any(x in parameter for x in grain_size): 
            plt.imshow(data_plane*1e6, origin = 'lower', norm = matplotlib.colors.LogNorm())
        
        elif 'urad' in self.parameter:
            plt.imshow(data_plane, origin = 'lower', norm = matplotlib.colors.LogNorm())
        
        elif 'cos(psi)' in self.parameter:
            plt.imshow(np.arccos(data_plane)*180/math.pi, origin = 'lower')
        
        else:
            plt.imshow(data_plane, origin = 'lower')
            
        xaxis = np.array(dat['domain']['xmin'], dat['domain']['xmax'], 7)
        yaxis = np.array(dat['domain']['ymin'], dat['domain']['ymax'], 7)
        plt.xticks(np.linspace(0, dat['domain']['Nx'], 7), xaxis)
        plt.yticks(np.linspace(0, dat['domain']['Ny'], 7), yaxis)

        if 'plane_xy' in plane:
            plt.xlabel('x [au]')
            plt.ylabel('y [au]')
        else: 
            plt.xlabel('x [au]')
            plt.ylabel('z [au]')
            
        title_prop = np.array([r'$\sf{\rm T_{\rm d} (K)}$', r'$\sf{\rm T_{\rm g} (K)}$', r'$\sf{\rm a_{\rm align} (\mu m)}$', 
                             r'$\sf{\rm u_{\rm rad}/u_{\rm ISRF}}$', r'$\sf{\rm \psi}$', r'$\sf{\gamma_{\rm rad}}$',
                             r'$\sf{\rm a_{\rm min,aJ}^{\rm low-J} (\mu m)}$', r'$\sf{\rm a_{\rm max,aJ}^{\rm high-J} (\mu m)}$'])					
        plt.title(str(dat[parameter]) + '-' + str(title_local))
        plt.colorbar()
        plt.show()

class fits_file_pol_map:
	def __init__(self, filename):
		self.filename = filename

	def read_fits(self):
		"""
		Args:
			filename: directory to fits file
			index_wavelength: index of wavelengths in fits file. Default: 0: one single wavelengths in fits file

		Returns:
			dictionary of fits file
		"""

		hdulist = fits.open(self.filename)
		hdu = hdulist[0]

		dat = dict()
		dat_ = hdu.data
		head_ = hdu.header
		dat['Stoke_I'] = dat_[0, :, :, :]
		dat['Stoke_Q'] = dat_[1, :, :, :]
		dat['Stoke_U'] = dat_[2, :, :, :]
		dat['Stoke_V'] = dat_[3, :, :, :]
		dat['P'] = np.sqrt(pow(dat['Stoke_Q'], 2) + pow(dat['Stoke_U'], 2)) / dat['Stoke_I']
		dat['optical_depth'] = dat_[4, :, :, :]
		dat['column_density'] = dat_[5, :, :, :]
		dat['nr_wave'] = head_['NAXIS3']

		if dat['nr_wave'] == 1:
			dat['wavelength'] = head_['HIERARCH WAVELENGTH1']
		else:
			wave = np.zeros(dat['nr_wave'])
			for i in range (0, dat['nr_wave']):
				wave[i] =  head_['HIERARCH WAVELENGTH'+str(i+1)]
			dat['wavelength'] = wave


		dat['Nx'] = head_['NAXIS1']
		dat['Ny'] = head_['NAXIS2']
 
		dat['dx'] = head_['CDELT1B']
		dat['dy'] = head_['CDELT2B']
		dat['x'] = head_['CRVAL1B'] + np.arange(dat['Nx']) * dat['dx']
		dat['y'] = head_['CRVAL2B'] + np.arange(dat['Ny']) * dat['dy']
		dat['xmin'] = head_['CRVAL1B'] - 0.5 * dat['dx']
		dat['ymin'] = head_['CRVAL2B'] - 0.5 * dat['dy']
		dat['xmax'] = dat['x'].max() + 0.5 * dat['dx']
		dat['ymax'] = dat['y'].max() + 0.5 * dat['dy']

		return dat

	def optical_depth_distance(self, dat):
		optical_depth = dat['optical_depth']
		nr_radial = int(dat['Nx']/2)
		optical_depth_radial = optical_depth[nr_radial, nr_radial:]
		radial_direction = dat['x'][nr_radial:] 
		return radial_direction, optical_depth_radial

	def index_radius(self, dictionary, distance, nr_pixel_beam_size):  
		"""
		Args: 
			dictionary: dictionary of fits file
			distance: radius from the center region
			nr_pixel_beam_size: number of pixel in the beam size 

		Returns:
			coor_x: position of cell in the observed region on x direction
			coor_y: position of cell in the observed region on y direction
		"""
		nr_pixel = dictionary['Nx']
		x = np.arange(0, nr_pixel, 1)   
		X, Y = np.meshgrid(x, x)
		R = np.sqrt((X - nr_pixel//2)**2 + (Y - nr_pixel//2)**2)

		if ((distance != 0) and (distance - nr_pixel_beam_size > 0)): 
			index = np.where(R < distance)		# remove grid inside the distance
			coor_x = index[0]
			coor_y = index[1]
			R[coor_x, coor_y] = np.nan

		#if ((distance != (nr_pixel//2)-1) and (distance + nr_pixel_beam_size < nr_pixel)):
		index = np.where(R > distance+nr_pixel_beam_size)	#remove grid outside the distance
		coor_x = index[0]
		coor_y = index[1]
		R[coor_x, coor_y] = np.nan	

		index = np.where(np.isnan(R) == False)
		coor_x = index[0]
		coor_y = index[1]
		return coor_x, coor_y


	def P_I_matrix(self, dictionary, nr_pixel_beam_size = 1):           
		"""
		This function is to get the 
		Args:
			dictionary: dictionary of fits file
			nr_pixel_beam_size: number of pixel within the beam size, default is 1 pixel

		Returns:
			2D matrix of [nr_pixel_in_outer_boundary x nr_pixel_radial_direction]
			of I, P, tau, and NH
			number of value (# 0) on each column
		"""

		max_grid = 0
		nr_pixel = dictionary['Nx']//2
		for i in range(0, nr_pixel):
			coor_x, coor_y = self.index_radius(dictionary, i, nr_pixel_beam_size)
			if (len(coor_x) > max_grid):
				max_grid = len(coor_x)

		I_matrix = np.zeros([max_grid, nr_pixel])   
		P_matrix = I_matrix.copy()                  		                                           
		tau_matrix = I_matrix.copy()
		NH_matrix = I_matrix.copy()
		nr_point = np.zeros(nr_pixel)

		for i in range(0, nr_pixel):	
			coor_x, coor_y = self.index_radius(dictionary, i,  nr_pixel_beam_size)

			I = dictionary['Stoke_I'][0, coor_x, coor_y]
			P = dictionary['P'][0, coor_x, coor_y]
			tau = dictionary['optical_depth'][0, coor_x, coor_y]
			NH = dictionary['column_density'][0, coor_x, coor_y]*1e-4 # cm-2

			cout = 0
			for j in range(0, len(P)):
				I_matrix[j,i] = I[j]
				P_matrix[j,i] = P[j]				
				tau_matrix[j,i] = tau[j]
				NH_matrix[j,i] = NH[j]
				cout += 1
			nr_point[i] = cout
                    
		return I_matrix, P_matrix, tau_matrix, NH_matrix, nr_point 


	def mean_value(self, I_matrix, P_matrix, tau_matrix, NH_matrix, nr_point):
		"""
		Args: 
			2D matrix of [nr_pixel_in_outer_boundary x nr_pixel_radial_direction] of
			I, P, tau, NH
			nr of point (value # 0) in each column of the matrix

		Returns:
			array of the mean value over the circle at different position from the center 
			to the outer boundary
		"""
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
        

	def polarization_curve(self, dictionary, index_x, index_y):
		"""
		Args:
			dictionary: dictionary of fits file
			index_x: position of pixel in the observed region on x direction
			index_y: position of pixel in the observed region of y direction

		Returns:
			P: array of [nr_wavelength] of P in the observed region
			P_mean: array of length [nr_wavelength] of mean polarization degree in the observed region
		"""
		P = dictionary['P'][:, index_x, index_y]
		# shape of P: [wavelength x len(index_x)]


		nr_wave = dictionary['nr_wave']
		P_mean = np.zeros(nr_wave)
		for i in range(0, nr_wave):
			P_mean[i] = sum(P[i, :])/len(P[i, :])
		return P, P_mean




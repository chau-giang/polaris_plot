import os
import numpy as np
import matplotlib.pyplot as plt
from plot_polaris import plot

fig, ax = plt.subplots(figsize = (6, 4.5))
color_name = np.array(['red', 'orange', 'green', 'deepskyblue', 'blue', 'purple', 'black'])

wavelength = "2mm"
file_name = np.array(["para_grain_" + wavelength + ".fits", "Ncl_10_" + wavelength + ".fits", "Ncl_50_" + wavelength + ".fits", "Ncl_100_" + wavelength + ".fits", "Ncl_1000_" + wavelength + ".fits"])

label_name = np.array(['paramagnetic grains', r'$\sf{\rm N_{\rm cl} = 10}$',  r'$\sf{\rm N_{\rm cl} = 50}$',  r'$\sf{\rm N_{\rm cl} = 100}$',  r'$\sf{\rm N_{\rm cl} = 1000}$',r'$\sf{\rm N_{\rm cl} = 10^{4}}$'])

directory = "/home/chaugiang/Dropbox/POLARIS-/projects/Bok_globule/amax_10um/larmor/zoom_in_4000au/data"

fig, ax = plt.subplots(figsize = (6, 4.5))

for i in range(0, len(file_name)):
    pol_ave = plot.polarized_thermal(os.path.join(directory, file_name[i]), color_name[i], label_name[i])
    pol_ave.plot_P_I()


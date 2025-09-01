#-------------------------------------------------------------------------------
# Librairies
#-------------------------------------------------------------------------------

import shutil, os, pickle, skfmm, math, vtk, imageio, random
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
from PIL import Image
from vtk.util.numpy_support import vtk_to_numpy

# Own
from Read_PF_outputs import *
from Estimate_ElasticParameters import *

#-------------------------------------------------------------------------------

def create_folder(name):
    '''
    Create a new folder. If it already exists, it is erased.
    '''
    if Path(name).exists():
        shutil.rmtree(name)
    os.mkdir(name)

#-------------------------------------------------------------------------------

def index_to_str_3(j):
  '''
  An integer is converted to a float with 3 components
  '''
  if j < 10:
      j_str = '00'+str(j)
  elif 10 <= j and j < 100:
      j_str = '0'+str(j)
  else :
      j_str = str(j)
  return j_str

#-------------------------------------------------------------------------------

def index_to_str_4(j):
  '''
  An integer is converted to a float with 4 components
  '''
  if j < 10:
      j_str = '000'+str(j)
  elif 10 <= j and j < 100:
      j_str = '00'+str(j)
  elif 100 <= j and j < 1000:
      j_str = '0'+str(j)
  else :
      j_str = str(j)
  return j_str

#-------------------------------------------------------------------------------
# Plots
#-------------------------------------------------------------------------------

plot_data_extracted = False
plot_maps_bin = False
plot_maps_pf = False
plot_maps_bin_output = False

#-------------------------------------------------------------------------------
# Preparation
#-------------------------------------------------------------------------------

create_folder('output')
if plot_data_extracted:
    create_folder('output/data_extracted')
if plot_maps_bin:
    create_folder('output/maps_bin')
if plot_maps_pf:
    create_folder('output/maps_pf')
create_folder('data')
create_folder('e')
create_folder('i')
create_folder('vtk')
create_folder('csv')

# dictionnaries
dict_pf = {
    'plot_data_extracted': plot_data_extracted,
    'plot_maps_bin': plot_maps_bin,
    'plot_maps_pf': plot_maps_pf,
    'plot_maps_bin_output': plot_maps_bin_output
}
dict_loading = {}

#-------------------------------------------------------------------------------
# Read data
#-------------------------------------------------------------------------------
print('Extract data from CTscans')

# determine the size (> 150 for REV)
size = 20

# Image are 1341 x 1200 x 1200
# definition of the extraction zone
i_x_min = random.randint(450, 970-size)
i_x_max = i_x_min + size
i_y_min = random.randint(590, 970-size)
i_y_max = i_y_min + size
i_z_min = random.randint(200, 1480-size)
i_z_max = i_z_min + size

# extract data
data_extracted = np.zeros((i_x_max-i_x_min+1, i_y_max-i_y_min+1, i_z_max-i_z_min+1))
# iterate on the slices
if plot_data_extracted:
    images = []
for i_z in range(i_z_min, i_z_max+1):
    # read the image
    with Image.open('Tengattini2023/CGB29AT'+ index_to_str_4(i_z) +'.png') as im:
        # convert PIL into numpy
        data = np.array(im)
        # extract data
        for i_y in range(i_y_min, i_y_max+1):
            for i_x in range(i_x_min, i_x_max+1):
                data_extracted[i_x-i_x_min, i_y-i_y_min, i_z-i_z_min] = data[i_x, i_y].copy() 

        # plot data extracted
        if plot_data_extracted : 
            fig, (ax1) = plt.subplots(nrows=1,ncols=1,figsize=(16,9))
            ax1.imshow(data_extracted[:,:, i_z-i_z_min], extent=(0, i_x_max-i_x_min+1, 0, i_y_max-i_y_min+1))
            ax1.set_xlabel('axis x (pixels)')
            ax1.set_ylabel('axis y (pixels)')
            fig.tight_layout()
            fig.savefig('output/data_extracted/'+index_to_str_3(i_z-i_z_min)+'.png')
            plt.close(fig)

            # save for the gif 
            images.append(imageio.imread('output/data_extracted/'+index_to_str_3(i_z-i_z_min)+'.png'))
if plot_data_extracted:
    # generate gif
    imageio.mimsave('output/data_extracted/ctscans.gif', images)

# compute the histogram of the pixel values
L_p, L_values = np.histogram(data_extracted, bins=np.arange(256))
# pp L_values
L_values_pp = []
for i in range(len(L_values)-1):
    L_values_pp.append(0.5*L_values[i+1] + 0.5*L_values[i])
# define threshold
pore_cement = 50
cement_grain = 200
# plot
fig, (ax1) = plt.subplots(nrows=1,ncols=1,figsize=(16,9))
ax1.plot(L_values_pp, L_p)
ax1.plot([pore_cement, pore_cement], [min(L_p), max(L_p)], color='k', linestyle='dotted')
ax1.plot([cement_grain, cement_grain], [min(L_p), max(L_p)], color='k', linestyle='dotted')
ax1.set_xlabel('pixel value (-)')
ax1.set_ylabel('number (-)')
fig.tight_layout()
fig.savefig('output/histo_pixel_value.png')
plt.close(fig)

# extract grain and cement maps
data_grain = np.zeros((data_extracted.shape[0], data_extracted.shape[1], data_extracted.shape[2]))
data_cement = np.zeros((data_extracted.shape[0], data_extracted.shape[1], data_extracted.shape[2]))
data_matter = np.zeros((data_extracted.shape[0], data_extracted.shape[1], data_extracted.shape[2]))
for i_z in range(data_extracted.shape[2]):
    for i_x in range(data_extracted.shape[0]):
        for i_y in range(data_extracted.shape[1]):
            if cement_grain <= data_extracted[i_x, i_y, i_z]:
                data_grain[i_x, i_y, i_z] = 1
                data_matter[i_x, i_y, i_z] = 1
            if pore_cement <= data_extracted[i_x, i_y, i_z] and data_extracted[i_x, i_y, i_z] < cement_grain:
                data_cement[i_x, i_y, i_z] = 1
                data_matter[i_x, i_y, i_z] = 1
 
    # plot
    if plot_maps_bin:
        fig, (ax1, ax2) = plt.subplots(nrows=1,ncols=2,figsize=(16,9))
        ax1.imshow(data_grain[:, :, i_z], cmap='binary')
        ax1.set_title('grain')
        ax2.imshow(data_cement[:, :, i_z], cmap='binary')
        ax2.set_title('cement')
        fig.tight_layout()
        fig.savefig('output/maps_bin/'+index_to_str_3(i_z)+'.png')
        plt.close(fig)

# save
dict_pf['i_x_min'] = i_x_min
dict_pf['i_x_max'] = i_x_max
dict_pf['i_y_min'] = i_y_min
dict_pf['i_y_max'] = i_y_max
dict_pf['i_z_min'] = i_z_min
dict_pf['i_z_max'] = i_z_max
dict_pf['M_grain_0'] = data_grain
dict_pf['M_cement_0'] = data_cement

#-------------------------------------------------------------------------------
# Compute mesh
#-------------------------------------------------------------------------------

m_size = 1

x_L = np.arange(-m_size*(data_extracted.shape[0]-1)/2, m_size*(data_extracted.shape[0]-1)/2+0.1*m_size, m_size)
y_L = np.arange(-m_size*(data_extracted.shape[1]-1)/2, m_size*(data_extracted.shape[1]-1)/2+0.1*m_size, m_size)
z_L = np.arange(-m_size*(data_extracted.shape[2]-1)/2, m_size*(data_extracted.shape[2]-1)/2+0.1*m_size, m_size)

# save
dict_pf['m_size'] = m_size
dict_pf['x_L'] = x_L
dict_pf['y_L'] = y_L
dict_pf['z_L'] = z_L

#-------------------------------------------------------------------------------
# Compute sdfs
#-------------------------------------------------------------------------------
print('Compute sdfs')

bin_map_matter = -np.ones(data_matter.shape)
# iterate on x
for i_x in range(data_grain.shape[0]):
    # iterate on y
    for i_y in range(data_grain.shape[1]):
        # iterate on z  
        for i_z in range(data_grain.shape[2]):
            if data_matter[i_x, i_y, i_z] == 1:
                bin_map_matter[i_x, i_y, i_z] = 1
                
# compute sdf
sdf_matter = skfmm.distance(bin_map_matter, dx=np.array([m_size, m_size, m_size]))
        
#-------------------------------------------------------------------------------
# Compute phase variables
#-------------------------------------------------------------------------------
print('Compute pfs')

# pf parameters
w_int_pf = m_size*6/1

# init
pf_map_matter = np.zeros(sdf_matter.shape)

# compute the phase field variables
for i_z in range(len(z_L)):
    for i_x in range(len(x_L)):
        for i_y in range(len(y_L)):
            # matter
            if sdf_matter[i_x, i_y, i_z] > w_int_pf/2: # inside the cement
                pf_map_matter[i_x, i_y, i_z] = 1
            elif sdf_matter[i_x, i_y, i_z] < -w_int_pf/2: # outside the cement
                pf_map_matter[i_x, i_y, i_z] = 0
            else : # in the interface
                pf_map_matter[i_x, i_y, i_z] = 0.5*(1+math.cos(math.pi*(-sdf_matter[i_x, i_y, i_z]+w_int_pf/2)/w_int_pf))
                
    # plot
    if plot_maps_pf:
        fig, (ax1) = plt.subplots(nrows=1,ncols=1,figsize=(16,9))
        ax1.imshow(pf_map_matter[:, :, i_z], cmap='binary', vmin=0, vmax=1)
        ax1.set_title('matter = grain+cement')
        fig.tight_layout()
        fig.savefig('output/maps_pf/'+index_to_str_3(i_z)+'.png')
        plt.close(fig)

dict_pf['w_int_pf'] = w_int_pf

#-------------------------------------------------------------------------------
# Write phase variables
#-------------------------------------------------------------------------------

file_to_write_matter = open('data/matter.txt','w')
# x
file_to_write_matter.write('AXIS X\n')
line = ''
for x in x_L:
    line = line + str(x)+ ' '
line = line + '\n'
file_to_write_matter.write(line)
# y
file_to_write_matter.write('AXIS Y\n')
line = ''
for y in y_L:
    line = line + str(y)+ ' '
line = line + '\n'
file_to_write_matter.write(line)
# z
file_to_write_matter.write('AXIS Z\n')
line = ''
for z in z_L:
    line = line + str(z)+ ' '
line = line + '\n'
file_to_write_matter.write(line)
# data
file_to_write_matter.write('DATA\n')
for i_z in range(len(z_L)):
    for i_y in range(len(y_L)):
        for i_x in range(len(x_L)):
            file_to_write_matter.write(str(pf_map_matter[i_x, i_y, i_z])+'\n')
# close
file_to_write_matter.close()

#-------------------------------------------------------------------------------
# Define pf parameter
#-------------------------------------------------------------------------------

# mesh pf
f_mesh_pf = 1

# free energy
W_pf = 1
ed_pf = 0.4*W_pf
kappa_pf = W_pf*w_int_pf*w_int_pf/9.86

# resolution
crit_res_pf = 1e-3
n_ite_pf_max = 50
dt_pf = 0.02
n_proc = 4

# save
dict_pf['f_mesh_pf'] = f_mesh_pf
dict_pf['W_pf'] = W_pf
dict_pf['ed_pf'] = ed_pf
dict_pf['kappa_pf'] = kappa_pf
dict_pf['crit_res_pf'] = crit_res_pf
dict_pf['n_ite_pf_max'] = n_ite_pf_max
dict_pf['dt_pf'] = dt_pf
dict_pf['n_proc'] = n_proc

#-------------------------------------------------------------------------------
# Compute ed
#-------------------------------------------------------------------------------

ed_map = np.zeros(data_grain.shape)
# iterate on x
for i_x in range(data_grain.shape[0]):
    # iterate on y
    for i_y in range(data_grain.shape[1]):
        # iterate on z  
        for i_z in range(data_grain.shape[2]):
            ed_map[i_x, i_y, i_z] = ed_pf

#-------------------------------------------------------------------------------
# Write ed
#-------------------------------------------------------------------------------

file_to_write_ed = open('data/ed.txt','w')
# x
file_to_write_ed.write('AXIS X\n')
line = ''
for x in x_L:
    line = line + str(x)+ ' '
line = line + '\n'
file_to_write_ed.write(line)
# y
file_to_write_ed.write('AXIS Y\n')
line = ''
for y in y_L:
    line = line + str(y)+ ' '
line = line + '\n'
file_to_write_ed.write(line)
# z
file_to_write_ed.write('AXIS Z\n')
line = ''
for z in z_L:
    line = line + str(z)+ ' '
line = line + '\n'
file_to_write_ed.write(line)
# data
file_to_write_ed.write('DATA\n')
for i_z in range(len(z_L)):
    for i_y in range(len(y_L)):
        for i_x in range(len(x_L)):
            file_to_write_ed.write(str(ed_map[i_x, i_y, i_z])+'\n')
# close
file_to_write_ed.close()

#-------------------------------------------------------------------------------
# Write pf input
#-------------------------------------------------------------------------------

file_to_write = open('PF_Dissolution.i','w')
file_to_read = open('PF_Dissolution_Matter_template.i','r')
lines = file_to_read.readlines()
file_to_read.close()

# one pf for the grains and the cement combined
j = 0
for line in lines :
    j = j + 1
    if j == 4:
        line = line[:-1] + ' ' + str(int(len(x_L)*f_mesh_pf)) + '\n'
    if j == 5:
        line = line[:-1] + ' ' + str(int(len(y_L)*f_mesh_pf)) + '\n'
    if j == 6:
        line = line[:-1] + ' ' + str(int(len(z_L)*f_mesh_pf)) + '\n'
    if j == 7:
        line = line[:-1] + ' ' + str(min(x_L)) + '\n'
    if j == 8:
        line = line[:-1] + ' ' + str(max(x_L)) + '\n'
    if j == 9:
        line = line[:-1] + ' ' + str(min(y_L)) + '\n'
    if j == 10:
        line = line[:-1] + ' ' + str(max(y_L)) + '\n'
    if j == 11:
        line = line[:-1] + ' ' + str(min(z_L)) + '\n'
    if j == 12:
        line = line[:-1] + ' ' + str(max(z_L)) + '\n'
    if j == 69:
        line = line[:-1] + " '" + str(1) + ' ' + str(kappa_pf) + "'\n"
    if j == 86:
        line = line[:-1] + " '" + str(W_pf) + "'\n"
    if j == 124 or j == 125 or j == 128 or j == 129:
        line = line[:-1] + " " + str(crit_res_pf) + "\n"
    if j == 132:
        line = line[:-1] + " " + str(10*n_ite_pf_max) + "\n"
    if j == 133:
        line = line[:-1] + " " + str(dt_pf*n_ite_pf_max) + "\n"
    if j == 137:
        line = line[:-1] + " " + str(dt_pf) + "\n"
    file_to_write.write(line)

file_to_write.close()

#-------------------------------------------------------------------------------
# Run pf and sort outputs
#-------------------------------------------------------------------------------

# MOOSE simulation
os.system('mpiexec -n '+str(n_proc)+' ~/projects/moose/modules/phase_field/phase_field-opt -i PF_Dissolution.i')

# sort .i and .e files
os.rename('PF_Dissolution.i','i/PF_Dissolution.i')
os.rename('PF_Dissolution_out.e','e/PF_Dissolution_out.e')

# sort .vtk
j = 0
j_str = index_to_str_3(j)
filepath = Path('PF_Dissolution_other_'+j_str+'.pvtu')
while filepath.exists():
    for i_proc in range(n_proc):
        os.rename('PF_Dissolution_other_'+j_str+'_'+str(i_proc)+'.vtu','vtk/PF_Dissolution_other_'+j_str+'_'+str(i_proc)+'.vtu')
    os.rename('PF_Dissolution_other_'+j_str+'.pvtu','vtk/PF_Dissolution_other_'+j_str+'.pvtu')
    j = j + 1
    j_str = index_to_str_3(j)
    filepath = Path('PF_Dissolution_other_'+j_str+'.pvtu')
# save last indice
last_j = j-1
last_j_str = index_to_str_3(j-1)

# sort .csv
os.rename('PF_Dissolution_csv.csv','csv/PF_Dissolution_csv.csv')

# delete files
shutil.rmtree('e')
shutil.rmtree('i')
shutil.rmtree('data')
shutil.rmtree('csv')

# user print
print('\nEnd of the simulation')

#-------------------------------------------------------------------------------
# Read csv output 
#-------------------------------------------------------------------------------

#Read_PF_csv(crit_res_pf)

#-------------------------------------------------------------------------------
# Define parameters loadings
#-------------------------------------------------------------------------------

# grain
young_grain = 1
poisson_grain = 0.3

# cement
young_cement = 1
poisson_cement = 0.3

# pore
young_pore = 0.01
poisson_pore = 0.3

# resolution
crit_res_fem = 1e-3
dt_fem = 0.05
n_proc = n_proc

# plan the loadings
# compression, shearing, triaxial, isotropic
loading = ['triaxial', 'isotropic']

# define compression
if 'compression' in loading:
    dict_loading['compression_strain'] = -0.1
# define shearing
if 'shearing' in loading:
    pass
# define triaxial
if 'triaxial' in loading:
    dict_loading['triaxial_strain'] = -0.1
# define compression
if 'isotropic' in loading:
    dict_loading['isotropic_strain'] = -0.1/3

# save
dict_loading['young_grain'] = young_grain
dict_loading['poisson_grain'] = poisson_grain
dict_loading['young_cement'] = young_cement
dict_loading['poisson_cement'] = poisson_cement
dict_loading['young_pore'] = young_pore
dict_loading['poisson_pore'] = poisson_pore
dict_loading['crit_res_fem'] = crit_res_fem
dict_loading['dt_fem'] = dt_fem
dict_loading['n_proc'] = n_proc
dict_loading['loading'] = loading

#-------------------------------------------------------------------------------
# Estimate the evolution of the elastic parameters
#-------------------------------------------------------------------------------

# initialization
L_young = []
L_poisson = []
L_shear = []
L_bulk = []

# iteration on time
for iteration in range(last_j+1):
    print(iteration+1,'/', last_j+1)

    # Read vtk output 
    M_grain, M_cement = Read_PF_vtk(index_to_str_3(iteration), dict_loading, dict_pf)

    # Prepare the MOOSE simulation
    create_folder('data')
    create_folder('i')
    create_folder('csv')
    create_folder('e')
    if plot_maps_bin_output:
        create_folder('output/microstructure_'+index_to_str_3(iteration))

    # Generate the png file for Moose FEM simulation
    Generate_png(M_grain, M_cement, plot_maps_bin_output, index_to_str_3(iteration))

    if 'compression' in dict_loading['loading']:
        # Write the compression .i file
        Write_compression_i(x_L, y_L, z_L, dict_loading['compression_strain'], young_pore, poisson_pore, young_grain, poisson_grain, young_cement, poisson_cement, crit_res_fem, dt_fem)
        # Run fem MOOSE simulation
        os.system('mpiexec -n '+str(n_proc)+' ~/projects/moose/modules/solid_mechanics/solid_mechanics-opt -i FEM_Loading_Compression.i')
        
        # Read the csv output
        L_strain_xx, L_strain_xy, L_strain_xz, L_strain_yy, L_strain_yz, L_strain_zz,\
        L_stress_xx, L_stress_xy, L_stress_xz, L_stress_yy, L_stress_yz, L_stress_zz = Read_FEM_csv('FEM_Loading_Compression_csv.csv', M_grain, M_cement)
        # compute the strain
        L_strain = []
        for i_strain in range(len(L_stress_zz)):
            L_strain.append(i_strain/(len(L_stress_zz)-1)*dict_loading['compression_strain']*(max(z_L)-min(z_L)))
        # sort .i, .csv, .e files
        os.rename('FEM_Loading_Compression.i','i/FEM_Loading_Compression.i')
        os.rename('FEM_Loading_Compression_csv.csv','csv/FEM_Loading_Compression_csv.csv')
        os.rename('FEM_Loading_Compression_out.e','e/FEM_Loading_Compression_out.e')
        # interpolate elastic parameters
        YoungModulusSample = Interpolate_compression_prop(L_strain, L_stress_zz)
        # save
        if not 'triaxial' in dict_loading['loading']: # triaxial test has the priority
            L_young.append(YoungModulusSample)

    if 'triaxial' in dict_loading['loading']:
        # Write the triaxial .i file
        Write_triaxial_i(x_L, y_L, z_L, dict_loading['triaxial_strain'], young_pore, poisson_pore, young_grain, poisson_grain, young_cement, poisson_cement, crit_res_fem, dt_fem)
        # Run fem MOOSE simulation
        os.system('mpiexec -n '+str(n_proc)+' ~/projects/moose/modules/solid_mechanics/solid_mechanics-opt -i FEM_Loading_Triaxial.i')
        
        # Read the csv output
        L_disp_x, L_disp_y, L_disp_z,\
        L_strain_xx, L_strain_xy, L_strain_xz, L_strain_yy, L_strain_yz, L_strain_zz,\
        L_stress_xx, L_stress_xy, L_stress_xz, L_stress_yy, L_stress_yz, L_stress_zz = Read_FEM_csv('FEM_Loading_Triaxial_csv.csv', M_grain, M_cement)
        # pp data
        L_strain_x = []
        L_strain_y = []
        L_strain_z = []
        for i_disp in range(len(L_disp_z)):
            L_strain_x.append(L_disp_x[i_disp]/(max(x_L)-min(x_L)))
            L_strain_y.append(L_disp_y[i_disp]/(max(y_L)-min(y_L)))
            L_strain_z.append(L_disp_z[i_disp]/(max(z_L)-min(z_L)))
        # sort .i, .csv, .e files
        os.rename('FEM_Loading_Triaxial.i','i/FEM_Loading_Triaxial.i')
        os.rename('FEM_Loading_Triaxial_csv.csv','csv/FEM_Loading_Triaxial_csv.csv')
        os.rename('FEM_Loading_Triaxial_out.e','e/FEM_Loading_Triaxial_out.e')
        # interpolate elastic parameters
        YoungModulusSample, ShearModulusSample, PoissonRatioSample = Interpolate_triaxial_props(L_strain_x, L_strain_y, L_strain_z, L_stress_zz)
        # save
        L_young.append(YoungModulusSample)
        L_poisson.append(PoissonRatioSample)

    # TO DO same for shearing and isotropic
    if 'shearing' in dict_loading['loading']:
        pass

    if 'isotropic' in dict_loading['loading']:
        # Write the isotropic .i file
        Write_isotropic_i(x_L, y_L, z_L, dict_loading['isotropic_strain'], young_pore, poisson_pore, young_grain, poisson_grain, young_cement, poisson_cement, crit_res_fem, dt_fem)
        # Run fem MOOSE simulation
        os.system('mpiexec -n '+str(n_proc)+' ~/projects/moose/modules/solid_mechanics/solid_mechanics-opt -i FEM_Loading_Isotropic.i')
        
        # Read the csv output
        L_disp_x, L_disp_y, L_disp_z,\
        L_strain_xx, L_strain_xy, L_strain_xz, L_strain_yy, L_strain_yz, L_strain_zz,\
        L_stress_xx, L_stress_xy, L_stress_xz, L_stress_yy, L_stress_yz, L_stress_zz = Read_FEM_csv('FEM_Loading_Isotropic_csv.csv', M_grain, M_cement)
        # compute the stress
        # pp data
        L_vol_strain = []
        for i_disp in range(len(L_disp_z)):
            strain_x_i = L_disp_x[i_disp]/(max(x_L)-min(x_L))
            strain_y_i = L_disp_y[i_disp]/(max(y_L)-min(y_L))
            strain_z_i = L_disp_z[i_disp]/(max(z_L)-min(z_L))
            L_vol_strain.append(strain_x_i + strain_y_i + strain_z_i)
        # sort .i, .csv, .e files
        os.rename('FEM_Loading_Isotropic.i','i/FEM_Loading_Isotropic.i')
        os.rename('FEM_Loading_Isotropic_csv.csv','csv/FEM_Loading_Isotropic_csv.csv')
        os.rename('FEM_Loading_Isotropic_out.e','e/FEM_Loading_Isotropic_out.e')
        # interpolate elastic parameters
        BulkModulusSample = Interpolate_isotropic_prop(L_stress_xx, L_stress_yy, L_stress_zz, L_vol_strain)
        # save
        L_bulk.append(BulkModulusSample)

# save
dict_loading['L_young'] = L_young
dict_loading['L_poisson'] = L_poisson
dict_loading['L_shear'] = L_shear
dict_loading['L_bulk'] = L_bulk

#-------------------------------------------------------------------------------
# Plot evolution elasic parameters
#-------------------------------------------------------------------------------

fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(nrows=2,ncols=2,figsize=(16,9))
ax1.plot(L_young)
ax2.plot(L_poisson)
ax3.plot(L_shear)
ax4.plot(L_bulk)
ax1.set_title('Young modulus (-)')
ax2.set_title('Poisson ratio (-)')
ax3.set_title('Shear modulus(-)')
ax4.set_title('Bulk modulus (-)')
fig.tight_layout()
fig.savefig('output/evol_time_elastic_parameters.png')
plt.close(fig)

#-------------------------------------------------------------------------------
# Close
#-------------------------------------------------------------------------------

# delete files
shutil.rmtree('e')
shutil.rmtree('i')
shutil.rmtree('data')
shutil.rmtree('csv')

# user print
print('\nEnd of the simulation')

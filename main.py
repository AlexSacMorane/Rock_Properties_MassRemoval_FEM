#-------------------------------------------------------------------------------
# Librairies
#-------------------------------------------------------------------------------

import shutil, os, pickle, skfmm, math, vtk
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
# Preparation
#-------------------------------------------------------------------------------

create_folder('output')
create_folder('data')
create_folder('e')
create_folder('i')
create_folder('vtk')
create_folder('csv')

#-------------------------------------------------------------------------------
# Read data
#-------------------------------------------------------------------------------

# Image are 1341 x 1200 x 1200
# definition of the extraction zone
i_x_min = 500
i_x_max = 600
i_y_min = 500
i_y_max = 600
i_z_min = 550
i_z_max = 600

# extract data
data_extracted = np.zeros((i_x_max-i_x_min+1, i_y_max-i_y_min+1, i_z_max-i_z_min+1))
# iterate on the slices
for i_z in range(i_z_min, i_z_max):
    # read the image
    with Image.open('Tengattini2023/CGB29AT'+ index_to_str_4(i_z_str) +'.png') as im:
        # convert PIL into numpy
        data = np.array(im)
        # extract data
        for i_y in range(i_y_min, i_y_max+1):
            for i_x in range(i_x_min, i_x_max+1):
                data_extracted[i_x-i_x_min, i_y-i_y_min, i_z] = data[i_x, i_y].copy() 

        # plot data extracted
        fig, (ax1) = plt.subplots(nrows=1,ncols=1,figsize=(16,9))
        ax1.imshow(data_extracted, extent=(0, i_y_max-i_y_min+1, 0, i_z_max-i_z_min+1))
        ax1.set_xlabel('axis y (pixels)')
        ax1.set_ylabel('axis z (pixels)')
        fig.tight_layout()
        fig.savefig('output/data_extracted_'+index_to_str_3(i_z-i_z_min)+'.png')
        plt.close(fig)

raise ValueError('stop')

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
data_grain = np.zeros((data_extracted.shape[0], data_extracted.shape[1]))
data_cement = np.zeros((data_extracted.shape[0], data_extracted.shape[1]))
data_matter = np.zeros((data_extracted.shape[0], data_extracted.shape[1]))
for i_y in range(data_extracted.shape[0]):
    for i_z in range(data_extracted.shape[1]):
        if cement_grain <= data_extracted[i_y, i_z]:
            data_grain[i_y, i_z] = 1
            data_matter[i_y, i_z] = 1
        if pore_cement <= data_extracted[i_y, i_z] and data_extracted[i_y, i_z] < cement_grain:
            data_cement[i_y, i_z] = 1
            data_matter[i_y, i_z] = 1
 
# plot
fig, (ax1, ax2) = plt.subplots(nrows=1,ncols=2,figsize=(16,9))
ax1.imshow(data_grain, cmap='binary')
ax1.set_title('grain')
ax2.imshow(data_cement, cmap='binary')
ax2.set_title('cement')
fig.tight_layout()
fig.savefig('output/maps_bin_cement_grain.png')
plt.close(fig)

#-------------------------------------------------------------------------------
# Compute mesh
#-------------------------------------------------------------------------------

m_size = 1

#x_L = np.arange(-m_size*(map_data.shape[0]-1)/2, m_size*(map_data.shape[0]-1)/2+0.1*m_size, m_size)
y_L = np.arange(-m_size*(data_extracted.shape[0]-1)/2, m_size*(data_extracted.shape[0]-1)/2+0.1*m_size, m_size)
z_L = np.arange(-m_size*(data_extracted.shape[1]-1)/2, m_size*(data_extracted.shape[1]-1)/2+0.1*m_size, m_size)

#-------------------------------------------------------------------------------
# Compute sdfs
#-------------------------------------------------------------------------------

bin_map_grain = -np.ones(data_grain.shape)
bin_map_cement = -np.ones(data_cement.shape)
bin_map_matter = -np.ones(data_matter.shape)
# iterate on y
for i_y in range(data_grain.shape[0]):
    # iterate on z  
    for i_z in range(data_grain.shape[1]):
        if data_grain[i_y, i_z] == 1:
            bin_map_grain[i_y, i_z] = 1
        if data_cement[i_y, i_z] == 1:
            bin_map_cement[i_y, i_z] = 1
        if data_matter[i_y, i_z] == 1:
            bin_map_matter[i_y, i_z] = 1
                
# compute sdf
sdf_grain = skfmm.distance(bin_map_grain, dx=np.array([m_size, m_size]))
sdf_cement = skfmm.distance(bin_map_cement, dx=np.array([m_size, m_size]))
sdf_matter = skfmm.distance(bin_map_matter, dx=np.array([m_size, m_size]))
        
# plot
fig, (ax1, ax2, ax3) = plt.subplots(nrows=1,ncols=3,figsize=(16,9))
ax1.imshow(sdf_grain, extent=(min(y_L), max(y_L), min(z_L), max(z_L)), cmap='bwr')
ax1.set_title('grain')
ax2.imshow(sdf_cement, extent=(min(y_L), max(y_L), min(z_L), max(z_L)), cmap='bwr')
ax2.set_title('cement')
ax3.imshow(sdf_matter, extent=(min(y_L), max(y_L), min(z_L), max(z_L)), cmap='bwr')
ax3.set_title('matter')
fig.tight_layout()
fig.savefig('output/maps_sdf_cement_grain_matter.png')
plt.close(fig)

#-------------------------------------------------------------------------------
# Compute phase variables
#-------------------------------------------------------------------------------

# pf parameters
w_int_pf = m_size*6/1

# init
pf_map_grain = np.zeros(sdf_grain.shape)
pf_map_cement = np.zeros(sdf_cement.shape)
pf_map_matter = np.zeros(sdf_matter.shape)

# compute the phase field variables
for i_y in range(len(y_L)):
    for i_z in range(len(z_L)):
        # grain
        if sdf_grain[i_y, i_z] > w_int_pf/2: # inside the grain
            pf_map_grain[i_y, i_z] = 1
        elif sdf_grain[i_y, i_z] < -w_int_pf/2: # outside the grain
            pf_map_grain[i_y, i_z] = 0
        else : # in the interface
            pf_map_grain[i_y, i_z] = 0.5*(1+math.cos(math.pi*(-sdf_grain[i_y, i_z]+w_int_pf/2)/w_int_pf))
        # cement
        if sdf_cement[i_y, i_z] > w_int_pf/2: # inside the cement
            pf_map_cement[i_y, i_z] = 1
        elif sdf_cement[i_y, i_z] < -w_int_pf/2: # outside the cement
            pf_map_cement[i_y, i_z] = 0
        else : # in the interface
            pf_map_cement[i_y, i_z] = 0.5*(1+math.cos(math.pi*(-sdf_cement[i_y, i_z]+w_int_pf/2)/w_int_pf))
        # matter
        if sdf_matter[i_y, i_z] > w_int_pf/2: # inside the cement
            pf_map_matter[i_y, i_z] = 1
        elif sdf_matter[i_y, i_z] < -w_int_pf/2: # outside the cement
            pf_map_matter[i_y, i_z] = 0
        else : # in the interface
            pf_map_matter[i_y, i_z] = 0.5*(1+math.cos(math.pi*(-sdf_matter[i_y, i_z]+w_int_pf/2)/w_int_pf))

# plot
fig, (ax1, ax2, ax3) = plt.subplots(nrows=1,ncols=3,figsize=(16,9))
ax1.imshow(pf_map_grain, cmap='bwr')
ax1.set_title('grain')
ax2.imshow(pf_map_cement, cmap='bwr')
ax2.set_title('cement')
ax3.imshow(pf_map_matter, cmap='bwr')
ax3.set_title('matter')
fig.tight_layout()
fig.savefig('output/maps_pf_cement_grain.png')
plt.close(fig)

#-------------------------------------------------------------------------------
# Write phase variables
#-------------------------------------------------------------------------------

file_to_write_grain = open('data/grain.txt','w')
file_to_write_cement = open('data/cement.txt','w')
file_to_write_matter = open('data/matter.txt','w')
# x
file_to_write_grain.write('AXIS X\n')
file_to_write_cement.write('AXIS X\n')
file_to_write_matter.write('AXIS X\n')
line = ''
for y in y_L:
    line = line + str(y)+ ' '
line = line + '\n'
file_to_write_grain.write(line)
file_to_write_cement.write(line)
file_to_write_matter.write(line)
# y
file_to_write_grain.write('AXIS Y\n')
file_to_write_cement.write('AXIS Y\n')
file_to_write_matter.write('AXIS Y\n')
line = ''
for z in z_L:
    line = line + str(z)+ ' '
line = line + '\n'
file_to_write_grain.write(line)
file_to_write_cement.write(line)
file_to_write_matter.write(line)
# data
file_to_write_grain.write('DATA\n')
file_to_write_cement.write('DATA\n')
file_to_write_matter.write('DATA\n')
for l in range(len(z_L)):
    for c in range(len(y_L)):
        file_to_write_grain.write(str(pf_map_grain[-1-l][c])+'\n')
        file_to_write_cement.write(str(pf_map_cement[-1-l][c])+'\n')
        file_to_write_matter.write(str(pf_map_matter[-1-l][c])+'\n')
# close
file_to_write_grain.close()
file_to_write_cement.close()
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

#-------------------------------------------------------------------------------
# Compute ed
#-------------------------------------------------------------------------------

ed_map = np.zeros(data_grain.shape)
# iterate on y
for i_y in range(data_grain.shape[0]):
    # iterate on z  
    for i_z in range(data_grain.shape[1]):
        if data_grain[i_y, i_z] == 1 or data_cement[i_y, i_z] == 1:
            ed_map[i_y, i_z] = ed_pf

#-------------------------------------------------------------------------------
# Write ed
#-------------------------------------------------------------------------------

file_to_write_ed = open('data/ed.txt','w')
# x
file_to_write_ed.write('AXIS X\n')
line = ''
for y in y_L:
    line = line + str(y)+ ' '
line = line + '\n'
file_to_write_ed.write(line)
# y
file_to_write_ed.write('AXIS Y\n')
line = ''
for z in z_L:
    line = line + str(z)+ ' '
line = line + '\n'
file_to_write_ed.write(line)
# data
file_to_write_ed.write('DATA\n')
for l in range(len(z_L)):
    for c in range(len(y_L)):
        file_to_write_ed.write(str(ed_map[-1-l][c])+'\n')
# close
file_to_write_ed.close()

#-------------------------------------------------------------------------------
# Write pf input
#-------------------------------------------------------------------------------

# PF_Dissolution_GrainCement_template.i or PF_Dissolution_Matter_template.i
name_template = 'PF_Dissolution_Matter_template.i'

file_to_write = open('PF_Dissolution.i','w')
file_to_read = open(name_template,'r')
lines = file_to_read.readlines()
file_to_read.close()

# one pf for the grains and one pf for the cement
if name_template == 'PF_Dissolution_GrainCement_template.i':
    j = 0
    for line in lines :
        j = j + 1
        if j == 4:
            line = line[:-1] + ' ' + str(int(len(y_L)*f_mesh_pf)) + '\n'
        if j == 5:
            line = line[:-1] + ' ' + str(int(len(z_L)*f_mesh_pf)) + '\n'
        if j == 7:
            line = line[:-1] + ' ' + str(min(y_L)) + '\n'
        if j == 8:
            line = line[:-1] + ' ' + str(max(y_L)) + '\n'
        if j == 9:
            line = line[:-1] + ' ' + str(min(z_L)) + '\n'
        if j == 10:
            line = line[:-1] + ' ' + str(max(z_L)) + '\n'
        if j == 93:
            line = line[:-1] + " '" + str(1) + ' ' + str(kappa_pf) + ' '\
                                    + str(1) + ' ' + str(kappa_pf) + "'\n"
        if j == 110:
            line = line[:-1] + " '" + str(W_pf) + "'\n"
        if j == 123:
            line = line[:-1] + " '" + str(W_pf) + "'\n"
        if j == 165 or j == 166 or j == 169 or j == 170:
            line = line[:-1] + " " + str(crit_res_pf) + "\n"
        if j == 173:
            line = line[:-1] + " " + str(10*n_ite_pf_max) + "\n"
        if j == 174:
            line = line[:-1] + " " + str(dt_pf*n_ite_pf_max) + "\n"
        if j == 178:
            line = line[:-1] + " " + str(dt_pf) + "\n"
        file_to_write.write(line)

# one pf for the grains and the cement combined
if name_template == 'PF_Dissolution_Matter_template.i':
    j = 0
    for line in lines :
        j = j + 1
        if j == 4:
            line = line[:-1] + ' ' + str(int(len(y_L)*f_mesh_pf)) + '\n'
        if j == 5:
            line = line[:-1] + ' ' + str(int(len(z_L)*f_mesh_pf)) + '\n'
        if j == 7:
            line = line[:-1] + ' ' + str(min(y_L)) + '\n'
        if j == 8:
            line = line[:-1] + ' ' + str(max(y_L)) + '\n'
        if j == 9:
            line = line[:-1] + ' ' + str(min(z_L)) + '\n'
        if j == 10:
            line = line[:-1] + ' ' + str(max(z_L)) + '\n'
        if j == 61:
            line = line[:-1] + " '" + str(1) + ' ' + str(kappa_pf) + "'\n"
        if j == 78:
            line = line[:-1] + " '" + str(W_pf) + "'\n"
        if j == 116 or j == 117 or j == 120 or j == 121:
            line = line[:-1] + " " + str(crit_res_pf) + "\n"
        if j == 124:
            line = line[:-1] + " " + str(10*n_ite_pf_max) + "\n"
        if j == 125:
            line = line[:-1] + " " + str(dt_pf*n_ite_pf_max) + "\n"
        if j == 129:
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
j_str = index_to_str(j)
filepath = Path('PF_Dissolution_other_'+j_str+'.pvtu')
while filepath.exists():
    for i_proc in range(n_proc):
        os.rename('PF_Dissolution_other_'+j_str+'_'+str(i_proc)+'.vtu','vtk/PF_Dissolution_other_'+j_str+'_'+str(i_proc)+'.vtu')
    os.rename('PF_Dissolution_other_'+j_str+'.pvtu','vtk/PF_Dissolution_other_'+j_str+'.pvtu')
    j = j + 1
    j_str = index_to_str(j)
    filepath = Path('PF_Dissolution_other_'+j_str+'.pvtu')
# save last indice
last_j = j-1
last_j_str = index_to_str(j-1)

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

#Read_PF_csv(name_template, crit_res_pf)

#-------------------------------------------------------------------------------
# Define parameters loadings
#-------------------------------------------------------------------------------

# grain
young_grain = 1
poisson_grain = 0.3

# cement
young_cement = young_grain/5
poisson_cement = 0.3

# pore
young_pore = young_cement/100
poisson_pore = 0.3

# resolution
crit_res_fem = 1e-3
dt_fem = 0.05
n_proc = 4

#-------------------------------------------------------------------------------
# Estimate the evolution of the elastic parameters
#-------------------------------------------------------------------------------

# preparation
create_folder('output/maps_bin_cement_grain_output')

# initialization
L_young = []
L_poisson = []
L_shear = []
L_bulk = []

# check if the map is known (TO DO: create a database)
L_L_i_XYZ_not_used = []
L_XYZ_used = None

# iteration on time
for iteration in range(last_j+1):
    print(iteration+1,'/', last_j)

    # Read vtk output 
    L_L_i_XYZ_not_used, L_XYZ_used, M_grain, M_cement = Read_PF_vtk(name_template, index_to_str(iteration), L_L_i_XYZ_not_used, L_XYZ_used, n_proc, pf_map_matter, y_L, z_L, data_grain, data_cement)

    # Prepare the MOOSE simulation
    create_folder('data')
    create_folder('i')
    create_folder('csv')
    create_folder('e')

    # Generate the png file for Moose FEM simulation
    Generate_png(M_grain, M_cement)

    # Write the compression .i file
    Write_compression_i(y_L, z_L, young_pore, poisson_pore, young_grain, poisson_grain, young_cement, poisson_cement, crit_res_fem, dt_fem)
    # Run fem MOOSE simulation
    os.system('mpiexec -n '+str(n_proc)+' ~/projects/moose/modules/tensor_mechanics/tensor_mechanics-opt -i FEM_Loading_Compression.i')
    # Read the csv output
    L_strain, L_stress_xx, L_stress_xy, L_stress_yy = Read_FEM_csv('FEM_Loading_Compression_csv.csv', M_grain, M_cement, z_L)
    # sort .i, .csv, .e files
    os.rename('FEM_Loading_Compression.i','i/FEM_Loading_Compression.i')
    os.rename('FEM_Loading_Compression_csv.csv','csv/FEM_Loading_Compression_csv.csv')
    os.rename('FEM_Loading_Compression_out.e','e/FEM_Loading_Compression_out.e')
    # interpolate elastic parameters
    YoungModulusSample, PoissonRatioSample = Interpolate_compression_props(L_strain, L_stress_xx, L_stress_yy)
    # save
    L_young.append(YoungModulusSample)
    L_poisson.append(PoissonRatioSample)

    # TO DO same for shearing and isotropic

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

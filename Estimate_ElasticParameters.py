#-------------------------------------------------------------------------------
# Librairies
#-------------------------------------------------------------------------------

import shutil, os, pickle, skfmm, math, vtk
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
from PIL import Image
from vtk.util.numpy_support import vtk_to_numpy

#-------------------------------------------------------------------------------
# Convert a index to str
#-------------------------------------------------------------------------------

def index_to_str(j):
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
# Generate png for the loading of the microstructure
#-------------------------------------------------------------------------------

def Generate_png(M_grain, M_cement, plot_maps_bin_output, iteration_str):
    '''
    Generate a png file (input of the FEM Moose simulation) from the numpy arrays.
    '''
    # initialize the png
    data_png = np.array(np.zeros((M_grain.shape[0], M_grain.shape[1], M_grain.shape[2], 3)))
    if plot_maps_bin_output:
        M_matter = np.array(np.zeros(M_grain.shape))

    # iterate on the mesh
    for i_x in range(M_grain.shape[0]):
        for i_y in range(M_grain.shape[1]):
            for i_z in range(M_grain.shape[2]):
                # create image
                if M_grain[i_x, i_y, i_z] == 0 and M_cement[i_x, i_y, i_z] == 0 :
                    # pore
                    data_png[i_x, i_y, i_z, :] = [0, 0, 0]
                elif M_grain[i_x, i_y, i_z] == 1:
                    # grain
                    data_png[i_x, i_y, i_z, :] = [1/255, 125/255, 125/255]
                    if plot_maps_bin_output:
                        M_matter[i_x, i_y, i_z] = 1
                else:
                    # cement
                    data_png[i_x, i_y, i_z, :] = [2/255, 250/255, 250/255]
                    if plot_maps_bin_output:
                        M_matter[i_x, i_y, i_z] = 0.5
    # iterate on z
    for i_z in range(M_grain.shape[2]):
        # generate the .png file
        plt.imsave('data/microstructure'+index_to_str(i_z+1)+'.png', data_png[:, :, i_z, :])
        # save the .png file
        if plot_maps_bin_output:
            fig, (ax1) = plt.subplots(nrows=1,ncols=1,figsize=(16,9))
            ax1.imshow(M_matter[:,:,i_z], extent=(0, M_matter.shape[0], 0, M_matter.shape[1]))
            ax1.set_xlabel('axis x (pixels)')
            ax1.set_ylabel('axis y (pixels)')
            fig.tight_layout()
            fig.savefig('output/microstructure_'+iteration_str+'/'+index_to_str(i_z)+'.png')
            plt.close(fig)
            
#-------------------------------------------------------------------------------
# Write FEM inputs
#-------------------------------------------------------------------------------

def Write_compression_i(x_L, y_L, z_L, compression_strain, young_pore, poisson_pore, young_grain, poisson_grain, young_cement, poisson_cement, crit_res_fem, dt_fem):
    '''
    Generate from a template the input file for Moose simulation.

    The sample is under compression.
    '''
    file_to_write = open('FEM_Loading_Compression.i','w')
    file_to_read = open('FEM_Loading_Compression_template.i','r')
    lines = file_to_read.readlines()
    file_to_read.close()
    j = 0
    for line in lines :
        j = j + 1
        if j == 5:
            line = line[:-1] + ' ' + str(int(len(x_L))) + '\n'
        if j == 6:
            line = line[:-1] + ' ' + str(int(len(y_L))) + '\n'
        if j == 7:
            line = line[:-1] + ' ' + str(int(len(z_L))) + '\n'
        if j == 8:
            line = line[:-1] + ' ' + str(min(x_L)) + '\n'
        if j == 9:
            line = line[:-1] + ' ' + str(max(x_L)) + '\n'
        if j == 10:
            line = line[:-1] + ' ' + str(min(y_L)) + '\n'
        if j == 11:
            line = line[:-1] + ' ' + str(max(y_L)) + '\n'
        if j == 12:
            line = line[:-1] + ' ' + str(min(y_L)) + '\n'
        if j == 13:
            line = line[:-1] + ' ' + str(max(y_L)) + '\n'
        if j == 48:
            line = line[:-1] + ' ' + str(compression_strain*(max(z_L)-min(z_L))) + '*t\n'
        if j == 94:
            line = line[:-1] + ' ' + str(young_pore) + '\n'
        if j == 95:
            line = line[:-1] + ' ' + str(poisson_pore) + '\n'
        if j == 105:
            line = line[:-1] + ' ' + str(young_grain) + '\n'
        if j == 106:
            line = line[:-1] + ' ' + str(poisson_grain) + '\n'
        if j == 116:
            line = line[:-1] + ' ' + str(young_cement) + '\n'
        if j == 117:
            line = line[:-1] + ' ' + str(poisson_cement) + '\n'
        if j == 137 or j == 139 or j == 140:
            line = line[:-1] + ' ' + str(crit_res_fem) + '\n'
        if j == 146:
            line = line[:-1] + ' ' + str(dt_fem) + '\n'
        file_to_write.write(line)
    file_to_write.close()

#-------------------------------------------------------------------------------

def Write_triaxial_i(x_L, y_L, z_L, triaxial_stress, young_pore, poisson_pore, young_grain, poisson_grain, young_cement, poisson_cement, crit_res_fem, dt_fem):
    '''
    Generate from a template the input file for Moose simulation.

    The sample is under triaxial condition.
    '''
    file_to_write = open('FEM_Loading_Triaxial.i','w')
    file_to_read = open('FEM_Loading_Triaxial_template.i','r')
    lines = file_to_read.readlines()
    file_to_read.close()
    j = 0
    for line in lines :
        j = j + 1
        if j == 5:
            line = line[:-1] + ' ' + str(int(len(x_L))) + '\n'
        if j == 6:
            line = line[:-1] + ' ' + str(int(len(y_L))) + '\n'
        if j == 7:
            line = line[:-1] + ' ' + str(int(len(z_L))) + '\n'
        if j == 8:
            line = line[:-1] + ' ' + str(min(x_L)) + '\n'
        if j == 9:
            line = line[:-1] + ' ' + str(max(x_L)) + '\n'
        if j == 10:
            line = line[:-1] + ' ' + str(min(y_L)) + '\n'
        if j == 11:
            line = line[:-1] + ' ' + str(max(y_L)) + '\n'
        if j == 12:
            line = line[:-1] + ' ' + str(min(y_L)) + '\n'
        if j == 13:
            line = line[:-1] + ' ' + str(max(y_L)) + '\n'
        if j == 48:
            line = line[:-1] + " '" + str(triaxial_stress) + "*t'\n"
        if j == 68:
            line = line[:-1] + ' ' + str(young_pore) + '\n'
        if j == 69:
            line = line[:-1] + ' ' + str(poisson_pore) + '\n'
        if j == 79:
            line = line[:-1] + ' ' + str(young_grain) + '\n'
        if j == 80:
            line = line[:-1] + ' ' + str(poisson_grain) + '\n'
        if j == 90:
            line = line[:-1] + ' ' + str(young_cement) + '\n'
        if j == 91:
            line = line[:-1] + ' ' + str(poisson_cement) + '\n'
        if j == 111 or j == 113 or j == 114:
            line = line[:-1] + ' ' + str(crit_res_fem) + '\n'
        if j == 120:
            line = line[:-1] + ' ' + str(dt_fem) + '\n'
        file_to_write.write(line)
    file_to_write.close()

#-------------------------------------------------------------------------------

def Write_shearing_i():
    '''
    Generate from a template the input file for Moose simulation.

    The sample is under shearing.
    '''
    # TO DO
    pass

#-------------------------------------------------------------------------------

def Write_isotropic_i(x_L, y_L, z_L, isotropic_stress, young_pore, poisson_pore, young_grain, poisson_grain, young_cement, poisson_cement, crit_res_fem, dt_fem):
    '''
    Generate from a template the input file for Moose simulation.

    The sample is under isotropic loading.
    '''
    file_to_write = open('FEM_Loading_Isotropic.i','w')
    file_to_read = open('FEM_Loading_Isotropic_template.i','r')
    lines = file_to_read.readlines()
    file_to_read.close()
    j = 0
    for line in lines :
        j = j + 1
        if j == 5:
            line = line[:-1] + ' ' + str(int(len(x_L))) + '\n'
        if j == 6:
            line = line[:-1] + ' ' + str(int(len(y_L))) + '\n'
        if j == 7:
            line = line[:-1] + ' ' + str(int(len(z_L))) + '\n'
        if j == 8:
            line = line[:-1] + ' ' + str(min(x_L)) + '\n'
        if j == 9:
            line = line[:-1] + ' ' + str(max(x_L)) + '\n'
        if j == 10:
            line = line[:-1] + ' ' + str(min(y_L)) + '\n'
        if j == 11:
            line = line[:-1] + ' ' + str(max(y_L)) + '\n'
        if j == 12:
            line = line[:-1] + ' ' + str(min(y_L)) + '\n'
        if j == 13:
            line = line[:-1] + ' ' + str(max(y_L)) + '\n'
        if j == 48 or j == 60 or j == 72:
            line = line[:-1] + " '" + str(isotropic_stress) + "*t'\n"
        if j == 80:
            line = line[:-1] + ' ' + str(young_pore) + '\n'
        if j == 81:
            line = line[:-1] + ' ' + str(poisson_pore) + '\n'
        if j == 91:
            line = line[:-1] + ' ' + str(young_grain) + '\n'
        if j == 92:
            line = line[:-1] + ' ' + str(poisson_grain) + '\n'
        if j == 102:
            line = line[:-1] + ' ' + str(young_cement) + '\n'
        if j == 103:
            line = line[:-1] + ' ' + str(poisson_cement) + '\n'
        if j == 123 or j == 125 or j == 126:
            line = line[:-1] + ' ' + str(crit_res_fem) + '\n'
        if j == 132:
            line = line[:-1] + ' ' + str(dt_fem) + '\n'
        file_to_write.write(line)
    file_to_write.close()

#-------------------------------------------------------------------------------

def Read_FEM_csv(namefile, M_grain, M_cement):
    '''
    Read the csv file after a Moose FEM simulation.
    '''
    f = open(namefile, "r")
    lines = f.readlines()
    f.close()
    # init data
    L_time = []
    L_stress_xx_grain = []
    L_stress_xy_grain = []
    L_stress_xz_grain = []
    L_stress_yy_grain = []
    L_stress_yz_grain = []
    L_stress_zz_grain = []
    
    L_strain_xx_grain = []
    L_strain_xy_grain = []
    L_strain_xz_grain = []
    L_strain_yy_grain = []
    L_strain_yz_grain = []
    L_strain_zz_grain = []
    
    L_stress_xx_cement = []
    L_stress_xy_cement = []
    L_stress_xz_cement = []
    L_stress_yy_cement = []
    L_stress_yz_cement = []
    L_stress_zz_cement = []
    
    L_strain_xx_cement = []
    L_strain_xy_cement = []
    L_strain_xz_cement = []
    L_strain_yy_cement = []
    L_strain_yz_cement = []
    L_strain_zz_cement = []

    # prepare homogenization
    L_stress_xx = []
    L_stress_xy = []
    L_stress_xz = []
    L_stress_yy = []
    L_stress_yz = []
    L_stress_zz = []
    
    L_strain_xx = []
    L_strain_xy = []
    L_strain_xz = []
    L_strain_yy = []
    L_strain_yz = []
    L_strain_zz = []
    
    s_grain = np.sum(M_grain)
    s_cement = np.sum(M_cement)

    # iterate on lines
    for line in lines[1:]:
        line = line.replace("\n", "")
        data = line.split(',')
        # read data
        L_time.append(float(data[0]))

        L_stress_xx_grain.append(float(data[14]))
        L_stress_xy_grain.append(float(data[16]))
        L_stress_xz_grain.append(float(data[18]))
        L_stress_yy_grain.append(float(data[20]))
        L_stress_yz_grain.append(float(data[22]))
        L_stress_zz_grain.append(float(data[24]))

        L_strain_xx_grain.append(float(data[2]))
        L_strain_xy_grain.append(float(data[4]))
        L_strain_xz_grain.append(float(data[6]))
        L_strain_yy_grain.append(float(data[8]))
        L_strain_yz_grain.append(float(data[10]))
        L_strain_zz_grain.append(float(data[12]))

        L_stress_xx_cement.append(float(data[13]))
        L_stress_xy_cement.append(float(data[15]))
        L_stress_xz_cement.append(float(data[17]))
        L_stress_yy_cement.append(float(data[19]))
        L_stress_yz_cement.append(float(data[21]))
        L_stress_zz_cement.append(float(data[23]))

        L_strain_xx_cement.append(float(data[1]))
        L_strain_xy_cement.append(float(data[3]))
        L_strain_xz_cement.append(float(data[5]))
        L_strain_yy_cement.append(float(data[7]))
        L_strain_yz_cement.append(float(data[9]))
        L_strain_zz_cement.append(float(data[11]))

        # compute homogenization
        L_stress_xx.append((L_stress_xx_grain[-1]*s_grain + L_stress_xx_cement[-1]*s_cement)/(s_grain+s_cement))
        L_stress_xy.append((L_stress_xy_grain[-1]*s_grain + L_stress_xy_cement[-1]*s_cement)/(s_grain+s_cement))
        L_stress_xz.append((L_stress_xz_grain[-1]*s_grain + L_stress_xz_cement[-1]*s_cement)/(s_grain+s_cement))
        L_stress_yy.append((L_stress_yy_grain[-1]*s_grain + L_stress_yy_cement[-1]*s_cement)/(s_grain+s_cement))
        L_stress_yz.append((L_stress_yz_grain[-1]*s_grain + L_stress_yz_cement[-1]*s_cement)/(s_grain+s_cement))
        L_stress_zz.append((L_stress_zz_grain[-1]*s_grain + L_stress_zz_cement[-1]*s_cement)/(s_grain+s_cement))

        L_strain_xx.append((L_strain_xx_grain[-1]*s_grain + L_strain_xx_cement[-1]*s_cement)/(s_grain+s_cement))
        L_strain_xy.append((L_strain_xy_grain[-1]*s_grain + L_strain_xy_cement[-1]*s_cement)/(s_grain+s_cement))
        L_strain_xz.append((L_strain_xz_grain[-1]*s_grain + L_strain_xz_cement[-1]*s_cement)/(s_grain+s_cement))
        L_strain_yy.append((L_strain_yy_grain[-1]*s_grain + L_strain_yy_cement[-1]*s_cement)/(s_grain+s_cement))
        L_strain_yz.append((L_strain_yz_grain[-1]*s_grain + L_strain_yz_cement[-1]*s_cement)/(s_grain+s_cement))
        L_strain_zz.append((L_strain_zz_grain[-1]*s_grain + L_strain_zz_cement[-1]*s_cement)/(s_grain+s_cement))

    return L_stress_xx, L_stress_xy, L_stress_xz, L_stress_yy, L_stress_yz, L_stress_zz, \
           L_strain_xx, L_strain_xy, L_strain_xz, L_strain_yy, L_strain_yz, L_strain_zz

#-------------------------------------------------------------------------------

def lsm_linear(L_y, L_x):
    '''
    Least square method to determine y = ax + b.
    '''
    # compute sums
    s_1 = 0
    s_2 = 0
    s_3 = 0
    s_4 = 0
    s_5 = 0
    for i in range(len(L_y)):
        s_1 = s_1 + 1*L_x[i]*L_x[i]
        s_2 = s_2 + 1*L_x[i]
        s_3 = s_3 + 1
        s_4 = s_4 + 1*L_x[i]*L_y[i]
        s_5 = s_5 + 1*L_y[i]
    # solve problem
    M = np.array([[s_1, s_2],[s_2, s_3]])
    V = np.array([s_4, s_5])
    result = np.linalg.solve(M, V)
    a = result[0]
    b = result[1]
    # correlation linear
    cov = 0
    vx = 0
    vy = 0
    for i in range(len(L_y)):
        cov = cov + (L_x[i]-np.mean(L_x))*(L_y[i]-np.mean(L_y))
        vx = vx + (L_x[i]-np.mean(L_x))*(L_x[i]-np.mean(L_x))
        vy = vy + (L_y[i]-np.mean(L_y))*(L_y[i]-np.mean(L_y))
    corr = cov/(math.sqrt(vx*vy))
    return a, b, corr

#-------------------------------------------------------------------------------

def Interpolate_compression_props(L_strain_zz, L_stress_zz):
    '''
    Interpolate the mechanical properties from a compression test:
        - Young Modulus Y    
    '''
    # interpolate function
    a, b, corr = lsm_linear(L_stress_zz, L_strain_zz)
    # print result
    #print('\nYoung Modulus interpolation (y=ax+b):')
    #print('a:', a, 'b:', b, 'cor:', corr)
    # save parameter
    YoungModulusSample = a

    return YoungModulusSample

#-------------------------------------------------------------------------------

def Interpolate_triaxial_props(L_strain_xx, L_strain_yy, L_strain_zz, L_stress_zz):
    '''
    Interpolate the mechanical properties from a compression test:
        - Young Modulus Y
        - Shear Modulus G     
        - Poisson ratio v     
    '''
    # interpolate function
    a, b, corr = lsm_linear(L_stress_zz, L_strain_zz)
    # print result
    #print('\nYoung Modulus interpolation (y=ax+b):')
    #print('a:', a, 'b:', b, 'cor:', corr)
    # save parameter
    YoungModulusSample = a

    # interpolate function
    a, b, corr = lsm_linear(L_strain_xx, L_strain_zz)
    # print result
    #print('\nYoung Modulus interpolation (y=ax+b):')
    #print('a:', a, 'b:', b, 'cor:', corr)
    # save parameter
    PoissonRatioSample_xz = -a
    # interpolate function
    a, b, corr = lsm_linear(L_strain_yy, L_strain_zz)
    # print result
    #print('\nYoung Modulus interpolation (y=ax+b):')
    #print('a:', a, 'b:', b, 'cor:', corr)
    # save parameter
    PoissonRatioSample_yz = -a
    # compute the mean
    PoissonRatioSample = PoissonRatioSample_xz*0.5 + PoissonRatioSample_yz*0.5

    return YoungModulusSample, PoissonRatioSample
    
#-------------------------------------------------------------------------------

def Interpolate_shearing_prop(L_strain, L_stress_xy):
    '''
    Interpolate the mechanical propertie from a shear test:
        - Shear Modulus G     
    ''' 
    # TO DO
    pass

#-------------------------------------------------------------------------------

def Interpolate_isotropic_prop(L_strain):
    '''
    Interpolate the mechanical propertie from an isotropic test:
        - Bulk Modulus K     
    ''' 
    # pp data
    L_volumetric_strain = []
    for i_strain in range(len(L_strain_xx)):
        L_volumetric_strain.append(L_strain_xx[i_strain]+L_strain_yy[i_strain]+L_strain_zz[i_strain])    
    # interpolate function
    a, b, corr = lsm_linear(L_mean_stress, L_volumetric_strain)
    # print result
    #print('\nBulk Modulus interpolation (y=ax+b):')
    #print('a:', a, 'b:', b, 'cor:', corr)
    # save parameter
    BulkModulusSample = a
    
    return BulkModulusSample

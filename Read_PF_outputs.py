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

def Read_PF_csv(crit_res_pf):
    '''
    Read the csv output of the PF simulation.
    '''
    print('Read csv')

    # read file
    f = open('csv/PF_Dissolution_csv.csv', "r")
    lines = f.readlines()
    f.close()
    # init data
    time_pp = []
    matter_pp = []
    grain_pp = []
    cement_pp = []

    # iterate on lines
    for line in lines[1:]:
        line = line.replace("\n", "")
        data = line.split(',')
        # read data
        time_pp.append(float(data[0]))
        matter_pp.append(float(data[1]))
    
    # save
    dict_save = {
        'crit_res_pf': crit_res_pf,
        'time': time_pp,
        'matter': matter_pp
    }
    pickle.dump(dict_save, open('pp/hydration_csv/'+str(crit_res_pf)+'.dict','wb'))

#-------------------------------------------------------------------------------

def Read_PF_vtk(iteration_str, L_L_i_XYZ_not_used, L_XYZ_used, n_proc, pf_map_matter, x_L, y_L, z_L, data_grain, data_cement):
    '''
    Read the vtk outputs of the PF simulation.
    '''
    if iteration_str == '000':
        print('Read vtk')

    # template of the files read
    template_file = 'vtk/PF_Dissolution_other_'

    # initialization
    L_matter = []

    # Help the algorithm to know which node to used
    if L_L_i_XYZ_not_used == []:
        L_XYZ_used = []
        know_map = False
    else : 
        know_map = True

    # iterate on the proccessors used
    for i_proc in range(n_proc):
        if iteration_str == '000':
            print('proc', i_proc+1, '/', n_proc)

        # name of the file to load
        namefile = template_file+iteration_str+'_'+str(i_proc)+'.vtu'

        # load a vtk file as input
        reader = vtk.vtkXMLUnstructuredGridReader()
        reader.SetFileName(namefile)
        reader.Update()

        # Grab a scalar from the vtk file
        nodes_vtk_array= reader.GetOutput().GetPoints().GetData()
        matter_vtk_array = reader.GetOutput().GetPointData().GetArray("matter")
        
        #Get the coordinates of the nodes and the scalar values
        nodes_array = vtk_to_numpy(nodes_vtk_array)
        matter_array = vtk_to_numpy(matter_vtk_array)
        
        # Help the algorithm to know which nodes to use
        if not know_map:
            L_i_XYZ_not_used = []
        
        # First iteration must detect common zones between processors
        if not know_map:
            # iterate on the nodes
            for i_XYZ in range(len(nodes_array)) :
                XYZ = nodes_array[i_XYZ]
                # Do not consider twice a point
                if list(XYZ) not in L_XYZ_used :
                    L_XYZ_used.append(list(XYZ))
                    L_matter.append(matter_array[i_XYZ])
                # Help the algorithm to know which nodes to used
                else :
                    L_i_XYZ_not_used.append(i_XYZ)
            # Help the algorithm to know which nodes to used
            L_L_i_XYZ_not_used.append(L_i_XYZ_not_used)
        # Algorithm knows already where to look
        else :
            L_i_XYZ_not_used = L_L_i_XYZ_not_used[i_proc]
            # all data are considered
            if L_i_XYZ_not_used == []:
                L_matter = L_matter + list(matter_array)
            # not all data are considered
            else :
                L_matter = L_matter + list(matter_array[:L_i_XYZ_not_used[0]])
                
    # Rebuild maps     
    M_grain, M_cement = Map_from_list(pf_map_matter, L_XYZ_used, L_matter, x_L, y_L, z_L, data_grain, data_cement, iteration_str)

    return L_L_i_XYZ_not_used, L_XYZ_used, M_grain, M_cement

#-------------------------------------------------------------------------------

def Map_from_list(pf_map_matter, L_XYZ_used, L_matter, x_L, y_L, z_L, data_grain, data_cement, iteration_str):
    '''
    Rebuild numpy array from a list of values.
    '''
    # initialization
    M_matter = np.zeros(pf_map_matter.shape)

    # iterate on the domain
    for i in range(len(L_XYZ_used)):
        # interpolate meshes
        find_ix = abs(np.array(x_L)-L_XYZ_used[i][0])
        find_iy = abs(np.array(y_L)-L_XYZ_used[i][1])
        find_iz = abs(np.array(z_L)-L_XYZ_used[i][2])
        i_x = list(find_ix).index(min(find_ix))
        i_y = list(find_iy).index(min(find_iy))
        i_z = list(find_iz).index(min(find_iz))
        # rebuild
        M_matter[i_x, i_y, i_z] = L_matter[i]

    # Apply masks to retrieve grain and cement
    M_grain, M_cement = Grain_Cement_from_Matter(M_matter, data_grain, data_cement, iteration_str)

    return M_grain, M_cement

#-------------------------------------------------------------------------------

def Grain_Cement_from_Matter(M_matter, data_grain, data_cement, iteration_str):
    '''
    Apply the initial mask to retrieve grain and cement from matter.
    '''
    # initialization
    M_cement = np.zeros(M_matter.shape)
    M_grain = np.zeros(M_matter.shape)    
    M_matter_pp = np.zeros(M_matter.shape)
    # iterate on the mesh
    for i_x in range(M_matter.shape[0]):
        for i_y in range(M_matter.shape[1]):
            for i_z in range(M_matter.shape[1]):
                if data_grain[i_x, i_y, i_z] == 1 and M_matter[i_x, i_y, i_z] > 0.5:
                    M_grain[i_x, i_y, i_z] = 1
                    M_matter_pp[i_x, i_y, i_z] = 1
                if data_cement[i_x, i_y, i_z] == 1 and M_matter[i_x, i_y, i_z] > 0.5:
                    M_cement[i_x, i_y, i_z] = 1   
                    M_matter_pp[i_x, i_y, i_z] = 0.5
    # plot
    #fig, (ax1, ax2) = plt.subplots(nrows=1,ncols=2,figsize=(16,9))
    #ax1.imshow(M_grain, cmap='binary')
    #ax1.set_title('grain')
    #ax2.imshow(M_cement, cmap='binary')
    #ax2.set_title('cement')
    #fig, (ax1) = plt.subplots(nrows=1,ncols=1,figsize=(16,9))
    #ax1.imshow(M_matter_pp, cmap='binary')
    #ax1.set_title('grain (black) and cement (grey)')
    #fig.suptitle(iteration_str)
    #fig.tight_layout()
    #fig.savefig('output/maps_bin_cement_grain_output/'+iteration_str+'.png')
    #plt.close(fig)

    return M_grain, M_cement




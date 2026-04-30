# -*- coding: utf-8 -*-
"""
Evan Marschall

Script to convert gmsh file to RSQSim fault file, this will use tets

Good luck
"""
import numpy as np
import pandas as pd
#%% 
"""
Define the paths and some parameters first, names of faults, fault indices from gmsh, 
number of faults, sliprates,rakes, etc. 
"""
indir = "/Users/emarschall/Desktop/RSQSim_stuff/Open_Source_Stuff/Getting_started_example/meshing_gmsh_example/" #directory where msh file is
gfile = "step_over.msh" # mesh file
outdir = indir # Change if desired 
outfile = "Step_over.flt"

g_surfaces = [102,103] #surface numbers in g mesh
section_number = [102,103]
rsqsim_names = ["seg1","seg2"] #names for RSQSim segments
sliprates = [1e-09,1e-09] #sliprates for segments (m/s)
rakes = [180,180]

#%%
# Get the total number of nodes 
file = open(indir+gfile)
num_nodes = int(file.readlines()[4])
file.close()

# Panda data frame with all nodes 
df_nodes = pd.read_csv(indir+gfile, sep = " ", header =None, skiprows = 5, nrows = num_nodes)

#%%
# Get all the element data
len2elements = 4+ num_nodes+ 3 # number of lines to get to # of elements 

file = open(indir+gfile)
num_elements = int(file.readlines()[len2elements])
file.close()

df_elements = pd.read_csv(indir+gfile, sep = " ", header =None, skiprows = len2elements+1, nrows = num_elements, usecols= [0,1,2,3,4,5,6,7])

#Filter the data fram to only get values for the fault surfaces list in the (g_surfaces list at the start opf script)
filtered_df = df_elements[df_elements[3].isin(g_surfaces)]
filtered_df = filtered_df.drop([0,1,2],axis=1) 
#replace g_surface ewith RSQSim seg name
#%%

for i in range(len(g_surfaces)):
    filtered_df[3] = filtered_df[3].replace(g_surfaces[i],rsqsim_names[i] )

name_list = filtered_df[3].tolist()
num_list = filtered_df[4].tolist()   

rake_list = name_list
sliprate_list = name_list
for i in range(len(g_surfaces)):
    value_to_replace = rsqsim_names[i]
    new_value_rake = rakes[i]
    new_value_sliprate = sliprates[i]
    rake_list = [new_value_rake if item == value_to_replace else item for item in rake_list]
    sliprate_list = [new_value_sliprate if item == value_to_replace else item for item in sliprate_list]
    




#%% Convert to RSQSim style 

vertex_array = np.full((len(filtered_df),9),np.nan)
## x1 y1 z1 x2 y2 z2 x3 y3 z3 rake sliprate section_number section_name (no white spaces in the name)

#have to subtract 1 becasue python is 0 indexed and gmsh file is not
for i in range(len(filtered_df)):
    index1 = filtered_df[5].values[i] -1
    index2 = filtered_df[6].values[i] -1
    index3 = filtered_df[7].values[i] -1
    
    vertex_array[i,0] = df_nodes[1].values[index1]
    vertex_array[i,1] = df_nodes[2].values[index1]
    vertex_array[i,2] = df_nodes[3].values[index1]
    vertex_array[i,3] = df_nodes[1].values[index2]
    vertex_array[i,4] = df_nodes[2].values[index2]
    vertex_array[i,5] = df_nodes[3].values[index2]
    vertex_array[i,6] = df_nodes[1].values[index3]
    vertex_array[i,7] = df_nodes[2].values[index3]
    vertex_array[i,8] = df_nodes[3].values[index3]
    
    

df_rsqsim = pd.DataFrame(vertex_array)
df_rsqsim['rake'] = rake_list
df_rsqsim['sliprate'] = sliprate_list
df_rsqsim['segnum'] = num_list
df_rsqsim['segname'] = name_list

#%%
#Write file
df_rsqsim.to_csv(outdir+outfile,sep =" ",index = False, header = False)











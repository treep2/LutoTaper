# -*- coding: utf-8 -*-
"""
Created on Sat Nov 19 03:57:11 2022

@author: treep
"""




import numpy as np
from scipy.constants import c

import sys

#sys.path.append(r'P:\1People files\Tom Reep\Measurement files\SA measurements 1005')

#import utilities as ut

import sys, os
sys.path.append(r"C:\Program Files\Lumerical\v231\api\python") #Default windows lumapi path

import os, sys
import numpy as np
import pandas as pd
import scipy as sc

import lumapi
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import pickle

from scipy.interpolate import interp1d


scolours = ['tab:blue','tab:red','tab:orange','tab:green','tab:purple','tab:brown','tab:pink','tab:grey','tab:olive','tab:cyan','blue','red','purple','navy','sienna', 'maroon','violet','teal', 'yellowgreen']
scolours = scolours*500 #lifehack

saving_directory = r'C:\Users\treep\OneDrive - UGent\Documents\UGent PhD\Side projects\Lionix MLL\Taper optimisation\Data'
#%%
mode = lumapi.MODE()

mset = mode.set
#%% Define parameters 

SiN_width = 3e-6                # Width of SiN WG
SiN_height = 300e-9             # Height of SiN
SiN_depth = 0e-9                # How far is the SiN buried in oxide

#Silicon taper parameters
taper_length_ini=80e-6          # First time simulation of taper length (doesnt matter)

Si_height = 400e-9              # Height of the silicon waveguide
Si_width_tip = 100e-9           # Silicon tip width | if 0 simulator will find the tip width
Si_width_WG = 2e-6              # Final silicon waveguide width

#Thickness of BCB
BCB_inter = 50e-9               # thickness BCB between source target
AlOx_thickness = 35e-9
BCB_clad = 0e-6                 # If BCB cladding is used (usually the case...)

# Number of cells
Ebeam_wtolerance = 30e-9        # For the tolerant design, find out how wide it can go

CellN = round((Si_width_WG-Si_width_tip)/Ebeam_wtolerance) #Calculate the number of cells required
Number_cells_perGroup = 3       #For each cell, add 3 intermediate cells

MaxLength= 500e-6               # Not critical, just make sure its big

## EME stuff
n_background = 1                # Usually no oxide background for EME simulations
wavelength = 1550e-9            # Wavelength for the simulation
Simwidth = 10e-6                # Width of the simulation region
SimHeight = 4.5e-6              # Height of the simulation region
num_modes = 15                  # Number of simulation modes foreach cell

#%% Building the simulation region
## MATERIAL
try:
    opt_material=mode.addmaterial('Dielectric');
    mode.setmaterial(opt_material,'name','SiN');
    mode.setmaterial('SiN','Refractive Index',2);
    mode.setmaterial('SiN',"color", np.array([166/255, 41/255, 225/255, 220/255]));
    
    opt_material=mode.addmaterial('Dielectric');
    mode.setmaterial(opt_material,'name','SiOx');
    mode.setmaterial('SiOx','Refractive Index',1.445);
    mode.setmaterial('SiOx',"color", np.array([153/255,255/255, 255/255, 255/255]));
    
    opt_material=mode.addmaterial('Dielectric');
    mode.setmaterial(opt_material,'name','BCB');
    mode.setmaterial('BCB','Refractive Index',1.54);
    mode.setmaterial('BCB',"color", np.array([153/255, 0, 0,120/255]));
    
    opt_material=mode.addmaterial('Dielectric');
    mode.setmaterial(opt_material,'name','air');
    mode.setmaterial('air','Refractive Index',1);
    mode.setmaterial('air',"color", np.array([0, 0, 0, 0]));
    
except:
    print('Materials already added')

## CLEAR SESSION
mode.switchtolayout()
mode.selectall()
mode.delete()


### Setting up the simulation
mode.addrect()
mset('name','BOX')
mset('x min',-10e-6)
mset('x max',MaxLength + 10e-6)
mset('y min',-100e-6)
mset('y max',100e-6)
mset('z min',-3e-6)
mset('z max',-SiN_height-SiN_depth)
mset('material','SiOx')
mset('alpha',0.3)
mset("override mesh order from material database",1)
mset("mesh order",3)


mode.addrect()
mset('name','AlOx')
mset('x min',-10e-6)
mset('x max',MaxLength + 10e-6)
mset('y min',-SiN_width/2)
mset('y max',SiN_width/2)
mset('z min',0)
mset('z max',AlOx_thickness)
mset('index',1.7)
mset('alpha',0.3)
mset("override mesh order from material database",1)
mset("mesh order",3)


#INPUT WAVEGUIDE
mode.addrect()
mset('name','wg SiN')
mset('x min',-10e-6)
mset('x max',MaxLength + 10e-6)
mset('y min',-SiN_width/2)
mset('y max',SiN_width/2)
mset('z min',-SiN_height-SiN_depth)
mset('z max',-SiN_depth)
mset('material','SiN')

#BCB
mode.addrect()
mset('name','BCB layer')
mset('x min',-10e-6)
mset('x max',MaxLength + 10e-6)
mset('y min',-100e-6)
mset('y max',100e-6)
mset('z min',AlOx_thickness)
mset('z max',AlOx_thickness+BCB_inter+BCB_clad)
mset('material','BCB')
mset("override mesh order from material database",1)
mset("mesh order",4)
mset('alpha',0.3)

#OUTPUT WAVEGUIDE
mode.addrect()      
mset("name", "Silicon output")
mset("x min", taper_length_ini)
mset("x max", MaxLength+10e-6)
mset("y", 0)
mset("y span", Si_width_WG)
mset("z min", AlOx_thickness+BCB_inter)
mset("z max", AlOx_thickness+BCB_inter+Si_height)
mset("material","Si (Silicon) - Palik")
mset("override mesh order from material database",1)
mset("mesh order",1)
mset("alpha",0.7)

dmesh = 15e-9 #Meshsize
mode.addmesh()
mset("name","mesh fine")
mset("y min",-max(SiN_width,Si_width_WG)/2-40e-9)
mset("y max",max(SiN_width,Si_width_WG)/2+40e-9)
mset("z min",-SiN_height-SiN_depth-40e-9)
mset("z max",AlOx_thickness+BCB_inter+Si_height+40e-9)
mset('x min',-10e-6)
mset('x max',MaxLength + 10e-6)
mset("dx",dmesh)
mset("dy",dmesh)   # changed
mset("dz",dmesh)   # changed

#%% Add the adiabatic mode converted defined in c-Si (linear taper)
taper_x = np.append(np.linspace(0,taper_length_ini,CellN+1),np.linspace(taper_length_ini,0,CellN+1)) #Actually the CellN is a bit overkill
taper_w = np.append(np.linspace(Si_width_tip/2,Si_width_WG/2,CellN+1),np.linspace(-Si_width_WG/2,-Si_width_tip/2,CellN+1)) #Just a linear taper

mode.addpoly()
mset('vertices',np.array([list(x) for x in list(zip(taper_x,taper_w))]))
mset("name", "Taper silicon")
mset("x",0)
mset("y",0)
mset("z min", BCB_inter+AlOx_thickness)
mset("z max", BCB_inter+Si_height+AlOx_thickness)
mset("material","Si (Silicon) - Palik")
mset("override mesh order from material database",1)
mset("mesh order",1)
mset("alpha",0.7)

#%% Add the EME solver

mode.addeme()
mset("solver type", "3D: X Prop")
mset("index", 1)
mset("wavelength", wavelength)
mset("y min", -Simwidth/2)
mset("y max", Simwidth/2)
mset("z max",SimHeight/2)
mset("z min",-SimHeight/2)
mset("x min", -5e-6)
mset("number of cell groups", CellN+2) #Input -> taper = 2
mset("display cells", 1)
mset("number of modes for all cell groups", num_modes)
mset("number of periodic groups", 1)
mset("energy conservation", "make passive") # or "none", "conserve energy"
mset("subcell method", np.array([0]+[0]*CellN+[0]))                             #No CVCS as it seems to do strange stuff in this optimisation approach
mset("z min bc", "Metal")
mset("z max bc", "Metal")
mset("y min bc", "Metal")
mset("y max bc", "Metal")
mset("cells", np.array([1]+[Number_cells_perGroup]*(CellN)+[1]))                #Number of cells per section, can also be optimised for
mset("group spans",np.array([5e-6]+[taper_length_ini/(CellN)]*CellN + [5e-6]))  #Lengths of each section

#Set coarse mesh settings
mode.setnamed("EME","define y mesh by",1);
mode.setnamed("EME","define z mesh by",1);
mode.setnamed("EME","dy",200e-9); 
mode.setnamed("EME","dz",200e-9);


#update port configuration (input and output port)
mode.setnamed("EME::Ports::port_1", "y", 0);
mode.setnamed("EME::Ports::port_1", "y span", 8e-6);
mode.setnamed("EME::Ports::port_1", "z", 0);
mode.setnamed("EME::Ports::port_1", "z span",8e-6);
mode.select("EME::Ports::port_1");
mode.setnamed("EME::Ports::port_1", "mode selection", "fundamental TE mode");
mode.updateportmodes();


mode.setnamed("EME::Ports::port_2", "y", 0);
mode.setnamed("EME::Ports::port_2", "y span", 8e-6);
mode.setnamed("EME::Ports::port_2", "z", 0);
mode.setnamed("EME::Ports::port_2", "z span", 8e-6);
mode.select("EME::Ports::port_2");
mode.setnamed("EME::Ports::port_2", "mode selection", "fundamental TE mode");
mode.updateportmodes();
#%% Run simulation

mode.save(r'C:\\Temp_lumerical\\Running_taper.lms')
mode.run()

#%% Optimising the simulation | Tricky part since you need to find best transmission vs length
#mode.load(saving_directory+'\\EME_opt1')

maxCell1 = 3e-6                                                                     # First sweep max length per cell
Sweeppoint1 = 101                                                                   # Number of sweep points in first cell
maxCell2 = 70e-6                                                                    # If first sweep is not long enough go to long sweep
Sweeppoint2 = 301                                                                   # Number of points in second sweep

mode.setemeanalysis("subcell method",np.array([0]+[0]*CellN+[0]))                   # Set subcell to no CVCS, doesn't work well currently and Lumerical doesn't explain how they approximate intermediate sections using this method

mode.setemeanalysis("Propagation sweep",1)                                          # Enable propagation sweep
df_final = pd.DataFrame(data={'Widths':taper_w[:len(taper_w)//2]*2})                # Define a df to save the results
df_final.loc[0,'x_opt'] = 0                                                         # Used to calculate the taper length

list_sweeps = []
mode.setemeanalysis("group spans",np.array([5e-6]+[25e-6/(CellN)]*CellN+[5e-6]))    # Start with real shitty taper

def getEMEsweep(Groupname = '_', stop = 10e-6,points=10):
    mode.setemeanalysis("parameter",Groupname) 
    mode.setemeanalysis("start",0)
    mode.setemeanalysis("stop",stop)
    mode.setemeanalysis("number of points",points)
    mode.emesweep()
    
    sweep_results = mode.getemesweep("S")                                           # Get results from first EME sweep
    return pd.DataFrame(                                                              # Store these results in a dataframe
            data={'x':[x[0] for x in sweep_results[Groupname.replace(' ','_')]],\
            'S21_abs2':abs(sweep_results['s21'])**2})
        

import scipy.signal as scs
def findPeaks(df_EMEsweep,maxdiff = 0.01):
    if maxdiff == None:
        height=0        
    else:
        height=df_EMEsweep['S21_abs2'].max()-maxdiff
        
    peakidxs = scs.find_peaks(df_EMEsweep['S21_abs2'],height=height)
    df_temp = pd.DataFrame({'idx':peakidxs[0],'Ts':peakidxs[1]['peak_heights']}).set_index('idx')
    if(len(peakidxs[0]) > 0):
        if(maxdiff==None):
            returnidx = df_temp['Ts'].idxmax()
        else:
            returnidx = min(peakidxs[0])
    else:
        returnidx = 1 #sometimes the peak is at 0 length, but don't do steps so just do the minimal sweep length

    return  df_EMEsweep.loc[returnidx,'x']#,peakidxs[1]['peak_heights'][0]

mincounter = 0
for CN in range(1,CellN+1):
    Groupname = "group span "+str(CN+1)                                             # +1 because first cell is not sweepable
    df_EMEsweep = getEMEsweep(Groupname, maxCell1, Sweeppoint1)                     # Perform the EME sweep for cell CN
    normal_max = df_EMEsweep.loc[df_EMEsweep['S21_abs2'].idxmax(),'x']              # Calculated to find where the current optimum is. If it is larger than maxCell1/2 the longer sweep is performed
    if(normal_max >maxCell1/2 ):                                                    # Perform a longer sweep if the optimal length seems to be longer (mainly time saving measure)
        df_EMEsweep = getEMEsweep(Groupname, maxCell2, Sweeppoint2)                 # Perform longer EME sweep

    # Not so nice way of selecting the optimal length of the section
    try:
        df_final.loc[CN,'dx_opt']  = findPeaks(df_EMEsweep,maxdiff=0.01)
    except:        
        try:
            df_final.loc[CN,'dx_opt']  = findPeaks(df_EMEsweep,maxdiff=0.005)
        except:
            df_final.loc[CN,'dx_opt']  = findPeaks(df_EMEsweep,maxdiff=None)
    
    
    #Store the data in df_final:
    df_final.loc[CN,'x_opt'] = df_final.loc[CN-1,'x_opt'] + df_final.loc[CN,'dx_opt']

   #Update EME with optimised data before next sweep
    EME_spans = [x[0] for x in mode.getemeanalysis("group spans")]
    EME_spans[CN] = df_final.loc[CN,'dx_opt']
    mode.setemeanalysis("group spans",np.array(EME_spans))
    

    dx_opt = df_final.loc[CN,'dx_opt']
    x_opt= df_final.loc[CN,'x_opt']
      
        
    print(f'{round(CN)}/max(CN) | dx_opt = {round(dx_opt*1e6,2)}, length={round(x_opt*1e6,2)}')

#%% Making plot for paper
fig, ax = plt.subplots()
mode.setemeanalysis("group spans",np.array([5e-6]+[15e-6/(CellN)]*CellN+[5e-6]))    # Start with real shitty taper
mode.setemeanalysis("subcell method",np.array([0]+[0]*CellN+[0]))                   # Set subcell to no CVCS, doesn't work well currently and Lumerical doesn't explain how they approximate intermediate sections using this method

mincounter = 0
dfs = []
for CN in range(1,7):
    Groupname = "group span "+str(CN+1)                                             # +1 because first cell is not sweepable
    df_EMEsweep = getEMEsweep(Groupname, maxCell2, Sweeppoint2)                     # Perform the EME sweep for cell CN
    normal_max = df_EMEsweep.loc[df_EMEsweep['S21_abs2'].idxmax(),'x']              # Calculated to find where the current optimum is. If it is larger than maxCell1/2 the longer sweep is performed
    # if(normal_max >maxCell1/2 ):                                                    # Perform a longer sweep if the optimal length seems to be longer (mainly time saving measure)
    #     df_EMEsweep = getEMEsweep(Groupname, maxCell2, Sweeppoint2)                 # Perform longer EME sweep

    # Not so nice way of selecting the optimal length of the section
    try:
        df_final.loc[CN,'dx_opt']  = findPeaks(df_EMEsweep,maxdiff=0.01)
    except:        
        try:
            df_final.loc[CN,'dx_opt']  = findPeaks(df_EMEsweep,maxdiff=0.005)
        except:
            df_final.loc[CN,'dx_opt']  = findPeaks(df_EMEsweep,maxdiff=None)
    
    
    #Store the data in df_final:
    df_final.loc[CN,'x_opt'] = df_final.loc[CN-1,'x_opt'] + df_final.loc[CN,'dx_opt']

   #Update EME with optimised data before next sweep
    EME_spans = [x[0] for x in mode.getemeanalysis("group spans")]
    EME_spans[CN] = df_final.loc[CN,'dx_opt']
    mode.setemeanalysis("group spans",np.array(EME_spans))
    

    dx_opt = df_final.loc[CN,'dx_opt']
    x_opt= df_final.loc[CN,'x_opt']
      
    dfs+=[df_EMEsweep]
    print(f'{round(CN)}/max(CN) | dx_opt = {round(dx_opt*1e6,2)}, length={round(x_opt*1e6,2)}, T:{normal_max}')
    df_EMEsweep.plot(x='x',y='S21_abs2',ax=ax,color=scolours[CN])
    

#%% Check simulation and make the taper more resilient

def returnemeprop(df):
    mode.setemeanalysis("group spans",np.array([5e-6]+list(df)+[5e-6]))    # Start with really shitty taper
    mode.emepropagate()
    res = mode.getresult("EME",'user s matrix')
    return abs(res[0][1])**2

mode.setemeanalysis("subcell method", np.array([0]+[1]*(len(df_final)-1)+[0])) #turn on CVCS now, but don't think its necessary

T_opt = returnemeprop(df_final.loc[1:,'dx_opt'])

print("Optimal taper transmission is {T_opt} for taper length of {round(df_final['x_opt'].max()*1e6,2)} um")

#Calculate the resilient lengths
Ncell_look = 1 # Look one cell to either side, or two sides. Normally do 1
Weights = [1,1,1] # Take the maximum of 1 to left, 1 of center or 1 of right. Or alternatively do a weight. Sometimes we got better results with [0.75,1,0.75]
df_final['resilient_x'] = 0

for i in range(1,len(df_final)):
    minval = int(max(i-Ncell_look,0))
    maxval = int(min(i+Ncell_look, len(df_final)))
    if(i<1):
        length_next = max((df_final.loc[minval:maxval,'dx_opt']).fillna(0)*Weights[1:])
    elif(i<len(df_final)-1):
        length_next = max((df_final.loc[minval:maxval,'dx_opt']).fillna(0)*Weights)
    else:
        length_next = max((df_final.loc[minval:maxval,'dx_opt']).fillna(0)*Weights[:2])
    df_final.loc[i,'dx_opt_resilient'] = length_next
    df_final.loc[i,'resilient_x'] =  df_final.loc[i-1,'resilient_x']+length_next
    


fig,ax = plt.subplots()
ax.plot(df_final['x_opt']*1e6,df_final['Widths']*1e6,linestyle='-',marker='.',color='tab:blue',label='_')
ax.plot(df_final['resilient_x']*1e6,df_final['Widths']*1e6,linestyle='-',marker='.',color='tab:red',label='_resilient')

ax.set_xlabel('Taper length [µm]')
ax.set_ylabel('Taper width [µm]')
ax.grid()
ax.legend()    


## Save all the data
datastuff = {'df':df_final, 'maxx base':max(df_final['x_opt']),'maxx resilient':max(df_final['resilient_x'])}
file = open(r'C:\Users\treep\OneDrive - UGent\Documents\UGent PhD\Side projects\Lionix MLL\Taper optimisation\Data\taper_120_2um_075.pik','wb')
pickle.dump( datastuff,file)
file.close()

#Python 2 compatible (in case you are using ancient ipkiss like me)
file = open(r'C:\Users\treep\OneDrive - UGent\Documents\UGent PhD\Side projects\Lionix MLL\Taper optimisation\Data\taper_120_2um_p2_075.pik','wb')
pickle.dump( datastuff,file,protocol=2)
file.close()

df_final['dx'] = df_final['dx_opt_resilient']
df_final.to_csv(r'C:\Users\treep\OneDrive - UGent\Documents\UGent PhD\Side projects\Lionix MLL\Taper optimisation\Data\dataframe_120_2um_075.csv')
#%% Simulate EME offset_x and width variation (this will take hours... probably best to do it overnight)

Resilient = False

import pickle
file = open(r'C:\Users\treep\OneDrive - UGent\Documents\UGent PhD\Side projects\Lionix MLL\Taper optimisation\Data\taper_120_2um_075.pik','rb')

datastuff = pickle.load(file)
file.close()
df_taper = datastuff['df']
    
fig,ax = plt.subplots()
df_taper_test = pd.DataFrame(data={'x':df_taper['resilient_x'] if Resilient else df_taper['x_opt'],'Widths':df_taper['Widths']})
df_taper_test['xum'] = df_taper_test['x']*1e6
df_taper_test.plot(x='xum',y='Widths',marker='.',ax=ax)   

#Create optimised taper shape
def updateEME_offsetsimulation(width_offset = 0e-9, ymisalignment=0e-9, BCBthickness=BCB_inter,CellsEME = 150):
    taper_x = np.concatenate((df_taper_test['x'],\
                              np.flip(df_taper_test['x'])) )
    taper_w = np.concatenate(((df_taper_test['Widths']+width_offset)/2,np.flip((df_taper_test['Widths']+width_offset)*-1/2)))
    
    WGheight = mode.getnamed("Taper silicon","z max")-mode.getnamed("Taper silicon","z min")
    mode.switchtolayout()
    mode.setnamed('Taper silicon','vertices',np.array([list(x) for x in list(zip(taper_x,taper_w))]))#%%
    mode.setnamed('Taper silicon','y',ymisalignment)#%%

    mode.setnamed('Silicon output','y span',max(taper_w)*2)#%%
    mode.setnamed('Silicon output','y',ymisalignment)#%%
    mode.setnamed('Silicon output','x min',max(taper_x))
    
    mode.setnamed("Taper silicon","z min", AlOx_thickness+BCBthickness)
    mode.setnamed("Taper silicon","z max",AlOx_thickness+BCBthickness+WGheight)
    mode.setnamed("Silicon output","z min", AlOx_thickness+BCBthickness)
    mode.setnamed("Silicon output","z max",AlOx_thickness+BCBthickness+WGheight)
    
    mode.setnamed("mesh fine","y min", min([min(taper_w)+ymisalignment,-3.5e-6/2]))
    mode.setnamed("mesh fine","y max", max([max(taper_w)+ymisalignment,3.5e-6/2]))
    mode.setnamed("mesh fine","z max", AlOx_thickness+BCBthickness+WGheight+10e-9)
    
    mode.setnamed("EME","number of cell groups",3) #Input ->taper-> output
    mode.setnamed("EME","group spans",np.array([5e-6,max(taper_x),5e-6])) #Lengths of each section
    mode.setnamed("EME","subcell method", np.array([0,1,0]))
    mode.setnamed("EME","cells", np.array([1,CellsEME,1])) #Number of cells per section, can also be optimised for
    
#df_misalignment = pd.DataFrame(data={'y misalignment':[0e-9,100e-9,200e-9]})

df_wvar = pd.DataFrame(data={'w variation':[-60e-9,-45e-9,-30e-9,0,30e-9,45e-9,60e-9]})



# for misidx, misalignment in enumerate(df_misalignment['y misalignment']):
#     updateEME_offsetsimulation(ymisalignment = misalignment)
#     mode.save(r'C:\\Temp_lumerical\\Test')
#     mode.run()
    
#     mode.emepropagate()
    
#     Result = mode.getresult("EME","user s matrix")
    
#     df_misalignment.loc[misidx,'Transmission'] = np.abs(Result[0][1])**2
    

# datastuff = {'df_misalignment':df_misalignment,'df_wvar':df_wvar}
# ## Save intermediate stuff in case things crash
# import pickle
# file = open(r'C:\Users\treep\OneDrive - UGent\Documents\UGent PhD\Simulation\Lumerical\Tapers\Si for Lionix\Results\Analysis\taper_120_2um_Stijn_Lionix_300nmthick_onlyAlOx_base.pik','wb')
# pickle.dump(datastuff,file)
# file.close()

for widx, woffset in enumerate(df_wvar['w variation']):
    print(f'w = {woffset*1e9}')
    updateEME_offsetsimulation(width_offset = woffset)
    mode.save(r'C:\\Temp_lumerical\\Test')
    mode.run()
    
    mode.emepropagate()
    
    Result = mode.getresult("EME","user s matrix")
    
    df_wvar.loc[widx,'Transmission'] = np.abs(Result[0][1])**2
    print(Result)

df_wvar.to_csv('Temp_075_nonresilient.csv')


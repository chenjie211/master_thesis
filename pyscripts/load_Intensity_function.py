#!/usr/bin/env python
# coding: utf-8

# In[18]:


def LoadIntensity(file=None, N_training_data=100, test_simulation=False):
    
    import numpy as np
    
    data_file = "/Users/chenjie/pyQEDA/pyQEDA_simulations/" + file + "/"
    experimental_data_file =  data_file + file + "_experimental_data/"
    training_data_file = data_file + file + "_training_data"
    simulation_data_file = data_file + file + "/"
    
    # read the number of beams amd tilts
    tilt_num = len(open(experimental_data_file + "_Tilts.txt", "r").readlines()) - 1
    beam_num = len(open(experimental_data_file + "_BeamList.txt", "r").readlines()) - 1
    print("Number of Tilts: ", tilt_num)
    print("Number of Beams: ", beam_num)
            
    
    # load the Intensities of experimental pattern
    exp_Intensity = np.zeros([tilt_num,beam_num])
    for row in range(tilt_num):
        Ar = open(experimental_data_file + "Ar_" + str(row) + ".txt","r").readlines()
        excited_list = []
        for line in range(len(Ar)):
            a = Ar[line].split()
            excited = int(a[0].strip(' '))
            excited_list.append(excited)
        Intensity = open(experimental_data_file + "Intensities_" + str(row) + ".txt","r").readlines()
        Intensity_list = []
        for line in range(len(Intensity)):
            b = Intensity[line].split()
            intensity = float(b[0].strip(' '))
            Intensity_list.append(intensity)
        for num in range(len(Intensity_list)):
            beam = excited_list[num]
            exp_Intensity[row,beam]=Intensity_list[num]
    #exp_Intensity = exp_Intensity/exp_Intensity.max()
       
     
    # load the Intensities of random pattern        
    random_Intensity = []
    for training_data in range(N_training_data):
        pattern_Intensity = np.zeros([tilt_num,beam_num])
        for row in range(tilt_num):
            Ar = open(training_data_file + "_" + str(training_data)+ "/" + "Ar_" + str(row) + ".txt","r").readlines()
            excited_list = []
            for line in range(len(Ar)):
                a = Ar[line].split()
                excited = int(a[0].strip(' '))
                excited_list.append(excited)
            Intensity = open(training_data_file + "_" + str(training_data)+ "/" + "Intensities_" + str(row) + ".txt","r").readlines()
            Intensity_list = []
            for line in range(len(Intensity)):
                b = Intensity[line].split()
                intensity = float(b[0].strip(' '))
                Intensity_list.append(intensity)
            for num in range(len(Intensity_list)):
                beam = excited_list[num]
                pattern_Intensity[row,beam]=Intensity_list[num]
        #pattern_Intensity = pattern_Intensity/pattern_Intensity.max()
        random_Intensity.append(pattern_Intensity)
    random_Intensity = np.array(random_Intensity)
    
    
    # load the Intensities of simulation pattern
    if test_simulation:
        sim_Intensity = np.zeros([tilt_num,beam_num])
        for row in range(tilt_num):
            Ar = open(simulation_data_file + "Ar_" + str(row) + ".txt","r").readlines()
            excited_list = []
            for line in range(len(Ar)):
                a = Ar[line].split()
                excited = int(a[0].strip(' '))
                excited_list.append(excited)
            Intensity = open(simulation_data_file + "Intensities_" + str(row) + ".txt","r").readlines()
            Intensity_list = []
            for line in range(len(Intensity)):
                b = Intensity[line].split()
                intensity = float(b[0].strip(' '))
                Intensity_list.append(intensity)
            for num in range(len(Intensity_list)):
                beam = excited_list[num]
                sim_Intensity[row,beam]=Intensity_list[num]
    #sim_Intensity = sim_Intensity/sim_Intensity.max()
    
    
    if test_simulation:
        print("shape of sim_Intensity: ",sim_Intensity.shape)
        return sim_Intensity
    else:
        print("shape of exp_Intensity: ",exp_Intensity.shape)
        print("shape of random_Intensity: ",random_Intensity.shape)
        return tilt_num, beam_num, exp_Intensity, random_Intensity


# -------------------------------------------------------------------------------------------------------------------

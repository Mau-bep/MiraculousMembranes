import numpy as np
import matplotlib.pyplot as plt 
import os 
import json




def Read_data(folderpath):
    # THis are 2 bead simulations
    dirs = os.listdir(folderpath)

    f = open(folderpath + "Bending_data_phase.csv","w+")
    for dir in dirs:
        if os.path.isdir(folderpath + dir):
            print("Reading data from {}".format(dir))
            # I want the input file
            json_file = folderpath + dir +"/Input_file.json"
            with open(json_file) as fjson:
                d = json.load(fjson)

            # So now d has the whole information 
            for Energy in d["Energy"]:
                if(Energy["Name"] == "Coverage"):
                    covStrength = Energy["constants"][0]
                    coverage = Energy["constants"][2]
            distance = d["theta"]

            # That is the distance and the target cov
            # I would love to also get the final coverage,
            # I now need to read the output fil
            output_file = folderpath + dir + "Output_file.txt"
            Sim_data = np.loadtxt(output_file,skiprows=1)

            E_Bend = Sim_data[:,4]
            E_Cov = Sim_data[:,7]

            f.write("{} {} {} {} {} {}\n".format(dir,distance,coverage,E_Bend[-1],E_Cov[-1],covStrength))
        
    f.close()

Read_data("../Results/Wrapping_two_fix_cov/")

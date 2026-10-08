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
            for Energy in d["Energies"]:
                if(Energy["Name"] == "Coverage"):
                    covStrength = Energy["constants"][0]
                    coverage = Energy["constants"][2]
            distance = d["theta"]

            # That is the distance and the target cov
            # I would love to also get the final coverage,
            # I now need to read the output fil
            output_file = folderpath + dir + "/Output_data.txt"
            Sim_data = np.loadtxt(output_file,skiprows=1)

            E_Bend = Sim_data[:,4]
            E_Cov = Sim_data[:,7]

            f.write("{} {} {} {} {} {}\n".format(dir,distance,coverage,E_Bend[-1],E_Cov[-1],covStrength))
        
    f.close()


def Read_final_energies(folderpath, outname="Final_energies.txt"):
    # Two bead simulations with the energy terms of the TwoBeadsFast/Cluster runs:
    # Bending_tan, Volume_constraint, Surface_tension, Bead (x2), Total_E.
    # Writes: directory, distance between beads (theta in the input file),
    # the final value of each energy term and the final total energy.
    expected = ["Bending_tan", "Volume_constraint", "Surface_tension", "Bead", "Bead", "Total_E"]
    dirs = sorted(d for d in os.listdir(folderpath) if os.path.isdir(folderpath + d))

    f = open(folderpath + outname, "w")
    f.write("# directory distance Bending_tan Volume_constraint Surface_tension Bead_0 Bead_1 Total_E\n")
    for dir in dirs:
        json_file = folderpath + dir + "/Input_file.json"
        output_file = folderpath + dir + "/Output_data.txt"
        if not (os.path.isfile(json_file) and os.path.isfile(output_file)):
            print("Skipping {} (no Input_file.json / Output_data.txt)".format(dir))
            continue

        with open(json_file) as fjson:
            d = json.load(fjson)
        distance = d["theta"]

        with open(output_file) as fout:
            header = fout.readline().split()
        # Columns: time step Volume Area <energies...> Total_E ...
        if header[4:10] != expected:
            print("Skipping {} (different energy terms: {})".format(dir, header[4:10]))
            continue

        Sim_data = np.loadtxt(output_file, skiprows=1)
        final = Sim_data[-1] if Sim_data.ndim > 1 else Sim_data

        print("Reading data from {}".format(dir))
        f.write("{} {} {}\n".format(dir, distance, " ".join("{:.8g}".format(e) for e in final[4:10])))
    f.close()


# Read_data("../Results/Wrapping_two_coverage_HighBend/")

Read_final_energies("../Results/TwoBeadsFast/")

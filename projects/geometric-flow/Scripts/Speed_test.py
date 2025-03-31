import numpy as np 
import matplotlib.pyplot as plt



def main_2():

    # Lets load the data 

    # There are a lot of different directories

    Bend = "../Results/Timings/Bending_1.0000_Nsim_2"
    Tension = "../Results/Timings/Surface_tension_1.0000_Nsim_2"
    Bend_H1 = "../Results/Timings/H1_Bending_1.0000_1.0000_1.0000_0.0000_Nsim_2"
    Tension_H1 = "../Results/Timings/H2_Surface_tension_1.0000_1.0000_1.0000_0.0000_Nsim_2"
    Bend_H2 = "../Results/Timings/H2_Bending_1.0000_1.0000_1.0000_0.0000_Nsim_2"
    Folder_list = [Bend, Tension, Bend_H1, Tension_H1, Bend_H2]
    label_list = ["Bending", "Surface tension", "H1 Bending", "H1 Surface tension","H2 Bend"]
    # Lets load the data from sim 2
    

    # Now there will be two differnt plots
    Energy_evol = []
    Grad_evol = []
    Sim_times = []



    for i in range(len(Folder_list)):
        # Ok lets load the data  
        folder = Folder_list[i]
        label = label_list[i]
        Output = np.loadtxt("{}/Output_data.txt".format(folder), skiprows=1)
        # Now that i have the output i care about the evolution of the ENergy and the gradient evol

        Energy_evol.append(Output[:,6])
        Grad_evol.append(Output[:,7])

        Timings = np.loadtxt("{}/Simulation_timings.txt".format(folder), skiprows=1)

        # sum = np.sum(Timings[:-1,0] )+np.sum(Timings[:-1,1] )+np.sum(Timings[:-1,2] )

        # print("The sum is {} and it should be {} ".format(sum, Timings[-1,0]))
        remeshing = np.sum(Timings[:,0])
        gradients = np.sum(Timings[:,1])
        backtracking = np.sum(Timings[:,2])
        construction = np.sum(Timings[:,3])
        compute = np.sum(Timings[:,4])
        solve = np.sum(Timings[:,5])

        gradients = gradients -construction - compute - solve
        
        Sim_times.append(Timings[-1,0])
        plt.title("Timings for {}".format(label))
        if(i < 2):
            plt.pie([remeshing, backtracking, gradients ], labels=["Remeshing", "Backtracking", "Gradients"])
        
        else:
            plt.pie([remeshing, backtracking, construction, compute, solve], labels=["Remeshing", "Backtracking", "Construction", "Compute", "Solve"])
        
        plt.savefig("{}/Pie_chart.png".format(folder))
        # plt.show()
        plt.clf()


    # Ok i have the pie charts
    # I need to plot the energy evolutions


    steps = [i for i in range(0,len(Energy_evol[0])*100,100)]   
    plt.plot(steps,Energy_evol[0],color = "purple", label="{}, {:.2f} minutes".format(label_list[0], Sim_times[0]/(60*1000)))
    plt.axhline(y=Energy_evol[2][-1], color='black', ls='dashed')
    
    steps = [i for i in range(0,len(Energy_evol[2])*100,100)]
    plt.plot(steps,Energy_evol[2],ls='dashed',color = "purple", label="{}, {:.2f} minutes".format(label_list[2], Sim_times[2]/(60*1000)))
    plt.plot(Energy_evol[4], ls ='dotted',color="black",label="{}, {:.2f} minutes".format(label_list[4], Sim_times[4]/(60*1000)))
    print(Energy_evol[4])
    plt.legend()
    plt.savefig("../Results/Timings_progress/Energy_evolution_bending.png")
    plt.show()


    plt.plot(Energy_evol[1],color = "purple", label="{}, {} miliseconds".format(label_list[1], Sim_times[1]))
    plt.axhline(y=Energy_evol[3][-1], color='black', ls='dashed')

    plt.plot(Energy_evol[3],ls='dashed',color = "purple", label="{}, {} miliseconds".format(label_list[3], Sim_times[3]))
    
    plt.legend()
    plt.savefig("../Results/Timings_progress/Energy_evolution_surface_tension.png")
    plt.show()






def main():

    # Lets load the data 
    basic_folder = "../Results/Timings_progress/"
    # There are a lot of different directories

    Bend_tension = basic_folder+"Surface_tension_1.0000_Bending_1.0000_Nsim_1"
    # Bend_tension_H1 = "../Results/Timings/Surface_tension_1.0000_H1_Bending_1.0000_1.0000_1.0000_0.0000_Nsim_2"
    # Bend_H1_tension = "../Results/Timings/H1_Surface_tension_1.0000_1.0000_1.0000_0.0000_Bending_1.0000_Nsim_2"
    Bend_H1_tension_H1 = basic_folder+"H1_Surface_tension_1.0000_1.0000_1.0000_0.0000_H1_Bending_1.0000_1.0000_1.0000_0.0000_Nsim_1"
    Bend_H2_tension_H2 = basic_folder+"H2_Surface_tension_1.0000_1.0000_1.0000_0.0000_H2_Bending_1.0000_1.0000_1.0000_0.0000_Nsim_1"

    # Folder_list = [Bend_tension, Bend_tension_H1, Bend_H1_tension, Bend_H1_tension_H1,Bend_H2_tension_H2]
    Folder_list = [Bend_tension,Bend_H1_tension_H1,Bend_H2_tension_H2]

    # label_list = ["Bend Tension", "Bend Tension H1", "Bend H1 Tension", "Bend H1 Tension H1", "Bend H2 Tension H2"]
    label_list = ["Bend Tension","Bend H1 Tension H1", "Bend H2 Tension H2"]

    # Lets load the data from sim 2
    

    # Now there will be two differnt plots
    Energy_evol = []
    Grad_evol = []
    Sim_times = []

    bar_remeshing = []
    bar_gradients = []
    bar_backtracking = []
    bar_construction = []
    bar_compute = []
    bar_solve = []



    for i in range(len(Folder_list)):
        # Ok lets load the data  
        folder = Folder_list[i]
        label = label_list[i]
        Output = np.loadtxt("{}/Output_data.txt".format(folder), skiprows=1)
        # Now that i have the output i care about the evolution of the ENergy and the gradient evol

        Energy_evol.append(Output[:,6])
        Grad_evol.append(Output[:,7])
        print("The folder is {}".format(folder) )
        Timings = np.loadtxt("{}/Simulation_timings.txt".format(folder), skiprows=1)

        # sum = np.sum(Timings[:-1,0] )+np.sum(Timings[:-1,1] )+np.sum(Timings[:-1,2] )

        # print("The sum is {} and it should be {} ".format(sum, Timings[-1,0]))
        remeshing = np.sum(Timings[:,0])
        gradients = np.sum(Timings[:,1])
        backtracking = np.sum(Timings[:,2])
        construction = np.sum(Timings[:,3])
        compute = np.sum(Timings[:,4])
        solve = np.sum(Timings[:,5])

        gradients = gradients -construction - compute - solve
        
        bar_remeshing.append(remeshing)
        bar_gradients.append(gradients)
        bar_backtracking.append(backtracking)
        bar_construction.append(construction)
        bar_compute.append(compute)
        bar_solve.append(solve)

        Sim_times.append(Timings[-1,0])
        plt.title("Timings for {}".format(label))
        if(i < 1):
            plt.pie([remeshing, backtracking, gradients ], labels=["Remeshing", "Backtracking", "Gradients"])
        
        else:
            plt.pie([remeshing, backtracking, construction, compute, solve], labels=["Remeshing", "Backtracking", "Construction", "Compute", "Solve"])
        
        plt.savefig("{}/Pie_chart_{}.png".format(basic_folder,label_list[i]))
        plt.show()
        plt.clf()


    # Ok i have the pie charts
    # I need to plot the energy evolutions

    # STEPPING

    steps = [i for i in range(0,len(Energy_evol[0]),1)]

    plt.plot(steps,Energy_evol[0],color = "purple", label="{}, {} miliseconds".format(label_list[0], Sim_times[0]))
    steps = [i for i in range(0,len(Energy_evol[1]),1)]
    plt.plot(steps,Energy_evol[1],color = "purple",ls = 'dashed', label="{}, {} miliseconds".format(label_list[1], Sim_times[1]))
    steps = [i for i in range(0,len(Energy_evol[2]),1)]
    plt.plot(steps,Energy_evol[2],ls='dashdot',color = "purple", label="{}, {} miliseconds".format(label_list[2], Sim_times[2]))
    # steps = [i for i in range(0,len(Energy_evol[3])*100,100)]
    # plt.plot(steps,Energy_evol[3],ls='dotted',color = "purple", label="{}, {} miliseconds".format(label_list[3], Sim_times[3]))
    # steps = [i for i in range(0,len(Energy_evol[4]),1)]
    # plt.plot(steps,Energy_evol[4],ls='dotted',color = "black", label="{}, {} miliseconds".format(label_list[4], Sim_times[4]))
    plt.xlabel(r"""Timestep """,fontsize=15.0)
    plt.ylabel(r"E",fontsize=15.0)

    plt.axhline(y=Energy_evol[2][-1], color='black', ls='dashed')

    
    plt.legend()
    plt.savefig("../Results/Timings_progress/Energy_evolution_double_2.png")
    plt.show()

    species = (
    "L2 gradient ",
    "H1 gradient",
    "H2 gradient$",
    )
    weight_counts = {
        "Remeshing": bar_remeshing,
        "Gradients": bar_gradients,
        "Backtracking": bar_backtracking,
        "Construction": bar_construction,
        "Compute": bar_compute,
        "Solve": bar_solve
    }

    width = 0.5

    fig, ax = plt.subplots()
    bottom = np.zeros(3)

    for boolean, weight_count in weight_counts.items():
        p = ax.bar(species, weight_count, width, label=boolean, bottom=bottom)
        bottom += weight_count

    ax.set_title("Time distribution per sim")
    ax.set_ylabel("time [ms]")
    ax.legend(loc="upper right")
    plt.savefig("../Results/Timings_progress/Bar_plot_sphere.png",bbox_inches='tight')
    plt.show()



def main_beads():
    folder_basic = "../Results/Timings_progress/"
    # Lets load the data 

    # There are a lot of different directories

    # Bend_Bead = folder_basic+"Bending_1.0000_Bead_radius_0.3000_str_400.0000_Nsim_2"
    # Bend_H1_Bead_H1 = folder_basic+"H1_Bending_1.0000_1.0000_1.0000_0.0000_H1_Bead_radius_0.3000_str_400.0000_Nsim_2"
    # Bend_H2_Bead_H2 = folder_basic+"H2_Surface_tension_1.0000_1.0000_1.0000_0.0000_H2_Bending_1.0000_1.0000_1.0000_0.0000_Nsim_2"

    Bend_Bead = folder_basic+"Bending_1.0000_Bead_radius_0.3000_str_400.0000_Nsim_1"
    Bend_H1_Bead_H1 = folder_basic+"H1_Bending_1.0000_1.0000_1.0000_0.0000_H1_Bead_radius_0.3000_str_400.0000_Nsim_2"
    Bend_H2_Bead_H2 = folder_basic+"H2_Bending_1.0000_1.0000_1.0000_0.0000_H2_Bead_radius_0.3000_str_400.0000_Nsim_2"

    Bend_Bead = folder_basic+"Bending_10.0000_Surface_tension_20.0000_Bead_radius_0.3000_str_400.0000_Nsim_3"
    Bend_H1_Bead_H1 = folder_basic+"H1_Bending_10.0000_1.0000_1.0000_0.0000_H1_Surface_tension_20.0000_1.0000_1.0000_0.0000_H1_Bead_radius_0.3000_str_400.0000_Nsim_3"
    Bend_H2_Bead_H2 = folder_basic+"H2_Bending_10.0000_1.0000_1.0000_0.0000_H2_Surface_tension_20.0000_1.0000_1.0000_0.0000_H2_Bead_radius_0.3000_str_400.0000_Nsim_3"



    Folder_list = [Bend_Bead,Bend_H1_Bead_H1,Bend_H2_Bead_H2]

    label_list = ["Bend Bead", "Bend Bead H1", "Bend Bead H2"]
    # Lets load the data from sim 2
    

    # Now there will be two differnt plots
    Energy_evol = []
    Grad_evol = []
    Sim_times = []

    bar_remeshing = []
    bar_gradients = []
    bar_backtracking = []
    bar_construction = []
    bar_compute = []
    bar_solve = []


    for i in range(len(Folder_list)):
        # Ok lets load the data  
        folder = Folder_list[i]
        label = label_list[i]
        Output = np.loadtxt("{}/Output_data.txt".format(folder), skiprows=1)
        # Now that i have the output i care about the evolution of the ENergy and the gradient evol

        Energy_evol.append(Output[:,8])
        Grad_evol.append(Output[:,7])
        print("The folder is {}".format(folder) )
        Timings = np.loadtxt("{}/Simulation_timings.txt".format(folder), skiprows=1)

        # sum = np.sum(Timings[:-1,0] )+np.sum(Timings[:-1,1] )+np.sum(Timings[:-1,2] )

        # print("The sum is {} and it should be {} ".format(sum, Timings[-1,0]))
        remeshing = np.sum(Timings[:,0])
        gradients = np.sum(Timings[:,1])
        backtracking = np.sum(Timings[:,2])
        construction = np.sum(Timings[:,3])
        compute = np.sum(Timings[:,4])
        solve = np.sum(Timings[:,5])

        gradients = gradients -construction - compute - solve
        
        bar_remeshing.append(remeshing)
        bar_gradients.append(gradients)
        bar_backtracking.append(backtracking)
        bar_construction.append(construction)
        bar_compute.append(compute)
        bar_solve.append(solve)


        time = np.sum(Timings)-Timings[-1,0]
        # print("Ttal time is  {}".format(time- Timings[-1,0]))
        Sim_times.append(time)
        plt.title("Timings for {}".format(label))
        if(i < 1):
            plt.pie([remeshing, backtracking, gradients ], labels=["Remeshing", "Backtracking", "Gradients"])
        
        else:
            plt.pie([remeshing, backtracking, construction, compute, solve], labels=["Remeshing", "Backtracking", "Construction", "Compute", "Solve"])
        
        plt.savefig("{}/Pie_chart_{}.png".format(folder_basic,i))
        # plt.show()
        plt.clf()


    # Ok i have the pie charts
    # I need to plot the energy evolutions

    # STEPPING

    steps = [i for i in range(0,len(Energy_evol[0])*100,100)]
    print(len(steps))
    print(len(Energy_evol[0]))

    plt.plot(steps,Energy_evol[0],color = "purple", label="{}, {:.2f} minutes".format(label_list[0], Sim_times[0]/(60*1000)))
    steps = [i for i in range(0,len(Energy_evol[1])*100,100)]
    plt.plot(steps,Energy_evol[1],color = "purple",ls = 'dashed', label="{}, {:.2f} minutes".format(label_list[1], Sim_times[1]/(60*1000)))
    steps = [i for i in range(0,len(Energy_evol[2])*100,100)]
    plt.plot(steps,Energy_evol[2],ls='dashdot',color = "purple", label="{}, {:.2f} minutes".format(label_list[2], Sim_times[2]/(60*1000)))
    plt.xlabel(r"""Timestep """,fontsize=15.0)
    plt.ylabel(r"E",fontsize=15.0)
    # plt.axhline(y=np.min( [Energy_evol[1][-1],Energy_evol[0][-1],Energy_evol[2][-1]]), color='black', ls='dashed')
    plt.axhline(y = np.min( [Energy_evol[i][-1] for i in range(len(Energy_evol))]), color='black', ls='dashed' )
    # plt.yscale("log")
    # plt.xlim(10000)
    plt.legend()
    plt.savefig("../Results/Timings_progress/Energy_evolution_bead_2.png")
    plt.show()
    plt.clf()


    # I want to do a barplot now

    species = (
    "L2 gradient ",
    "H1 gradient",
    "H2 gradient$",
    )
    weight_counts = {
        "Remeshing": bar_remeshing,
        "Gradients": bar_gradients,
        "Backtracking": bar_backtracking,
        "Construction": bar_construction,
        "Compute": bar_compute,
        "Solve": bar_solve
    }

    width = 0.5

    fig, ax = plt.subplots()
    bottom = np.zeros(3)

    for boolean, weight_count in weight_counts.items():
        p = ax.bar(species, weight_count, width, label=boolean, bottom=bottom)
        bottom += weight_count

    ax.set_title("Time distribution per sim")
    ax.set_ylabel("time [ms]")
    ax.legend(loc="upper right")
    plt.savefig("../Results/Timings_progress/Bar_plot_beads_2.png",bbox_inches='tight')

    plt.show()




# main_2()
main_beads()

main()





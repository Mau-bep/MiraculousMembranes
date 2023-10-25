import numpy as np
import matplotlib.pyplot as plt 
import matplotlib.animation as animation
import os
# import imageio
import cv2
import sys
import ffmpeg
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
import subprocess

#print(cv2.getBuildInformation())

nu=float(sys.argv[1])
KB=float(sys.argv[2])



def Get_img(file_path,img_path):
    X_coord=[]
    Y_coord=[]
    Z_coord=[]
    triangles=[]

    fig = plt.figure(figsize=(3*4.11,3*3.94))
    ax = fig.add_subplot(projection='3d')
    ax.elev=0
    ax.azim =0
    counter=1000
    with open(file_path,"r") as file:
        
        line = file.readline()
        while(line):
            
            line=file.readline()
            Arr=np.array( line.split(' '))
            if(Arr[0])=='v':
                X_coord.append(float(Arr[1]))
                Y_coord.append(float(Arr[2]))
                Z_coord.append(float(Arr[3]))
            
            if(Arr[0]=='f'):
                # if(counter==1000):
                #     print(line)
                #     print('got to faces, yei\n')
                #     counter=102
                ind1=int(Arr[1])-1
                ind2=int(Arr[2])-1
                ind3=int(Arr[3])-1

                triangles.append(((X_coord[ind1],Y_coord[ind1],Z_coord[ind1]),(X_coord[ind2],Y_coord[ind2],Z_coord[ind2]),((X_coord[ind3],Y_coord[ind3],Z_coord[ind3]))))
                

            if(Arr[0]=='#'):
                print('Got here no end')
                if(Arr[1]=='End'):
                    break
        
 


                


    
    ax.add_collection(Poly3DCollection(triangles,alpha=0.6,edgecolors='purple',facecolors='pink',fc='pink',lw=0.1))
    ax.set_xlim([-1,1])
    ax.set_ylim([-0.6,0.6])
    ax.set_zlim([-1.2,1.2])
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_zlabel('z')
    ax.set_aspect('equal')
    plt.savefig(img_path,bbox_inches='tight')
    
    #fig.canvas.draw()
    
    plt.close()



def Get_img_2(file_path,img_path):
    X_coord=[]
    Y_coord=[]
    Z_coord=[]
    triangles=[]
    plt.clf()
    # fig = plt.figure(figsize=(3*4.11,3*3.94))
    ax = fig.add_subplot(projection='3d')
    ax.elev=0
    ax.azim =0
    counter=1000
    with open(file_path,"r") as file:
        
        line = file.readline()
        while(line):
            
            line=file.readline()
            Arr=np.array( line.split(' '))
            if(Arr[0])=='v':
                X_coord.append(float(Arr[1]))
                Y_coord.append(float(Arr[2]))
                Z_coord.append(float(Arr[3]))
            
            if(Arr[0]=='f'):
                # if(counter==1000):
                #     print(line)
                #     print('got to faces, yei\n')
                #     counter=102
                ind1=int(Arr[1])-1
                ind2=int(Arr[2])-1
                ind3=int(Arr[3])-1

                triangles.append(((X_coord[ind1],Y_coord[ind1],Z_coord[ind1]),(X_coord[ind2],Y_coord[ind2],Z_coord[ind2]),((X_coord[ind3],Y_coord[ind3],Z_coord[ind3]))))
                

            if(Arr[0]=='#'):
                print('Got here no end')
                if(Arr[1]=='End'):
                    break
        
 


                


    
    ax.add_collection(Poly3DCollection(triangles,alpha=0.6,edgecolors='purple',facecolors='pink',fc='pink',lw=0.1))
    ax.set_xlim([-1,1])
    ax.set_ylim([-0.6,0.6])
    ax.set_zlim([-1.2,1.2])
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_zlabel('z')
    ax.set_aspect('equal')
    plt.savefig(img_path,bbox_inches='tight')
    
    #fig.canvas.draw()
    
    # plt.close()

def create_gif(directory_path, gif_name, duration):
    # Initialize variables
    images = []
    image_idx = 0
    file_path = os.path.join(directory_path, f"{image_idx}.obj")
    print(file_path)
    image_path= os.path.join(directory_path,"image_side.jpg")
    # Read the images one by one until there are no more
    counter=0
    
    while os.path.exists(file_path):
        
        # print('Getting an image\n')
        Get_img(file_path,image_path)

        # print('DOne with img\n')
        counter+=1
        
        frame = cv2.imread(image_path)
        
        if frame is not None:
            frame = cv2.resize(frame, (frame_width, frame_height))
            video_writer.write(frame)
        image_idx += 10000
        file_path = os.path.join(directory_path, f"{image_idx}.obj")
        # print(file_path)






pre_folder='Mem3DG_IMG_correct/Curv_adap_0.10Min_rel_length_0.50/'
#pre_folder='Mem3DG_IMG_blood_correct_ini/Curv_adap_0.10Min_rel_length_0.50/'
#pre_folder='Tests_cil_regular/'
dir='nu_{:.3f}_c0_0.000_KA_10.000_KB_{:.6f}'.format(nu,KB)
folder_path= pre_folder+dir+'/'
# filename='100000.obj'
output_video = pre_folder+dir+".mp4"
output_file = pre_folder+dir+"_2.mp4"





fps = 3
frame_width = 640
frame_height = 480
fourcc = cv2.VideoWriter_fourcc(*'avc1')
video_writer = cv2.VideoWriter(output_video, fourcc, fps, (frame_width, frame_height))
print('Now creating gif')

# Uncomment when ready
create_gif(folder_path,'Animation_overmeshed_short.mp4',0.1)
video_writer.release()



#directory_path=folder_path

def Animate(frame):
    

    # Initialize variables
    images = []
    image_idx = frame
    file_path = os.path.join(directory_path, f"{image_idx}.obj")
    # print(file_path)
    image_path= os.path.join(directory_path,"image_side.jpg")
    # Read the images one by one until there are no more
    counter=0
    Get_img_2(file_path,image_path)




#fig = plt.figure(figsize=(3*4.11,3*3.94))

#step_ini=0
#dt=10000

#counter=0
#file_path = os.path.join(directory_path, f"{step_ini}.obj")
#print(file_path)
#while os.path.exists(file_path):
#    print('here yei')
#    step_fin=step_ini+dt*counter
#    file_path = os.path.join(directory_path, f"{step_fin+dt}.obj")
#    counter+=1

#interval=50
#step_fin=20000

#ani=animation.FuncAnimation(fig,Animate,np.arange(step_ini,step_fin,dt),interval=interval)
#f=r"/nfs/scistore16/wojtgrp/mrojasve/Cluster_Folders/projects/geometric-flow/build/"+output_file
#ani.save(f,writer="ffmpeg")













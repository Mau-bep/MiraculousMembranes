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
strg=float(sys.argv[3])




def Get_img(file_path,img_path):
    X_coord=[]
    Y_coord=[]
    Z_coord=[]
    triangles=[]

    fig = plt.figure(figsize=(3*4.11,3*3.94))
    ax = fig.add_subplot(projection='3d')
    ax.elev=-90
    ax.azim =90
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
        
 

    splitted_line=Bead_file.readline().split(' ')
    Positions=np.array( [splitted_line[0],splitted_line[1],splitted_line[3]],dtype=float)
    # print(Positions)


    # x = r * np.outer(np.cos(u), np.sin(v)) + Positions[0]
    # y = r * np.outer(np.sin(u), np.sin(v)) + Positions[1]
    # z = r * np.outer(np.ones(np.size(u)), np.cos(v)) + Positions[2]
    
    # ax.plot_surface(x,y,z,color='blue',alpha=0.5)
    
    ax.add_collection(Poly3DCollection(triangles,alpha=1.0,edgecolors='purple',facecolors='pink',fc='pink',lw=0.1))
    ax.set_xlim([-6,6])
    ax.set_ylim([-6,6])
    ax.set_zlim([-6,6])
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_zlabel('z')
    ax.set_aspect('equal')
    plt.savefig(img_path,bbox_inches='tight')
    
    #fig.canvas.draw()
    
    plt.close()

    #img = np.frombuffer(fig.canvas.tostring_rgb(),dtype=np.unit8)
    #img = img.reshape(fig.canvas.get_idth_height()[::-1]+(3,))
    # plt.show()
    #process.communicate(input=img.tobytes())

# Get_img(file_path,'./Mesh_screenshot.jpg')



# def Get_img_2(file_path,img_path):
#     X_coord=[]
#     Y_coord=[]
#     Z_coord=[]
#     triangles=[]
#     plt.clf()
#     # fig = plt.figure(figsize=(3*4.11,3*3.94))
#     ax = fig.add_subplot(projection='3d')
#     ax.elev=0
#     ax.azim =0
#     counter=1000
#     with open(file_path,"r") as file:
        
#         line = file.readline()
#         while(line):
            
#             line=file.readline()
#             Arr=np.array( line.split(' '))
#             if(Arr[0])=='v':
#                 X_coord.append(float(Arr[1]))
#                 Y_coord.append(float(Arr[2]))
#                 Z_coord.append(float(Arr[3]))
            
#             if(Arr[0]=='f'):
#                 # if(counter==1000):
#                 #     print(line)
#                 #     print('got to faces, yei\n')
#                 #     counter=102
#                 ind1=int(Arr[1])-1
#                 ind2=int(Arr[2])-1
#                 ind3=int(Arr[3])-1

#                 triangles.append(((X_coord[ind1],Y_coord[ind1],Z_coord[ind1]),(X_coord[ind2],Y_coord[ind2],Z_coord[ind2]),((X_coord[ind3],Y_coord[ind3],Z_coord[ind3]))))
                

#             if(Arr[0]=='#'):
#                 print('Got here no end')
#                 if(Arr[1]=='End'):
#                     break
        
 



    
#     ax.add_collection(Poly3DCollection(triangles,alpha=0.6,edgecolors='purple',facecolors='pink',fc='pink',lw=0.1))
#     ax.set_xlim([-1,1])
#     ax.set_ylim([-0.6,0.6])
#     ax.set_zlim([-1.2,1.2])
#     ax.set_xlabel('x')
#     ax.set_ylabel('y')
#     ax.set_zlabel('z')
#     ax.set_aspect('equal')
#     plt.savefig(img_path,bbox_inches='tight')
    
#     #fig.canvas.draw()
    
#     # plt.close()

def create_gif(directory_path, duration):
    # Initialize variables
    images = []
    image_idx = 0
    file_path = os.path.join(directory_path, f"{image_idx}.obj")
    # print(file_path)
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






pre_folder='../Results/Mem3DG_Beads/Curv_adap_0.10Min_rel_length_0.50/'
#pre_folder='Mem3DG_IMG_blood_correct_ini/Curv_adap_0.10Min_rel_length_0.50/'
#pre_folder='Tests_cil_regular/'
dir='nu_{:.3f}_c0_0.000_KA_10.000_KB_{:.6f}_Inter_{:.6f}'.format(nu,KB,strg)
folder_path= pre_folder+dir+'/'



Bead_path=folder_path+"Bead_data.txt"


Bead_file = open(Bead_path)
line = Bead_file.readline()
# print('THe first line is {}'.format(line) )
r = 2.0
u = np.linspace(0, 2 * np.pi, 100)
v = np.linspace(0, np.pi, 100)


mit_bead=False
if(mit_bead):
    output_video = pre_folder+dir+".mp4"
else:
    output_video = pre_folder+dir+"nobead_2.mp4"




fps = 3
frame_width = 640
frame_height = 480
fourcc = cv2.VideoWriter_fourcc(*'avc1')
video_writer = cv2.VideoWriter(output_video, fourcc, fps, (frame_width, frame_height))
print('Now creating gif')

create_gif(folder_path,0.1)
video_writer.release()














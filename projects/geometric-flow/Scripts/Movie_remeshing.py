import numpy as np
import os 
import matplotlib.pyplot as plt 
from matplotlib import rc
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
## for Palatino and other serif fonts use:
rc('font',**{'family':'serif','serif':['Palatino'], 'size':14})
rc('text', usetex=True)
# import sys
import cv2
from PIL import Image





def generate_video():
    image_folder = "../Results/Debug_remesh_trial/Bending_1.0000_Bead_radius_0.4000_str_400.0000_Nsim_1/Snapshots/"

    image_folder = "../Results/Timings_progress/Snap_sphere/"

    video_name = '../Results/Timings_progress/Sphere_comparison.avi'

    images = [img for img in os.listdir(image_folder) if img.endswith(".png")]
    print("Images:",images)

    frame = cv2.imread(os.path.join(image_folder, images[0]))
    height, width, layers = frame.shape

    video = cv2.VideoWriter(video_name,cv2.VideoWriter_fourcc(*'DIVX'), 1, (width,height))
    for i in range(1,len(images)+1):
        video.write(cv2.imread(os.path.join(image_folder, "{}.png".format(i))))
    # for image in images:
        # video.write(cv2.imread(os.path.join(image_folder, image)))

    video.release()
    cv2.destroyAllWindows()
    print("Video generated")


generate_video()
                    
                


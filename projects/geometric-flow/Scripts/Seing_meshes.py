import pyvista
import matplotlib.pyplot as plt
import os
import imageio
import cv2
import sys

nu=float(sys.argv[1])
KB=float(sys.argv[2])

def create_gif(directory_path, gif_name, duration):
    # Initialize variables
    images = []
    image_idx = 0
    file_path = os.path.join(directory_path, f"{image_idx}.obj")

    # Read the images one by one until there are no more
    counter=0
    while os.path.exists(file_path):
        plotter = pyvista.Plotter(off_screen=True)
        reader = pyvista.get_reader(file_path)
        mesh = reader.read()
        
        plotter.add_mesh(mesh,color='pink',show_edges=True,edge_color='black')
        plotter.camera_position = 'yz'
        plotter.camera.azimuth= 45
        plotter.camera.elevation = 45

        image_path=os.path.join(directory_path, "image_side.jpg")
        plotter.show(screenshot=image_path)
        # plt.imshow(plotter.image)
        # plt.show()


        # Load the image and append it to the list
        # images.append(imageio.imread(image_path))
        counter+=1
        # Increment the image index and update the file path
        
        
        frame = cv2.imread(image_path)
        
        if frame is not None:
            frame = cv2.resize(frame, (frame_width, frame_height))
            video_writer.write(frame)
        image_idx += 10000
        file_path = os.path.join(directory_path, f"{image_idx}.obj")
        



# nu=0.8
pre_folder='projects/geometric-flow/build/Mem3DG_IMG_parallel/Curv_adap_0.10Min_rel_length_0.50/'
dir='nu_{:.3f}_c0_0.000_KA_11.000_KB_{:.6f}'.format(nu,KB)
folder_path= pre_folder+dir+'/'
# filename='100000.obj'




output_video = pre_folder+dir+".mp4"
fps = 3
frame_width = 1702
frame_height = 835
fourcc = cv2.VideoWriter_fourcc(*"mp4v")
video_writer = cv2.VideoWriter(output_video, fourcc, fps, (frame_width, frame_height))

create_gif(folder_path,'Animation_overmeshed_short.mp4',0.1)


video_writer.release()
cv2.destroyAllWindows()
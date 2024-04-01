import numpy as np
import os 
import matplotlib.pyplot as plt 
from matplotlib import rc
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
## for Palatino and other serif fonts use:
rc('font',**{'family':'serif','serif':['Palatino'], 'size':14})
rc('text', usetex=True)
import sys

import movis as mv


Movie_dir="../Results/Videos_poster/"



title = mv.layer.Composition(size=(1920, 1080), duration=2.5)
title.add_layer(mv.layer.Rectangle(title.size, color="black"))  # Set background
pos = title.size[0] // 2, title.size[1] // 2
title.add_layer(
    mv.layer.Text(r'Cell shape through the eyes of geometry processing', font_size=50, font_family='sans-serif', color='#ffffff'),
    name='text',  # The layer item can be accessed by name
    offset=0.5,  # Show the text after one second
    position=pos,  # The layer is centered by default, but it can also be specified explicitly
    anchor_point=(0.0, 0.0),
    opacity=1.0, scale=1.0, rotation=0.0,  # anchor point, opacity, scale, and rotation are also supported
    blending_mode='normal')  # Blending mode can be specified for each layer.


canham = mv.layer.Composition(size=(1920, 1080), duration=2.5)
canham.add_layer(mv.layer.Rectangle(canham.size, color="black"))  # Set background
pos = canham.size[0] // 2, canham.size[1] // 2
canham.add_layer(
    mv.layer.Text(r'Canham-Helfrich Hamiltonian', font_size=50, font_family='sans-serif', color='#ffffff'),
    name='text',  # The layer item can be accessed by name
    offset=0.5,  # Show the text after one second
    position=pos,  # The layer is centered by default, but it can also be specified explicitly
    anchor_point=(0.0, 0.0),
    opacity=1.0, scale=1.0, rotation=0.0,  # anchor point, opacity, scale, and rotation are also supported
    blending_mode='normal')  # Blending mode can be specified for each layer.


scene = mv.layer.Composition(size=(1920, 1080), duration=2.0)
scene.add_layer(mv.layer.Rectangle(scene.size, color="black"))  # Set background
pos = scene.size[0] // 2, scene.size[1] // 2
scene.add_layer(
    mv.layer.Text(r'Prolate case v=0.8', font_size=50, font_family='sans-serif', color='#ffffff'),
    name='text',  # The layer item can be accessed by name
    offset=0.0,  # Show the text after one second
    position=pos,  # The layer is centered by default, but it can also be specified explicitly
    anchor_point=(0.0, 0.0),
    opacity=1.0, scale=1.0, rotation=0.0,  # anchor point, opacity, scale, and rotation are also supported
    blending_mode='normal')  # Blending mode can be specified for each layer.



scene2 = mv.layer.Composition(size=(1920, 1080), duration=2.0)
scene2.add_layer(mv.layer.Rectangle(scene.size, color="black"))  # Set background
pos = scene2.size[0] // 2, scene.size[1] // 2
scene2.add_layer(
    mv.layer.Text(r'Oblate case v=0.625', font_size=50, font_family='sans-serif', color='#ffffff'),
    name='text',  # The layer item can be accessed by name
    offset=0.0,  # Show the text after one second
    position=pos,  # The layer is centered by default, but it can also be specified explicitly
    anchor_point=(0.0, 0.0),
    opacity=1.0, scale=1.0, rotation=0.0,  # anchor point, opacity, scale, and rotation are also supported
    blending_mode='normal')  # Blending mode can be specified for each layer.



scene3 = mv.layer.Composition(size=(1920, 1080), duration=2.0)
scene3.add_layer(mv.layer.Rectangle(scene.size, color="black"))  # Set background
pos = scene.size[0] // 2, scene.size[1] // 2
scene3.add_layer(
    mv.layer.Text(r'Stomatocyte case v=0.58', font_size=50, font_family='sans-serif', color='#ffffff'),
    name='text',  # The layer item can be accessed by name
    offset=0.0,  # Show the text after one second
    position=pos,  # The layer is centered by default, but it can also be specified explicitly
    anchor_point=(0.0, 0.0),
    opacity=1.0, scale=1.0, rotation=0.0,  # anchor point, opacity, scale, and rotation are also supported
    blending_mode='normal')  # Blending mode can be specified for each layer.


scene4 = mv.layer.Composition(size=(1920, 1080), duration=2.0)
scene4.add_layer(mv.layer.Rectangle(scene.size, color="black"))  # Set background
pos = scene4.size[0] // 2, scene.size[1] // 2
scene4.add_layer(
    mv.layer.Text(r'Membrane wrapping', font_size=50, font_family='sans-serif', color='#ffffff'),
    name='text',  # The layer item can be accessed by name
    offset=0.0,  # Show the text after one second
    position=pos,  # The layer is centered by default, but it can also be specified explicitly
    anchor_point=(0.0, 0.0),
    opacity=1.0, scale=1.0, rotation=0.0,  # anchor point, opacity, scale, and rotation are also supported
    blending_mode='normal')  # Blending mode can be specified for each layer.


scene4_5 = mv.layer.Composition(size=(1920, 1080), duration=2.0)
scene4_5.add_layer(mv.layer.Rectangle(scene.size, color="black"))  # Set background
pos = scene.size[0] // 2, scene.size[1] // 2
scene4_5.add_layer(
    mv.layer.Text(r'No wrapping', font_size=50, font_family='sans-serif', color='#ffffff'),
    name='text',  # The layer item can be accessed by name
    offset=0.0,  # Show the text after one second
    position=pos,  # The layer is centered by default, but it can also be specified explicitly
    anchor_point=(0.0, 0.0),
    opacity=1.0, scale=1.0, rotation=0.0,  # anchor point, opacity, scale, and rotation are also supported
    blending_mode='normal')  # Blending mode can be specified for each layer.
 
scene5 = mv.layer.Composition(size=(1920, 1080), duration=2.0)
scene5.add_layer(mv.layer.Rectangle(scene.size, color="black"))  # Set background
pos = scene.size[0] // 2, scene.size[1] // 2
scene5.add_layer(
    mv.layer.Text(r'Partial wrapping', font_size=50, font_family='sans-serif', color='#ffffff'),
    name='text',  # The layer item can be accessed by name
    offset=0.0,  # Show the text after one second
    position=pos,  # The layer is centered by default, but it can also be specified explicitly
    anchor_point=(0.0, 0.0),
    opacity=1.0, scale=1.0, rotation=0.0,  # anchor point, opacity, scale, and rotation are also supported
    blending_mode='normal')  # Blending mode can be specified for each layer.


scene6 = mv.layer.Composition(size=(1920, 1080), duration=2.0)
scene6.add_layer(mv.layer.Rectangle(scene.size, color="black"))  # Set background
pos = scene.size[0] // 2, scene.size[1] // 2
scene6.add_layer(
    mv.layer.Text(r'Complete wrapping', font_size=50, font_family='sans-serif', color='#ffffff'),
    name='text',  # The layer item can be accessed by name
    offset=0.0,  # Show the text after one second
    position=pos,  # The layer is centered by default, but it can also be specified explicitly
    anchor_point=(0.0, 0.0),
    opacity=1.0, scale=1.0, rotation=0.0,  # anchor point, opacity, scale, and rotation are also supported
    blending_mode='normal')  # Blending mode can be specified for each layer.

prolate = mv.layer.Video(Movie_dir+"Prolate.avi")
oblate = mv.layer.Video(Movie_dir+"Redbloodcel.avi")
stomatocyte = mv.layer.Video(Movie_dir+"STOMATOCYTE.avi" )
partial_wrapping = mv.layer.Video(Movie_dir+"Wrapping_KA_100_1_twotry0001-0500.mp4")
complete_wrapping = mv.layer.Video(Movie_dir+"Wrapping_KA_100_4_twotry0001-0500.mp4")
no_wrapping = mv.layer.Video(Movie_dir+"Wrapping_KA_50_0.01_20001-0226.mp4")

scene = mv.concatenate([title,canham,scene, prolate,scene2,oblate,scene3,stomatocyte,scene4,scene4_5,no_wrapping,scene5,partial_wrapping,scene6,complete_wrapping])




scene.write_video(Movie_dir+'Poster_vid.mp4')

# from moviepy.editor import *




# clip1 = VideoFileClip(Movie_dir+"Prolate.avi")
# clip0 = VideoFileClip(Movie_dir+"Letters.avi")
# txt_clip1 = TextClip(r"Prolate case $(v=0.8)$",fontsize=40,font='EB-Garamond-12-Bold',color='white')
# txt_clip1= txt_clip1.set_pos('center').set_duration(2)
# txt_clip1.fps=24
# video = concatenate_videoclips([txt_clip1])
# # print(TextClip.list('font'))

# video.write_videofile(Movie_dir+"lettrs.avi",codec='png')
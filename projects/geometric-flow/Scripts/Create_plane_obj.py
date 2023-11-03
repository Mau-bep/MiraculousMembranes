import numpy as np

import sys 
import math 
import os 
import matplotlib.pyplot as plt


from matplotlib.collections import PatchCollection
from matplotlib.patches import Circle, Polygon, Wedge


N_base_triangles=int(sys.argv[1])

Lx=10
Ly=10*np.sqrt(3)/2



os.makedirs('../Meshes/',exist_ok=True)

# I need to calculate a few things
N_vertices_x=N_base_triangles+1
N_vertices_y=N_vertices_x+1


triangle_edge=Lx/N_base_triangles
print('Triangle sides are {}'.format(triangle_edge))
N_vertices= int(N_vertices_x+N_vertices_y)* math.ceil(N_base_triangles/2)
# For the ith vertex the coordinates are 
Vert_per_full_row=N_vertices_x+N_vertices_y
X_points=[]
Y_points=[]
Z_points=[-2.0]+N_vertices*[0.0]
X_points=[5.0]
Y_points=[5.0]


for i in range(N_vertices):
    # I have the number of vertices
    # I want their position now 
    Full_row_index=int((i)/Vert_per_full_row) # We add the -1 cause the indices will be from 0 to vertex_per_full_row-1
    Up= ((i%Vert_per_full_row)-N_vertices_x)>=0
     #Here we are determining if the vertex is in the bottom or upper row 
    # print('Full row index is {} and Up is {}'.format(Full_row_index,Up))
    Down =not(Up)
    Up_col_index=(i%Vert_per_full_row)-N_vertices_x
    x_coord=triangle_edge*Up*(    ((Up_col_index-1)>0)*(Up_col_index-1)  + (Up_col_index%(N_vertices_y-1)!=0)*0.5         )                                        
    x_coord+=triangle_edge*Down*(  (i%Vert_per_full_row)  )
    y_coord=2*triangle_edge*(Full_row_index)*np.sqrt(3)/2+Up*np.sqrt(3)/2*triangle_edge
    X_points.append(x_coord)
    Y_points.append(y_coord)

plot=False

if(plot):
    plt.axhline(0,color='purple')
    plt.axhline(10,color='purple')
    plt.axvline(0,color='purple')
    plt.axvline(10,color='purple')
    # plt.sc('equal')
    plt.scatter(X_points,Y_points,color='black')
    plt.gca().set_aspect('equal')
    plt.show()
    plt.clf()

        # If up is correct is in the upper part of the high row





# Lets fo for the base triangles

Triangles=[]


vert_index=0
# I need to decide if i add 6 or seven
add_six=True
for j in range(N_base_triangles):

    # print(vert_index)
    if(add_six):
        # print("This one ")
        # print([0,1+vert_index,1+vert_index+N_vertices_x])
        Triangles.append([0,1+vert_index,1+vert_index+N_vertices_x])
        # print([0,1+vert_index+N_vertices_x-1,1+vert_index+N_vertices_x+N_vertices_x-1])
        Triangles.append([0,1+vert_index+N_vertices_x-1,1+vert_index+N_vertices_x+N_vertices_x])
    else:
        # print("Now this one")
        # print([0,1+vert_index,1+vert_index+N_vertices_y])
        Triangles.append([0,1+vert_index,1+vert_index+N_vertices_y])
        # print([0,1+vert_index+N_vertices_y,1+vert_index+N_vertices_y+N_vertices_y])
        Triangles.append([0,1+vert_index+N_vertices_y-1,1+vert_index+N_vertices_y+N_vertices_y-2])

    
    if((j%2)==0):
        add_six=False
        vert_index+=N_vertices_x
    else:
        add_six=True
        vert_index+=N_vertices_y
    # print((add_six)*6)W

# Now i need to add the verticals
# Triangles=[]
for j in range(N_base_triangles):
    # print("Cuenteme")
    Triangles.append([0,1+j,1+j+1])


for j in range(N_base_triangles+1):
    print("Cuenteme")
    Triangles.append([0,j+N_vertices-N_vertices_y+1,j+N_vertices-N_vertices_y+2])



print(len(Triangles))



print(len(X_points))
print(len(Triangles))
print(max(Triangles))

Collection=[]
for Triangle in Triangles:
    # print(len(Triangle))
    # print('Here we go')
    polygon = Polygon([ [X_points[Triangle[0]],Y_points[Triangle[0]]],[X_points[Triangle[1]],Y_points[Triangle[1]]],[X_points[Triangle[2]],Y_points[Triangle[2]]]], closed=True)
    Collection.append(polygon)

# plot=False
if(plot):
    fig, ax = plt.subplots()
    colors = 100 * np.random.rand(len(Collection))
    p = PatchCollection(Collection, alpha=0.4)
    p.set_array(colors)
    ax.add_collection(p)
    fig.colorbar(p, ax=ax)
    ax.set_xlim(0,10)
    ax.set_ylim(0,10)


    plt.show()




for i in range(N_vertices):
    if(i==N_vertices-N_vertices_y):
        break 
    Up= ((i%Vert_per_full_row)-N_vertices_x)>=0
    Down=not(Up)
    if(Up):
        # print('WHen its up then its up')
        if(i%Vert_per_full_row==N_vertices_x):
            # This is for the first triangle
            Triangles.append([1+i,1+i+1,1+i+N_vertices_y])
            continue
        if((i+1)%Vert_per_full_row==0):
            continue
        
        Triangles.append([1+i,1+i+N_vertices_x,1+i+N_vertices_x+1])
        Triangles.append([1+i,1+i+1,1+i+N_vertices_y])
    else:
        # I add the triangle that is up
        Triangles.append([1+i,1+i+N_vertices_x,1+i+N_vertices_x+1])
        if( not(i%Vert_per_full_row==N_vertices_x-1)):
            # print('Add the triangle to the right')
            Triangles.append([1+i,1+i+1,1+i+1+N_vertices_x])
        
Collection=[]
for Triangle in Triangles:
    if(Triangle[0]==0):
        continue
    # print('Here we go')
    polygon = Polygon([ [X_points[Triangle[0]],Y_points[Triangle[0]]],[X_points[Triangle[1]],Y_points[Triangle[1]]],[X_points[Triangle[2]],Y_points[Triangle[2]]]], closed=True)
    Collection.append(polygon)



if(plot):

    fig, ax = plt.subplots()
    colors = 100 * np.random.rand(len(Collection))
    p = PatchCollection(Collection, alpha=0.4)
    p.set_array(colors)
    ax.add_collection(p)
    fig.colorbar(p, ax=ax)
    ax.set_xlim(0,10)
    ax.set_ylim(0,10)


    plt.show()


file=open('../Meshes/Plane_{}_base_triangles.obj'.format(N_base_triangles),'w')


file.write("# This is a mesh file with a {}x{} plane with {} base triangles\n".format(Lx,Ly,N_base_triangles))



for i in range(len(X_points)):
    #
    file.write("v {} {} {}\n".format(X_points[i],Y_points[i],Z_points[i]))


for triangle in Triangles:
    file.write("f {} {} {}\n".format(triangle[0]+1,triangle[1]+1,triangle[2]+1))






















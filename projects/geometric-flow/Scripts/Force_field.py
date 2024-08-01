import numpy as np 
import matplotlib.pyplot as plt 







# plt.xlim(-2,2)
# plt.ylim(-2,2)
# plt.arrow(0,0,1,1, head_width=0.1,length_includes_head=True,color='black')
# plt.show()
def Field_1():
    R = 4.0 
    theta = np.linspace(0,2*np.pi ,50)
    X = R*np.cos(theta)
    Y = R*np.sin(theta)

    factors = []
    Bead = np.array([6,0])
    for i in range(len(X)):

        # I want to draw an arrow for each of these.
        dist_x = Bead[0]-X[i]
        dist_y = Bead[1]-Y[i]
        norm = np.sqrt(dist_x**2+dist_y**2)
        dist_x = dist_x/norm
        dist_y = dist_y/norm
        
        # print("La distancia es {} luego la normal apunta hacia {} {} y el producto punto es {}".format(norm,dist_x,dist_y,))
        # Ok i have the distance
        # I want to evaluate the dot product in the normal direction first.
        factor = (dist_x*X[i]/(R*norm**2) + dist_y*Y[i]/(R*norm**2))
        print("La distancia es {} luego la normal apunta hacia {} {} y el producto punto es {}".format(norm,dist_x,dist_y,factor))
        
        plt.arrow( X[i],Y[i],factor*X[i]/R,factor*Y[i]/R, head_width=0.1,length_includes_head=True,color='black')

    plt.axis('equal')
    plt.show()


def Force_field_2():


    d = np.array([3.0,0.0])
    X, Y = np.meshgrid(np.linspace(-5, 5, 1024), np.linspace(-5, 5, 1024))
    Z = ( (X*d[0]+ Y*d[1]) /( np.sqrt((X-d[0])**2+(Y-d[1])**2  ))) 
    levels = np.geomspace(np.min(Z)-np.min(Z)+1, np.max(Z)-np.min(Z)+1, 100)

    print(levels)
    print(np.min(Z))
    print(np.max(Z))

    fig, ax = plt.subplots()

    ax.contourf(X, Y, Z, levels=levels)

    plt.show()

Force_field_2()
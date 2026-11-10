"""Reduced area difference da = M / (4 pi R0) of closed triangle meshes (OBJ), the quantity the Bilayer energy
works with, computed independently of the C++ code (numpy only).

    python3 area_difference.py membrane_0.obj Final_state.obj

Prints the area, the reduced volume v = V / (4 pi/3 R0^3), da and the sorted extents of every file. M = 1/2 sum_e l_e theta_e,
R0 = sqrt(A / 4 pi); outward normals give da = 1 for a sphere."""
import sys, numpy as np
def load(p):
    V=[];F=[]
    for l in open(p):
        t=l.split()
        if t and t[0]=='v': V.append([float(x) for x in t[1:4]])
        if t and t[0]=='f': F.append([int(x.split('/')[0])-1 for x in t[1:4]])
    return np.array(V),np.array(F)
def quantities(V,F):
    P=V[F]; n=np.cross(P[:,1]-P[:,0],P[:,2]-P[:,0]); a=0.5*np.linalg.norm(n,axis=1); n=n/(2*a[:,None])
    A=a.sum(); Vol=np.einsum('ij,ij->',P[:,0],np.cross(P[:,1],P[:,2]))/6
    # edges -> adjacent faces
    ed={}
    for fi,f in enumerate(F):
        for k in range(3):
            i,j=f[k],f[(k+1)%3]; ed.setdefault((min(i,j),max(i,j)),[]).append((fi,i,j))
    M=0
    for (i,j),fs in ed.items():
        if len(fs)!=2: continue
        (f1,a1,b1),(f2,_,_)=fs
        e=V[b1]-V[a1]; l=np.linalg.norm(e); e/=l
        th=np.arctan2(e@np.cross(n[f1],n[f2]),n[f1]@n[f2])
        M+=0.5*l*th
    R0=np.sqrt(A/4/np.pi)
    return A,Vol,M/(4*np.pi*R0),Vol/(4*np.pi/3*R0**3)
for p in sys.argv[1:]:
    V,F=load(p); A,Vol,da,v=quantities(V,F)
    ext=V.max(0)-V.min(0)
    print(p.split('/')[-3],"A=%.4f v=%.4f da=%.4f extents=%s"%(A,v,da,np.round(np.sort(ext),3)))

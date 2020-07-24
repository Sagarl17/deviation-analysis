import os
import sys
import json
import math
import laspy
import numpy as np
from scipy.spatial import ConvexHull,KDTree
from shapely.geometry import Polygon,Point,LineString


def distance(a,b):
    return math.sqrt((a[0]-b[0])*(a[0]-b[0])+(a[1]-b[1])*(a[1]-b[1]))



def pnt_to_coords(proj_mat,x,y,z):
    """ 
    Converts local coordinates to global coordinates with the help of 
    transformation matrix and gives the coordinates as the result.This equation is 
    just simple form of matrix multiplication of Proj matrix * inverse([x,y,z,1])
    """

    w=1/(proj_mat[12]*x+proj_mat[13]*y+proj_mat[14]*z+proj_mat[15])

    X=(proj_mat[0]*x+proj_mat[1]*y+proj_mat[2]*z+proj_mat[3])*w
    Y=(proj_mat[4]*x+proj_mat[5]*y+proj_mat[6]*z+proj_mat[7])*w
    Z=(proj_mat[8]*x+proj_mat[9]*y+proj_mat[10]*z+proj_mat[11])*w

    return X,Y,Z



def representative_point(points):
    representative_array=[]
    point_3d=np.vstack([points.x,points.y,points.z]).T
    clusters=points.cluster_id
    max_clusterid=int(np.amax(points.cluster_id))
    for c in range(1,max_clusterid+1):
        cand = [i for i in range(len(point_3d)) if clusters[i]==c]
        x = np.take(points.x,cand)
        y = np.take(points.y,cand)
        z = np.take(points.z,cand)
        x=np.sum(x)/len(x)
        y=np.sum(y)/len(y)
        z=np.sum(z)/len(z)
        representative_array.append([x,y,z])
        x,y,z=0,0,0
    return representative_array,clusters

slab_points=laspy.file.File('slabclustered.las')
beam_points=laspy.file.File('beamclustered.las')
column_points=laspy.file.File('columnclustered.las')
wall_points=laspy.file.File('wallclustered.las')


slab_reps,slab_clusters=representative_point(slab_points)
beam_reps,beam_clusters=representative_point(beam_points)
column_reps,column_clusters=representative_point(column_points)
wall_reps,wall_clusters=representative_point(wall_points)


with open('tm.json') as f:
     proj_mat = json.load(f)
proj_mat=proj_mat['tm']

ifc=json.load(open('ifc_flat.json'))

bim=open('Basement1.obj', 'r')                                                                                                   #Opening BIM
bim=bim.read()                                                                                                                  #Reading BIM                
bim=bim.split('g ')                                                                                                             #Separating BIM to array of objects
obj_names,obj_vertices,obj_polys,z_bounds=[],[],[],[]                                                                           #Initializing arrays to extract the data to
for i in bim:                                                                                                                   #Reading each object from BIM
    k=i.splitlines()                                                                                                            #Separating each object to array of lines

    vertices=[]                                                                                                                 #Initialize array to add vertices of BIM object to
    for j in k:                                                                                                                 #Reading each line from object                    
        try:
            if j[0]=='v' and j[1]!='n':                                                                                         #Checking if line contains vertices
                v=j.split(' ')                                                                                                  #Converting to array of coordinates                                                                                   
                x,y,z=pnt_to_coords(proj_mat,float(v[1]),float(v[2]),float(v[3]))                                               #Convertng local coordinates to BIM coordinates
                vertices.append([x,y,z])                                                                                        #Appending the converted coordinates to vertices array
        except:
            pass
    obj_names.append(k[0])
    obj_vertices.append(vertices)

obj_names.pop(0)
obj_vertices.pop(0)
infile=laspy.file.File('Basement.las')
point_3d=np.vstack([infile.x,infile.y,infile.z]).T
tree = KDTree(point_3d)

for i in range(len(obj_names)):  
    set_points_xy=[]                           
    maxx,minx,maxy,miny,maxz,minz=0,1000000000,0,10000000000,0,10000000000                                                          #Initialize max,min variables
    for j in obj_vertices[i]:                                                                                                       #Extracting x,y coordinates and max,min variables
        set_points_xy.append([j[0],j[1]])
        if j[2]>maxz:
            maxz=j[2]
        if j[2]<minz:
            minz=j[2]
        if j[0]<minx:
            minx=j[0]
        if j[0]>maxx:
            maxx=j[0]
        if j[1]<miny:
            miny=j[1]
        if j[1]>maxy:
            maxy=j[1]
    
    set_points_xy=list(set([tuple(t) for t in set_points_xy]))                                                              #Set of points to eliminate duplicates
    sp_xy=[]
    hull=ConvexHull(set_points_xy)                                                                                                  #Applying convex hull to set of points
    for h_v in hull.vertices:                                                                                                       #Extracting coordinates of hull to an array
        sp_xy.append(set_points_xy[h_v])
    poly_xy=Polygon(sp_xy)                                                                                                          #Create polygon of extracted coordinates
    obj_polys.append(poly_xy)
    z_bounds.append([minz,maxz])


out_x,out_y,out_z=np.array([1]),np.array([1]),np.array([1])


point_3d=np.vstack([beam_points.x,beam_points.y,beam_points.z]).T
polygon=0


print('beam')
for i in range(len(beam_reps)):

    point=Point((beam_reps[i][0],beam_reps[i][1]))
    for k in range(len(obj_polys)):
        if point.within(obj_polys[k]) and beam_reps[i][2]<=z_bounds[k][1] and beam_reps[i][2]>=z_bounds[k][0] and ifc[obj_names[k]]['type']=='IfcBeam':
            
            polygon=obj_polys[k]
            height_bounds=z_bounds[k]
            name=obj_names[k]
    
    d=1000
    if polygon==0:
        for k in range(len(obj_polys)):
            if ifc[obj_names[k]]['type']=='IfcBeam':
                ext_coords=list(obj_polys[k].exterior.coords)
                for ec in ext_coords:
                    if distance(beam_reps[i],ec)<d:
                        d=distance(beam_reps[i],ec)
                        polygon=obj_polys[k]
                        height_bounds=z_bounds[k]
                        name=obj_names[k]


    
    bounds=polygon.bounds
    indices = [c for c in range(len(point_3d)) if beam_clusters[c]==i]
    x = np.take(beam_points.x,indices)
    y = np.take(beam_points.y,indices)
    z = np.take(beam_points.z,indices)
    struct_3d=np.vstack([x,y,z]).T
    minx,miny,minz,maxx,maxy,maxz=np.amin(x),np.amin(y),np.amin(z),np.amax(x),np.amax(y),np.amax(z)

    if abs(minx-bounds[0])>0.15 or abs(miny-bounds[1])>0.15 or abs(minz-height_bounds[0])>0.15 or abs(maxx-bounds[2])>0.15 or abs(maxy-bounds[3])>0.15 or abs(maxz-height_bounds[1])>0.15:
        print(name)
        out_x=np.append(out_x,x)
        out_y=np.append(out_y,y)
        out_z=np.append(out_z,z)

    polygon=0


point_3d=np.vstack([column_points.x,column_points.y,column_points.z]).T
polygon=0


print('column')

for i in range(len(column_reps)):

    point=Point((column_reps[i][0],column_reps[i][1]))
    for k in range(len(obj_polys)):
        if point.within(obj_polys[k]) and column_reps[i][2]<=z_bounds[k][1] and column_reps[i][2]>=z_bounds[k][0] and ifc[obj_names[k]]['type']=='IfcColumn':
            
            polygon=obj_polys[k]
            height_bounds=z_bounds[k]
            name=obj_names[k]

    d=1000
    if polygon==0:
        for k in range(len(obj_polys)):
            if ifc[obj_names[k]]['type']=='IfcColumn':
                ext_coords=list(obj_polys[k].exterior.coords)
                for ec in ext_coords:
                    if distance(column_reps[i],ec)<d:
                        d=distance(column_reps[i],ec)
                        polygon=obj_polys[k]
                        height_bounds=z_bounds[k]
                        name=obj_names[k]

    
    bounds=polygon.bounds
    indices = [c for c in range(len(point_3d)) if column_clusters[c]==i]
    x = np.take(column_points.x,indices)
    y = np.take(column_points.y,indices)
    z = np.take(column_points.z,indices)
    struct_3d=np.vstack([x,y,z]).T
    minx,miny,minz,maxx,maxy,maxz=np.amin(x),np.amin(y),np.amin(z),np.amax(x),np.amax(y),np.amax(z)

    if abs(minx-bounds[0])>0.15 or abs(miny-bounds[1])>0.15 or abs(minz-height_bounds[0])>0.15 or abs(maxx-bounds[2])>0.15 or abs(maxy-bounds[3])>0.15 or abs(maxz-height_bounds[1])>0.15:
        print(name)
        out_x=np.append(out_x,x)
        out_y=np.append(out_y,y)
        out_z=np.append(out_z,z)

    polygon=0




point_3d=np.vstack([wall_points.x,wall_points.y,wall_points.z]).T
polygon=0

print('wall')

for i in range(len(wall_reps)):

    point=Point((wall_reps[i][0],wall_reps[i][1]))
    for k in range(len(obj_polys)):
        if point.within(obj_polys[k]) and wall_reps[i][2]<=z_bounds[k][1] and wall_reps[i][2]>=z_bounds[k][0] and ifc[obj_names[k]]['type']=='IfcWallStandardCase':
            
            polygon=obj_polys[k]
            height_bounds=z_bounds[k]
            name=obj_names[k]

    d=1000
    if polygon==0:
        for k in range(len(obj_polys)):
            if ifc[obj_names[k]]['type']=='IfcWallStandardCase':
                ext_coords=list(obj_polys[k].exterior.coords)
                for ec in ext_coords:
                    if distance(wall_reps[i],ec)<d:
                        d=distance(wall_reps[i],ec)
                        polygon=obj_polys[k]
                        height_bounds=z_bounds[k]
                        name=obj_names[k]
    
    bounds=polygon.bounds
    indices = [c for c in range(len(point_3d)) if wall_clusters[c]==i]
    x = np.take(wall_points.x,indices)
    y = np.take(wall_points.y,indices)
    z = np.take(wall_points.z,indices)
    struct_3d=np.vstack([x,y,z]).T
    minx,miny,minz,maxx,maxy,maxz=np.amin(x),np.amin(y),np.amin(z),np.amax(x),np.amax(y),np.amax(z)

    if abs(minx-bounds[0])>0.15 or abs(miny-bounds[1])>0.15 or abs(minz-height_bounds[0])>0.15 or abs(maxx-bounds[2])>0.15 or abs(maxy-bounds[3])>0.15 or abs(maxz-height_bounds[1])>0.15:
        print(name)
        out_x=np.append(out_x,x)
        out_y=np.append(out_y,y)
        out_z=np.append(out_z,z)
    polygon=0



point_3d=np.vstack([slab_points.x,slab_points.y,slab_points.z]).T
polygon=0


print('slab')
for i in range(len(slab_reps)):

    point=Point((slab_reps[i][0],slab_reps[i][1]))
    for k in range(len(obj_polys)):
        if point.within(obj_polys[k]) and slab_reps[i][2]<=z_bounds[k][1] and slab_reps[i][2]>=z_bounds[k][0] and ifc[obj_names[k]]['type']=='IfcSlab':
            
            polygon=obj_polys[k]
            height_bounds=z_bounds[k]
            name=obj_names[k]
    
    d=1000
    if polygon==0:
        for k in range(len(obj_polys)):
            if ifc[obj_names[k]]['type']=='IfcSlab':
                ext_coords=list(obj_polys[k].exterior.coords)
                for ec in ext_coords:
                    if distance(slab_reps[i],ec)<d:
                        d=distance(slab_reps[i],ec)
                        polygon=obj_polys[k]
                        height_bounds=z_bounds[k]
                        name=obj_names[k]


    
    bounds=polygon.bounds
    indices = [c for c in range(len(point_3d)) if slab_clusters[c]==i]
    x = np.take(slab_points.x,indices)
    y = np.take(slab_points.y,indices)
    z = np.take(slab_points.z,indices)
    struct_3d=np.vstack([x,y,z]).T
    minx,miny,minz,maxx,maxy,maxz=np.amin(x),np.amin(y),np.amin(z),np.amax(x),np.amax(y),np.amax(z)

    if abs(minx-bounds[0])>0.15 or abs(miny-bounds[1])>0.15 or abs(minz-height_bounds[0])>0.15 or abs(maxx-bounds[2])>0.15 or abs(maxy-bounds[3])>0.15 or abs(maxz-height_bounds[1])>0.15:
        print(name)
        out_x=np.append(out_x,x)
        out_y=np.append(out_y,y)
        out_z=np.append(out_z,z)

    polygon=0





out_x = np.delete(out_x, (0), axis=0)
out_y = np.delete(out_y, (0), axis=0)
out_z = np.delete(out_z, (0), axis=0)
red=np.ones(len(out_x))*255
green=np.zeros(len(out_x))
blue=np.zeros(len(out_x))

outfile=laspy.file.File('Test.las',mode='w',header=infile.header)
outfile.x=out_x
outfile.y=out_y
outfile.z=out_z
outfile.red=red
outfile.green=green
outfile.blue=blue
outfile.close()
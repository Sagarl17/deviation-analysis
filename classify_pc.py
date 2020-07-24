import os
import sys
import json
import math
import laspy
import numpy as np
import multiprocessing
from scipy.spatial import ConvexHull,KDTree
from shapely.geometry import Polygon,Point,LineString

sys.setrecursionlimit(100000000)


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

def get_points(index,d=0.5):
    print(index,len(obj_names))
    polygon=obj_polys[index]
    height_bounds=z_bounds[index]
    ext_coords=list(polygon.exterior.coords)
    bounds=polygon.bounds

    point_3ds=point_3d[(point_3d[:,2]>= height_bounds[0]-0.5) & (point_3d[:,2]<=height_bounds[1]+0.5)]                                                              
    point_3ds=point_3ds[(point_3ds[:,0]>= bounds[0]-0.5) & (point_3ds[:,0]<=bounds[2]+0.5)]                                                           
    point_3ds=point_3ds[(point_3ds[:,1]>= bounds[1]-0.5) & (point_3ds[:,1]<=bounds[3]+0.5)]

    lines,indices,nearest_structures=[],[],[]
    for i in range(0,len(ext_coords)-1):
        lines.append([ext_coords[i],ext_coords[i+1]])
    for line in lines:
        if line[1][0]==line[0][0]:
            x1,x2,x3,x4=line[1][0],line[1][0],line[0][0],line[0][0]
            y1,y2,y3,y4=line[1][1]-0.5,line[1][1]+0.5,line[0][1]-0.5,line[0][1]+0.5
        else:
            if line[1][1]==line[0][1]:
                y1,y2,y3,y4=line[1][1],line[1][1],line[0][1],line[0][1]
                x1,x2,x3,x4=line[1][0]-0.5,line[1][0]+0.5,line[0][0]-0.5,line[0][0]+0.5
            else:
                slope=(line[1][1]-line[0][1])/(line[1][0]-line[0][0])
                slope=-1/slope
                x1=line[0][0]+(d/(math.sqrt(1+slope*slope)))
                x2=line[0][0]-(d/(math.sqrt(1+slope*slope)))
                x3=line[1][0]+(d/(math.sqrt(1+slope*slope)))
                x4=line[0][0]-(d/(math.sqrt(1+slope*slope)))

                y1=line[0][1]+((d*slope)/(math.sqrt(1+slope*slope)))
                y2=line[0][1]-((d*slope)/(math.sqrt(1+slope*slope)))
                y3=line[1][1]+((d*slope)/(math.sqrt(1+slope*slope)))
                y4=line[1][1]-((d*slope)/(math.sqrt(1+slope*slope)))
        poly_line=Polygon([(x1,y1),(x2,y2),(x3,y3),(x4,y4)])
        for point in point_3ds:
            point_geom=Point((point[0],point[1]))
            if point_geom.within(poly_line) or point_geom.within(polygon):
                ind=np.where(np.all(point_3d==point,axis=1))
                indices.append(int(ind[0][0]))
    
    for point in ext_coords:
        query_points=tree.query_ball_point([point[0],point[1],height_bounds[0]],r = .5)
        indices=indices+query_points
        query_points=tree.query_ball_point([point[0],point[1],height_bounds[1]],r = .5)
        indices=indices+query_points
    
    indices=list(set(indices))
    
    for i in range(len(obj_polys)):
        if i!=index:
            poly_ext_coords=list(obj_polys[i].exterior.coords)
            for j in ext_coords:
                for k in poly_ext_coords:
                    if distance(j,k)<=0.5:
                        nearest_structures.append(i)
                        break
    removal_list=[]
    for p in nearest_structures:
        for i in indices:
            point_geom=Point((point_3d[i][0],point_3d[i][1]))
            if point_geom.within(obj_polys[p]) and not point_geom.within(polygon):
                removal_list.append(i)

    removal_list=list(set(removal_list))
    
    for i in removal_list:
        indices.remove(i)

    return indices

with open('tm.json') as f:
     proj_mat = json.load(f)
proj_mat=proj_mat['tm']


ifc=json.load(open('ifc_flat.json'))
structure=[]
for i in ifc:
    structure.append(ifc[i]['type'])
structure=list(set(structure))
structure_cluster=[0]*len(structure)

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

pool = multiprocessing.Pool(processes=multiprocessing.cpu_count() - 1)
indices=pool.map(get_points, range(len(obj_names)),chunksize=1)

structure_points=[]
for i in range(len(obj_names)):
    if len(indices[i])>0:
        for j in range(len(structure)):
            if ifc[obj_names[i]]['type']==structure[j]:
                s=np.ones(len(indices[i]))*structure_cluster[j]
                c=np.ones(len(indices[i]))*j
                indices[i]=np.take(infile.points,indices[i])
                if len(structure_points)==0:
                    structure_points=indices[i]
                    structure_clusters=s
                    structure_class=c
                else:
                    structure_points=np.append(structure_points,indices[i])
                    structure_clusters=np.append(structure_clusters,s)
                    structure_class=np.append(structure_class,c)
                structure_cluster[j]=structure_cluster[j]+1


structure_class = structure_class.astype(np.uint8, copy=False)
outfile=laspy.file.File('Test.las',mode='w',header=infile.header)
outfile.points=structure_points
outfile.close()
infile=laspy.file.File('Test.las',mode='r')

header = infile.header
outfile=laspy.file.File('classified_pc.las',mode='w',header=header)
outfile.define_new_dimension(name = 'cluster_id',data_type = 5, description = 'Cluster_id')
outfile.cluster_id=structure_clusters
outfile.x = infile.x
outfile.y = infile.y
outfile.z = infile.z
outfile.red = infile.red
outfile.green = infile.green
outfile.blue = infile.blue
outfile.classification=structure_class
outfile.close()
print(structure)

infile=laspy.file.File('classified_pc.las',mode='w',header=header)
classess=infile.classification
cand = [i for i in range(len(point_3d)) if classess[i]==0]
points = np.take(infile.points,cand)
outfile=laspy.file.File('slabclustered.las',mode='w',header=header)
outfile.points=points
outfile.close()
cand = [i for i in range(len(point_3d)) if classess[i]==1]
points = np.take(infile.points,cand)
outfile=laspy.file.File('columnclustered.las',mode='w',header=header)
outfile.points=points
outfile.close()
cand = [i for i in range(len(point_3d)) if classess[i]==2]
points = np.take(infile.points,cand)
outfile=laspy.file.File('beamclustered.las',mode='w',header=header)
outfile.points=points
outfile.close()
cand = [i for i in range(len(point_3d)) if classess[i]==4]
points = np.take(infile.points,cand)
outfile=laspy.file.File('wallclustered.las',mode='w',header=header)
outfile.points=points
outfile.close()


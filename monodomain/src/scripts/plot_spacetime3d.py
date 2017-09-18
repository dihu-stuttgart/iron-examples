import math
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from mpl_toolkits.mplot3d import Axes3D
import exnode_reader
import sys

# get folder name from command line argument
foldername="./"
if len(sys.argv) >1:
   foldername=sys.argv[1]

# append slash if necessary
if foldername[-1] != "/":
   foldername=foldername + "/" 

print "foldername: ", foldername

# creat image  name
imagename="{}".format(foldername)+"image.png"

# get all exnode files from the folder
ls = exnode_reader.get_exnode_files(foldername)
condition_file= lambda name: "Time" in name
ls_selected=list(np.extract(map(condition_file,ls),ls))
files=ls_selected
#files=sorted(ls_selected)
#print 'exnode files',files

# get all parts in parallel runs from the folder
condition_part= lambda name: "_0.part" in name
ls_parts=list(np.extract(map(condition_part,ls_selected),ls_selected))
n_parts=len(ls_parts) #No. of cores
print "n_parts: ",n_parts

n_files=len(files)
print "n_files: ", n_files
n_timestep=n_files/n_parts
print "No. of time steps: ", n_timestep

#get the number of data in each file
x=exnode_reader.parse_file(foldername+files[0],[["Coordinate",1]]) #field_name and component_no
n_data=len(x)
print "No. of data point: ",n_data
#print x 

# get the data
Vm=[0.0 for i in range(n_timestep)]
for i_time in range (n_timestep):
    Vm[i_time]=exnode_reader.parse_file(foldername+files[i_time*n_parts+1],[["Vm",1]])
    #print "time step: ", i_time
    #print Vm[i_time]  
#print "Vm: ",Vm

def makeplot(x,y,z):

    global foldername, imagename

    Y ,X=np.meshgrid(y,x) #different ordering, columnss are elements of y. There are len(x) rows.
#    print "lengthof x: ",len(x),"X: ", X
#    print "length of y:", len(y),"Y: ", Y
#    print "z: ", z
    Z = X.copy()
    for i in range(len(x)):
     for j in range(len(y)):
         Z[i,j]=z[i][j]
#         print "i: ",i,"j: ",j,"Z: ",Z[i,j] 
        
    fig=plt.figure(1)
    ax1=fig.gca(projection='3d')
    plt.tight_layout()
    ax1.plot_surface(X,Y,Z,cmap=cm.coolwarm, linewidth=0,rstride=1,cstride=1)
    #ax1.gridlines.set_visible(False)
    ax1.set_zlim(-80,50)
    plt.title("Vm")
    plt.show()    
    plt.savefig(imagename)
    print "Figure saved to {}".format(imagename)

#print foldername
j_e=len(foldername)-1
for j in range(j_e):
    if foldername[j]=='t':
       j_s=j+1

time_step=float("0."+foldername[j_s:j_e])
print "time_step",time_step
time=np.arange(n_timestep)*time_step
end_time=10.0
time=time*end_time/time[n_timestep-1]
print "time: ",time
makeplot(time,x,Vm)



#!/usr/bin/python 
# -*- coding: utf-8 -*-

########################################################################
#  File Name	: 'plot_v_a.py'
#  Author	: Steve NGUYEN
#  Contact      : steve.nguyen@college-de-france.fr
#  Created	: mercredi, septembre 14 2011
#  Revised	: 
#  Version	: 
#  Target MCU	: 
#
#  This code is distributed under the GNU Public License
# 		which can be found at http://www.gnu.org/licenses/gpl.txt
#
#
#  Notes:	notes
########################################################################


import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
#import numpy as np
import sys
from scipy import *


def load_file(filename,sizex,sizey):
    print "loading %s..." %(filename)
    f=open(filename,'r')
    #data=f.read()
    #x=[]
    #print data.split(' ')
    #for i in data.split(' '):
    #    if i!='':
    #        x.append(float(i))
    #print x
    
    #data=data.split(' ')
    #data=data[34:]
    #data=data.split
    #print data
    #res=map(float,data)
   # print res


    data=[]
    it=0
    for l in f:
        if l[0]!='#':
            l=l.replace(' \n','')
            l=l.split(' ')
            li=[]
            
            for ll in l:
                if ll!='':
                    # if float(ll)<50.0:
                    #     li.append(float(ll))
                    # else:
                    #     li.append(50.0)

                    li.append(float(ll))
            data.append(li)
    a=array(data)
    a=a.reshape(sizex,sizey)
    print a
    return a

if 'control' in sys.argv[1]:
	x=3 #51
	y=200 #51
else:
	x=51
	y=51

v=load_file(sys.argv[1],x,y)

fname=sys.argv[1].split('.')[0]
fname+='.pdf'

if 'control' in sys.argv[1]:
	plt.plot(range(y),v.transpose())
	#plt.show()
	plt.savefig(fname)
	exit()

o=ones(x*y)
o=o.reshape(x,y)

#v=log(o/v+o)



fig = plt.figure()
ax = fig.gca(projection='3d')
X = arange(0, 51, 1)
Y = arange(0, 51, 1)
X, Y = meshgrid(X, Y)

Z = v.transpose()

surf = ax.plot_surface(X, Y, Z, rstride=1, cstride=1, cmap=cm.jet,
        linewidth=0, antialiased=False)
#ax.set_zlim(-1.01, 1.01)

fig.colorbar(surf, shrink=0.5, aspect=5)

ax.set_xlabel(r'Position')
ax.set_ylabel(r'Velocity')
ax.set_zlabel(r'Quasi-distance')

ax.view_init(45, 0)



for i in range(0,360,10):
    ax.view_init(45, i)
    s='quasidist-%.2d.png'%((i/10))
    print s
    plt.savefig(s,bbox_inches='tight')





#plt.show()



#loc,lab=plt.xticks()
#print loc

# fig = plt.figure()
# ax = fig.add_subplot(111)
# #im=plt.imshow(Dist.transpose(),origin='lower',vmin=0.0,vmax=100.0)

# im=plt.imshow(v.transpose(),origin='lower')
# #im=plt.contourf(v.transpose(),50,origin='lower')

# #im=plt.imshow(Dist.transpose(),origin='lower',norm=LogNorm(vmin=0.0,vmax=10e10))
# im.set_interpolation('nearest')
# cb=plt.colorbar()
# #cb.ax.set_ylabel(t)

# #plt.gca().set_zscale('log')

# #plt.show()
# plt.savefig(fname)
                          

#!/usr/bin/python 
# -*- coding: utf-8 -*-

########################################################################
#  File Name   : 'plot_v_a.py'
#  Author   : Steve NGUYEN
#  Contact      : steve.nguyen@college-de-france.fr
#  Created  : mercredi, septembre 14 2011
#  Revised  : 
#  Version  : 
#  Target MCU  : 
#
#  This code is distributed under the GNU Public License
#     which can be found at http://www.gnu.org/licenses/gpl.txt
#
#
#  Notes:   notes
########################################################################


from scipy import *

from matplotlib import rc
import matplotlib.pyplot as plt
import sys


X=51
Y=51
THETA=51

def load_file(filename,sizex,sizey,sizetheta):
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
                    li.append(float(ll))
            data.append(li)
    a=array(data)
    a=a.reshape(sizex,sizey,sizetheta)
    #print a
    return a

if 'control' in sys.argv[1]:
   x=3 #51
   y=200 #51
else:
   x=X
   y=Y
   theta=THETA

v=load_file(sys.argv[1],x,y,theta)



fname=sys.argv[1].split('.')[0]
fname+='.pdf'

# if 'control' in sys.argv[1]:
#  plt.plot(range(y),v.transpose())
#  #plt.show()
#  plt.savefig(fname)
#  exit()

o=zeros(x*y)
o=o.reshape(x,y)

for i in range(x):
    for j in range(y):
        o[i][j]=v[i][j][theta/2]
        #o[i][j]=v[i][j][theta-1]
        #o[i][j]=v[i][j][0]
        #o[i][j]=v[i][j][theta/2+4]

        #o[i][j]=v[i][j][4]
        #o[i][j]=v[i][j][8]

        #o[i][j]=min(v[i][j])

#v=log(o/v+o)


fig = plt.figure()
ax = fig.add_subplot(111)
#im=plt.imshow(Dist.transpose(),origin='lower',vmin=0.0,vmax=100.0)

#im=plt.imshow(v.transpose(),origin='lower',vmin=9.7,vmax=10.0)
im=plt.imshow(o.transpose(),origin='lower',vmin=0,vmax=10.0)
#im=plt.imshow(o.transpose(),origin='lower',vmin=0,vmax=100.0)
#im=plt.imshow(o.transpose(),origin='lower',vmin=0,vmax=500.0)
#im=plt.imshow(o.transpose(),origin='lower')
#im=plt.contourf(o.transpose(),100,origin='lower')

plt.axis('equal')

#im=plt.contourf(v.transpose(),50,origin='lower')

#im=plt.imshow(Dist.transpose(),origin='lower',norm=LogNorm(vmin=0.0,vmax=10e10))
im.set_interpolation('nearest')
cb=plt.colorbar()
#cb.ax.set_ylabel(t)


#cb.set_ticks([0,10,20])
#cb.set_ticklabels(['-0.25','0','0.25'])


#plt.gca().set_zscale('log')


#ax.set_xlabel(r'Position (rad)')
#ax.set_ylabel(r'Velocity (rad/s)')
rc('text', usetex=True)

#loc,lab= plt.xticks()
#plt.xticks( (0,6.25,12.75,19.125,25.5,31.875,38.25,44.625,50), ('-180', '-135', '-90', '-45', '0','45', '90', '135', '180') )
#plt.xticks( (0,6.25,12.75,19.125,25.5,31.875,38.25,44.625,50), (r'$-\pi$', r'-$\frac{3\pi}{4}$', r'$-\frac{\pi}{2}$', r'$-\frac{\pi}{4}$', '$0$',r'$\frac{\pi}{4}$', r'$\frac{\pi}{2}$', r'$\frac{3\pi}{4}$', r'$\pi$') )
#plt.yticks( (0,8.5,17,25.5,34,42.5,50), ('$-1.5$', '$-1.0$', '$-0.5$', '$0.0$', '$0.5$','$1.0$', '$1.5$') )
#plt.yticks( (0,6.25,12.75,19.125,25.5,31.875,38.25,44.625,50), (r'$-2.5$', r'$-1.875$', r'$-1.25$', r'$-0.625$', '$0$',r'$0.625$', r'$1.25$', r'$1.875$', r'$2.5$') )


plt.show()
#plt.savefig('QD2D.pdf',bbox_inches='tight')
#plt.savefig('value2D_torqueangle.png',bbox_inches='tight')
#plt.savefig('qdpolicy.pdf',bbox_inches='tight')
#plt.savefig('bellpolicy_torqueangle.png',bbox_inches='tight')
                          

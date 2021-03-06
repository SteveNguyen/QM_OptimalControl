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


from scipy import *

import matplotlib.pyplot as plt
import sys



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
                    li.append(float(ll))
            data.append(li)
    a=array(data)
    a=a.reshape(sizex,sizey)
    print a
    return a





#ctrl=loadtxt('control.dat')
ctrl=loadtxt('simu_bel.dat')
#print ctrl


#print ctrl

# fig = plt.figure()
# ax = fig.add_subplot(111)
# #im=plt.imshow(Dist.transpose(),origin='lower',vmin=0.0,vmax=100.0)

# ax.set_aspect(1)
# plt.grid(True)

# plt.xlim(-1.2,1.2)
# plt.ylim(-1.2,1.2)


# plt.plot(0,0, 'bo')

for i in range(70):

    fig = plt.figure()
    ax = fig.add_subplot(111)

    #plt.clf()
    ax.set_aspect(1)
    plt.grid(True)

    plt.xlim(-1.2,1.2)
    plt.ylim(-1.2,1.2)
    plt.plot(0,0, 'bo')
    
    l = plt.plot(cos(ctrl[i][0]+pi/2.0),sin(ctrl[i][0]+pi/2.0), 'ro')           
    plt.plot([0,cos(ctrl[i][0]+pi/2.0)],[0,sin(ctrl[i][0]+pi/2.0)], 'b')

    plt.setp(l, 'markersize', 20)
    plt.setp(l, 'markerfacecolor', 'b')

    plt.draw()
    s='bell-penduleanim-%.2d.png'%(i)
    
    plt.savefig(s,bbox_inches='tight')


# im=plt.imshow(v.transpose(),origin='lower')

# im.set_interpolation('nearest')
# cb=plt.colorbar()

# ax.set_xlabel(r'Position')
# ax.set_ylabel(r'Velocity')


# plt.xticks( (0,6.25,12.75,19.125,25.5,31.875,38.25,44.625,50), ('-180', '-135', '-90', '-45', '0','45', '90', '135', '180') )
# plt.yticks( (0,8.5,17,25.5,34,42.5,50), ('-1.5', '-1.0', '-0.5', '0.0', '0.5','1.0', '1.5') )

#plt.show()

# for i in range(70):


#     #plt.plot((ctrl[:i+1])[i][0],(ctrl[:i+1])[i][1])

#     for j in range(i):
#         if (abs((ctrl[:i+1])[j+1][0]-(ctrl[:i+1])[j][0])<10) and (abs((ctrl[:i+1])[j+1][1]-(ctrl[:i+1])[j][1])<10):
#             plt.plot([(ctrl[:i+1])[j][0],(ctrl[:i+1])[j+1][0]],[(ctrl[:i+1])[j][1],(ctrl[:i+1])[j+1][1]],'k',lw=2)
#         else:
#             plt.plot((ctrl[:i+1])[j][0],(ctrl[:i+1])[j][1],'k.')

            
#     s='qdtrj-%.2d.png'%(i)
    
#     plt.savefig(s,bbox_inches='tight')

                          

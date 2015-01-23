#!/usr/bin/python 
# -*- coding: utf-8 -*-


from scipy import *

import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
import sys


meta_file=sys.argv[1]+"/meta.dat"
landscape_file=sys.argv[1]+"/landscape.dat"
#landscape_file=sys.argv[1]+"/learning.dat"
policy_file=sys.argv[1]+"/policy.dat"
control_file=sys.argv[1]+"/control.dat"
control_discr_file=sys.argv[1]+"/control_discr.dat"
        

def split_line(f):
    return map(int,f.readline().replace('\n','').split(' '))
# state dim cardinality
#l=f.readline()
#l=l.replace('\n','')
#l=l.split(' ')

# META
f=open(meta_file)
l=split_line(f)
stateCards=[]
for ll in l:
    stateCards.append(int(ll))

# landscape
l=split_line(f)
isProbs=l[0]

# control

#l=split_line(f)
#nb_watch=int(l[0])
#nb_t=int(l[1])

f.close()


def load_file(filename,sizex,sizey, nb_cont=1):
    print "loading %s..." %(filename)
    f=open(filename,'r')
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
    f.close()
    a=array(data)
    a=a.reshape(sizex,sizey)
    b=a
    for i in range(nb_cont-1):
        b=concatenate((b,a))
    return b

#CREATE FIG
fig = plt.figure()

#IMPORT LANDSCAPE
print stateCards
landscape=load_file(landscape_file, stateCards[0], stateCards[1],1)

if isProbs:
    o=ones(size(landscape))
    o=o.reshape(landscape.shape)
    landscape=log(o/landscape+o)

#IMPORT POLICY
policy=load_file(policy_file, stateCards[0], stateCards[1],1)

#IMPORT CONTROL
control=mlab.load(control_file)

nb_t=control.shape[0]
nb_watch=control.shape[1]

#IMPORT CONTROL_DISCR
control_discr=mlab.load(control_discr_file)
#control_discr[control_discr[:,0]<(3*(stateCards[0])/4),0]+=stateCards[0]

def plot_sub_map(data, sub):
    ax = fig.add_subplot(sub)
    im=plt.imshow(data.transpose(),origin='lower')
    im.set_interpolation('nearest')
    cb=plt.colorbar()
    # plt.plot(control_discr[:,0],control_discr[:,1])
#FILL FIG

sub=221

plot_sub_map(landscape, sub)

sub+=1
plot_sub_map(policy, sub)

sub+=1
ax = fig.add_subplot(sub)
plt.plot(range(nb_t),control)

plt.savefig(sys.argv[1]+"/plot.pdf")


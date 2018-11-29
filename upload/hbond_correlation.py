# ***************************************************************
# File : hbond_correlation.py 
# Used to calculate H-bond correlation functions S(t),C(t),Sd(t),Cd(t) 

# usage: python hbond_correltion.py [-h] f side atoms

# positional arguments:
#   f           Name of the trajectory file
#   side        Size of the box
#   atoms       Total number of atoms

# optional arguments:
#   -h, --help  show help message 
# ***************************************************************

from itertools import islice,izip
from hbond_pair import hbond,popvar
from multiprocessing import Pool
import matplotlib.pyplot as plt
import argparse,time
from math import log


def parse():                                        #function to parse command-line arguments   
    parser = argparse.ArgumentParser(description='Hydrogen Bond Analysis')
    parser.add_argument("f",    type=str,  help="Name of the trajectory file")
    parser.add_argument("side" ,type=float,help="Size of the box")
    parser.add_argument("atoms",type=int,  help="Total number of atoms")
    args=parser.parse_args() 
    return [args.f, args.side, args.atoms]


def get_h(tmp):                                     #process the frame and call fortran module                                
    tmp =[map(float,line.split()[3:6]) for line in tmp[:768]]
    return hbond(tmp,side)


def calculate():                                    #returns number of h-bond and their probability 
    hframe=[]
    pl=Pool()                                       #create a pool of worker processes
    with open(file) as f:
        while True:
            tmp=list(islice(f,9,9+atoms))           #read one time frame from file
            if not tmp: break
            hframe.append(pl.apply_async(get_h,[tmp])) #load the function call to a worker in parallel
    pl.close()
    pl.join()
    hframe=[i.get() for i in hframe]
    hframe,hdframe=zip(*hframe)
    shb,chb,shbd,chbd=popvar(hframe,hdframe)
    return shb,chb,shbd,chbd

file,side,atoms=parse()

def save(a,b,c):   									#save the data
    txt="\n".join("%s %s"%(i,j) for i,j in zip(a,b))
    with open(c+"_correlation.dat","w") as f:
        f.write(txt)

if __name__ == '__main__':
    shb,chb,shbd,chbd=calculate()
    x=[i/1000.0 for i in xrange(len(shb))]
    save(x,shb,"shb")
    save(x,chb,"chb")
    save(x,shbd,"shbd")
    save(x,chbd,"chbd")
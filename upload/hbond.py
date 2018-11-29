# ***************************************************************
# File : hbond.py 
# Used to calculate distribution of H-bond from trajectory file 

# usage: python hbond.py [-h] f side atoms

# positional arguments:
#   f           Name of the trajectory file
#   side        Size of the box
#   atoms       Total number of atoms

# optional arguments:
#   -h, --help  show help message 
# ***************************************************************

from itertools import islice
from hbond import hbond
from multiprocessing import Pool
import matplotlib.pyplot as plt
import argparse,time


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
    ll=[]
    pl=Pool()                                       #create a pool of worker processes
    with open(file) as f:
        while True:
            tmp=list(islice(f,9,9+atoms))           #read one time frame from file
            if not tmp: break
            ll.append(pl.apply_async(get_h,[tmp]))  #load the function call to a worker in parallel
    res=[i for j in ll for i in j.get()]            #extract result from the worker process
    return prob(res)

def prob(ll):                                       #calculate probabilities from raw data
    bins=(max(ll)+1)
    hist=[0]*bins
    for i in ll : hist[i]+=1
    hist=[i/float(len(ll)) for i in hist]
    return range(bins),hist

def save(a,b):                                      #save the data
    txt="Average= %s\n"%sum(i*j for i,j in zip(a,b))                                    
    txt+="\n".join("%s %s"%(i,j) for i,j in zip(a,b))
    with open("hbond_dist.dat","w") as f:
        f.write(txt)

file,side,atoms=parse()
if __name__ == '__main__':
    s=time.time()
    a,b=calculate()
    save(a,b)
    plt.plot(a,b,"ro-",markersize=1,lw=1)
    plt.title("H-bond distribution")
    plt.xlabel("Number of H-bond (n)", fontsize=10)
    plt.ylabel("Probability of getting a Water molecule with n H-bond ", fontsize=10)
    print "Average Number of H-bond", sum(i*j for i,j in zip(a,b))
    print "Time needed", time.time()-s
    # plt.show()
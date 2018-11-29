from itertools import islice,izip
from hbond import hbond,pop
from multiprocessing import Pool
import matplotlib.pyplot as plt
import argparse,time
from math import log


def get_h(tmp):                                     #process the frame and call fortran module                                
    tmp =[map(float,line.split()[3:6]) for line in tmp[:768]]
    return hbond(tmp,side)


def calculate():                                    #returns number of h-bond and their probability 
    hframe=[]
    pl=Pool()                                       #create a pool of worker processes
    s=time.time()
    with open(file) as f:
        while True:
            tmp=list(islice(f,9,9+atoms))           #read one time frame from file
            if not tmp: break
            hframe.append(pl.apply_async(get_h,[tmp]))  #load the function call to a worker in parallel
    pl.close()
    pl.join()
    hframe=[i.get() for i in hframe]
    print time.time()-s
    s=time.time()
    # ll=[pl.apply_async(get_sc,[tmp]) for tmp in izip(*hframe)]
    # ll=[i.get() for i in ll]
    shb,chb=pop(hframe)

    print time.time()-s
    return shb,chb




file,side,atoms=r"E:\Project\new_start\0\nve\for_shb\water_nve_u.traj",19.79338,768
def save(a,b,c):   
    txt="\n".join("%s %s"%(i,j) for i,j in zip(a,b))
    with open(c+"hbond_dist.dat","w") as f:
        f.write(txt)

if __name__ == '__main__':
    shb,chb=calculate()
    x=[i/1000.0 for i in xrange(len(shb))]
    save(x,shb,"shb")
    save(x,chb,"chb")
    plt.plot(x,shb,"r-")
    plt.plot(x,chb,"b-")
    plt.show()
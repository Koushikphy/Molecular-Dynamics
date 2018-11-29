#**************************************************************************
# File : rdf.py
# Used for obtaining RDF,RCN and KB integral from trajectory files
# Usage: python rdf.py [-h] [-dr DR] [-plot PLOT] f side atoms name1 name2
# positional arguments:
#   f           Name of the trajectory file
#   side        Size of the box
#   atoms       Total number of atoms
#   name1       Name of 1st element
#   name2       Name of 2nd element

# optional arguments:
#   -h, --help  show help message
#   -dr DR      Size of the bins
#   -plot PLOT  Plot the result
#
#****************************************************************************

from rdf import homo,hetero,normalise
from itertools import islice
import matplotlib.pyplot as plt
from multiprocessing import Pool
import argparse,time



def parse():                                        #function to parse command-line arguments   
    parser = argparse.ArgumentParser(description='Radial Distribution Function')
    parser.add_argument("f",    type=str,  help="Name of the trajectory file")
    parser.add_argument("side" ,type=float,help="Size of the box")
    parser.add_argument("atoms",type=int,  help="Total number of atoms")
    parser.add_argument("name1",type=str,  help="Name of 1st element")
    parser.add_argument("name2",type=str,  help="Name of 2nd element")
    parser.add_argument("-dr"  ,type=float,help="Size of the bins",default=.1)
    parser.add_argument("-plot",type=bool, help="Plot the result" ,default=False)
    args=parser.parse_args()  
    bins=int(args.side/(2.0*args.dr))
    return [args.f, args.side, args.atoms,args.name1,args.name2, args.dr, bins,args.plot]


def get(tm):                                        #returns non-normalised rdf data
    tm =[line.split()[2:6] for line in tm]
    tm1=[[map(float,i[1:])] for i in tm if i[0]==name1]
    if name1==name2:
        return homo(tm1,side,dr,bins)
    tm2=[[map(float,j[1:])] for j in tm if j[0]==name2]
    return hetero(tm1,tm2,side,dr,bins)

def get_num(tm):                                    #returns # of atoms of type 1 and 2
    tm=[line.split()[2] for line in tm]
    m=tm.count(name1)
    if name1==name2: return m,m
    return m,tm.count(name2)

def calculate():                                    #returns normalised rdf,rcn,kbi data
    p=Pool()                                        #create a pool of worker processes
    res=[]
    with open(file) as f:
        while True:
            tmp=list(islice(f,9,9+atoms))
            if not tmp: break
            res.append(p.apply_async(get,[tmp]))    #load the function call to a worker in parallel
            tm=tmp
    res=[r.get() for r in res]                      #extract result from the worker process
    m,n=get_num(tm)
    return normalise(res,side,dr,m,n)   			#normalise the data            

def save(r,rdf,rcn,kb):                             #saves the data
    txt="\n".join(["%s %s"%(i,j) for i,j in zip(*[r,rdf])])
    txt1="\n".join(["%s %s"%(i,j) for i,j in zip(*[r,rcn])])
    txt2="\n".join(["%s %s"%(i,j) for i,j in zip(*[r,kb])])
    with open("rdf_%s%s.dat"%(name1,name2),"w") as f:
        f.write(txt+"\n\n")
    with open("rcn_%s%s.dat"%(name1,name2),"w") as g:
        g.write(txt1+"\n\n")
    with open("kb_%s%s.dat"%(name1,name2),"w") as g:
        g.write(txt2+"\n\n")

def plot_result(r,rdf,rcn,kb):                         #plots the result
    plt.plot(r,rdf)
    plt.xlabel("$r(\AA)$", fontsize=18)
    plt.ylabel("$g_{{{}{}}}(r)$".format(name1,name2), fontsize=18)
    plt.figure()
    plt.xlabel("$r(\AA)$", fontsize=18)
    plt.ylabel("$N_{{{}{}}}(r)$".format(name1,name2), fontsize=18)
    plt.plot(r,rcn)
    plt.figure()
    plt.xlabel("$r(\AA)$", fontsize=18)
    plt.ylabel("$G_{{{}{}}}(r)$".format(name1,name2), fontsize=18)
    plt.plot(r,kb)
    plt.show()

file,side,atoms,name1,name2,dr,bins,plot=parse()  
if __name__ == '__main__':
    s=time.time()
    r,rdf,rcn,kb=calculate()
    save(r,rdf,rcn,kb)
    print "Time needed", time.time()-s
    if plot:
        plot_result(r,rdf,rcn,kb)

from rdft import homo,hetero,normalise
from itertools import islice
import matplotlib.pyplot as plt
from multiprocessing import Pool
import argparse


def parse():
    parser = argparse.ArgumentParser(description='Radial Distribution Function')
    parser.add_argument("f",  type=str,  help="Name of the trajectory file")
    parser.add_argument("name1",type=str,  help="Name of 1st element")
    parser.add_argument("name2",type=str,  help="Name of 2nd element")
    parser.add_argument("side" ,type=float,help="Size of the box")
    parser.add_argument("--dr" ,type=float,help="Size of the bins",default=.1)
    parser.add_argument("--plot",type=bool,help="Plot the result" ,default=False)
    args=parser.parse_args()  
    bins=int(args.side/(2.0*args.dr))
    return [args.f,args.name1,args.name2,args.side,args.dr,bins,args.plot]


def get_hist(tm):
    tm =[line.split()[2:6] for line in tm]
    tm1=[[map(float,i[1:])] for i in tm if i[0]==name1]
    if name1==name2:
        return homo(tm1,side,dr,bins)
    tm2=[[map(float,j[1:])] for j in tm if j[0]==name2]
    return hetero(tm1,tm2,side,dr,bins)

def get_num(tm):
    tm=[line.split()[2] for line in tm]
    m=tm.count(name1)
    if name1==name2: return m,m
    return m,tm.count(name2)

def calculate():
    p=Pool()
    result=[]
    with open(file) as f:
        while True:
            tmp=list(islice(f,9,531))
            if not tmp: break
            result.append(p.apply_async(get_hist,[tmp]))
            tm=tmp
    result=[r.get() for r in result]
    m,n=get_num(tm)
    return normalise(result,side,dr,m,n)

def save(r,rdf,rcn):
    txt="\n".join(["%s %s"%(i,j) for i,j in zip(*[r,rdf])])
    txt1="\n".join(["%s %s"%(i,j) for i,j in zip(*[r,rcn])])
    with open("rdf_%s%s.dat"%(name1,name2),"w") as f:
        f.write(txt+"\n\n")
    with open("rcn_%s%s.dat"%(name1,name2),"w") as g:
        g.write(txt1+"\n\n")

def plot_result(r,rdf,rcn):
    plt.plot(r,rdf)
    plt.xlabel("$r(\AA)$", fontsize=18)
    plt.ylabel("$g_{{{}{}}}(r)$".format(name1,name2), fontsize=18)
    plt.title("$Radial\ Distribution\ function$")
    plt.figure()
    plt.xlabel("$r(\AA)$", fontsize=18)
    plt.ylabel("$N_{{{}{}}}(r)$".format(name1,name2), fontsize=18)
    plt.title("$Running\ Coordination\ Number$")
    plt.plot(r,rcn)
    plt.show()

file,name1,name2,side,dr,bins,plot=parse()  
if __name__ == '__main__':
    r,rdf,rcn=calculate()
    save(r,rdf,rcn)
    if plot:
        plot_result(r,rdf,rcn)
from itertools import islice,izip
from msd import msd
import matplotlib.pyplot as plt
from multiprocessing import Pool
import argparse,time

def parse():                                        #function to parse command-line arguments   
    parser = argparse.ArgumentParser(description='Mean Square Displacement Calculation')
    parser.add_argument("f",    type=str,  help="Name of the trajectory file")
    parser.add_argument("atoms",type=int,  help="Total number of atoms")
    parser.add_argument("name", type=str,  help="Name of the element",nargs="?",default="all")
    args=parser.parse_args() 
    return [args.f,args.atoms,args.name]

def get(data,m):                                    #return [[t0,msd(t0)],[t1,msd(t1)]..tm,msd(tm)]
    return msd(data,m)    

def calculate(name): 
    data,n=[],0
    with open(file) as f:
        while True:
            tmp=list(islice(f,9,9+nat))             #read one time frame from file
            if not tmp: break
            # tmp =[map(float,line.split()[3:6]) for line in tmp]
            tmp =(line.split()[2:6] for line in tmp)
            tmp=((map(float,i[1:])) for i in tmp if i[0]==name)
            n+=1
            data.append(tmp)
    data=izip(*data)                                #get [r(t1),r(t2)..r(tn)] for each particle
    m=int(n/2.0)
    pl=Pool()                                       #create a pool of worker processes
    res=[pl.apply_async(get,[i,m]) for i in data]   #load the function call to a worker in parallel
    res=[r.get() for r in res]                      #[msd(t0),msd(t1)..msd(tm)] for each particle seperately
    x=[i*100 for i in xrange(m+1)]                  #list of [0,1,2...m]
    avg=[sum(i)/len(i) for i in izip(*res)]         #[msd(t0),msd(t1)..msd(tm)] averaged over all particle
    save(x,avg,name)
    return x,avg,res



    
def save(a,b,name):                                      #save the data
    txt="\n".join("%s %s"%(i,j) for i,j in zip(a,b))
    with open("msd_%s.dat"%name,"w") as f:
        f.write(txt)
    

file,nat,name=parse()
if __name__ == '__main__':
    s=time.time()
    handles=[]
    if name=="all":                                 # plot onle averages for each type of particle
        x,avg,_=calculate("H")
        handles.append(plt.plot(x,avg,"-",label="H")[0])    
        x,avg,_=calculate("O")
        handles.append(plt.plot(x,avg,"-",label="O")[0])    
        x,avg,_=calculate("Na")
        handles.append(plt.plot(x,avg,"-",label="Na")[0])    
        x,avg,_=calculate("Cl")
        handles.append(plt.plot(x,avg,"-",label="Cl")[0])
    else:                                           #plot each particle and their average for type "name"
        x,avg,res=calculate(name)
        for i in res:
            plt.plot(x,i,"r-",lw=.5)
        handles.append(plt.plot(x,avg,"r-",label=name)[0])
        handles.append(plt.plot(x,avg,"g-",label="Average")[0])
    plt.yscale('log')
    plt.xscale('log')
    plt.ylabel("$MSD(\\tau)\ (in\ \AA^2)$")
    plt.xlabel("$\\tau\ (in\ fs) $")
    plt.title("Mean Square Displacement Plot")
    plt.legend(handles=handles) 
    print "Time needed", time.time()-s
    # plt.show()
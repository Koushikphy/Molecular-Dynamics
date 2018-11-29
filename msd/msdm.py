from itertools import islice
from msd import msd,msd1,msd2
import matplotlib.pyplot as plt
from multiprocessing import Pool
import time
file=r"E:\Project\\new_sol\25\\nve\sol_nve_u.traj"
nat=818
side=20.3465
num_element=nat


def get_m(data,a):
    start,end,total=a
    return msd2(data,start,end,total)    #return t,msd(t)

def calculate1():
    data=[]
    with open(file) as f:
        for i in xrange(1000):
            tmp=list(islice(f,9,9+nat))
            if not tmp: break
            tmp =[map(float,line.split()[3:6]) for line in tmp]
            # tmp =[line.split()[2:6] for line in tmp]
            # tmp=[[map(float,i[1:])] for i in tmp if i[0]=="O"]
            data.append(tmp)
    n=len(data)
    m=int(n/2.0)
    result=msd1(data,m)      #return list of msd(t) for t=1,2,3...m
    return m,result

# def calculate1(name):
#     data,n=[],0
#     with open(file) as f:
#         while True:
#             tmp=list(islice(f,9,9+nat))
#             if not tmp: break
#             tmp =[map(float,line.split()[3:6]) for line in tmp]
#             # tmp =(line.split()[2:6] for line in tmp)
#             # tmp=((map(float,i[1:])) for i in tmp if i[0]==name)
#             n+=1
#             data.append(tmp)
#     data=izip(*data)
#     m=int(n/2.0)
#     result=[msd3(atom,m) for atom in data]
#     return m,result



def calculate():
    data=[]
    pl=Pool()
    with open(file) as f:
        for i in xrange(100):
            tmp=list(islice(f,9,9+nat))
            if not tmp: break
            tmp =[map(float,line.split()[3:6]) for line in tmp]
            data.append(tmp)
    m=int(len(data)/2.0)
    p=int(m/4)
    a=[1,p,p]
    b=[p+1,2*p,p]
    c=[2*p+1,3*p,p]
    d=[3*p+1,m,(m-3*p)]
    result=[]
    result.append(pl.apply_async(get_m,[data,a]))
    result.append(pl.apply_async(get_m,[data,b]))
    result.append(pl.apply_async(get_m,[data,c]))
    result.append(pl.apply_async(get_m,[data,d]))
    result=[j for i in result for j in i.get()]
    return result


if __name__ == '__main__':
    s=time.time()
    # res=calculate()
    m,y=calculate1()
    x=[i*100 for i in xrange(1,m+1)]
    for i in xrange(num_element):
        plt.plot(x,y[i],"r-",lw=.5)
    plt.yscale('log')
    plt.xscale('log')
    plt.ylabel("$MSD(\\tau)\ (in\ \AA^2)$")
    plt.xlabel("$\\tau(in\ ps) $")
    plt.show()

    # txt="\n".join(str(i) for i in res)
    # with open("result","w") as f:
    #     f.write(txt)

    # a,b=zip(*res)

    # plt.plot(a,b,"ro-",lw=1)
    print time.time()-s

    # plt.show()




# if __name__ == '__main__':
#     cpucount = cpu_count()
#     pool = Pool(processes=cpucount)  # start worker processes

#     # RDF setup
#     pdbdir = argv[1] # directory with PDB files (eg. "./pdb-ex-10k/")
#     boxsize = (51.16866708780796, 12.79216677195199, 12.79216677195199)
#     numbins = 300 # number of bins
#     outfile = pdbdir + "rdf.out"
    
#     # calculate RDF from center of mass of molecule
#     numAT = 3 # number of atoms in molecule
#     cm = False # default to false
#     #if len(argv) == 3:
#     #    cm = argv[2]
    
#     # prepare lists of pdb files
#     files = listdir(pdbdir) # list of files
#     nfiles = 0 # number of total pdb files
#     pdbs = [] # list of pdb files
#     for f in files:
#         if not f.endswith(".pdb"):
#             continue # skip other files
#         nfiles += 1
#         pdbs.append(f)
    
#     if cpucount > nfiles:
#         cpucount = nfiles
    
#     files_oncpu = int(floor(nfiles / cpucount))
#     dist_pdbs = [] # pdbs to distribute
#     p = 0
#     for i in range (cpucount):
#         tmp = []
#         for j in range(files_oncpu):
#             tmp.append(pdbs[p])
#             p += 1
#         dist_pdbs.append(tmp)
    
#     # distribute remaining files
#     rem = nfiles - (files_oncpu * cpucount)
#     if (rem > 0):
#         for i in range(rem):
#             dist_pdbs[i].append(pdbs[(files_oncpu * cpucount)+i])
    
#     # distribute lists to pools
#     results = [ pool.apply_async(makeRDF, [pdbdir,dist_pdbs[i],boxsize,numbins,cm,numAT])
#                 for i in range(len(dist_pdbs)) ]
    
#     # sum up results
#     hist = [0]*(numbins+1) # initialize list
#     cnt = 0;
#     for i in range(len(dist_pdbs)):
#         result = results[i].get()
#         cnt += result[1]
#         for j in range(numbins+1):
#             hist[j] += result[0][j]
    
#     #print hist

#     rdf = normalize(hist, boxsize, cnt)
    
#     print "Writing output file ... " +outfile+ " (" + time.strftime("%c") +")"
#     outfile = open(outfile, "w")
#     for r in sorted(rdf.iterkeys()): # sort before writing into file
#         outfile.write("%15.8g %15.8g\n"%(r, rdf[r]))
#     outfile.close()

from rdfs import homo,hetero
import matplotlib.pyplot as plt
import time,sys,argparse

s=time.time()
def parse():
    parser = argparse.ArgumentParser(description='Radial Distribution Function Calculator')
    parser.add_argument("--f",  type=str,  help="Name of the trajectory file")
    parser.add_argument("name1",type=str,  help="Name of 1st element")
    parser.add_argument("name2",type=str,  help="Name of 2nd element")
    parser.add_argument("side" ,type=float,help="Size of the box")
    parser.add_argument("--dr" ,type=float,help="Size of the bins",default=.1)
    args=parser.parse_args()  
    bins=int(args.side/(2.0*args.dr))
    return [args.f,args.name1,args.name2,args.side,args.dr,bins]



def file_read():
    with open(file) as f:
        d=filter(None,f.read().split("\n"))

    atoms,data=int(d[0]),d[2:]

    return [filter(None,i.split()) for i in data]

def calculate(data):
    data1=[map(float,i[1:4]) for i in data if i[0]==name1]
    if name1==name2:
        r,rdf=homo(data1,side,dr,int(side/(2.0*dr)))
    else:
        data2=[map(float,i[1:4]) for i in data if i[0]==name2]
        r,rdf=hetero(data1,data2,side,dr,int(side/(2.0*dr)))
    return r,rdf

def plot_result(dat):
    plt.plot(*dat)
    plt.xlabel("$r$", fontsize=18)
    plt.ylabel("$g_{{{}{}}}(r)$".format(name1,name2), fontsize=18)
    plt.title("$Radial\ Distribution\ function$")
    print time.time()-s
    plt.show()

if __name__ == '__main__':
    file,name1,name2,side,dr,bins=parse()
    plot_result(calculate(file_read()))
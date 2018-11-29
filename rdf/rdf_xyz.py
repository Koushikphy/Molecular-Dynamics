from math import pi,sqrt 
from rdfs import homo,hetero
import matplotlib.pyplot as plt
import time
s=time.time()

file="C:\Users\Koushik Naskar\Desktop\Project\New folder\water.xyz"

with open(file) as f:
    d=filter(None,f.read().split("\n"))

atoms,side,data=int(d[0]),float(d[1]),d[2:]


dr=.1
bins=int(side/(2.0*dr))
print bins
hist=[0]*bins
rdf=[]
rcn=[]
data =[filter(None,i.split()) for i in data]

data=[[i[0],float(i[1]),float(i[2]),float(i[3])] for i in data]
coord=[[x,y,z] for _,x,y,z in data]

datah=[[x,y,z] for n,x,y,z in data if n=="H"]
datao=[[x,y,z] for n,x,y,z in data if n=="O"]
count=0

#get g(r) between coord and coord1
coord=datah
coord1=datah
no=len(coord)
nh=len(coord1)
count=0
for i in range(len(coord)):
    for j in range(len(coord1)):
        if (coord==coord1) & (i==j): continue
        xx=coord[i][0]-coord1[j][0]
        yy=coord[i][1]-coord1[j][1]
        zz=coord[i][2]-coord1[j][2]
        xx-=side*round(xx/side)
        yy-=side*round(yy/side)
        zz-=side*round(zz/side)
        r=sqrt(xx*xx + yy*yy + zz*zz)
        bi=int(r/dr)
        if bi<bins:
            hist[bi]+=1
        else:
            count+=1




const1= 4*pi*dr**3
const2=const1/3.0
cn=0
for i in range(bins):
    rr=(i+.5)*dr
    vs=const2 + const1*i*(i+1)        #volume of shell
    rho=atoms/side**3                          #rho=N/V
    val1=hist[i]*side**3/(nh*no)          #atoms is actually number of samples
    val=val1/vs
    cn+=hist[i]/float(no)
    rdf.append([rr,val])
    rcn.append([rr,cn])



p,q=zip(*rdf)
x,y=zip(*rcn)
z=[1 for val in x]
# plt.plot(p,q)
plt.figure()
plt.plot(x,y)
print time.time()-s
plt.show()
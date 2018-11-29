from math import sqrt


M_Na=22.989999
M_Cl=35.450000 
M=1.5                     #concentration
bl=2					#minimum space
boxl=19.818  		#box length
con=6.023e-4*boxl**3


N_NaCl= int(M*con) # just an approximation actual is calculated from volume after npt
print N_NaCl


file="C:\Users\Koushik Naskar\Desktop\Project\water\\nve\system_nve.xyz"

with open(file) as f:
    d=filter(None,f.read().split("\n"))
atoms,unw,data=int(d[0]),d[1],d[2:]
mol=atoms/3
data =[i.split() for i in data]
data=[[i[0],float(i[1]),float(i[2]),float(i[3])] for i in data] #list of all atoms
wdata=[[data[3*i],data[3*i+1],data[3*i+2]]for i in xrange(mol)] #list of water molecule




#checks and returns water molecule lying within radius "bl" from cord
def check(cord):		
	ch=[]
	for i in wdata:
		for j in i:
			dis=sqrt(sum((a-b)**2 for a,b in zip(cord,j[1:])))
			if dis<bl: 
				ch.append(i)
				break
	return ch




#following loop takes 2 water molecule place Na and Cl atom in their oxygen positions and 
# remove the water molecules, then checks for any close water molecule and removes them

for i in xrange(N_NaCl):
    cord=wdata[0][0][1:]
    cord1=wdata[1][0][1:]
    nal=["Na"]+cord
    cll=["Cl"]+cord1
    wdata.pop(0)
    wdata.pop(0)
    c_na=check(cord)
    for mol in check(cord1):
    	wdata.remove(mol)
    for mol in check(cord):
    	wdata.remove(mol)
    wdata.append([nal,cll])



all_atoms=[p for i in wdata for p in i]  #list of all atoms altogether
txt="\n".join("{}{: 12.6f}{: 12.6f}{: 12.6f}".format(*i) for i in all_atoms)
txt="{}\n{}\n".format(len(all_atoms),unw)+txt
with open("sol.xyz","w") as f:
    f.write(txt)

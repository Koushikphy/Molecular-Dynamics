from random import uniform
M=#provide a concentration
boxl=19.79388                #box length
N_mol= int(M*6.023e-4*boxl**3)
#actual concentration due to integer conversion
M=N_mol/(6.023e-4*boxl**3)

#open the LAMMPS dump data file
with open("water.data") as f:
    dat=f.read()

#generate coordinate for Na and Cl atom
def gen_cord():
    x=uniform(xlo,xhi)
    y=uniform(ylo,yhi)
    z=uniform(zlo,zhi)
    co=[x,y,z]
    if co in coord:
        gen_cord()
    return co

def get_line(i,x,y):
    nal=[768+2*i+1,1.0,3.0,1.0]+x+[0,0,0]
    cll=[768+2*i+2,1.0,4.0,-1.0]+y+[0,0,0]
    return nal,cll

_,dat=dat.split("angle types")
data,dat=[i.strip() for i in dat.split("Masses")]

in_data="""
LAMMPS data file ,NaCl molecule {}, concentration {} M 

{} atoms
4 atom types
512 bonds
1 bond types
256 angles
1 angle types


""".format(N_mol,M,768+2*N_mol)
in_data+=data

#extract coordinate,bonds,angle etc. information form dump file
lim=[map(float,i.split()[:2]) for i in data.split("\n")]
xlo,xhi,ylo,yhi,zlo,zhi=[j for i in lim for j in i]

masses,dat=dat.split("Atoms # full")
masses=masses.strip()+"\n3 22.990000\n"+"4 35.450001\n"
atoms,dat=dat.split("Velocities")
_,dat=dat.split("Bonds")
bonds,angles=dat.split("Angles")      

atoms=atoms.strip().split("\n")
atoms=[i.split() for i in atoms]

atoms=[map(int,i[:3])+map(float,i[3:7]) + map(int,i[7:]) for i in atoms]
coord=[i[4:7] for i in atoms]

#generate N_mol number of Na and Cl atom
for i in xrange(N_mol):
    a=gen_cord()
    b=gen_cord()
    coord+=[a,b]
    atoms+=get_line(i,a,b)

atoms="\n".join(" ".join(str(j) for j in i) for i in atoms)

#create the data file for solution
txt="""{}

Masses

{}
Atoms # full

{}

Bonds

{}

Angles

{}
""".format(in_data.strip(),masses,atoms,bonds.strip(),angles.strip())

with open("sys.dat","w") as f:
    f.write(txt)
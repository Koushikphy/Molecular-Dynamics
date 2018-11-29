from random import uniform
from math import sqrt,cos,sin,pi

N_molecule=256  #provide number of water molecule
rho=1.0
D_oh=1.0
ANG_hoh=109.47*pi/180
RCUT_oo=2.7

H_mass=1.008 
O_mass=15.9994 
Total_mass=N_molecule*(H_mass*2+O_mass)/6.023e23
box_l=1e8*(Total_mass/rho)**(1/3.0)

#generate 2 hydrogen atom position for a water molecule, assuming 
#a oxygen at origin and plane of molecule in x-y plane for simplicity 
def gen_cord():
    phi=uniform(0,2*pi)
    z1=0                
    x1=D_oh*cos(phi)
    y1=D_oh*sin(phi)
    x2=cos(ANG_hoh)*x1-sin(ANG_hoh)*y1 #rotate x1,y1 by ANG_hoh
    y2=sin(ANG_hoh)*x1+cos(ANG_hoh)*y1
    z2=0
    return [x1,y1,z1],[x2,y2,z2]

#check for any oxygen atom within cut_off RCUT_oo
def check_no_ox(tmp):
    for j in cord:
        if sqrt(sum([(i-j)**2 for i,j in zip(j,tmp)]))<=RCUT_oo:
            return False
    return True

cord=[]
atom_list="%s\n%s\n"%(3*N_molecule,box_l)

#generate N_molecule number water molecule
for i in xrange(N_molecule):
    while True:
        tmpo=[uniform(-box_l/2,box_l/2) for _ in xrange(3)]
        if check_no_ox(tmpo): break 
    cord.append(tmpo)
    atom_list+="{} {} {} {} {}\n".format(*(["O"]+tmpo+[O_mass]))

    tmp,tmp1=gen_cord()

    tmph1=[i+j for i,j in zip(tmpo,tmp)]
    atom_list+="{} {} {} {} {}\n".format(*(["H"]+tmph1+[H_mass]))

    tmph2=[i+j for i,j in zip(tmpo,tmp1)]
    atom_list+="{} {} {} {} {}\n".format(*(["H"]+tmph2+[H_mass]))


with open("watertest.xyz","w") as f:
    f.write(atom_list)

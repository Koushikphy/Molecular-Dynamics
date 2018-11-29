file="sol.xyz"   #get it from nacl_build_inplace_water.py
with open(file) as f:
    d=filter(None,f.read().split("\n"))

atoms,_,data=int(d[0]),d[1],d[2:]
side=19.818
sid=side/2.0
mol=atoms/3
data =[filter(None,i.split()) for i in data]

data=[[i[0],float(i[1]),float(i[2]),float(i[3])] for i in data]

w_mol=[i[0] for i in data].count("O")

coord=[i[1:] for i in data]
wdata=[[coord[3*i],coord[3*i+1],coord[3*i+2]]for i in xrange(w_mol)]

new_wdata=[]


for o,h1,h2 in wdata:
    ox,oy,oz=o
    h1x,h1y,h1z=h1
    h2x,h2y,h2z=h2
    if ox-h1x>sid:
        h1x+=side
    elif h1x-ox>sid:
        h1x-=side
    if oy-h1y>sid:
        h1y+=side
    elif h1y-oy>sid:
        h1y-=side
    if oz-h1z>sid:
        h1z+=side
    elif h1z-oz>sid:
        h1z-=side

    if ox-h2x>sid:
        h2x+=side
    elif h2x-ox>sid:
        h2x-=side
    if oy-h2y>sid:
        h2y+=side
    elif h2y-oy>sid:
        h2y-=side
    if oz-h2z>sid:
        h2z+=side
    elif h2z-oz>sid:
        h2z-=side

    h1=["H",h1x,h1y,h1z]
    h2=["H",h2x,h2y,h2z]
    o=["O",ox,oy,oz]
    new_wdata+=[o,h1,h2]




new_wdata+=data[w_mol*3:]  #list of all atoms altogether

txt="\n".join("{}{: 12.6f}{: 12.6f}{: 12.6f}".format(*i) for i in new_wdata)
txt="{}\n{}\n".format(len(new_wdata),side)+txt
with open("sol1_u.xyz","w") as f:
    f.write(txt)

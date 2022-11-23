from vpython import *
import math as m
import matplotlib.pyplot as plt

# Stars interacting gravitationally
Vlist=[]
Rlist=[]
Tlist=[]
# Bruce Sherwood
scene.width = 1400
scene.height =700

# Display text below the 3D graphics:
scene.title = "Gravity assist simulation"

scene.caption = """    """

Nstars = 5  # change this to have more or fewer objects

G = 6.67e-11 # Universal gravitational constant

# Constant values
Msun = 1.99E30 #Mass of sun
Rsun = 5e10 #Radius of sun in our model (It's chosen in such a way that we can see it in our simulation)
L = 1.5e11 #Astronomical Unit
q= m.pi/180 # Degree to radian

scene.range = 8.5*L
#scene.forward = vec(-1,-1,-1)

#We define axes in here
xaxis = curve(color=color.gray(0.5), radius=3e8)
xaxis.append(vec(0,0,0))
xaxis.append(vec(L,0,0))
yaxis = curve(color=color.gray(0.5), radius=3e8)
yaxis.append(vec(0,0,0))
yaxis.append(vec(0,L,0))
zaxis = curve(color=color.gray(0.5), radius=3e8)
zaxis.append(vec(0,0,0))
zaxis.append(vec(0,0,L))
# We choose color of every object here
Stars = []
star_colors = [color.red, color.yellow, color.blue,
              color.orange, color.white, color.cyan]
#Momentum
psum = vec(0,0,0)

#Now we define our objects(their shape, initial position, mass, initial momentum and color)

#sun
star= sphere(pos=L*vector(0,0,0), make_trail=True, retain=500, trail_radius=Rsun/5)
star.radius = Rsun
star.mass = Msun
star.momentum=vector(0,0,0)*star.mass
star.color = star.trail_color = star_colors[1 % 6]
Stars.append( star )
psum = psum + star.momentum

#earth
star= sphere(pos=L*vector(0,1,0), make_trail=True, retain=500, trail_radius=Rsun/6)
star.radius = 4*Rsun/7
star.mass = 5.97e24
star.momentum=vector(sqrt(G*Msun/L),0,0)*star.mass
star.color = star.trail_color = star_colors[2 % 6]
Stars.append( star)
psum = psum + star.momentum  

#jupiter
tetaj=-13.8
star= sphere(pos=5.2*L*vector(m.cos(tetaj*q),m.sin(tetaj*q),0), make_trail=True, retain=5000, trail_radius=3*Rsun/8)
#star= sphere(pos=5.2*L*vector(0,1,0), make_trail=True, retain=5000, trail_radius=Rsun/5)
star.radius =3*Rsun/4
star.mass = 1.9e27
star.momentum=sqrt(G*Msun/(5.2*L))*vector(m.sin(tetaj*q),-m.cos(tetaj*q),0)*star.mass
#star.momentum=sqrt(G*Msun/(5.2*L))*vector(1,0,0)*star.mass
star.color = star.trail_color = star_colors[3 % 6]
Stars.append( star )
psum = psum + star.momentum  



#saturn
tetas=-53.402
star= sphere(pos=9.5*L*vector(m.cos(tetas*q),m.sin(tetas*q),0), make_trail=True, retain=5000, trail_radius=Rsun/3)
star.radius = 2*Rsun/3
star.mass = 5.68e26
star.momentum=sqrt(G*Msun/(9.5*L))*vector(m.sin(tetas*q),-m.cos(tetas*q),0)*star.mass
star.color = star.trail_color = star_colors[4 % 6]
Stars.append( star)
psum = psum + star.momentum  

#uranus
tetau=-109.566
star= sphere(pos=19.8*L*vector(m.cos(tetau*q),m.sin(tetau*q),0), make_trail=True, retain=5000, trail_radius=Rsun/4)
star.radius = Rsun/2
star.mass = 8.68e25
star.momentum=sqrt(G*Msun/(19.8*L))*vector(m.sin(tetau*q),-m.cos(tetau*q),0)*star.mass
star.color = star.trail_color = star_colors[5 % 6]
Stars.append( star )
psum = psum + star.momentum  

#satellite
star= sphere(pos=L*vector(0,1+1e-2,0), make_trail=True, retain=5000, trail_radius=Rsun/12)
star.radius = 2*Rsun/5
star.mass = 1e4
star.momentum=vector(39000,0,0)*star.mass
star.color = star.trail_color = star_colors[6 % 6]
Stars.append( star )
psum = psum + star.momentum  


dt = 100000 #runtime duration
T=0

#Now we compute forces on each object in following function and then update their momentum
def computeForces():
    global Stars
    N = len(Stars)
    for i in range(N):
        si = Stars[i]
        if si is None: continue
        F = vec(0,0,0)
        pos1 = si.pos
        m1 = si.mass
        for j in range(N):
            if i == j: continue
            sj = Stars[j]
            if sj is None: continue
            r = sj.pos - pos1
            F = F + (G*m1*sj.mass*r.hat / mag(r)**2)
        si.momentum = si.momentum + F*dt

i=0
while True:

    rate(100)
    # Calculate all forces on all stars
    computeForces()
    T=T+dt
    # Having updated all momenta, now we update all positions
    for star in Stars:
        star.pos = star.pos + star.momentum*(dt/star.mass)
        if star.mass==1e4:
            Vlist.append(mag(star.momentum)/(1e3*star.mass))
            Rlist.append(mag(star.pos)/1.5e11)
            Tlist.append(T/(86400*365.25))
    if Rlist[i]>25:
        break
    i=i+1
    
    
# Now we plot velocity-time graph of our satellite
plt.scatter(Rlist,Vlist)    
plt.ylabel("Velocity (km/s)")  
plt.xlabel("Distance from sun (AU)")      
plt.show()
plt.scatter(Tlist,Vlist)
plt.ylabel("Velocity (km/s)")  
plt.xlabel("Time (yr)") 
plt.show()
  
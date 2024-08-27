# user input ----------------------
iarea   = 1
im0     = 10
nmax    = 12000
ft      = 0.3
rfd     = 0.15
pout    = 80000.0  
rfpout  = 1.0      
pstart  = 90000.0
vxstart = 100.0
rhostart = 1.0
pt1     = 100000.0
tt1     = 300.0                  
cpzucv  = 1.4
rg      = 287.06
ermax   = 0.000001
anival  = 300
anistep = 200
# end user input ------------------

import matplotlib.pyplot as plt
import numpy as np
import matplotlib.animation as animation


# constant thermodynamic values
cpzurg = cpzucv/(cpzucv-1.0)
cvzurg = 1.0/(cpzucv-1.0)
rgzucp = 1.0/cpzurg
rgzucv = 1.0/cvzurg  
  
x      = []
r      = []
r1     = []    
vx     = []
p      = []  
qe     = [] 
qm     = []   
rho    = [] 
moment = []
energy = []  
flux   = []  
fblc   = []
fblp   = []        
rla    = []    
step   = []  
aix    = []     
a      = []  
rla    = []  
mach   = []  
maca1  = []
vxold  = []   

#---------------------------------------------------------------
if iarea == 1:
# nozzle contour
    with open('./nozzle_coord.txt', 'r') as reader:
        lines = reader.readlines()
        im = -1
        for line in lines:
    #      x.append(float(line.split()[1]))
           r1.append(float(line.split()[0]))
           maca1.append(float(line.split()[2]))            
           im = im+1
    reader.close() 

    i2 = -1
    for i in range(0,im+1,2):
        i2 = i2+1
        r.append(r1[i])

    im = i2    
#--------------------------------------------------------------- 
else:       
    
#  constant area contour ==> 0D problem
   im = im0
   for i in range(0,im+1):
       r.append(1.0) 

#---------------------------------------------------------------        

print('im = ',im)

# 2D array for movies
t = (nmax+1,im+1)
y2 = np.zeros(t)
y1 = np.zeros(im+1)    


for i in range(0,im+1):
    x.append(float(i)/float(im))     

# expand arrays with initial conditions  
for i in range(0,im+1):  
    vx.append(vxstart) 
    vxold.append(vxstart)
    p.append(pstart)    
    rho.append(rhostart)   
    moment.append(1.0)      
    energy.append(1.0) 
    qe.append(1.0)
    qm.append(1.0)
    aix.append(1.0)                             
    flux.append(1.0)
    fblc.append(1.0)
    fblp.append(1.0)
    step.append(1.0)
    flux.append(1.0)
    a.append(1.0)  
    rla.append(1.0)  
    mach.append(0.0) 
                                             

# areas and conservative variables
for i in range(0,im+1):  
    aix[i]    = 2.0*r[i]
    moment[i] = rho[i]*vx[i]
    energy[i] = (p[i]/rgzucv+0.5*rho[i]*vx[i]**2) 

# initialization massratio
massratio = 100.0 

# begin iteration loop
        
for n in range(0,nmax+1):   

    for i in range(0,im+1):
        # calculate speed of sound
        aq       = max(100.0,cpzucv*p[i]/rho[i])
        a[i]     = (aq)**0.5
        vxold[i] = vx[i]   


    # calculation of local timesteps     

    for i in range(0,im):
         
        # volume mean values for speed of sound and velocity 
        rav = 0.5*(a[i]+a[i+1])
        vxv = 0.5*(vx[i]+vx[i+1])
        # area components at center of single cell volume
        rix = 0.5*(aix[i]+aix[i+1])
        # product at center of volume of
        # (speed of sound + velocity) times area-components
        ctx = rav*rix
        # scalar product x-area*velocity
        cai = rix*vxv
        # delta(time)/delta(volume) for i-direction
        rdti   = abs(cai) + ctx
        rla[i] = rdti      

    # inlet and outlet weighting
    rla[0]    = 2.0*rla[0]
    rla[im-1] = 2.0*rla[im-1]

    # sum up for nodes    
    for i in range(1,im):
        step[i] = rla[i-1]+rla[i]
    # boundary values
    step[0]  = rla[0]   
    step[im] = rla[im-1]

    # invert to correct step values delta(time)/delta(volume)
    for i in range(0,im+1):
        step[i] = 1.0/step[i]


    #  update rho
    #  **********    

    # fluxes rho
    for i in range(0,im+1):  
        flux[i] = rho[i]*vx[i]*aix[i] 

    # fluxbalance rho single elements
    for i in range(0,im):  
        fblc[i] = flux[i]-flux[i+1]

    # sum up for nodes    
    for i in range(1,im):
        fblp[i] = fblc[i-1]+fblc[i]
    # boundary values
    fblp[0]  = fblc[0]   
    fblp[im] = fblc[im-1]

    # average for single element             
    for i in range(0,im):
        qe[i] = 0.5*(rho[i+1]+rho[i])

    # sum up for nodes    
    for i in range(1,im):
        qm[i] = 0.5*(qe[i-1]+qe[i])
    # boundary values
    qm[0]  = qe[0]   
    qm[im] = qe[im-1]  

    # update rho for new timestep
    for i in range(0,im+1):  
        rho[i]  = (1.0-rfd)*rho[i]+rfd*qm[i]+ft*fblp[i]*step[i]

 
    #  update moment
    #  *************    

    # fluxes moment
    for i in range(0,im+1):  
        flux[i] = (moment[i]*vx[i]+p[i])*aix[i] 

    # fluxbalance moment single elements
    for i in range(0,im):  
        fblc[i] = flux[i]-flux[i+1]+0.5*(p[i+1]+p[i])*(aix[i+1]-aix[i])

    # sum up for nodes    
    for i in range(1,im):
        fblp[i] = fblc[i-1]+fblc[i]
    # boundary values
    fblp[0]  = fblc[0]   
    fblp[im] = fblc[im-1]

    # average for single element             
    for i in range(0,im):
        qe[i] = 0.5*(moment[i+1]+moment[i])

    # sum up for nodes    
    for i in range(1,im):
        qm[i] = 0.5*(qe[i-1]+qe[i])
    # boundary values
    qm[0]  = qe[0]   
    qm[im] = qe[im-1]  

    # update rho for new timestep
    for i in range(0,im+1):  
        moment[i]  = (1.0-rfd)*moment[i]+rfd*qm[i]+ft*fblp[i]*step[i]


    #  update energy
    #  *************    

    # fluxes energy
    for i in range(0,im+1):  
        flux[i] = (energy[i]+p[i])*vx[i]*aix[i] 

    # fluxbalance energy single elements
    for i in range(0,im):  
        fblc[i] = flux[i]-flux[i+1]

    # sum up for nodes    
    for i in range(1,im):
        fblp[i] = fblc[i-1]+fblc[i]
    # boundary values
    fblp[0]  = fblc[0]   
    fblp[im] = fblc[im-1]

    # average for single element             
    for i in range(0,im):
        qe[i] = 0.5*(energy[i+1]+energy[i])

    # sum up for nodes    
    for i in range(1,im):
        qm[i] = 0.5*(qe[i-1]+qe[i])
    # boundary values
    qm[0]  = qe[0]   
    qm[im] = qe[im-1]  

    # update rho for new timestep
    for i in range(0,im+1):  
        energy[i]  = (1.0-rfd)*energy[i]+rfd*qm[i]+ft*fblp[i]*step[i]


    #  update primaries
    #  ****************    

    for i in range(0,im+1):  
        vx[i]    = moment[i]/rho[i]
        p [i]    = (energy[i]-0.5*rho[i]*vx[i]**2)*rgzucv
        mach[i]  = vx[i]/(cpzucv*p[i]/rho[i])**0.5            
    
            
    #  inlet boundary
    #  **************

    p[0]      = p[1]
    ts1       = tt1*(p[0]/pt1)**rgzucp
    rma       = (2.0*cvzurg*(tt1/ts1-1.0))**0.5
    rho[0]    = p[0]/(rg*ts1)
    vx[0]     = rma*(cpzucv*rg*ts1)**0.5
    moment[0] = rho[0]*vx[0]
    energy[0] = p[0]/rgzucv+0.5*rho[0]*vx[0]**2 


    #  outlet boundary
    #  ***************

    p[im]   = (1.0-rfpout)*p[im]+rfpout*pout


    # maximum vx-changes
    vxc  = 0.0
    ivxc = 0
    for i in range(0,im+1):  
        if abs(vx[i]-vxold[i]) > vxc:
            vxc = abs(vx[i]-vxold[i])
            ivxc = i

    if vxc > 100.0:
        print('claculation diverges!')                   
        break  
                  
    # convergence check
    massflow_in  = rho[0]*vx[0]*aix[0] 
    massflow_out = rho[im]*vx[im]*aix[im] 
    massratio    = massflow_out/massflow_in 
    if abs(1.0-massratio) < ermax:
        print('convergence criterium met')
        print('Vx Eintritt: ',"%.3f" % vx[0],) 
        print('Vx Austritt: ',"%.3f" % vx[im],) 
        print('p  Eintritt: ',"%.3f" % p[0],) 
        print('p  Austritt: ',"%.3f" % p[im],)        
        break   

    print('it ',n,' mratio ',"%.6f" % massratio,' vxc ',"%.6f" % vxc,' i ',ivxc, 'rfd ',rfd)
       

    # for movie

    nstep = n  

    for i in range(0,im+1): 
        y2[n][i] = mach[i]
   
    # end iteration loop
              
        
# result values
        
print('  ')          
print('Results: ')    
print('  ')  
massflow_in  = rho[0]*vx[0]*aix[0] 
massflow_out = rho[im]*vx[im]*aix[im]       
print('Inlet  Massflow: ',massflow_in)    
print('Outlet Massflow: ',massflow_out)    
print('Massflow Ratio Out/in: ',massflow_out/massflow_in)   
print('  ')      


# plot
   
plt.plot(x,mach,color='red',linestyle='',marker='x')
plt.grid(True)
plt.axis([0.0,1.0,0.0,1.80])
plt.title('Vx with pout='+str(pout)+' after '+str(n)+' timesteps')
plt.xlabel('x-Koordinate')
plt.ylabel('Vx')    
plt.show()


# calculate total pressure ratio

# total pressure and temperature at inlet

# inlet
ts        = p[0]/(rho[0]*rg)
aq        = cpzucv*rg*ts
rmaq      = vx[0]**2/aq
tt        = ts*(1.0+0.5*(cpzucv-1.0)*rmaq)
pt        = p[0]/((ts/tt)**cpzurg)
ptinlt    = pt
ttinlt    = tt

# outlet
ts    = p[im]/(rho[im]*rg)
aq    = cpzucv*rg*ts
rmaq  = vx[im]**2/aq
tt    = ts*(1.0+0.5*(cpzucv-1.0)*rmaq)
pt    = p[im]/((ts/tt)**cpzurg)
ptout = pt
ttout = tt

ptv   = ptout/ptinlt
ttv   = ttout/ttinlt
             

result_dat = open('Demo1D.output','w')  
result_dat.write("results at inlet and outlet\n") 
result_dat.write("***************************\n")   
result_dat.write("   \n")  
result_dat.write("Input nmax                         = %d\n" % nmax)    
result_dat.write("Input pt1                          = %.2f [Pa]  \n" % pt1)    
result_dat.write("Input tt1                          = %.2f [K]   \n" % tt1)       
result_dat.write("Input pout                         = %.2f [Pa]  \n" % pout)  
result_dat.write("Input ft                           = %.2f       \n" % ft)    
result_dat.write("Input rfd                          = %.2f       \n" % rfd)         
result_dat.write("Input ermax                        = %.6f   \n" % ermax)       
 
result_dat.write("   \n")        
result_dat.write("Count of nodes                     = %d\n" % im)    
result_dat.write("Massflow           inlet           = %.2f [Kg/s]\n" % massflow_in)    
result_dat.write("Massflow          outlet           = %.2f [Kg/s]\n" % massflow_out)       
result_dat.write("Total Temperature  inlet           = %.2f [K]   \n" % ptinlt) 
result_dat.write("Total Temperature  inlet           = %.2f [K]   \n" % ttinlt)
result_dat.write("Total Pressure    outlet           = %.2f [Pa]  \n" % ptout)
result_dat.write("Total Temperature outlet           = %.2f [Pa]  \n" % ttout)
result_dat.write("Total Pressure Ratio    outlet/inlet  = %.6f    \n" % ptv)
result_dat.write("Total Temperature Ratio outlet/inlet  = %.6f    \n" % ttv)  
result_dat.write("  \n")      
result_dat.write("i     vx        p         rho \n")

for i in range(0,im+1):    
    result_dat.write("%d\t%.2f\t%.2f\t%.3f\n" % (i,vx[i],p[i],rho[i]))

# thrust and momentum check
if iarea == 1:    
    fpw  = 0.0
    for i in range(0,im):
        fpw = fpw + 0.5*(p[i+1]+p[1])*(aix[i+1]-aix[i])
            
    fpw = fpw*0.0001
        
    # momentum 
    momentry =  (rho[1] *vx[1]**2 +p[1]) *aix[1] *0.0001
    momexit  = -(rho[im]*vx[im]**2+p[im])*aix[im]*0.0001  

    #momentum check
    fcheck   = (momentry+momexit+fpw)/fpw

    thrust   = 0.0001*rho[im]*vx[im]**2    
     
    result_dat.write("Momentum Entry = %.5f [t]   \n" % momentry) 
    result_dat.write("Momentum Exit  = %.5f [t]   \n" % momexit)    
    result_dat.write("P-Force  Wall  = %.5f [t]   \n" % fpw)     
    result_dat.write("Momentum-Check = %.5f [-]   \n" % fcheck) 
    result_dat.write("Thrust         = %.5f [t]   \n" % thrust)              
   
result_dat.close()


mach_dat = open('nozzle_plot.dat','w')    
for i in range(0,im+1):  
    mach[i] = vx[i]/(cpzucv*p[i]/rho[i])**0.5 
    mach_dat.write("%d\t%.2f\t%.2f\t%.2f\n" % (i,x[i],p[i],mach[i]))        
mach_dat.close()


#movie

fig, ax = plt.subplots()

line, = ax.plot(x,y1,marker="x", linestyle=" ", color="red")

def init():
    ax.set_xlim(0.0, 1.0)
    ax.set_ylim(0.0, 1.80)    
    plt.grid(True)    
    return line,

def animate(frame):
    for i in range(im+1):
        y1[i] = y2[frame][i]
    #print(frame,y1)

    line.set_ydata(y1)  # update the data.
    return line,  

print('nstep = ',nstep)

ani = animation.FuncAnimation(
    fig, animate, interval=anival, frames=range(0,nstep,anistep), init_func=init, blit=True, repeat=False)

# save the animation using Pillow as a gif
writer = animation.PillowWriter(fps=15,
                                 metadata=dict(artist='Me'),
                                 bitrate=1800)
ani.save('Demo1D.gif', writer=writer)


plt.show()


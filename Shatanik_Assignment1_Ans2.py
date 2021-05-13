'''
   Author: Shatanik Bhattacharya
   Department: Department of Astronomy and Astrophyics
'''

import numpy as np
from scipy.integrate import odeint, solve_ivp
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d
import time


#FUNCTION TO BE INTEGRATED
def D_vec_r_Dt_D_vec_v_Dt( t , Config1, we,dummy  ):

    Config0=np.reshape(Config1, (2, 3))
    if (np.linalg.norm(Config0[0])<R or np.linalg.norm(Config0[0])>2.*R):
        return None
    ret    = [[0,0,0],[0,0,0]]
    Omega  = np.array([0,0,we])
    OmegaXR= np.cross(Omega,np.array(Config0[0]))
    
    ret[0] = Config0[1]
    # ACCELERATION IN THE ROTATING FRAME
    ret[1] = np.ndarray.tolist(-G*M/( np.linalg.norm(Config0[0])**3 ) * np.array(Config0[0]) - 2*np.cross(Omega,np.array(Config0[1])) - np.cross(Omega,OmegaXR))
    return np.array(ret).flatten()
    
# THE FUNCTION THAT DOES OUR JOB FOR DIFFERENT VALUES OF THE PARAMETERS USED IN THE CODE
def do_job(R,G,M,w1,  lat0,long0,  alpha,v):
    
    theta0, phi0=90.-lat0, long0
    x0= R * np.sin((theta0)*np.pi/180.) * np.cos(phi0*np.pi/180.)
    y0= R * np.sin((theta0)*np.pi/180.) * np.sin(phi0*np.pi/180.)
    z0= R * np.cos((theta0)*np.pi/180.)
        
    Vx0 = v*((np.cos(alpha))*(np.sin(phi0*np.pi/180.))+(np.sin(alpha))*(np.sin(theta0*np.pi/180.))*(np.cos(phi0*np.pi/180.)))
    
    Vy0 = v*( -(np.cos(alpha))*(np.cos(phi0*np.pi/180.))+(np.sin(alpha))*(np.sin(theta0*np.pi/180.))*(np.sin(phi0*np.pi/180.)))
    
    Vz0 = v*(np.sin(alpha))*(np.sin((90.-theta0)*np.pi/180.))


    r0=[x0,y0,z0]
    v0=[Vx0,Vy0,Vz0]
    
    Config0=[r0,v0]
    Config1=np.array(Config0).flatten()
    Config0=np.reshape(Config1, (2, 3))
    
    #INTEGRATING
    start=time.time()
    tfin=8640.
    ts = np.linspace(0.,tfin,10000)
    Xs = solve_ivp(D_vec_r_Dt_D_vec_v_Dt , [0.,tfin],np.array([x0,y0,z0,Vx0,Vy0,Vz0]), args=(w1,0), method='RK45',    atol=1e-6,max_step=0.1,first_step=0.1,dense_output=False)
    #print(time.time()-start)
    ts=Xs.t
    ys = Xs.y[1]
    
    Longi = (np.arctan(Xs.y[1][-1]/Xs.y[0][-1]))*(180./np.pi)
    Lati = 90.-(np.arccos(Xs.y[2][-1]/R))*(180./np.pi)
    if Longi>=0:
        if Lati>=0:
            print('The missile will land at', Lati ,'N', Longi,'E'," in ",ts[-1]," seconds")
        else:
            print('The missile will land at', -Lati ,'S', Longi,'E'," in ",ts[-1]," seconds")
    else:
        if Lati>=0:
            print('The missile will land at', Lati ,'N', -Longi,'W'," in ",ts[-1]," seconds")
        else:
            print('The missile will land at', -Lati ,'S', -Longi,'W'," in ",ts[-1]," seconds")


    theta, phi = np.linspace(0, 2*np.pi, 50), np.linspace(0, np.pi, 50)
    THETA, PHI = np.meshgrid(theta, phi)
    Ri = 1.0
    Xi = Ri*np.sin(PHI)*np.cos(THETA)
    Yi = Ri*np.sin(PHI)*np.sin(THETA)
    Zi = Ri*np.cos(PHI)
    fig = plt.figure()
    ax = fig.add_subplot(1,1,1, projection='3d')
    ax.plot(Xs.y[0]/R, Xs.y[1]/R, Xs.y[2]/R, 'blue')
    ax.plot_wireframe(Xi, Yi, Zi,lw = 0.1 , color = 'red',alpha =0.5)
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    plt.show()
    
    Want=input("Want to plot graphs of the different parameters? (not necessary) Then press Y or y: ")
    if (Want=='Y' or Want=='y'):
        plt.plot(Xs.t,Xs.y[0],'b-')
        plt.title(r"x")
        plt.show()
    
        plt.plot(Xs.t,Xs.y[3],'b-')
        plt.title(r"$v_{x}$")
        plt.show()
    
        plt.plot(Xs.t,Xs.y[1],'b-')
        plt.title(r"y")
        plt.show()
    
        plt.plot(Xs.t,Xs.y[4],'b-')
        plt.title(r"$v_{y}$")
        plt.show()
    
        plt.plot(Xs.t,Xs.y[2],'b-')
        plt.title(r"z")
        plt.show()
    
        plt.plot(Xs.t,Xs.y[5],'b-')
        plt.title(r"$v_{z}$")
        plt.show()
    
    
# INITIALIZATION OF VARIABLES
R=6378000
G=6.67430*10**(-11)
M=5.97219*10**24
w0=2*np.pi/86400


print("\n\n When rocket is released at 45 degrees from 8.5 degrees N 77 degrees E with 8000 kmph; Earth is rotating with angular frequency ",2*np.pi/86400," Hz :::==> ")
do_job(R,G,M,w0,  8.5,77.,  45.*(np.pi/180.),8000.*5./18.)


print("\n\n When rocket is released at 45 degrees from 30 degrees N 77 degrees E with 8000 kmph; Earth is rotating with angular frequency ",2*np.pi/86400," Hz :::==> ")
do_job(R,G,M,w0,  30.,77.,  45.*(np.pi/180.),8000.*5./18.)


print("\n\n When rocket is released at 30 degrees from 8.5 degrees N 77 degrees E with 8000 kmph; Earth is rotating with angular frequency ",2*np.pi/86400," Hz :::==> ")
do_job(R,G,M,w0,  8.5,77.,  30.*(np.pi/180.),8000.*5./18.)


print("\n\n When rocket is released at 60 degrees from 8.5 degrees N 77 degrees E with 8000 kmph; Earth is rotating with angular frequency ",2*np.pi/86400," Hz :::==> ")
do_job(R,G,M,w0,  8.5,77.,  60.*(np.pi/180.),8000.*5./18.)

print("\n\n When rocket is released at 45 degrees from 8.5 degrees N 77 degrees E with 8000 kmph; Earth is rotating with angular frequency ",2*np.pi/86400," Hz :::==> ")
do_job(R,G,M,w0,  8.5,77.,  45.*(np.pi/180.),8000.*5./18.)


print("\n\n When rocket is released at 45 degrees from 8.5 degrees N 77 degrees E with 8000 kmph; Earth is not rotating :::==> ")
do_job(R,G,M,0.,  8.5,77.,  45.*(np.pi/180.),8000.*5./18.)


print("\n\n When rocket is released at 45 degrees from 8.5 degrees N 77 degrees E with 8000 kmph; Earth is rotating with angular frequency ",4*np.pi/86400," Hz (double of normal) :::==> ")
do_job(R,G,M,2*w0,  8.5,77.,  45.*(np.pi/180.),8000.*5./18.)


AbsDeltaPhi=list()
for wi in np.linspace(0.,4*w0,20):

    
    theta0, phi0=90.-0., 0.
    x0= R * np.sin((theta0)*np.pi/180.) * np.cos(phi0*np.pi/180.)
    y0= R * np.sin((theta0)*np.pi/180.) * np.sin(phi0*np.pi/180.)
    z0= R * np.cos((theta0)*np.pi/180.)
    
    alpha=45.*np.pi/180.
    v=8000.*5./18.
     
    Vx0 = v*((np.cos(alpha))*(np.sin(phi0*np.pi/180.))+(np.sin(alpha))*(np.sin(theta0*np.pi/180.))*(np.cos(phi0*np.pi/180.)))
    
    Vy0 = v*( -(np.cos(alpha))*(np.cos(phi0*np.pi/180.))+(np.sin(alpha))*(np.sin(theta0*np.pi/180.))*(np.sin(phi0*np.pi/180.)))
    
    Vz0 = v*(np.sin(alpha))*(np.sin((90.-theta0)*np.pi/180.))


    r0=[x0,y0,z0]
    v0=[Vx0,Vy0,Vz0]
    
    Config0=[r0,v0]
    Config1=np.array(Config0).flatten()
    
    tfin=86400.
    ts = np.linspace(0.,tfin,1000)
    Xs = solve_ivp(D_vec_r_Dt_D_vec_v_Dt , [0.,tfin],Config1, args=(wi,0), method='RK45',    atol=1e-6,max_step=0.1,first_step=0.1,dense_output=False)
    
    AbsDelPhi = abs((np.arctan(Xs.y[1][-1]/Xs.y[0][-1]))*(180./np.pi))
    AbsDeltaPhi.append(AbsDelPhi)
      

plt.plot(np.linspace(0.,4.,20),AbsDeltaPhi,label=r"Variation of the magnitude of longitude of impact w.r.t. varying $\omega$ of earth when released from 0N 0E towards west with $\alpha=45^{o}$ and $v=8000$ kmph")
plt.xlabel(r"$\dfrac{\omega}{\omega_{0}}$")
plt.ylabel(r"$|\Delta\Phi|$")
plt.legend()
plt.show()    
 
    

# CUSTOM INPUT
print("\n\n Custom input variables (for fun): ")

lat0,long0=float(input("Latitude: ")),float(input("Longitude: "))
alpha = float(input("Enter the angle at which the missile is lauched with respect to west direction: "))*(np.pi/180.)
v = float(input("Enter the speed (in kmph) at which the missile is launched: "))*5./18.

do_job(R,G,M,w0,  lat0,long0,  alpha,v)

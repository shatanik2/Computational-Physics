'''
   Author: Shatanik Bhattacharya
   Department: Department of Astronomy and Astrophyics
'''

import numpy as np
import scipy
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp


mass  = 0.5
R     = 0.1
I     = (mass * R**2)/2
w0    = np.pi

m     = 0.02
I_add = I + m * R**2
wdl   = w0 /  2.
wdh   = 2  * w0

A=0.002
wd=wdl
b=1e-2

g=9.8


def dU_dtI(t, U, b,I,g,w0,A,wd,check_added_mass):
    dU=np.zeros_like(U)
    dU[0]=U[1]
    dU[1]=(-b*U[1]/I) - (w0**2 * U[0]) + A*np.sin(wd*t)/I +  check_added_mass * m * g * R * np.sin(U[0])/I
    return dU

       
U0  = [45.*np.pi/180., 0.]
ts  = np.linspace(0, 10, 10000)


'''
#WEAKLY DAMPED CASE (UNCOMMENT THIS SECTION TO VIEW THE PHASE SPACE TRAJECTORY)
sol  =    solve_ivp(dU_dtI,[0.,10.],  np.array(U0), args=(0.001,I,g,w0,0,wd,0) , method='RK45', atol=1e-6,t_eval=ts)

plt.plot(sol.y[0] , I*sol.y[1],'b-',label=r'b='+str(0.001)+", mass added="+str(0)+", $\omega_{0}=\pi\:$ Hz, Forcing="+str(0)+", $\omega_{d}=\dfrac{\omega_{0}}{2}$")
plt.xlabel(r"Angle $\phi$ in radians")
plt.ylabel(r"Angular momentum $I_{tot}\:\dot{\phi}$ in $kg\:m^{2}\:s^{-1}$")
plt.title("Phase space trajectory")
plt.legend()
plt.show()

plt.xlabel(r"Time $t$ in seconds")
plt.ylabel(r"Angle $\phi$ in degrees")
plt.title("Pohl's oscillator")
plt.plot(sol.t,sol.y[0]*180./np.pi,'b-')

plt.show()
'''

'''
#UNDAMPED CASE (UNCOMMENT THIS SECTION TO VIEW THE PHASE SPACE TRAJECTORY)
sol  =    solve_ivp(dU_dtI,[0.,10.],  np.array(U0), args=(0.0,I,g,w0,0,wd,0) , method='RK45', atol=1e-6,t_eval=ts)

plt.plot(sol.y[0] , I*sol.y[1],'b-',label=r'b='+str(0.0)+", mass added="+str(0)+", $\omega_{0}=\pi\:$ Hz, Forcing="+str(0)+", $\omega_{d}=\dfrac{\omega_{0}}{2}$")
plt.xlabel(r"Angle $\phi$ in radians")
plt.ylabel(r"Angular momentum $I_{tot}\:\dot{\phi}$ in $kg\:m^{2}\:s^{-1}$")
plt.title("Phase space trajectory")
plt.legend()
plt.show()

plt.xlabel(r"Time $t$ in seconds")
plt.ylabel(r"Angle $\phi$ in degrees")
plt.title("Pohl's oscillator")
plt.plot(sol.t,sol.y[0]*180./np.pi,'b-')

plt.show()
'''


'''
#HIGHLY DAMPED CASE (UNCOMMENT THIS SECTION TO VIEW THE PHASE SPACE TRAJECTORY)
sol  =    solve_ivp(dU_dtI,[0.,10.],  np.array(U0), args=(0.05,I,g,w0,0,wd,0) , method='RK45', atol=1e-6,t_eval=ts)

plt.plot(sol.y[0] , I*sol.y[1],'b-',label=r'b='+str(0.05)+", mass added="+str(0)+", $\omega_{0}=\pi\:$ Hz, Forcing="+str(0)+", $\omega_{d}=\dfrac{\omega_{0}}{2}$")
plt.xlabel(r"Angle $\phi$ in radians")
plt.ylabel(r"Angular momentum $I_{tot}\:\dot{\phi}$ in $kg\:m^{2}\:s^{-1}$")
plt.title("Phase space trajectory")
plt.legend()
plt.show()

plt.xlabel(r"Time $t$ in seconds")
plt.ylabel(r"Angle $\phi$ in degrees")
plt.title("Pohl's oscillator")
plt.plot(sol.t,sol.y[0]*180./np.pi,'b-')

plt.show()
'''


'''
#WEAKLY DAMPED CASE with ADDED MASS (UNCOMMENT THIS SECTION TO VIEW THE PHASE SPACE TRAJECTORY)
sol  =    solve_ivp(dU_dtI,[0.,10.],  np.array(U0), args=(0.001,I_add,g,w0,0,wd,1) , method='RK45', atol=1e-6,t_eval=ts)

plt.plot(sol.y[0] , I*sol.y[1],'g-',label=r'b='+str(0.001)+", mass added="+str(0.02)+"kg, $\omega_{0}=\pi\:$ Hz, Forcing="+str(0)+", $\omega_{d}=\dfrac{\omega_{0}}{2}$")
plt.xlabel(r"Angle $\phi$ in radians")
plt.ylabel(r"Angular momentum $I_{tot}\:\dot{\phi}$ in $kg\:m^{2}\:s^{-1}$")
plt.title("Phase space trajectory")
plt.legend()
plt.show()

plt.xlabel(r"Time $t$ in seconds")
plt.ylabel(r"Angle $\phi$ in degrees")
plt.title("Pohl's oscillator")
plt.plot(sol.t,sol.y[0]*180./np.pi,'g-')

plt.show()
'''



'''
#WEAKLY DAMPED CASE WITH DRIVE OF FREQUENCY W0/2 (UNCOMMENT THIS SECTION TO VIEW THE PHASE SPACE TRAJECTORY)
sol  =    solve_ivp(dU_dtI,[0.,10.],  np.array(U0), args=(0.001,I,g,w0,0.002,w0/2.,0) , method='RK45', atol=1e-6,t_eval=ts)

plt.plot(sol.y[0] , I*sol.y[1],'b-',label=r'b='+str(0.001)+", mass added="+str(0)+", $\omega_{0}=\pi\:$ Hz, Forcing="+str(0.002)+", $\omega_{d}=\dfrac{\omega_{0}}{2}$")
plt.xlabel(r"Angle $\phi$ in radians")
plt.ylabel(r"Angular momentum $I_{tot}\:\dot{\phi}$ in $kg\:m^{2}\:s^{-1}$")
plt.title("Phase space trajectory")
plt.legend()
plt.show()

plt.xlabel(r"Time $t$ in seconds")
plt.ylabel(r"Angle $\phi$ in degrees")
plt.title("Pohl's oscillator")
plt.plot(sol.t,sol.y[0]*180./np.pi,'b-')

plt.show()
'''

'''
#WEAKLY DAMPED CASE WITH DRIVE OF FREQUENCY W0 (UNCOMMENT THIS SECTION TO VIEW THE PHASE SPACE TRAJECTORY)
sol  =    solve_ivp(dU_dtI,[0.,10.],  np.array(U0), args=(0.001,I,g,w0,0.002,w0/1.,0) , method='RK45', atol=1e-6,t_eval=ts)

plt.plot(sol.y[0] , I*sol.y[1],'b-',label=r'b='+str(0.001)+", mass added="+str(0)+", $\omega_{0}=\pi\:$ Hz, Forcing="+str(0.002)+", $\omega_{d}=\omega_{0}$")
plt.xlabel(r"Angle $\phi$ in radians")
plt.ylabel(r"Angular momentum $I_{tot}\:\dot{\phi}$ in $kg\:m^{2}\:s^{-1}$")
plt.title("Phase space trajectory")
plt.legend()
plt.show()

plt.xlabel(r"Time $t$ in seconds")
plt.ylabel(r"Angle $\phi$ in degrees")
plt.title("Pohl's oscillator")
plt.plot(sol.t,sol.y[0]*180./np.pi,'b-')

plt.show()
'''

'''
#WEAKLY DAMPED CASE WITH DRIVE OF FREQUENCY 2*W0 (UNCOMMENT THIS SECTION TO VIEW THE PHASE SPACE TRAJECTORY)
sol  =    solve_ivp(dU_dtI,[0.,10.],  np.array(U0), args=(0.001,I,g,w0,0.002,2.*w0,0) , method='RK45', atol=1e-6,t_eval=ts)

plt.plot(sol.y[0] , I*sol.y[1],'b-',label=r'b='+str(0.001)+", mass added="+str(0)+", $\omega_{0}=\pi\:$ Hz, Forcing="+str(0.002)+", $\omega_{d}=2\:\omega_{0}$")
plt.xlabel(r"Angle $\phi$ in radians")
plt.ylabel(r"Angular momentum $I_{tot}\:\dot{\phi}$ in $kg\:m^{2}\:s^{-1}$")
plt.title("Phase space trajectory")
plt.legend()
plt.show()

plt.xlabel(r"Time $t$ in seconds")
plt.ylabel(r"Angle $\phi$ in degrees")
plt.title("Pohl's oscillator")
plt.plot(sol.t,sol.y[0]*180./np.pi,'b-')

plt.show()
'''

'''
#WEAKLY DAMPED CASE WITH DRIVE OF FREQUENCY W0/2 with ADDED MASS (UNCOMMENT THIS SECTION TO VIEW THE PHASE SPACE TRAJECTORY)
sol  =    solve_ivp(dU_dtI,[0.,10.],  np.array(U0), args=(0.001,I_add,g,w0,0.002,w0/2.,1) , method='RK45', atol=1e-6,t_eval=ts)

plt.plot(sol.y[0] , I*sol.y[1],'g-',label=r'b='+str(0.001)+", mass added="+str(0.02)+"kg, $\omega_{0}=\pi\:$ Hz, Forcing="+str(0.002)+", $\omega_{d}=\dfrac{\omega_{0}}{2}$")
plt.xlabel(r"Angle $\phi$ in radians")
plt.ylabel(r"Angular momentum $I_{tot}\:\dot{\phi}$ in $kg\:m^{2}\:s^{-1}$")
plt.title("Phase space trajectory")
plt.legend()
plt.show()

plt.xlabel(r"Time $t$ in seconds")
plt.ylabel(r"Angle $\phi$ in degrees")
plt.title("Pohl's oscillator")
plt.plot(sol.t,sol.y[0]*180./np.pi,'g-')

plt.show()
'''

'''
#WEAKLY DAMPED CASE WITH DRIVE OF FREQUENCY W0 with ADDED MASS (UNCOMMENT THIS SECTION TO VIEW THE PHASE SPACE TRAJECTORY)
sol  =    solve_ivp(dU_dtI,[0.,10.],  np.array(U0), args=(0.001,I_add,g,w0,0.002,w0,1) , method='RK45', atol=1e-6,t_eval=ts)

plt.plot(sol.y[0] , I*sol.y[1],'g-',label=r'b='+str(0.001)+", mass added="+str(0.02)+"kg, $\omega_{0}=\pi\:$ Hz, Forcing="+str(0.002)+", $\omega_{d}=\omega_{0}$")
plt.xlabel(r"Angle $\phi$ in radians")
plt.ylabel(r"Angular momentum $I_{tot}\:\dot{\phi}$ in $kg\:m^{2}\:s^{-1}$")
plt.title("Phase space trajectory")
plt.legend()
plt.show()

plt.xlabel(r"Time $t$ in seconds")
plt.ylabel(r"Angle $\phi$ in degrees")
plt.title("Pohl's oscillator")
plt.plot(sol.t,sol.y[0]*180./np.pi,'g-')

plt.show()
'''


#WEAKLY DAMPED CASE WITH DRIVE OF FREQUENCY 2*W0 with ADDED MASS (UNCOMMENT THIS SECTION TO VIEW THE PHASE SPACE TRAJECTORY)
sol  =    solve_ivp(dU_dtI,[0.,10.],  np.array(U0), args=(0.001,I_add,g,w0,0.002,2.*w0,1) , method='RK45', atol=1e-6,t_eval=ts)

plt.plot(sol.y[0] , I*sol.y[1],'g-',label=r'b='+str(0.001)+", mass added="+str(0.02)+"kg, $\omega_{0}=\pi\:$ Hz, Forcing="+str(0.002)+", $\omega_{d}=2\:\omega_{0}$")
plt.xlabel(r"Angle $\phi$ in radians")
plt.ylabel(r"Angular momentum $I_{tot}\:\dot{\phi}$ in $kg\:m^{2}\:s^{-1}$")
plt.title("Phase space trajectory")
plt.legend()
plt.show()

plt.xlabel(r"Time $t$ in seconds")
plt.ylabel(r"Angle $\phi$ in degrees")
plt.title("Pohl's oscillator")
plt.plot(sol.t,sol.y[0]*180./np.pi,'g-')

plt.show()











'''
   Author: Shatanik Bhattacharya
   Department: Department of Astronomy and Astrophyics
'''

import numpy as np
from scipy import constants as const
from scipy import sparse as sparse
from scipy.sparse.linalg import eigsh
 
hbar,e,c = const.hbar,const.e,const.c
alpha=0.471*hbar*c
sigma=0.2025*10**18*(e**2/(hbar*c))

m=1.32*10**9
mu=m*e/(2*c**2)

# Potential matrix
def V(r):
    return sparse.diags((- alpha/r + sigma*r))

# Angular degree matrix
def l_term(r,l):
    return sparse.diags((l * (l + 1) / r**2))

# Laplacian's Matrix Form
def Laplacian_matrix_rep(r):
    h = r[1] - r[0]     
    diag = -490.0 / (180*h**2) * np.ones(N)     
    off_diag1  =  270.0 / (180*h**2) * np.ones(N - 1)
    off_diag2  =  -27.0 / (180*h**2) * np.ones(N - 2)
    off_diag3  =  2.0 / (180*h**2) * np.ones(N - 3)
    Laplacian_matrix_rep = sparse.diags([diag, off_diag3, off_diag3, off_diag2, off_diag2, off_diag1, off_diag1], (0, -3, 3, -2, 2, -1, 1))
    return Laplacian_matrix_rep
   
# Hamiltonian whose eigenvalues and eigenvectors are to be found  
def Hamil(r,l):     
    return (-hbar**2 / (2.0 * mu) * (Laplacian_matrix_rep(r) - l_term(r,l)) + V(r))

print("=======================================================================")
print("    Value of l \t\tEnergy (in eV) \t\t  Mass (in GeV/c^2) ")
print("=======================================================================")

N = 2000
l = 0
r = np.linspace(2e-14, 0.0, N, endpoint=False)
H = Hamil(r,l)

spectrum=list()

#Solving the eigenvalue problem and sorting eigenvalue
num_eigv = 2
vals, vecs = eigsh(H, k=num_eigv, which='SM')

vecs = np.array(vecs.T)
vals = np.sort(vals)

energy0 = [(vals[i]/e) for i in range(num_eigv)]


mass0 = [((vals[i]/e)+2*m)*10**(-9) for i in range(num_eigv)]
for i in range(len(energy0)):
   print("\t",l,"\tE",i,"= ",energy0[i],"\t",mass0[i])
#print("\t",0,"\tE1=",energy0[1],"\t",mass0[1])


N = 2000
l = 1
r = np.linspace(2e-14, 0.0, N, endpoint=False)
H = Hamil(r,l)

#Solving the eigenvalue problem and sorting eigenvalue
num_eigv = 1
vals, vecs = eigsh(H, k=num_eigv, which='SM')
 
vecs = np.array(vecs.T)
vals = np.sort(vals)

energy1 = [(vals[i]/e) for i in range(num_eigv)]

mass1 = [((vals[i]/e)+2*m)*10**(-9) for i in range(num_eigv)]
for i in range(len(energy1)):
   print("\t",l,"\tE",len(energy0)+i,"= ",energy1[i],"\t",mass1[i])


energy=energy0+energy1
mass=mass0+mass1
print("\n\nSpectrum energy values (in MeV) ==>")
for i in range(len(energy)):
    for j in range(i+1,len(energy)):
         spectrum.append(abs(energy[i]-energy[j])/10**6)

print(sorted(spectrum))


'''
print("\n")
print("PDG values for reference===>")
print("Mass of lightest S wave charmonium= " ,3.096," GeV/c^2")
print("Mass of lightest P wave charmonium= " ,3.525," GeV/c^2")
'''


import numpy as np

C = 2.9979e10     # speed of light in cgs
M_sol = 1.98892e30  # kg

## divide the quantity in cgs over these to get in ETK units

cgs_to_ETK_Mass   = M_sol*1e3            # code unit mass to g 
cgs_to_ETK_Length = 1.477*100*1000       # code unit length to cm 
cgs_to_ETK_Time   = cgs_to_ETK_Length/C  # code unit time to s 

cgs_to_ETK_density  = cgs_to_ETK_Mass/(cgs_to_ETK_Length)**3  # code unit desnsity to g/cm^3
cgs_to_ETK_pressure = cgs_to_ETK_Mass/(cgs_to_ETK_Length*cgs_to_ETK_Time**2)  # code unit pressure to g/(cm s^2)

cgs_to_ETK_K       = lambda gamma: cgs_to_ETK_pressure/np.power(cgs_to_ETK_density,gamma) # code unit Kappa to g/(cm s^2)



rho_min = 5e6  # min value of rest mass density
rho_max = 3e15  # max value of rest mass density

print("Min and max densities in ETK units:")
print(rho_min/cgs_to_ETK_density)
print(rho_max/cgs_to_ETK_density)


# Create SLy crust

rho1 = pow(10,14.7)/cgs_to_ETK_density 
rho2 = pow(10,15.0)/cgs_to_ETK_density

print(75*"*")
print("Dividing densities: ")
print(75*"*")
print(f"{rho1:.20e}\n",f"{rho2:.20e}\n")

# Set SLy crust

rhoL_1 = 2.62789e12/cgs_to_ETK_density
rhoL_2 = 3.78358e11/cgs_to_ETK_density
rhoL_3 = 2.44034e7/cgs_to_ETK_density
rhoL_4 = 0.0

GammaL_1 = 1.35692
GammaL_2 = 0.62223
GammaL_3 = 1.28733
GammaL_4 = 1.58425

KL_1 = 3.99874e-8 * pow(cgs_to_ETK_density, GammaL_1-1)  # notice a missing c^2 in Ki values in Table II of Read et al. 2009
KL_2 = 5.32697e+1 * pow(cgs_to_ETK_density, GammaL_2-1) 
KL_3 = 1.06186e-6 * pow(cgs_to_ETK_density, GammaL_3-1)  
KL_4 = 6.80110e-9 * pow(cgs_to_ETK_density, GammaL_4-1)  

print(75*"*")
print("Low density region dividing densities: ")
print(75*"*")
print(f"{rhoL_1:.20e}\n",f"{rhoL_2:.20e}\n",f"{rhoL_3:.20e}\n",f"{rhoL_4:.20e}\n")
print(75*"*")
print("Low density region Ks:")
print(75*"*")
print(f"{KL_1:.20e}\n", f"{KL_2:.20e}\n", f"{KL_3:.20e}\n", f"{KL_4:.20e}\n")



## Choose EOS 

EOS = "APR4" 
perturb_collapse = True
delta_K = 1.01
set_crust = True 

if(EOS=="SLy"):
    p1 = pow(10.0,34.384)/cgs_to_ETK_density/C**2
    Gamma1 = 3.005
    Gamma2 = 2.988
    Gamma3 = 2.851

if(EOS=="H4"):
    p1 = pow(10.0,34.669)/cgs_to_ETK_density/C**2
    Gamma1 = 2.909
    Gamma2 = 2.246
    Gamma3 = 2.144

if(EOS=="APR4"):
    p1 = pow(10.0,34.269)/cgs_to_ETK_density/C**2
    Gamma1 = 2.830
    Gamma2 = 3.445
    Gamma3 = 3.348



print(45*"*",EOS,45*"*","\n p1",p1)
print(75*"*")
print("\n Gammas:")
print(75*"*")
print(f"{Gamma1}\n{Gamma2}\n{Gamma3}\n")


K1 = p1 / pow(rho1,Gamma1)
K2 = K1 * pow( rho1, Gamma1-Gamma2)
K3 = K2 * pow( rho2, Gamma2-Gamma3)

if perturb_collapse:
    K1 *= delta_K
    K2 *= delta_K
    K3 *= delta_K

print(75*"*")
print("\n K_i:")
print(75*"*")
print(f"K1\n{K2}\n{K3}\n")


# Low density region (crust)

# calculate eps and alpha

epsL_4 = 0.0
alphaL_4 = 0.0
epsL_3 = (1+alphaL_4)*rhoL_3 + KL_4/(GammaL_4 - 1)*pow(rhoL_3, GammaL_4)
alphaL_3 = epsL_3/rhoL_3 - 1 - KL_3/(GammaL_3 - 1)*pow(rhoL_3, GammaL_3 -1)
epsL_2 = (1+alphaL_3)*rhoL_2 + KL_3/(GammaL_3 - 1)*pow(rhoL_2, GammaL_3)
alphaL_2 = epsL_2/rhoL_2 - 1 - KL_2/(GammaL_2 - 1)*pow(rhoL_2, GammaL_2 -1)
epsL_1 = (1+alphaL_2)*rhoL_1 + KL_2/(GammaL_2 - 1)*pow(rhoL_1, GammaL_2)
alphaL_1 = epsL_1/rhoL_1 - 1 - KL_1/(GammaL_1 - 1)*pow(rhoL_1, GammaL_1 -1)


# calculate pressure in crust
pL_3 = KL_3*pow(rhoL_3,GammaL_3)
pL_2 = KL_2*pow(rhoL_2,GammaL_2)
pL_1 = KL_1*pow(rhoL_1,GammaL_1)

print(75*"*")
print("Eps for the low density region (crust): ")
print(75*"*")
print(f"{epsL_4}\n{epsL_3}\n{epsL_2}\n{epsL_1}")


# Interface

if set_crust:
    rho0 = pow(KL_1/K1,1.0/(Gamma1-GammaL_1))    # Through continuity of the pressure at the dividing density
    eps0 = (1.0+alphaL_1)*rho0 + KL_1/(GammaL_1-1.0)*pow(rho0,GammaL_1)
    p0 = KL_1*pow(rho0,GammaL_1)
else: 
    rho0 = 5e6/cgs_to_ETK_density 
    p0 = K1 * pow(rho0,Gamma1)
    eps0 = p0/(Gamma1-1)+rho0
    h0 = (eps0 + p0)/rho0
print(75*"*")
print("Dividing density between low and high density regions:")
print(75*"*")
print(f"rho: {rho0:.20e}\n",f"eps0: {eps0}\n",f"p0:  {p0}\n",f"e {eps0/rho0-1}\n")

# High density 

alpha1 = eps0/rho0 - 1 - K1/(Gamma1 - 1)*pow(rho0, Gamma1 -1)
eps_int_1 = alpha1 + K1/(Gamma1-1)*pow(rho1,Gamma1-1)
eps1 = (1+alpha1)*rho1 + K1/(Gamma1 - 1)*pow(rho1, Gamma1)

alpha2 = eps1/rho1 - 1 - K2/(Gamma2 - 1)*pow(rho1, Gamma2 -1)
eps_int_2 = alpha2 + K2/(Gamma2-1)*pow(rho2,Gamma2-1)
eps2 = (1+alpha2)*rho2 + K2/(Gamma2 - 1)*pow(rho2, Gamma2)

alpha3 = eps2/rho2 - 1 - K3/(Gamma3 - 1)*pow(rho2, Gamma3 -1)

p2 = K2*pow(rho2,Gamma2)

print(75*"*")
print("Proper energy densitty for the high density region: ")
print(75*"*")
print(f"{eps1}\n{eps2}\n")

print(75*"*")
print("Specific internal energy for the high density region: ")
print(75*"*")
print(f"{eps_int_1}\n{eps_int_2}\n")


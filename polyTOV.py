import numpy as np

def EOS_p2erc(p, K=100., Gamma=2.):
  """
  Equation of state (EOS)
  Given pressure, return energy density, rest-mass density and sound speed
  """      
  ene = (p/K)**(1./Gamma) + p/(Gamma-1.)
  rho = (p/K)**(1./Gamma)                       
  cs2 = K*Gamma*(Gamma-1)/(Gamma-1 + K*Gamma*rho**(Gamma-1))*rho**(Gamma-1)
  return ene, rho, cs2

def EOS_r2pe(rho, K=100., Gamma=2.):
  """
  Equation of state (EOS)
  Given rest-mass density, return energy density and pressure
  Polytropic EOS: P = k rho^Gamma
  """
  p = K*rho**Gamma
  e = rho + p/(Gamma-1.);
  return p, e


def make_TOV(K, Gamma):
    """
    Returns a TOV function with specific EOS
    """
    def TOV(t, y):
        """
        Tolmann-Oppenheimer-Volkhoff equations
        d/dt y(t) = R.H.S. 
        """
        r = t
        # Unpack current state array
        m_b = y[0] # baryon mass
        m   = y[1] # mass of a sphere of radius r
        p   = y[2] # pressure
        nu  = y[3] # potential 

        # Call the EOS
        ene, rho, _ = EOS_p2erc(p,K1,Gamma1) 

        # Set the RHS
        dy = np.empty_like(y)
        dy[0] = 4*np.pi*rho*r**2/np.sqrt(1-2*m/r)          # eq for baryon mass
        dy[1] = 4*np.pi*ene*r**2                           # eq for grav mass                          
        dy[2] = -(ene+p)*(m + 4*np.pi*r**3*p)/(r*(r-2*m))  # eq for pressure 
        dy[3] = (m + 4*np.pi*r**3*p)/(r*(r-2*m))           # eq for potential 
        return dy
    return TOV

found_radius = lambda t,y ,pfloor=1e-10: ((y[2]-pfloor)<=0.)
# finds pressure=0. integration stops when this function returns True
# should add another factor here to make it more robust, maybe the density or enthalpy would work, pressure changes too fast

def inital_TOV(rho_0,K, Gamma,verbose=True):
    """
    initialize TOV
    """
    p0,e0 = EOS_r2pe(rho_0,K,Gamma)
    if verbose: 
        print("Cold energy density",e0)
        print("Cold Pressure",p0)
        print("rest-mass density",rho_0)
    m_b0  =  4./3.*np.pi*rho_0*rmin**3
    m0    = 4./3.*np.pi*e0*rmin**3  
    nu0   = 0 # this is probably not zero 
    sol0 = [m_b0, m0, p0,nu0]
    
    return sol0 

def solve_ode_rk4(t_array, y0, dydt_fun, stop_when=None, verbose=False):
    """
     RK4 integrator
    Returns the solution at the last step or at the surface.
    """
    N = len(t_array)
    dt = np.diff(t_array)[0]
    y = y0
    for i in range(N):
        t = t_array[i]
        y_prev = np.copy(y)
        k1 = dydt_fun(t, y)
        k2 = dydt_fun(t + dt/2, y + dt/2 * k1)
        k3 = dydt_fun(t + dt/2, y + dt/2 * k2)
        k4 = dydt_fun(t + dt, y + dt * k3)
        y += dt/6 * (k1 + 2*k2 + 2*k3 + k4)
        if verbose:
            print(t, y)
        if stop_when:
            if bool(stop_when(t, y)):
                print(50*"*"+"\nSurface found.\n"+50*"*")
                return t_array[i-1], y_prev
    if stop_when:
        print("No surface reached")
    return t_array[-1], y

# get grid 
rmin, rmax = 1e-6, 20. # technically should start in zero, improve this
N = 10000 # number of points between rmin and rmax
rspan = np.linspace(rmin,rmax,N)

# set up eos  and initialize model 
K1 = 100
Gamma1 = 2
rho_0 =  1.28e-3   # Central (maximal) rest-mass density


# calculate initial conditions
sol0 = inital_TOV(rho_0,K1)  # includes m_b0, m0, p0,nu0
TOV = make_TOV(K1, Gamma1)

# Solve
t, sol =  solve_ode_rk4(rspan, sol0, TOV, stop_event=found_radius, verbose=False)
# Get mass and radius
R = t*1.477  # km 
M_bar = sol[0] 
M = sol[1] # Msun
pmin = sol[2]



print("************************************************")
print(f"Min pressure: {pmin}")  
print(f"Radius:       {R} km ")        
print(f"Mass:         {M} $M_\odot$")   
print(f"Baryon Mass:  {M_bar} $M_\odot$ ")
print("************************************************")





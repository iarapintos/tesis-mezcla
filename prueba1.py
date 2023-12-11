import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt

#mi ecuacion es dM/dt = dot M_b + dot M_ac

# function that returns dy/dt
def model(M,t):
    
    pi = 3.14
    A = 5/2
    G = 6.67
    c = 300
    rho = 7
    P = 9
    
    dot_M_ac = 4*pi*A*G**2*c**(-5)*M**2*(c**2*rho-P)
    
    dMdt = dot_M_ac
    return dMdt

# initial condition
y0 = 8

# time points
t = np.linspace(-100, 100, 1000)

# solve ODE
sol = solve_ivp(model, [t[0], t[-1]], [y0], method='RK45', t_eval=t)

t = sol.t
y = sol.y[0]

# plot results
plt.plot(t,y)
plt.xlabel('time')
plt.ylabel('y(t)')
plt.xscale('log')
#plt.yscale('log')
plt.show()

#%% 

import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt

#mi ecuacion es dM/dt = dot M_b + dot M_ac

# function that returns dy/dt
def model(M,t):
    
    pi = 3.14
    A = 5/2
    G = 6.67
    c = 300
    rho = 7
    P = 9
    
    dot_M_ac = 4*pi*A*G**2*c**(-5)*M**2*(c**2*rho-P)
    
    dMdt = dot_M_ac
    return dMdt

# initial condition
y0 = 8

# time points
t = np.linspace(-100, 100, 1000)

# solve ODE
sol = solve_ivp(model, [t[0], t[-1]], [y0], method='RK45', t_eval=t)

t = sol.t
y = sol.y[0]

# plot results
plt.plot(t,y)
plt.xlabel('time')
plt.ylabel('y(t)')
plt.xscale('log')
#plt.yscale('log')
plt.show()

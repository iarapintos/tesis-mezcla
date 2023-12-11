#ojo!!! quizas esta todo muy mal derivado!!!!!!

import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt

#mi ecuacion es dM/dt = dot M_b + dot M_ac

# function that returns dy/dt
def model(M,t):
    
    m_0 = 1
    a_b = 1
    w = 1/3
    T_b = 1
    num = a_b*2*(1-w)* t * (1+(t/T_b)*2)*(-2/3-w/3)
    den = 3*T_b**2
    
    dot_M_b = num/den 
    
    dMdt = dot_M_b 
    return dMdt

# initial condition
y0 = 20
# time points
t = np.linspace(-1000, 1000, 1000)


# solve ODE
sol = solve_ivp(model, [t[0], t[-1]], [y0], method='RK45', t_eval=t)

y = sol.y

# plot results
plt.plot(t,y)
plt.xlabel('time')
plt.ylabel('y(t)')

#plt.yscale('log')
#plt.xscale('log')
plt.show()

#%% 


import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt

t = np.linspace(-1000, 1000, 1000)

a_b = 1
T_b = 1
w = 1/3
m_0 = 1

a_t = a_b *(1+(t/T_b)**2)**((1/3)*(1-w))

M = m_0*a_t 




# plot results
plt.plot(t,M)
plt.xlabel('tiempo')
plt.ylabel('masa')

#plt.yscale('log')
#plt.xscale('log')
plt.show()

#%% 

import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt

#mi ecuacion es dM/dt = dot M_b + dot M_ac

# function that returns dy/dt
def model(t, M):
    
    m_0 = 1
    a_b = 1
    w = 1/3
    T_b = 1
    
    num = a_b*2*(1-w)* t * (1+(t/T_b)**2)**(-2/3-w/3)
    den = 3*T_b**2
    
    dot_M_b = num/den 
    
    
    dMdt = dot_M_b 
    return dMdt

# initial condition
y0 = 1

# time points
tiempo = np.linspace(0, 1000, 1000)

print(t[0])
print(t[-1])

# solve ODE
sol = solve_ivp(model, [tiempo[0], tiempo[-1]], [y0], method='RK45', t_eval=tiempo)

t = sol.t
y = sol.y[0]

# plot results
plt.plot(t,y)
plt.xlabel('time')
plt.ylabel('y(t)')

plt.yscale('log')
#plt.xscale('log')
plt.show()

#%% 

from scipy.integrate import quad

def integrando(x, a, b):
    m_0 = 1
    a_b = 1
    w = 1/3
    T_b = 1
    
    num = a_b*2*(1-w)* t * (1+(t/T_b)**2)**(-2/3-w/3)
    den = 3*T_b**2
    
    dot_M_b = num/den
    return dot_M_b 

integral = quad(integrando, -1000, 1000, )



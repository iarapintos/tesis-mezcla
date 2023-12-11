
import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt

#mi ecuacion es dM/dt = dot M_b + dot M_ac

# function that returns dy/dt

m_b = 2.5e20
w = 1/3

def model(t, M):
    a_b = 1
    T_b = 1
    
    num = a_b*2*(1-w)* t * (1+(t/T_b)**2)**((-2/3)-(w/3))
    den = 3*T_b**2
    
    dot_M_b = m_b*(num/den)
    
    #dot_M_b = t
    
    dMdt = dot_M_b 
    return dMdt

# initial condition
#y0 = (1+1000*2)**(2/9)

# time points
t = np.linspace(-1e6, 1e6, 1000)
t_inic = t[0]

#cond inicial 
y0 =  m_b *(1+t_inic**2)**(1/3*(1-w))

# solve ODE
sol = solve_ivp(model, [t[0], t[-1]], [y0], method='RK45', t_eval=t)

#t = sol.t
y = sol.y[0]
print(y)   


# plot results
plt.plot(t,y)
plt.xlabel('time')
plt.ylabel('y(t)')

#plt.xlim(-0.1, 0.1)

#plt.yscale('log')
#plt.xscale('log')
plt.show()

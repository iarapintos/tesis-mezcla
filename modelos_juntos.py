import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt

#mi ecuacion es dM/dt = dot M_b + dot M_rh

m_b = 2.5e20
w = 1/3

def model_hawking(t, M): 
    if M> 1e17:     
        A_M = 5.3e25
        #A_M = 2
        dMdt_h = -A_M / M**2
    else: 
        A_M = 7.8e26
        #A_M = 20
        dMdt_h = -A_M / M**2
    return dMdt_h

def model_bounce(t,M): 
    a_b = 1
    T_b = 1
   
    num = a_b*2*(1-w)* t * (1+(t/T_b)**2)**((-2/3)-(w/3))
    den = 3*T_b**2
   
    dot_M_b = m_b*(num/den) 
    
    return dot_M_b


# function that returns dy/dt
def model_gral(t, M):
    
    a_b = 1
    T_b = 1
   
   
    num = a_b*2*(1-w)* t * (1+(t/T_b)**2)**((-2/3)-(w/3))
    den = 3*T_b**2
    
    dot_M_b = m_b*(num/den)
    
    if M> 1e17:     
        A_M = 5.3e25
        dot_M_rh = -A_M / M**2
    else: 
        A_M = 7.8e26
        dot_M_rh = -A_M / M**2

    
    dMdt = dot_M_b + dot_M_rh
    return dMdt

# time points
t = np.linspace(-1e6, 1e6, 1000)
t_inic = t[0]

# initial condition de las masas
y0 =  m_b *(1+t_inic**2)**((1/3)*(1-w))



# solve ODE
sol_gral = solve_ivp(model_gral, [t[0], t[-1]], [y0], method='RK45', t_eval=t)

#sol_haw = solve_ivp(model_hawking, [t[0], t[-1]], [y0], method='RK45', t_eval=t)

#sol_bounce = solve_ivp(model_bounce, [t[0], t[-1]], [y0], method='RK45', t_eval=t)

#t = sol.t
y_gral = sol_gral.y[0] 
#y_haw = sol_haw.y[0]
#y_bounce = sol_bounce.y[0]

# plot results
plt.plot(t,y_gral, label= 'solucion general')
#plt.plot(t,y_haw, label= 'rad hawking')
#plt.plot(t, y_bounce, label='dinamica del espacio-tiempo')
plt.xlabel('time')
plt.ylabel('y(t)')

plt.legend()

#plt.xlim(-0.1, 0.1)

#plt.yscale('log')
#plt.xscale('log')
plt.show()

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import PillowWriter

# Parameters
L = 2  # Length of the domain
n = 20  # Number of spatial steps
x = np.linspace(0, L, n + 1)  # Spatial discretization
h = x[1] - x[0]  # Spatial step size
T = 1  # Total time
tau = 0.005  # Time step
m = int(T / tau)  # Number of time steps
t = np.linspace(0, T, m) 

gam = tau / h ** 2


# Initial Conditions
def u0(x):
    return np.zeros_like(x)


# Source term
def f(t, x):
    return (1 + np.pi * np.pi * t) * np.sin(np.pi * x)


# True solution
def u_true(t, x):
    return t * np.sin(np.pi * x)


# Initialize u and unew
u = u0(x)
unew = np.zeros_like(u)

# Set up the figure
fig, ax = plt.subplots()
line1, = ax.plot(x, u, 'b-', label='Numerical')  # Numerical solution
line2, = ax.plot(x, u_true(0, x), 'r--', label='True')  # True solution
ax.grid()
ax.legend()

frames = []


for frame in range(m):
    tt = t[frame]  

    
    unew[1:-1] = u[1:-1] + gam * (u[0:-2] - 2 * u[1:-1] + u[2:]) + tau * f(tt, x[1:-1])
    unew[0] = 0.0  # Dirichlet boundary condition
    unew[-1] = 0.0  # Dirichlet boundary condition
    u[:] = unew 

    
    abs_error = np.abs(u - u_true(tt, x))
    max_error = np.max(abs_error)
    print(f'Explicit Scheme: Time t = {tt:.4f}, Max Abs Error = {max_error:.4e}')
    
    line1.set_ydata(u)
    line2.set_ydata(u_true(tt, x))
    ax.set_ylim(np.min(u_true(tt, x)) - 0.1, np.max(u_true(tt, x)) + 0.1)
    ax.set_title(f'Explicit Scheme: Time t = {tt:.4f},  Error = {max_error:.6f}')  
    plt.draw()
    plt.pause(0.01) 
    frames.append(np.array(fig.canvas.buffer_rgba()))  


writer = PillowWriter(fps=20)  
with writer.saving(fig, "Explicit Scheme.gif", dpi=100):
    for frame in frames:
        writer.grab_frame()

plt.close(fig) 



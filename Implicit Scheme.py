import numpy as np
import matplotlib.pyplot as plt
from scipy.sparse import diags
from scipy.sparse.linalg import spsolve
from matplotlib.animation import PillowWriter

# Parameters
L = 2  # Length of the domain
n = 50  # Number of spatial intervals
x = np.linspace(0, L, n + 1)  # Spatial discretization
h = x[1] - x[0]  # Spatial step size
T = 1  # Total time
tau = 0.01  # Time step
m = int(T / tau)  # Number of time steps
t = np.linspace(tau, T, m)  # Time array

# Stability condition
gam = tau / h ** 2
print('gamma = %6.4f' % gam)


# Initial Conditions
def u0(x):
    return np.zeros_like(x)


# Source term
def f(t, x):
    return (1 + np.pi ** 2 * t) * np.sin(np.pi * x)


# True solution
def u_true(t, x):
    return t * np.sin(np.pi * x)


# Construct the system matrix as a sparse matrix
maindiag = (1 + 2 * gam) * np.ones(n + 1)
leftdiag = -gam * np.ones(n)
rghtdiag = -gam * np.ones(n)

# Dirichlet conditions at the boundaries
maindiag[0] = 1
maindiag[-1] = 1
leftdiag[-1] = 0
rghtdiag[0] = 0

# Create the sparse matrix A
A = diags([maindiag, leftdiag, rghtdiag], [0, -1, 1], format='csr')

# Initial condition
u = u0(x)

# Set up the figure for the animation
fig, ax = plt.subplots()
line_numerical, = ax.plot(x, u, 'b-', label='Numerical Solution')
line_true, = ax.plot(x, u_true(0, x), 'r--', label='True Solution')
ax.grid()
ax.set_ylim(-1.1, 1.1)  # Set fixed y limits for clarity
ax.set_title('Heat Equation Simulation: t = 0.00')
ax.set_xlabel('x')
ax.set_ylabel('u')
ax.legend()

# List to hold frames for the GIF
frames = []

# Time loop for the numerical solution
for tt in t:
    # Set up right-hand side
    b = u + tau * f(tt, x)

    # Apply Dirichlet conditions
    b[0] = 0.0
    b[n] = 0.0

    # Solve the linear system using the sparse solver
    u = spsolve(A, b)

    # Update plot lines
    line_numerical.set_ydata(u)
    line_true.set_ydata(u_true(tt, x))

    # Calculate the maximum absolute error
    error = np.max(np.abs(u - u_true(tt, x)))

    # Update the title to include the current time and error
    ax.set_title(f'Implicit Scheme: t = {tt:.2f},  Error = {error:.6f}')

    # Capture the current frame
    plt.draw()
    plt.pause(0.01)  # Allow for a brief pause to render the frame
    frames.append(np.array(fig.canvas.buffer_rgba()))  # Capture the frame

# Save the frames as a GIF using Pillow
writer = PillowWriter(fps=10)  # Set frames per second
with writer.saving(fig, "Implicit Scheme.gif", dpi=100):
    for frame in frames:
        writer.grab_frame()

plt.close(fig)  # Close the figure window



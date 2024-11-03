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
t = np.linspace(0, T, m)  # Time array from 0 to T

# Stability condition for explicit scheme
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

# List to hold frames for the GIF
frames = []

# Update function for each time step
for frame in range(m):
    tt = t[frame]  # Actual time value at the current frame

    # Update the numerical solution
    unew[1:-1] = u[1:-1] + gam * (u[0:-2] - 2 * u[1:-1] + u[2:]) + tau * f(tt, x[1:-1])
    unew[0] = 0.0  # Dirichlet boundary condition
    unew[-1] = 0.0  # Dirichlet boundary condition
    u[:] = unew  # Update u with the new values

    # Calculate the absolute error
    abs_error = np.abs(u - u_true(tt, x))
    max_error = np.max(abs_error)

    # Print the current time and maximum absolute error
    print(f'Explicit Scheme: Time t = {tt:.4f}, Max Abs Error = {max_error:.4e}')  # Print the actual time

    # Update the data for plotting
    line1.set_ydata(u)
    line2.set_ydata(u_true(tt, x))

    # Dynamically adjust the y-axis limits
    ax.set_ylim(np.min(u_true(tt, x)) - 0.1, np.max(u_true(tt, x)) + 0.1)

    # Update the title with the current time and max error
    ax.set_title(f'Explicit Scheme: Time t = {tt:.4f},  Error = {max_error:.6f}')  # Ensure title updates

    # Draw the current figure to capture the frame
    plt.draw()
    plt.pause(0.01)  # Pause to allow the figure to update

    # Save the current frame
    frames.append(np.array(fig.canvas.buffer_rgba()))  # Capture the current frame

# Save the frames as a GIF using Pillow
writer = PillowWriter(fps=20)  # Set frames per second
with writer.saving(fig, "Explicit Scheme.gif", dpi=100):
    for frame in frames:
        writer.grab_frame()

plt.close(fig)  # Close the figure window



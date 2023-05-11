import numpy as np
import matplotlib.pyplot as plt

def advection_simulation():
    # Select numerical parameters
    method = input('Choose a numerical method (FTCS, Lax, Lax-Wendroff): ')
    N = int(input('Enter number of grid points: '))
    L = 1.0  # System size
    h = L / N  # Grid spacing
    c = 1  # Wave speed
    print('Time for wave to move one grid spacing is', h / c)
    tau = float(input('Enter time step: '))
    coeff = -c * tau / (2.0 * h)  # Coefficient used by all schemes
    coefflw = 2 * coeff ** 2  # Coefficient used by L-W scheme
    print('Wave circles system in', L / (c * tau), 'steps')
    nStep = int(input('Enter number of steps: '))

    # Set initial and boundary conditions
    sigma = 0.1  # Width of the Gaussian pulse
    k_wave = np.pi / sigma  # Wave number of the cosine
    x = ((np.arange(N) - 1 / 2) * h) - L / 2  # Coordinates of grid points
    # Initial condition is a Gaussian-cosine pulse
    a = np.cos(k_wave * x) * np.exp(-x ** 2 / (2 * sigma ** 2))
    # Use periodic boundary conditions
    ip = np.roll(np.arange(N), -1)  # ip = i+1 with periodic b.c.
    im = np.roll(np.arange(N), 1)  # im = i-1 with periodic b.c.

    # Initialize plotting variables
    iplot = 1  # Plot counter
    aplot = np.zeros((N, nStep // plotStep + 1))  # Record the initial state
    tplot = np.zeros(nStep // plotStep + 1)  # Record the initial time (t=0)
    nplots = 50  # Desired number of plots
    plotStep = nStep // nplots  # Number of steps between plots

    # Loop over desired number of steps
    for iStep in range(1, nStep + 1):
        # Compute new values of wave amplitude using FTCS, Lax, or Lax-Wendroff method
        if method == 'FTCS':
            a = a + coeff * (a[ip] - a[im])
        elif method == 'Lax':
            a = (a[ip] - a[im]) / (2 * tau) + coeff * (a[ip] - a[im]) / (2 * h)
        elif method == 'Lax-Wendroff':
            a = a + coeff * (a[ip] - a[im]) + coefflw * (a[ip] + a[im] - 2 * a)
        else:
            print('Invalid numerical method:', method)
            return

        # Periodically record a(t) for plotting
        if iStep % plotStep == 0:
            iplot += 1
            aplot[:, iplot] = a  # Record a(i) for plotting
            tplot[iplot] = tau * iStep
            print(iStep, 'out of', nStep, 'steps completed')

    # Plotting
    plt.figure(1)
    plt.plot(x, aplot[:, 0], '-', x, a, '--')
    plt.legend(['Initial', 'Final'])
    plt

import numpy as np
import scipy.integrate as si

def EBOLA(x, T):
    S = x[0]  # Susceptible
    E = x[1]  # Exposed (Latent)
    I = x[2]  # Infected
    D = x[3]  # Dead
    B = x[4]  # Buried
    R = x[5]  # Recovered

    # Parameters for ebola disease model, these have been got from online sources and other works
    N = 60000  # Total population
    Bi = 0.3   # Effective contact rate between infected and susceptible
    Bd = 1.0   # Effective contact rate between dead and susceptible
    f = 1/6    # Rate of progression from latent to infectious
    m = 1/7.5  # Death rate due to Ebola
    mu = 1.0   # Burial rate of dead bodies
    r = 1/10   # Recovery rate in the community

    # The force of infection
    gamma = (Bi * I / N) + (Bd * D / N)

    # The differential equations
    Y = np.zeros(6)
    Y[0] = -gamma * S                          # dS/dt: Susceptible individuals
    Y[1] = gamma * S - f * E                   # dE/dt: Exposed individuals
    Y[2] = f * E - (m + r) * I                 # dI/dt: Infected individuals
    Y[3] = m * I                               # dD/dt: Dead individuals
    Y[4] = mu * D                              # dB/dt: Buried individuals
    Y[5] = r * I                               # dR/dt: Recovered individuals

    return Y

# Duration of simulation
start_time = 0    # Start time
end_time = 500    # End time
time_step = 1500  # Time steps
T = np.linspace(start_time, end_time, time_step)

# Initial conditions: [S, E, I, D, B, R]
solut0 = np.array([59990, 0, 10, 0, 0, 0]) 

# Solving the ODE
solut2 = si.odeint(EBOLA, solut0, T)

# Transpose the solution for easy plotting
final_solut = solut2.T

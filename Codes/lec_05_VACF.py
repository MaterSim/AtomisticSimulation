import numpy as np
import matplotlib.pyplot as plt
from scipy.fftpack import fft


def harmonic_force(r1, r2, k, r_eq):
    # Function to compute force due to harmonic potential
    r12 = np.linalg.norm(r2 - r1)
    force_mag = -k * (r12 - r_eq)
    force = -force_mag * (r2 - r1) / r12
    return force

def MD(r1, r2, v1, v2, N_steps):
    # MD simulation using the Velocity-Verlet algorithm
    velocities = np.zeros([N_steps, 6])
    F12 = harmonic_force(r1, r2, k, r_eq)

    for step in range(N_steps):
        # Compute the force on each atom
        #if step % 100 == 0: print(step, np.linalg.norm(r2-r1)*1e10, F12, v1, v2)

        # Update velocities and positions
        r1 += v1 * dt + 0.5 * F12 / mass * dt ** 2
        r2 += v2 * dt - 0.5 * F12 / mass * dt ** 2

        F12_new = harmonic_force(r1, r2, k, r_eq)

        v1 += 0.5 * (F12 + F12_new) * dt / mass
        v2 -= 0.5 * (F12 + F12_new) * dt / mass

        F12 = F12_new
        #velocities.append(v2 - v1)
        velocities[step][:3] = v1
        velocities[step][3:6] = v2
    # Convert velocities list to a numpy array
    #velocities = np.array(velocities)
    return velocities

if __name__ == "__main__":

    # Parameters for the simulation
    dt = 1e-15              # Time step in seconds (1 fs)

    data = [
            ('O$_2$ (vibration)', 1.16, 1180, 2.66e-26, 2000, False),
            ('O$_2$ (vibration + rotation)', 1.16, 1180, 2.66e-26, 5000, True),
            #('H2', 0.74,  510, 1.67e-27),
            #('N2', 1.10, 2294, 2.33e-26)
           ]
    fig, axs = plt.subplots(2, len(data), figsize=(12, 6))

    for i, (name, r, k, mass, N_steps, rotate) in enumerate(data):
        print(name, r, k, mass)
        # Initial positions of the diatomic molecules (in ang)
        r_eq = r * 1e-10
        r1 = np.array([0.0, 0.0, 0.0])
        r2 = np.array([0.0, 0.0, r_eq*1.2])

        # Initialize velocities
        v1 = np.zeros(3)
        v2 = np.zeros(3)
        if rotate:
            v1[0] += 50
            v2[0] -= 50

        # MD simulation
        velocities = MD(r1, r2, v1, v2, N_steps)

        # Plot VACF
        VACF = np.array([np.dot(velocities[0], velocities[t]) for t in range(N_steps)])
        #VACF = np.correlate(velocities, velocities, mode='full')
        #VACF = VACF[VACF.size // 2:]
        axs[0, i].plot(np.arange(N_steps)*dt*1e12, VACF)
        axs[0, i].set_title(name)
        axs[0, i].set_xlabel('Time (ps)')
        axs[0, i].set_ylabel('VACF')

        # Plot VDOS
        # Frequency axis for the VDOS
        # Fourier transform of the VACF to get VDOS
        VDOS = np.abs(fft(VACF))**2
        freqs = np.fft.fftfreq(N_steps, dt) / 1e12
        axs[1, i].plot(freqs[:N_steps//2], VDOS[:N_steps//2])
        axs[1, i].set_xlabel('Frequency (THz)')
        axs[1, i].set_ylabel('log-VDOS')
        axs[1, i].set_xlim([0, 60])
        axs[1, i].set_yscale('log')
    plt.tight_layout()
    plt.savefig('lec_05_VACF.png')

import numpy as np
import matplotlib.pyplot as plt

def lj_potential(r, epsilon=1.0, sigma1=1.0, sigma2=1.0, m=12, n=6):
    """
    Compute the Lennard-Jones potential for a given distance r.
    """
    return 4 * epsilon * ((sigma1 / r)**m - (sigma2 / r)**n)

# Generate a range of r values
r = np.linspace(0.8, 3, 500)  # Avoid r = 0 to prevent division by zero

# Plotting the LJ potential
plt.plot(r, lj_potential(r, m=12, n=6), label=r'$4 \left[ \left(\frac{1}{r}\right)^{12} - \left(\frac{1}{r}\right)^{6} \right]$')
plt.plot(r, lj_potential(r, m=9, n=3),  label=r'$4 \left[ \left(\frac{1}{r}\right)^{9} - \left(\frac{1}{r}\right)^{3} \right]$')
plt.plot(r, lj_potential(r, m=12, sigma2=0), ls='--', label=r'$4 \left(\frac{1}{r}\right)^{12}$')
plt.plot(r, lj_potential(r, m=6, sigma1=0),  ls='--', label=r'$-4 \left(\frac{1}{r}\right)^{6}$')
plt.axhline(0, color='gray', linestyle='--', linewidth=0.5)  # Horizontal line at V=0
plt.xlabel('Distance r')
plt.ylabel('Potential V(r)')
plt.title('Lennard-Jones Potential')
plt.ylim([-2.0, 10])
plt.legend()
plt.grid(True)
plt.savefig('lec_01_LJ.pdf')

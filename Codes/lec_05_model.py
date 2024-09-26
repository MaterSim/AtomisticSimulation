import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches

def plot_harmonic_model_with_arrows():
    fig, axs = plt.subplots(1, 2, figsize=(10, 5))

    # Define circles' positions
    pos1 = (1, 1)
    pos2 = (3, 1)
    radius = 0.2

    # Subplot 1: Horizontal arrows (indicating initial velocities)
    circle1 = plt.Circle(pos1, radius, color='b', fill=True)
    circle2 = plt.Circle(pos2, radius, color='b', fill=True)

    axs[0].add_patch(circle1)
    axs[0].add_patch(circle2)

    # Draw spring between the circles
    x_spring = np.linspace(1.2, 2.8, 100)
    y_spring = 0.2 * np.sin(10 * np.pi * x_spring)
    axs[0].plot(x_spring, 1 + y_spring, color='gray')

    # Add horizontal velocity arrows
    axs[0].arrow(1.0, 1, -0.3, 0, head_width=0.1, head_length=0.1, fc='r', ec='r')
    axs[0].arrow(3.0, 1, 0.3, 0, head_width=0.1, head_length=0.1, fc='r', ec='r')

    # Formatting subplot 1
    axs[0].set_aspect('equal')
    axs[0].set_xlim(0, 4)
    axs[0].set_ylim(0, 2)
    axs[0].set_title("A Pure Harmonic Oscillator")
    axs[0].axis('off')

    # Subplot 2: Non-horizontal arrows (indicating initial velocities)
    circle1 = plt.Circle(pos1, radius, color='b', fill=True)
    circle2 = plt.Circle(pos2, radius, color='b', fill=True)

    axs[1].add_patch(circle1)
    axs[1].add_patch(circle2)

    # Draw spring between the circles
    axs[1].plot(x_spring, 1 + y_spring, color='gray')

    # Add non-horizontal velocity arrows
    axs[1].arrow(1.0, 1, 0.3, 0.4, head_width=0.1, head_length=0.1, fc='r', ec='r')
    axs[1].arrow(3.0, 1, -0.3, -0.4, head_width=0.1, head_length=0.1, fc='r', ec='r')

    # Formatting subplot 2
    axs[1].set_aspect('equal')
    axs[1].set_xlim(0, 4)
    axs[1].set_ylim(0, 2)
    axs[1].set_title("A Harmonic Oscillator with Rotation")
    axs[1].axis('off')

    plt.tight_layout()
    plt.savefig('lec_05_diatomic.png')

plot_harmonic_model_with_arrows()

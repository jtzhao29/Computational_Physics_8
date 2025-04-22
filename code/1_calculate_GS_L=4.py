from find_ground_state import ising_monte_carlo, visualize_spin, calculate_energy_in_neighbors,energy_final
import numpy as np
import matplotlib.pyplot as plt
import random

if __name__ == "__main__":
    N = 4
    J = 1.0
    T=1
    steps = 1000
    spins = ising_monte_carlo(N, steps,T)
    visualize_spin(spins,N)
    print("energy of ground state:",energy_final(spins,N))
import numpy as np
import matplotlib.pyplot as plt
import os

L = 4
T = 1.0
beta = 1 / T
J = 1
n_steps = 550000
burn_in = 450000

def initial_config(L):
    return np.random.choice([-1, 1], size=(L, L))

def calc_energy(config):
    energy = 0
    for i in range(L):
        for j in range(L):
            S = config[i, j]
            neighbors = config[(i+1)%L, j] + config[i, (j+1)%L] + config[(i-1)%L, j] + config[i, (j-1)%L]
            energy -= J * S * neighbors / 2
    return energy

def metropolis_step(config, beta):
    for _ in range(L*L):
        i = np.random.randint(0, L)
        j = np.random.randint(0, L)
        S = config[i, j]
        neighbors = config[(i+1)%L, j] + config[i, (j+1)%L] + config[(i-1)%L, j] + config[i, (j-1)%L]
        dE = 2 * J * S * neighbors
        if dE <= 0 or np.random.rand() < np.exp(-beta * dE):
            config[i, j] *= -1
    return config

def run_simulation():
    config = initial_config(L)
    energies = []

    for step in range(n_steps):
        config = metropolis_step(config, beta)
        if step >= burn_in:
            E = calc_energy(config)
            energies.append(E)

    return np.array(energies)

if __name__ == "__main__":
    os.makedirs("./images", exist_ok=True)
    energies = run_simulation()
    avg_energy = np.mean(energies)
    print(f"Average Energy (L=4, T=1): {avg_energy:.4f}")

    plt.plot(energies)
    plt.xlabel("MC steps")
    plt.ylabel("Energy")
    plt.title("L=4, T=1 Ising energy change over time")
    plt.grid(True)
    plt.savefig("./images/energy_L4_T1.png")
    plt.show()

import numpy as np
import matplotlib.pyplot as plt
import os
from numba import njit

L = 4
T = 1.0
beta = 1 / T
J = 1
n_steps = 550000
burn_in = 450000
n_samples = 100  # 100个粒子（独立轨迹）

@njit
def initial_config(L):
    return np.random.choice(np.array([-1, 1]), size=(L, L))

@njit
def calc_energy(config, L):
    energy = 0.0
    for i in range(L):
        for j in range(L):
            S = config[i, j]
            neighbors = config[(i+1)%L, j] + config[i, (j+1)%L] + config[(i-1)%L, j] + config[i, (j-1)%L]
            energy -= J * S * neighbors / 2
    return energy

@njit
def metropolis_step(config, beta, L):
    for _ in range(L * L):
        i = np.random.randint(0, L)
        j = np.random.randint(0, L)
        S = config[i, j]
        neighbors = config[(i+1)%L, j] + config[i, (j+1)%L] + config[(i-1)%L, j] + config[i, (j-1)%L]
        dE = 2 * J * S * neighbors
        if dE <= 0 or np.random.rand() < np.exp(-beta * dE):
            config[i, j] *= -1
    return config

@njit
def run_single_simulation(L, beta, n_steps, burn_in):
    config = initial_config(L)
    energies = []

    for step in range(n_steps):
        config = metropolis_step(config, beta, L)
        if step >= burn_in:
            E = calc_energy(config, L)
            energies.append(E)

    return np.array(energies)

def run_ensemble(L, beta, n_steps, burn_in, n_samples):
    """进行n_samples个独立模拟，返回能量样本"""
    all_energies = []
    for _ in range(n_samples):
        energies = run_single_simulation(L, beta, n_steps, burn_in)
        print(f"Sample {_+1}/{n_samples} completed")
        all_energies.append(np.mean(energies))  # 对每个轨迹内部平均
    return np.array(all_energies)

if __name__ == "__main__":
    os.makedirs("./images", exist_ok=True)
    energies = run_ensemble(L, beta, n_steps, burn_in, n_samples)

    avg_energy = np.mean(energies)
    std_energy = np.std(energies)

    print(f"Ensemble Average Energy (L=4, T=1): {avg_energy:.4f} ± {std_energy:.4f}")

    plt.plot(energies, marker='o', linestyle='none')
    plt.xlabel("Sample index")
    plt.ylabel("Average Energy per sample")
    plt.title("Ensemble of energies (L=4, T=1)")
    plt.grid(True)
    plt.savefig("./images/ensemble_energy_L4_T1.png")
    plt.show()

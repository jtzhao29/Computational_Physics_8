import numpy as np
import matplotlib.pyplot as plt

L = 16
T = 1/(0.5 * np.log(1 + np.sqrt(2)) ) # 临界温度
n_steps = 2000000
J = 1

def init_lattice(L):
    return np.random.choice([-1, 1], size=(L, L))

def calc_energy(lattice):
    E = 0
    for i in range(L):
        for j in range(L):
            S = lattice[i, j]
            neighbors = lattice[(i+1)%L, j] + lattice[i, (j+1)%L] + \
                        lattice[(i-1)%L, j] + lattice[i, (j-1)%L]
            E -= J * S * neighbors
    return E / 2  # 避免重复计数

def metropolis_step(lattice, T):
    for _ in range(L * L):
        i, j = np.random.randint(0, L, size=2)
        S = lattice[i, j]
        neighbors = lattice[(i+1)%L, j] + lattice[i, (j+1)%L] + \
                    lattice[(i-1)%L, j] + lattice[i, (j-1)%L]
        dE = 2 * J * S * neighbors
        if dE <= 0 or np.random.rand() < np.exp(-dE / T):
            lattice[i, j] *= -1
    return lattice

energies = []
lattice = init_lattice(L)
for t in range(n_steps):
    lattice = metropolis_step(lattice, T)
    E = calc_energy(lattice)
    energies.append(E/L**2)  # 归一化能量

plt.figure(figsize=(6,4))
plt.plot(energies)
plt.xlabel("Monte Carlo Time Step $t$")
plt.ylabel("Average Energy per Spin $\\langle E(t) \\rangle$")
plt.title("Energy Relaxation ($L=16$, $T=T_c$)")
plt.grid()
plt.tight_layout()
plt.savefig("./images/energy_vs_time_L16.png")
plt.show()
import numpy as np
import matplotlib.pyplot as plt
from utils import initial_lattice, compute_energy, compute_magnetization, metropolis_step

def run_simulation(L, T, n_eq_steps=5000, n_meas_steps=10000):
    beta = 1.0 / T
    lattice = initial_lattice(L)

    # 平衡过程
    for _ in range(n_eq_steps):
        lattice = metropolis_step(lattice, beta)

    E_total, E2_total = 0, 0
    M_total, M2_total, Mabs_total = 0, 0, 0
    for _ in range(n_meas_steps):
        lattice = metropolis_step(lattice, beta)
        E = compute_energy(lattice)
        M = compute_magnetization(lattice)
        E_total += E
        E2_total += E**2
        M_total += M
        M2_total += M**2
        Mabs_total += abs(M)

    N = L * L
    avg_E = E_total / n_meas_steps
    avg_E2 = E2_total / n_meas_steps
    avg_M2 = M2_total / n_meas_steps
    avg_Mabs = Mabs_total / n_meas_steps

    c = beta**2 * (avg_E2 - avg_E**2) / N
    chi = beta * (avg_M2 - avg_Mabs**2) / N
    m2 = avg_M2 / (N**2)

    return avg_E / N, c, chi, m2

def plot_observables(T_list, L_list):
    for L in L_list:
        E_list, c_list, chi_list, m2_list = [], [], [], []
        for T in T_list:
            E, c, chi, m2 = run_simulation(L, T)
            E_list.append(E)
            c_list.append(c)
            chi_list.append(chi)
            m2_list.append(m2)

        plt.figure(1)
        plt.plot(T_list, m2_list, label=f"L={L}")
        plt.figure(2)
        plt.plot(T_list, c_list, label=f"L={L}")
        plt.figure(3)
        plt.plot(T_list, chi_list, label=f"L={L}")

    for i, name in enumerate(["m2", "c", "chi"], start=1):
        plt.figure(i)
        plt.xlabel("Temperature T")
        plt.ylabel(name)
        plt.title(f"{name} vs T")
        plt.legend()
        plt.grid(True)
        plt.savefig(f"images/{name}.png")
        plt.close()

if __name__ == '__main__':
    T_list = np.arange(1.5, 3.1, 0.1)
    L_list = [8, 16, 32]
    plot_observables(T_list, L_list)

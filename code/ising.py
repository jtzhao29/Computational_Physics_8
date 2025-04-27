import os
import numpy as np
import matplotlib.pyplot as plt
from utils import initialize_lattice, calculate_energy, calculate_magnetization, metropolis_step


def run_simulation(L, T, n_eq=600000, n_meas=1600000):
    """
    对 L x L 格子在温度 T 下进行 MCMC 模拟，返回能量和磁化强度统计
    """
    beta = 1.0 / T
    lattice = initialize_lattice(L)
    # 平衡热化
    for _ in range(n_eq):
        lattice = metropolis_step(lattice, beta)
    # 测量
    E_list, M_list = [], []
    for _ in range(n_meas):
        lattice = metropolis_step(lattice, beta)
        E_list.append(calculate_energy(lattice))
        M_list.append(calculate_magnetization(lattice))
    return np.array(E_list), np.array(M_list)

def plot_observables(T_list, L_list):
    os.makedirs("images", exist_ok=True)
    # Plot 1: Magnetization squared
    plt.figure(figsize=(6,4))
    for L in L_list:
        m2_all = []
        for T in T_list:
            E, M = run_simulation(L, T)
            N = L * L
            m2_all.append(np.mean(M**2) / N**2)
        plt.plot(T_list, m2_all, label=f"L={L}")
    plt.xlabel("Temperature T")
    plt.ylabel(r"$\langle m^2 \rangle$")
    plt.title(r"Magnetization squared $\langle m^2 \rangle$")
    plt.legend()
    plt.grid(True)
    plt.savefig("./images/magnetization_squared.png")
    plt.close()

    # Plot 2: Specific heat
    plt.figure(figsize=(6,4))
    for L in L_list:
        c_all = []
        for T in T_list:
            E, M = run_simulation(L, T)
            N = L * L
            beta = 1.0 / T
            c_all.append(beta**2 * (np.mean(E**2) - np.mean(E)**2) / N)
        plt.plot(T_list, c_all, label=f"L={L}")
    plt.xlabel("Temperature T")
    plt.ylabel("Specific heat c")
    plt.title("Specific heat as a function of temperature")
    plt.legend()
    plt.grid(True)
    plt.savefig("./images/specific_heat.png")
    plt.close()

    # Plot 3: Susceptibility
    plt.figure(figsize=(6,4))
    for L in L_list:
        chi_all = []
        for T in T_list:
            E, M = run_simulation(L, T)
            N = L * L
            beta = 1.0 / T
            chi_all.append(beta * (np.mean(M**2) - np.mean(np.abs(M))**2) / N)
        plt.plot(T_list, chi_all, label=f"L={L}")
    plt.xlabel("Temperature T")
    plt.ylabel(r"Susceptibility $\chi$")
    plt.title("Susceptibility as a function of temperature")
    plt.legend()
    plt.grid(True)
    plt.savefig("./images/susceptibility.png")
    plt.close()
    
if __name__ == '__main__':
    T_list = np.arange(1.5, 3.1, 0.1)
    L_list = [8, 16, 32]
    plot_observables(T_list, L_list)

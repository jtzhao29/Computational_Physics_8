import numpy as np
import matplotlib.pyplot as plt
from scipy.ndimage import uniform_filter1d
from scipy.optimize import curve_fit
from numba import njit

# 参数设置
L = 16
T = 0.5 * np.log(1 + np.sqrt(2))  # 临界温度
J = 1
n_steps = 2000
n_runs = 200

@njit
def init_lattice(L):
    return np.random.choice(np.array([-1, 1]), size=(L, L))

@njit
def calc_energy(lattice):
    L = lattice.shape[0]
    E = 0.0
    for i in range(L):
        for j in range(L):
            S = lattice[i, j]
            neighbors = lattice[(i+1)%L, j] + lattice[i, (j+1)%L] + \
                        lattice[(i-1)%L, j] + lattice[i, (j-1)%L]
            E -= J * S * neighbors
    return E / 2  # 避免重复计算

@njit
def metropolis_step(lattice, T):
    L = lattice.shape[0]
    for _ in range(L * L):
        i = np.random.randint(0, L)
        j = np.random.randint(0, L)
        S = lattice[i, j]
        neighbors = lattice[(i+1)%L, j] + lattice[i, (j+1)%L] + \
                    lattice[(i-1)%L, j] + lattice[i, (j-1)%L]
        dE = 2 * J * S * neighbors
        if dE <= 0 or np.random.rand() < np.exp(-dE / T):
            lattice[i, j] *= -1
    return lattice

@njit
def run_single_simulation(L, T, n_steps):
    lattice = init_lattice(L)
    energies = np.zeros(n_steps)
    for t in range(n_steps):
        lattice = metropolis_step(lattice, T)
        energies[t] = calc_energy(lattice) / L**2
    return energies

def run_multiple_simulations():
    all_energies = []
    for run in range(n_runs):
        print(f"Running simulation {run+1}/{n_runs}")
        energies = run_single_simulation(L, T, n_steps)
        all_energies.append(energies)
    all_energies = np.array(all_energies)
    mean_energies = np.mean(all_energies, axis=0)
    return mean_energies

def power_law(t, a, b):
    return a * t**b

def exp_decay(t, a, tau):
    return a * np.exp(-t / tau)

if __name__ == "__main__":
    energies = run_multiple_simulations()
    energies_smooth = uniform_filter1d(energies, size=1000)

    E_inf = np.mean(energies[-5000:])
    Delta = energies - E_inf

    plt.figure(figsize=(7,5))
    plt.plot(energies, label="Energy (smoothed)")
    plt.axhline(E_inf, color='r', linestyle='--', label=f"E(\u221E)={E_inf:.4f}")
    plt.xlabel("Monte Carlo Time Step")
    plt.ylabel("Average Energy per Spin")
    plt.title(f"Energy Relaxation (L={L}, T=T_c)")
    plt.legend()
    plt.grid()
    plt.tight_layout()
    plt.savefig("./images/energy_relaxation_smooth.png")
    plt.show()

    t_arr = np.arange(len(Delta))
    plt.figure(figsize=(7,5))
    plt.semilogy(t_arr, np.abs(Delta), label=r"$|\Delta(t)|$")
    plt.xlabel("Monte Carlo Time Step")
    plt.ylabel(r"$|\Delta(t)|$")
    plt.title("Relaxation Difference (log scale)")
    plt.grid()
    plt.legend()
    plt.tight_layout()
    plt.savefig("./images/energy_relaxation_logscale.png")
    plt.show()

    mask = (t_arr <1000)
    try:
        popt, pcov = curve_fit(power_law, t_arr[mask], np.abs(Delta[mask]), p0=(1.0, -0.5))
        a_fit, b_fit = popt
        print(f"Power-law fit: \u0394(t) ~ t^{b_fit:.3f}")

        plt.figure(figsize=(7,5))
        plt.semilogy(t_arr, np.abs(Delta), label="Data")
        plt.semilogy(t_arr, power_law(t_arr, *popt), '--', label=fr"Fit: $\Delta(t) \propto t^{{{b_fit:.2f}}}$")
        plt.xlabel("Monte Carlo Time Step")
        plt.ylabel(r"$|\Delta(t)|$")
        plt.title("Relaxation Fit (Power Law)")
        plt.legend()
        plt.grid()
        plt.tight_layout()
        plt.savefig("./images/energy_relaxation_fit.png")
        plt.show()
    except Exception as e:
        print(f"Fitting failed: {e}")

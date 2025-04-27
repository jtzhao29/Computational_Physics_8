import numpy as np
import matplotlib.pyplot as plt
from scipy.ndimage import uniform_filter1d
from scipy.optimize import curve_fit
from numba import njit

# 参数设置
L = 16
T = 1/(0.5 * np.log(1 + np.sqrt(2)) ) # 临界温度
J = 1
n_steps = 20000
n_runs = 1000
@njit
def init_lattice(L):
    return np.random.choice(np.array([-1, 1]), size=(L, L))

@njit
def calc_energy(lattice,L):
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
def metropolis_step(lattice, T,L):
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
        lattice = metropolis_step(lattice, T,L)
        energies[t] = calc_energy(lattice,L) / L**2
    return energies

def run_multiple_simulations(L):
    all_energies = []
    for run in range(n_runs):
        print(f"Running simulation {run+1}/{n_runs}")
        energies = run_single_simulation(L, T, n_steps)
        all_energies.append(energies)
    all_energies = np.array(all_energies)
    mean_energies = np.mean(all_energies, axis=0)
    return mean_energies

def power_law(t, a, b,c):
    return a * t**b+c

def exp_decay(t, a, tau):
    return a * np.exp(-t / tau)

def plot(energies):
    energies_smooth = uniform_filter1d(energies, size=1000)

    E_inf = np.mean(energies[-1000:])
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

    Delta_inf = np.mean(abs(Delta[-1000:]))
    t_arr = np.arange(len(Delta))
    plt.figure(figsize=(7,5))
    plt.semilogy(t_arr, np.abs(Delta), label=r"$|\Delta(t)|$")
    plt.axhline((Delta_inf), color='r', linestyle='--', label=f"$\\hat{{|\\Delta(t)| }}= {((Delta_inf)):.10f}$")
    x = 170
    plt.axvline(x, color='r', linestyle='--', label=f" t={x}")
    plt.xlabel("Monte Carlo Time Step")
    plt.ylabel(r"$|\Delta(t)|$")
    # plt.rc('text', usetex=True)  # Enable LaTeX rendering for labels
    plt.title("Relaxation Difference (log scale)")
    plt.grid()
    plt.legend()
    plt.tight_layout()
    plt.savefig("./images/energy_relaxation_logscale.png")
    plt.show()

    # mask = (t_arr < 200)
    # try:
    #     popt, pcov = curve_fit(power_law, t_arr[mask], np.abs(Delta[mask]), p0=(1.0, -0.5))
    #     a_fit, b_fit,c_fit = popt
    #     print(f"Power-law fit: \u0394(t) ~ t^{b_fit:.3f}")

    #     plt.figure(figsize=(7,5))
    #     plt.semilogy(t_arr, np.abs(Delta), label="Data")
    #     plt.semilogy(t_arr, power_law(t_arr, *popt), '--', label=fr"Fit: $\Delta(t) \propto t^{{{b_fit:.2f}}}$")
    #     plt.xlabel("Monte Carlo Time Step")
    #     plt.ylabel(r"$|\Delta(t)|$")
    #     plt.title("Relaxation Fit (Power Law)")
    #     plt.legend()
    #     plt.grid()
    #     plt.tight_layout()
    #     plt.savefig("./images/energy_relaxation_fit.png")
    #     plt.show()
    # except Exception as e:
    #     print(f"Fitting failed: {e}")
    # --- 新增：指数衰减拟合 Delta(t) ---

    # 选择拟合的时间区间
    fit_start = 50   # 跳过最开始涨落剧烈部分，比如50步之后
    fit_end = 1500   # 不一定要拟合到最后，可以自己调整
    mask = (t_arr >= fit_start) & (t_arr <= fit_end)

    # 使用指数衰减模型拟合
    try:
        popt, pcov = curve_fit(exp_decay, t_arr[mask], np.abs(Delta[mask]), p0=(Delta[fit_start], 500))
        a_fit, tau_fit = popt
        print(f"exp decay: Δ(t) ≈ {a_fit:.4e} * exp(-t/{tau_fit:.2f})")

        # 绘制拟合曲线
        plt.figure(figsize=(7,5))
        plt.semilogy(t_arr, np.abs(Delta), label=r"$|\Delta(t)|$ original data")
        plt.semilogy(t_arr, exp_decay(t_arr, *popt), '--r', label=fr"exp fit: $\tau = {tau_fit:.2f}$")
        plt.axhline((Delta_inf), color='g', linestyle='--', label=f"$\\hat{{|\\Delta(t)|}}={((Delta_inf)):.2e}$")
        plt.xlabel("Monte Carlo Time Step")
        plt.ylabel(r"$|\Delta(t)|$")
        plt.title("Relaxation Difference with Exponential Fit (log scale)")
        plt.legend()
        plt.grid()
        plt.tight_layout()
        plt.savefig("./images/energy_relaxation_expfit.png")
        plt.show()


    except Exception as e:
        print(f"指数拟合失败: {e}")
def plot_for_diff_L(L_list: list) -> None:
    energies_dict = {}
    Delta_dict = {}

    # Step 1: 跑模拟并保存数据
    for L in L_list:
        energies = run_multiple_simulations(L)
        # energies_smooth = uniform_filter1d(energies, size=1000)

        E_inf = np.mean(energies[-1000:])
        Delta = energies - E_inf

        energies_dict[L] = energies
        Delta_dict[L] = Delta

    # Step 2: 画能量随时间变化
    plt.figure(figsize=(7,5))
    for L in L_list:
        energies = energies_dict[L]
        E_inf = np.mean(energies[-1000:])
        plt.plot(energies, label=f"Energy when L={L}")
        plt.axhline(E_inf, linestyle='--', label=f"E(∞)={E_inf:.4f} when L={L}")

    plt.xlabel("Monte Carlo Time Step")
    plt.ylabel("Average Energy per Spin")
    plt.title(f"Energy Relaxation for Different L at T_c")
    plt.legend()
    plt.grid()
    plt.tight_layout()
    plt.savefig(f"./images/energy_relaxation_smooth_when_L={L_list}.png")
    plt.show()

    # Step 3: 画Delta随时间变化（对数坐标）
    plt.figure(figsize=(7,5))
    for L in L_list:
        Delta = Delta_dict[L]
        t_arr = np.arange(len(Delta))
        Delta_inf = np.mean(np.abs(Delta[-1000:]))

        plt.semilogy(t_arr, np.abs(Delta), label=rf"$|\Delta(t)|$, L={L}")
        plt.axhline(Delta_inf, linestyle='--', label=rf"$\hat{{|\Delta(t)|}}={Delta_inf:.2e}$, L={L}")

    plt.xlabel("Monte Carlo Time Step")
    plt.ylabel(r"$|\Delta(t)|$")
    plt.title("Relaxation Difference (log scale) for Different L")
    plt.grid()
    plt.legend()
    plt.tight_layout()
    plt.savefig(f"./images/energy_relaxation_logscale_when_L={L_list}.png")
    plt.show()

        
            

if __name__ == "__main__":
    # energies = run_multiple_simulations(L)
    # plot(energies)
    L_list = [8, 16, 32, 64]
    plot_for_diff_L(L_list)
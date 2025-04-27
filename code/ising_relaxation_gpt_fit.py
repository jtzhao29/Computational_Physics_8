import numpy as np
import matplotlib.pyplot as plt
from scipy.ndimage import uniform_filter1d
from scipy.optimize import curve_fit

# 参数设置
L = 16
T = 1/(0.5 * np.log(1 + np.sqrt(2)) ) # 临界温度
J = 1
n_steps = 200000  # 增加步数
n_runs = 5        # 多次独立实验

# 初始化格子
def init_lattice(L):
    return np.random.choice([-1, 1], size=(L, L))

# 计算能量
def calc_energy(lattice):
    E = 0
    for i in range(L):
        for j in range(L):
            S = lattice[i, j]
            neighbors = lattice[(i+1)%L, j] + lattice[i, (j+1)%L] + \
                        lattice[(i-1)%L, j] + lattice[i, (j-1)%L]
            E -= J * S * neighbors
    return E / 2  # 避免重复计算

# Metropolis单步

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

# 单次模拟，返回时间序列

def run_single_simulation():
    lattice = init_lattice(L)
    energies = []
    for t in range(n_steps):
        lattice = metropolis_step(lattice, T)
        E = calc_energy(lattice)
        energies.append(E / L**2)  # 归一化能量
    return np.array(energies)

# 多次模拟，取平均

def run_multiple_simulations():
    all_energies = []
    for run in range(n_runs):
        print(f"Running simulation {run+1}/{n_runs}")
        energies = run_single_simulation()
        all_energies.append(energies)
    all_energies = np.array(all_energies)
    mean_energies = np.mean(all_energies, axis=0)
    return mean_energies

# 矩形拟合函数： power law 或指数衰减
def power_law(t, a, b):
    return a * t**b

def exp_decay(t, a, tau):
    return a * np.exp(-t / tau)

if __name__ == "__main__":
    energies = run_multiple_simulations()
    
    # 平滑处理
    energies_smooth = uniform_filter1d(energies, size=1000)

    # 维持能量估计
    E_inf = np.mean(energies_smooth[-5000:])
    Delta = energies_smooth - E_inf

    # 绘制能量弥转曲线
    plt.figure(figsize=(7,5))
    plt.plot(energies_smooth, label="Energy (smoothed)")
    plt.axhline(E_inf, color='r', linestyle='--', label=f"E(\u221E)={E_inf:.4f}")
    plt.xlabel("Monte Carlo Time Step")
    plt.ylabel("Average Energy per Spin")
    plt.title(f"Energy Relaxation (L={L}, T=T_c)")
    plt.legend()
    plt.grid()
    plt.tight_layout()
    plt.savefig("./images/energy_relaxation_smooth.png")
    plt.show()

    # 绘制Delta(t)对数坐标图
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

    # 尝试用power law拟合
    mask = (t_arr > 1000)  # 克隆定前期
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
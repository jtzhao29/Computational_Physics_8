import numpy as np

def initialize_lattice(L):
    """
    初始化LxL的Ising模型晶格，随机选择+1或-1作为自旋状态
    """
    return 2 * np.random.randint(2, size=(L, L)) - 1

def calculate_energy(lattice, J=1):
    """
    计算当前状态下的能量
    """
    energy = 0
    L = lattice.shape[0]
    for i in range(L):
        for j in range(L):
            # 周期性边界条件
            neighbor_sum = lattice[(i + 1) % L, j] + lattice[i, (j + 1) % L] + \
                           lattice[(i - 1) % L, j] + lattice[i, (j - 1) % L]
            energy -= lattice[i, j] * neighbor_sum * J
    return energy / 2  # 因为每个相邻对会被重复计算

def calculate_magnetization(lattice):
    """
    计算磁化强度
    """
    return np.sum(lattice)

def metropolis_step(lattice, beta, J=1):
    """
    执行一次Metropolis步
    """
    L = lattice.shape[0]
    i, j = np.random.randint(L), np.random.randint(L)
    spin = lattice[i, j]

    # 周期性边界条件
    neighbor_sum = lattice[(i + 1) % L, j] + lattice[i, (j + 1) % L] + \
                   lattice[(i - 1) % L, j] + lattice[i, (j - 1) % L]
    
    delta_energy = 2 * J * spin * neighbor_sum

    # Metropolis判断准则
    if delta_energy < 0 or np.random.rand() < np.exp(-beta * delta_energy):
        lattice[i, j] = -spin  # 翻转自旋

def simulate_mcmc(L, T, num_steps=10000):
    """
    通过Monte Carlo方法模拟一个LxL的Ising模型
    """
    lattice = initialize_lattice(L)
    beta = 1 / T
    energies = []
    magnetizations = []
    
    for step in range(num_steps):
        metropolis_step(lattice, beta)
        
        # 记录当前能量与磁化强度
        energies.append(calculate_energy(lattice))
        magnetizations.append(calculate_magnetization(lattice))

    return np.array(energies), np.array(magnetizations)

import numpy as np

def initialize_lattice(L):
    """
    初始化 L x L 的 Ising 模型格子，自旋取值 +1 或 -1
    """
    return 2 * np.random.randint(2, size=(L, L)) - 1


def calculate_energy(lattice, J=1):
    """
    计算当前格子的总能量，周期性边界条件
    """
    L = lattice.shape[0]
    energy = 0
    for i in range(L):
        for j in range(L):
            S = lattice[i, j]
            neighbors = (
                lattice[(i + 1) % L, j] + lattice[i, (j + 1) % L]
                + lattice[(i - 1) % L, j] + lattice[i, (j - 1) % L]
            )
            energy -= J * S * neighbors
    return energy / 2  # 每对交互计算两次，故除以 2


def calculate_magnetization(lattice):
    """
    计算当前格子的总磁化强度
    """
    return np.sum(lattice)


def metropolis_step(lattice, beta, J=1):
    """
    在格子上执行一次 Metropolis 更新
    """
    L = lattice.shape[0]
    i, j = np.random.randint(L), np.random.randint(L)
    S = lattice[i, j]
    neighbors = (
        lattice[(i + 1) % L, j] + lattice[i, (j + 1) % L]
        + lattice[(i - 1) % L, j] + lattice[i, (j - 1) % L]
    )
    dE = 2 * J * S * neighbors
    if dE <= 0 or np.random.rand() < np.exp(-beta * dE):
        lattice[i, j] = -S
    return lattice

import numpy as np
from numba import njit

@njit
def initialize_lattice(L):
    """
    初始化 L x L 的 Ising 模型格子，自旋取值 +1 或 -1
    """
    lattice = np.empty((L, L), dtype=np.int32)
    for i in range(L):
        for j in range(L):
            if np.random.rand() < 0.5:
                lattice[i, j] = 1
            else:
                lattice[i, j] = -1
    return lattice

@njit
def calculate_energy(lattice, J=1):
    """
    计算当前格子的总能量，周期性边界条件
    """
    L = lattice.shape[0]
    energy = 0.0
    for i in range(L):
        for j in range(L):
            S = lattice[i, j]
            neighbors = (
                lattice[(i + 1) % L, j] +
                lattice[i, (j + 1) % L] +
                lattice[(i - 1) % L, j] +
                lattice[i, (j - 1) % L]
            )
            energy -= J * S * neighbors
    return energy / 2.0

@njit
def calculate_magnetization(lattice):
    """
    计算当前格子的总磁化强度
    """
    return np.sum(lattice)

@njit
def metropolis_step(lattice, beta, J=1):
    """
    在格子上执行一次 Metropolis 更新
    """
    L = lattice.shape[0]
    i = np.random.randint(0, L)
    j = np.random.randint(0, L)
    S = lattice[i, j]
    neighbors = (
        lattice[(i + 1) % L, j] +
        lattice[i, (j + 1) % L] +
        lattice[(i - 1) % L, j] +
        lattice[i, (j - 1) % L]
    )
    dE = 2 * J * S * neighbors
    if dE <= 0 or np.random.rand() < np.exp(-beta * dE):
        lattice[i, j] = -S
    return lattice

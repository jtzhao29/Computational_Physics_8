import numpy as np
import itertools

def compute_energy(config, L):
    """计算给定自旋配置的能量，周期性边界"""
    E = 0
    config = np.array(config).reshape((L, L))
    for i in range(L):
        for j in range(L):
            S = config[i, j]
            # 最近邻：右、下
            E -= S * (config[i, (j + 1) % L] + config[(i + 1) % L, j])
    return E

def exact_ising(L, T):
    """精确计算配分函数、平均能量和自由能"""
    beta = 1.0 / T
    states = itertools.product([-1, 1], repeat=L*L)
    
    Z = 0.0
    E_avg = 0.0
    
    for state in states:
        E = compute_energy(state, L)
        weight = np.exp(-beta * E)
        Z += weight
        E_avg += E * weight
    
    E_avg /= Z
    F = -np.log(Z) / beta
    return Z, E_avg, F

if __name__ == "__main__":
    L = 4
    T = 1.0
    Z, E_avg, F = exact_ising(L, T)
    print(f"Partition function Z = {Z:.5e}")
    print(f"Average energy <E> = {E_avg:.5f}")
    print(f"Free energy F = {F:.5f}")

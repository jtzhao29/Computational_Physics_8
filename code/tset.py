import numpy as np
from scipy.optimize import fsolve  # 修正导入


# $$s(m)=-\left[\frac{1+m}{2}\ln\left(\frac{1+m}{2}\right)+\frac{1-m}{2}\ln\left(\frac{1-m}{2}\right)\right]$$

def calculate_entropy(m):
    """
    计算给定磁化强度 m 的熵
    """
    if m == 1 or m == -1:
        return 0.0
    else:
        return -(0.5 * (1 + m) * np.log(0.5 * (1 + m)) + 0.5 * (1 - m) * np.log(0.5 * (1 - m)))
    

def calculate_energy(m)->float:
    #-\frac{1}{2} \cdot 4 \cdot 1 \cdot 16 \cdot (0.999329)^2 
    return -0.5 * 4 * 1 * 16 * (m)**2

def slove_func(m)->float:
    #m = \tanh(4m)
    f = m - np.tanh(4*m)
    return f




if __name__ == "__main__":
    m =0.9993256730151082
    # entropy = calculate_entropy(m)
    # print(f"Entropy for m={m}: {entropy:.4f}")
    # a =calculate_energy(m)
    # print(f"Energy: {a:.4f}")
    # m = fsolve(slove_func, 0.999329, xtol=1e-10)[0]
    # print(m)
    print(-31.9569 - 0.048 )
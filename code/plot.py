import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import linregress

L_list = [8, 16, 32, 64]
T_list = [100, 200, 700, 2400]

def plot_results(L, T):
    # 对数变换数据
    log_L = np.log(L)
    log_T = np.log(T)

    # 直线拟合
    slope, intercept, r_value, p_value, std_err = linregress(log_L, log_T)

    # 输出拟合表达式
    print(f"Fitted line: log(T) = {slope:.4f} * log(L) + {intercept:.4f}")
    print(f"Or equivalently: T = exp({intercept:.4f}) * L^{slope:.4f}")

    # 绘图
    plt.figure(figsize=(10, 6))
    plt.scatter(log_L, log_T, label="Data points", color="blue", marker='o')
    plt.plot(log_L, slope * log_L + intercept, label=f"Fit: log(T) = {slope:.4f} * log(L) + {intercept:.4f}", color="red")
    plt.xlabel("log(L)")
    plt.ylabel("log(T)")
    plt.title("Log-Log Plot with Linear Fit")
    plt.legend()
    plt.grid()
    plt.savefig("./images/linear_fit_loglog.png")
    plt.close()

if __name__ == "__main__":
    plot_results(L_list, T_list)
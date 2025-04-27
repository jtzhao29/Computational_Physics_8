import json
import numpy as np
import matplotlib.pyplot as plt

with open('data\energies_L16.json', 'r') as f:
    data = json.load(f)

print(len(data))

time = np.arange(len(data))
energy = np.array(data)

energy = abs(energy - np.mean(energy[-500:]) ) # 去除均值
energy = energy[:200]


def plot_Log( energy):
    time = np.arange(len(energy))
    log_energy = np.log(energy)
    k,b = np.polyfit(time, log_energy, 1) # 线性拟合
    print(f'k: {k}, b: {b}')
    plt.plot(time, energy)
    plt.plot(time, np.exp(k*time + b), label=f'Fitted line: log(y) = {k:.4f}*x+{b:.4f}', color='red')
    plt.yscale('log')
    plt.xlabel('Time',fontsize=15)
    plt.ylabel('Energy',fontsize=15)
    plt.title('Log of Energy vs Time',fontsize=15)
    plt.legend()
    path = 'images\log_energy_vs_time.png'
    plt.savefig(path)
    
    plt.show()

    
if __name__ == "__main__":
    plot_Log(energy)
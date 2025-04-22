import numpy as np
import matplotlib.pyplot as plt

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, Circle


# 取格点尺寸大小为N，计算基态
# 通过三维数组【N,N,2]来分别代表A,B两种格点

def ising_monte_carlo(N:int,step:int,T:float)->np.ndarray:
    """
    用于计算基态，使用蒙特卡洛方法，
    """
    system = np.random.choice([-1, 1], size=(2,N, N))
    # print(system)
    for step_num in range(step):
        type_gedian = np.random.randint(0, 2) #0代表A，1代表B
        # 随机选取一个格点
        i, j = np.random.randint(0, N, 2)
        energy_change = calculate_energy_change_in_neighbors(system, i, j, type_gedian,N)
        w = calculate_change_possibility(energy_change,T)
        if np.random.rand() < w:
            system[type_gedian, i, j] *= -1
    return system




def calculate_energy_change_in_neighbors(system:np.ndarray,x:int,y:int,type_gedian:int,N:int)->int:
    """
    计算邻居的能量变化
    """
    neighbor = find_neibor(x,y,type_gedian,N)
    energy_old = 0
    for type_gedian_nei,x_nei,y_nei in neighbor:
        energy_old += system[type_gedian_nei,x_nei,y_nei] * system[type_gedian,x,y]
    

    return  energy_old*-2


def calculate_energy_in_neighbors(system:np.ndarray,x:int,y:int,type_gedian:int,N:int)->int:
    """
    计算邻居的能量变化
    """
    neighbor = find_neibor(x,y,type_gedian,N)
    # print(f"for {type_gedian} in {x},{y},neighbor num is {len(neighbor)}")
    energy = 0
    for type_gedian_nei,x_nei,y_nei in neighbor:
        energy += system[type_gedian_nei,x_nei,y_nei] * system[type_gedian,x,y]
    
   
    return  energy


def calculate_change_possibility(energy_change:int,T=1e-30)->float:
    """
    计算能量变化的可能性
    """
    beta = 1/T
    if energy_change < 0:
        return 1.0
    else:
        return np.exp(-energy_change*beta)



def visualize_spin(spin: np.ndarray, N: int):
    """
    使用 FancyArrowPatch 可视化 AB 格点的自旋状态。
    
    参数:
        spin: np.ndarray, 形状为 (2, N, N)，第一维为 A/B 子格。
        N: int, 元胞的边长数量。
    """
    fig, ax = plt.subplots(figsize=(8, 8))

    for i in range(N):
        for j in range(N):
            # A 格点
            x_a, y_a = 2 * i, 2 * j
            s_a = spin[0, i, j]
            ax.add_patch(Circle((x_a, y_a), 0.3, color='pink', zorder=2))

            arrow_a = FancyArrowPatch(
                (x_a, y_a), (x_a, y_a + 0.2 * s_a),
                arrowstyle='-|>',
                color='white',
                mutation_scale=15,
                linewidth=1.5,
                zorder=3
            )
            ax.add_patch(arrow_a)

            # B 格点
            x_b, y_b = 2 * i + 1, 2 * j + 1
            s_b = spin[1, i, j]
            ax.add_patch(Circle((x_b, y_b), 0.3, color='blue', zorder=2))

            arrow_b = FancyArrowPatch(
                (x_b, y_b), (x_b, y_b + 0.2 * s_b),
                arrowstyle='-|>',
                color='white',
                mutation_scale=15,
                linewidth=1.5,
                zorder=3
            )
            ax.add_patch(arrow_b)

    # 添加 A-A 横向黑线
    for i in range(N - 1):
        for j in range(N):
            x1, y1 = 2 * i, 2 * j
            x2, y2 = 2 * (i + 1), 2 * j
            ax.plot([x1, x2], [y1, y2], color='black', linewidth=1, zorder=1)

    # 添加 B-B 纵向黑线
    for i in range(N):
        for j in range(N - 1):
            x1, y1 = 2 * i + 1, 2 * j + 1
            x2, y2 = 2 * i + 1, 2 * (j + 1) + 1
            ax.plot([x1, x2], [y1, y2], color='black', linewidth=1, zorder=1)

    # B 到周围四个 A 的连线
    for i in range(N):
        for j in range(N):
            x_b, y_b = 2 * i + 1, 2 * j + 1
            neighbors = [
                (x_b - 1, y_b + 1),  # 左上
                (x_b - 1, y_b - 1),  # 左下
                (x_b + 1, y_b + 1),  # 右上
                (x_b + 1, y_b - 1)   # 右下
            ]
            for x_a, y_a in neighbors:
                if 0 <= x_a < 2 * N and 0 <= y_a < 2 * N:
                    ax.plot([x_b, x_a], [y_b, y_a], color='black', linewidth=1, zorder=1)

    ax.set_xlim(-1, 2 * N + 1)
    ax.set_ylim(-1, 2 * N + 1)
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_aspect('equal')
    ax.grid(False)
    plt.title("ground state",fontsize=15)

    path = "./image/show_ground_state_fancyarrow_3.png"
    plt.savefig(path)
    plt.close()
def find_neibor(i:int, j:int, type_gedian:int, N:int) -> list:
    """
    找到格点的邻居（使用整除处理周期性边界条件）
    """
    neighbor = []

    if type_gedian == 0:
        # A 格点邻居：AA 横向 + AB 右上、右下、下方、左下三个
        neighbor.append((0, (i - 1) % N, j))      # A 左
        neighbor.append((1, (i - 1) % N, j))      # B 左下
        neighbor.append((0, (i + 1) % N, j))      # A 右
        neighbor.append((1, i, j))                # B 正下
        neighbor.append((1, i, (j - 1) % N))      # B 左
        neighbor.append((1, (i - 1) % N, (j - 1) % N))  # B 左下角

    elif type_gedian == 1:
        # B 格点邻居：BB 上下 + BA 左、右、右上、上
        neighbor.append((0, i, j))                # A 中
        neighbor.append((1, i, (j - 1) % N))      # B 下
        neighbor.append((1, i, (j + 1) % N))      # B 上
        neighbor.append((0, i, (j + 1) % N))      # A 上
        neighbor.append((0, (i + 1) % N, j))      # A 右
        neighbor.append((0, (i + 1) % N, (j + 1) % N))  # A 右上

    return neighbor


def energy_final(system:np.ndarray,N:int)->float:
    energy_total = 0
    for i in range(N):
        for j in range(N):
            e1 = calculate_energy_in_neighbors(system,i,j,0,N)
            # print(f"for A in {i},{j}, energy is {e1}")
            energy_total+=e1
            e2 = calculate_energy_in_neighbors(system,i,j,1,N)
            energy_total+=e2
            # print(f"for B in {i},{j}, energy is {e2}")
    return energy_total/2


    

if __name__ == "__main__":
    N = 8
    J = 1.0
    steps = 10000
    spins = ising_monte_carlo(N, steps)
    visualize_spin(spins,N)
    print("energy of ground state:",energy_final(spins,N))
# if __name__ == "__main__":
#     N = 8
#     spins = np.array([
#         [[1, 1, 1, 1, 1, 1, 1, 1], [1, 1, 1, 1, 1, 1, 1, 1], [1, 1, 1, 1, 1, 1, 1, 1], [1, 1, 1, 1, 1, 1, 1, 1],
#          [1, 1, 1, 1, 1, 1, 1, 1], [1, 1, 1, 1, 1, 1, 1, 1], [1, 1, 1, 1, 1, 1, 1, 1], [1, 1, 1, 1, 1, 1, 1, 1]],
#         [[-1, -1, -1, -1, -1, -1, -1, -1], [-1, -1, -1, -1, -1, -1, -1, -1], [-1, -1, -1, -1, -1, -1, -1, -1],
#          [-1, -1, -1, -1, -1, -1, -1, -1], [-1, -1, -1, -1, -1, -1, -1, -1], [-1, -1, -1, -1, -1, -1, -1, -1],
#          [-1, -1, -1, -1, -1, -1, -1, -1], [-1, -1, -1, -1, -1, -1, -1, -1]]
#     ])
#     print("energy of ground state:", energy_final(spins, N))
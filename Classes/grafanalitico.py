from PDES import PDES
import PDE
from Disc_tokenfix import df
import SERKF45
import SERKF45_corrigido as RK
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.animation import FuncAnimation
import time

disc_n = 11

resultado_analitico1 = []
for i in range(disc_n):
    for j in range(disc_n):
        x_ = i/(disc_n-1)
        y_ = j/(disc_n-1)
        t = 1
        F_analitico = t + x_ + y_
        # G_analitico = t + np.cosh(x_ - y_)
        resultado_analitico1.append(F_analitico)
        # resultado_analitico.append(G_analitico)


print("Resultado Analítico:")
print(resultado_analitico1)


resultado_analitico2 = []
for i in range(disc_n):
    for j in range(disc_n):
        x_ = i / (disc_n - 1)
        y_ = j / (disc_n - 1)
        t = 1
        F_analitico = float((t+1)*np.sin(x_) + (t+2)*np.cos(y_))
        # G_analitico = t + np.cosh(x_ - y_)
        resultado_analitico2.append(F_analitico)
        # resultado_analitico.append(G_analitico)

print("Resultado Analítico:")
print(resultado_analitico2)

disc_n = 21

tfinal = 1.0
resultado_analitico3 = []
for i in range(disc_n):
    for j in range(disc_n):
        x_ = i / (disc_n - 1)
        y_ = j / (disc_n - 1)
        F_analitico = float(np.tanh(tfinal) * (np.sin(np.pi * x_) + np.cos(np.pi * y_)))
        resultado_analitico3.append(F_analitico)

print("Resultado Analítico:")
print(resultado_analitico3)

def heatmap_generate(data, disc_n, title):
    data_reshaped = np.array(data).reshape((disc_n, disc_n), order='F')
    plt.imshow(data_reshaped, cmap='RdYlBu_r', interpolation='bilinear', origin='lower', extent=(0, 1, 0, 1))
    plt.colorbar(label='Value')
    plt.title(title)
    plt.xlabel('x')
    plt.ylabel('y')
    plt.show()
    
heatmap_generate(resultado_analitico1, 11, "t = 1")
heatmap_generate(resultado_analitico2, 11, "t = 1")
heatmap_generate(resultado_analitico3, 21, "t = 1")

def a3d_animation(data, disc_n, title):
    data_reshaped = np.array(data).reshape((disc_n, disc_n), order='F')
    x = np.linspace(0, 1, disc_n)
    y = np.linspace(0, 1, disc_n)
    X, Y = np.meshgrid(x, y)
    
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    
    def update(frame):
        ax.clear()
        ax.plot_surface(X, Y, data_reshaped, cmap='RdYlBu_r', edgecolor='none')
        ax.set_zlim(np.min(data_reshaped), np.max(data_reshaped))
        ax.set_title(f"{title}")
        ax.set_xlabel('X')
        ax.set_ylabel('Y')
        ax.set_zlabel('Valor')
    
    ani = FuncAnimation(fig, update, frames=1, repeat=False)
    plt.show()

a3d_animation(resultado_analitico1, 11, "t = 1")
a3d_animation(resultado_analitico2, 11, "t = 1")
a3d_animation(resultado_analitico3, 21, "t = 1")
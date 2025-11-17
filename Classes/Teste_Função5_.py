import time
from PDES import PDES 
import PDE
from Disc_tokenfix import df
import SERKF45_corrigido as RK
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.animation import FuncAnimation


#! Valores Iniciais

disc_n=21

PDE1 = PDE.PDE('dT/dt = -dT/dx - dT/dy + 10*sech(t)**2 * x *(sin(x)+cos(y))+10*tanh(t)*(sin(x)+cos(y)+x*cos(x)-x *sin(y))',
                ['T'], ['x','y'],[disc_n,disc_n], [(0,1),(0,1)], '0')



tfinal = 1
resultado_analitico = []
for i in range(disc_n):
    for j in range(disc_n):
        x_ = i/(disc_n-1)
        y_ = j/(disc_n-1)
        F_analitico = float(10*x_*np.tanh(tfinal)*(np.sin(x_)+np.cos(y_)))
        resultado_analitico.append(F_analitico)
        # resultado_analitico.append(G_analitico)


print("Resultado Analítico:")
print(resultado_analitico)

PDES1 = PDES([PDE1], ['x','y'], ['T'])

start_time_df = time.time()
resultado = df(
    PDES1, [disc_n, disc_n],
    west_bd="Dirichlet",
    method="backward",
    north_bd="Dirichlet", south_bd="Dirichlet", east_bd="Dirichlet",  # west via inlet
    north_func_bd='10*x*tanh(t)*(sin(x)+cos(1))',
    south_func_bd = '10*x*tanh(t)*(sin(x)+1)',
    east_func_bd = '10*tanh(t)*(sin(1)+cos(y))',
    west_func_bd = '0'
)
exec_time_df = time.time() - start_time_df

start_time_rk_cuda = time.time()

testar = RK.SERKF45_cuda(resultado[0], ['t'], resultado[1], PDES1.ic, 0, tfinal, 1000, 1, len(PDES1.sp_vars))

exec_time_rk_cuda = time.time() - start_time_rk_cuda

exec_time_total = time.time() - start_time_df



print("Resultado Numérico:")
print(testar[1][0][-1])

print("Erro Absoluto:")
erro_absoluto = np.abs(np.array(testar[1][0][-1]) - np.array(resultado_analitico))
print(erro_absoluto)

print("Erro Relativo (%):")
erro_relativo = (erro_absoluto / np.where(np.array(resultado_analitico) != 0, np.array(resultado_analitico), 1e-5)) * 100
print(erro_relativo)

print("Erro Médio Absoluto:")
erro_medio_absoluto = np.mean(erro_absoluto)
print(erro_medio_absoluto)

print("Erro Médio Relativo:")
erro_medio_relativo = np.mean(erro_relativo)
print(erro_medio_relativo)

print(f"Tempo de execução df: {exec_time_df:.2f} segundos")
print(f"Tempo de execução RK CUDA: {exec_time_rk_cuda:.2f} segundos")
print(f"Tempo de execução total: {exec_time_total:.2f} segundos")


x = np.linspace(0,1,disc_n)
y = np.linspace(0,1,disc_n)
X, Y = np.meshgrid(x, y, indexing='xy')


if len(PDES1.sp_vars) == 2:
    
    vetor = np.array(testar[1][0][-1], dtype=float)  # Convertendo para um array NumPy
    vetor = vetor.reshape((disc_n, disc_n), order='F')  # Reshape para matriz 2D (coluna maior que linha)
    plt.imshow(vetor, cmap='RdYlBu_r', interpolation='bilinear', extent=(0, 1, 0, 1), origin='lower')
    plt.show()
    
    lista_vetores = testar[1][0]  
    n_frames = len(lista_vetores)

    # 1) Empilha e reshape para (n_frames, disc_n, disc_n)
    data3d = np.array(lista_vetores, dtype=float) \
             .reshape((n_frames, disc_n, disc_n), order='F')
    
    # 2) Plot 3D
    # cria grelha de coordenadas X, Y no domínio [0,1]×[0,1]
    x = np.linspace(0, 1, disc_n)
    y = np.linspace(0, 1, disc_n)
    X, Y = np.meshgrid(x, y)

    # 2) Setup inicial da figura 3D
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.set_xlim(0,1)
    ax.set_ylim(0,1)
    ax.set_zlim(np.min(data3d), np.max(data3d))
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Valor')

    # plota a primeira fatia
    surf = ax.plot_surface(X, Y, data3d[0], cmap='RdYlBu_r', edgecolor='none')

    # 3) Função de update
    def update(frame):
        global surf
        # remove superfície antiga
        surf.remove()
        # desenha nova
        surf = ax.plot_surface(
            X, Y, data3d[frame],
            cmap='RdYlBu_r', edgecolor='none'
        )
        ax.set_title(f"Frame {frame+1}/{n_frames}")

    # 4) Cria animação
    ani = FuncAnimation(
        fig, update,
        frames=n_frames,
        interval=10,      # ms entre frames
        repeat=False,
        blit=False         # blit não funciona bem em 3D
    )

    plt.show()
    from matplotlib.animation import FFMpegWriter

    # configura o writer:
    writer = FFMpegWriter(
        fps=10,                # quadros por segundo
        metadata=dict(artist='Você'),
        bitrate=1800           # taxa de bits (quanto maior, melhor qualidade/tamanho)
    )




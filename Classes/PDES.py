import PDE
from Disc import df
import SERKF45
import RKF45_novo as RK
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from mpl_toolkits.mplot3d import Axes3D  # habilita o 3D em matplotlib
from matplotlib.animation import FuncAnimation


class PDES:
    eqs = []
    ivars = []
    sp_vars = []
    funcs = []
    ic = []
    n_sp = 1
    n_temp = 1
    def __init__(self, pdes, sp_vars, n_sp=1, n_temp=1):
        self.eqs = [pde.eq for pde in pdes]
        self.ivars = list([pde.ivar[i] for pde in pdes for i in range(len(pde.ivar))])
        self.funcs = list([pde.func[i] for pde in pdes for i in range(len(pde.func))])
        self.ic = list([pde.ic[i] for pde in pdes for i in range(len(pde.ic))])
        self.sp_vars = sp_vars
        self.n_sp = n_sp
        self.n_temp = n_temp
        
        pass
    
    
    def xs(self,vars):
        nvars = vars.copy()
        for i in range(len(nvars)): nvars[i] = f'XX{i}'
        return nvars


#! Valores Iniciais

# PDE2 = PDE.PDE('dT/dt = 0.3*d2T/dx2 + 0.5*d2T/dy2', ['t', 'x', 'y'], ['T'])
# PDE1 = PDE.PDE('dC/dt = 0.1*d2C/dx2', ['t', 'x'], ['C'])
# PDE3 = PDE.PDE('dD/dt = 0.1*d2D/dx2', ['t', 'x'], ['D'])
# **2 * x_ *(sin(x_)+cos(y_))+10*tanh(t)*(sin(x_)+cos(y_)+x_*cos(x_)-x_ *sin(y_))

disc_n=5

PDE1 = PDE.PDE('dF/dt = sinh(x+y) + d2F/dx2 + d2F/dy2 - 2*F - (d2G/dx2 + d2G/dy2) + 2*(G - t)',  
                ['F'], ['x','y'],[disc_n,disc_n], [(0,1),(0,1)], '0')
PDE2 = PDE.PDE('dG/dt = d2G/dx2 + d2G/dy2 - 2*(G - t) - (d2F/dx2 + d2F/dy2 - 2*F) + 1 ',  
                ['G'], ['x','y'],[disc_n,disc_n], [(0,1),(0,1)], 'cosh(x-y)')

resultado_analitico = []
for i in range(disc_n):
    for j in range(disc_n):
        x_ = i/(disc_n-1)
        y_ = j/(disc_n-1)
        t = 1
        F_analitico = t * np.sinh(x_ + y_)
        G_analitico = t + np.cosh(x_ - y_)
        resultado_analitico.append(F_analitico)
        resultado_analitico.append(G_analitico)

print("Resultado Analítico:")
print(resultado_analitico)

PDES1 = PDES([PDE1,PDE2], ['x','y'], ['F','G'])


print(PDES1.ic)
resultado = df(PDES1, [disc_n,disc_n], inlet=[[f't * sinh({i/(disc_n-1)})' for i in range(disc_n)],[f't + cosh({-i/(disc_n-1)})' for i in range(disc_n)]], method="backward") #! erro forward




# print(SERKF45.SERKF45(resultado[0], ['t'], resultado[1], initial_values, 0, 0.1, 10,2,len(PDES1.sp_vars)))
testar = RK.SERKF45_cuda(resultado[0], ['t'], resultado[1], PDES1.ic, 0, 1, 10, 1, len(PDES1.sp_vars))
#testar = SERKF45.SERKF45(resultado[0], ['t'], resultado[1], PDES1.ic, 0, 1, 10, 1, len(PDES1.sp_vars))
print(testar[1][0][-1])

df = pd.DataFrame(testar[1][0][-1], columns=["Valores"])

df.to_excel("saida.xlsx", index=False)

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
        interval=200,      # ms entre frames
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

    # salva a animação:
    ani.save('animacao3d.mp4', writer=writer)

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
import CN as CN

start_time_df = time.time()
# ! Valores Iniciais

disc_n = 100

PDE1 = PDE.PDE('dF/dt =  0.1*d2F/dx2',
               ['F'], ['x'], [disc_n], [(0, 1)], 'sin(pi * x)')

resultado_analitico = []
# for i in range(disc_n):
#     for j in range(disc_n):
#         x_ = i / (disc_n - 1)
#         y_ = j / (disc_n - 1)
#         t = 1
#         F_analitico = float(x_ + np.exp(t * x_))  # t=1
#         # G_analitico = t + np.cosh(x_ - y_)
#         resultado_analitico.append(F_analitico)
#         # resultado_analitico.append(G_analitico)

print("Resultado Analítico:")
print(resultado_analitico)

PDES1 = PDES([PDE1], ['x'], ['F'])

start_time_df = time.time()

resultado = df(
    PDES1, [disc_n],
    west_func_bd="0",  # x=0
    west_bd="Dirichlet",
    method="central",
    north_bd="Dirichlet", south_bd="Dirichlet", east_bd="Dirichlet",  
    # north_func_bd='(t+1)*sin(x)+(t+2)*cos(1)',   # y=1
    # south_func_bd='(t+1)*sin(x)+(t+2)*cos(0)',   # y=0          
    east_func_bd='0'     # x=1
)

exec_time_df = time.time() - start_time_df

exec_time_total = time.time() - start_time_df


print("Resultado Discretizado:")
print(resultado[0])

N = 50           # Malha 50x50
L = 1.0          # Dimensão da placa
dx = dy = L/(N-1)

start_time_cn = time.time()

resultado_final_cn = CN.cn_1d(
    resultado[0],
    resultado[1],
    nt=1000,
    dt=0.001,
    u_init_val=1
)

exec_time_cn = time.time() - start_time_cn
start_time_rk_cuda = time.time()

resultado_final_rk = RK.SERKF45_cuda(
    oldexpr=resultado[0],          # Lista de strings discretizadas
    ivar=['t'],                    # Variável independente
    funcs=resultado[1],            # LISTA DE VARIÁVEIS (XX0_0_0, etc.)
    yn=PDES1.ic,                   # Condição Inicial (Vetor de 100 posições)
    x0=0,                          # Tempo inicial
    xn=1,                        # Tempo final (IGUAL ao tempo simulado no CN)
    n=100,                         # Sugestão de divisões (o RK adaptativo ignora, mas use um valor base)
    n_funcs=1,                     # Número de PDEs (no seu caso, apenas a 'F')
    sp_vars=PDES1.sp_vars,         # Variáveis espaciais
    dt_init=0.001,                # PARÂMETRO CHAVE: Comece com um passo bem pequeno
    dt_max=0.1,                  # PARÂMETRO CHAVE: Limite o passo para não estourar a estabilidade
    tol=1e-6                       # Tolerância: 1e-6 é um ótimo equilíbrio
)
exec_time_rk_cuda = time.time() - start_time_rk_cuda

print("Resultado Final da Simulação cn 1D:")
print(resultado_final_cn)


print("Resultado Final da Simulação RK 1d:")
print(resultado_final_rk[1][0][-1])


exec_time_total = time.time() - start_time_df

print(f"Tempo de execução df: {exec_time_df:.2f} segundos")
print(f"Tempo de execução RK CUDA: {exec_time_rk_cuda:.2f} segundos")
print(f"Tempo de execução CN: {exec_time_cn:.2f} segundos")
print(f"Tempo de execução total: {exec_time_total:.2f} segundos")

x = np.linspace(0, 1, 100)

# 2. Cálculo da Solução Analítica (Série de Fourier)
# Ajustamos o tempo (t) para coincidir com o perfil do CN (aprox t=1.0 se alpha=0.1)
def analytical(x, t, alpha=0.1):
    u = np.zeros_like(x)
    for n in range(1, 50, 2):
        u += (4/(n*np.pi)) * np.sin(n*np.pi*x) * np.exp(-alpha*(n*np.pi)**2 * t)
    return u

res_analitico = analytical(x, t=1.0) # Alvo do CN

# # 3. Configuração do Gráfico 3D
# fig = plt.figure(figsize=(12, 8))
# ax = fig.add_subplot(111, projection='3d')

# # Plotando as linhas
# ax.plot(x, np.zeros_like(x), resultado_final_cn, label='Crank-Nicolson (Numérico)', color='blue', linewidth=2)
# ax.plot(x, np.ones_like(x)*0.1, resultado_final_rk[1][0][-1], label='Runge-Kutta (Numérico)', color='green', linewidth=2)
# ax.plot(x, np.zeros_like(x), res_analitico, '--', label='Analítico (Alvo CN)', color='red', alpha=0.7)

# # Adicionando superfícies translúcidas para preencher o volume (opcional, melhora a estética)
# ax.add_collection3d(plt.fill_between(x, 0, resultado_final_cn, color='blue', alpha=0.1), zs=0, zdir='y')
# ax.add_collection3d(plt.fill_between(x, 0, resultado_final_rk[1][0][-1], color='green', alpha=0.1), zs=0.1, zdir='y')

# # Customização
# ax.set_xlabel('Posição (x)')
# ax.set_ylabel('Eixo de Comparação')
# ax.set_zlabel('Temperatura (u)')
# ax.set_title('Validação da Biblioteca: EDP de Calor 1D')
# ax.legend()

# # Ajuste do ângulo de visão para melhor profundidade
# ax.view_init(elev=25, azim=-45)

# plt.show()

from sklearn.metrics import mean_absolute_error, mean_squared_error, r2_score

# --- 1. Sincronização de Dados ---
u_cn = np.array(resultado_final_cn)
u_rk = np.array(resultado_final_rk[1][0][-1])
u_ana = np.array(res_analitico)

# --- 2. Função para Erro Relativo Médio (%) ---
def mean_relative_error(y_true, y_pred):
    # Evitamos divisão por zero nas bordas (onde u=0) usando um pequeno epsilon
    mask = y_true > 1e-5 
    return np.mean(np.abs((y_true[mask] - y_pred[mask]) / y_true[mask])) * 100

# --- 3. Cálculo das Métricas ---

# Métricas Crank-Nicolson
mae_cn = mean_absolute_error(u_ana, u_cn)
rmse_cn = np.sqrt(mean_squared_error(u_ana, u_cn))
r2_cn = r2_score(u_ana, u_cn)
mre_cn = mean_relative_error(u_ana, u_cn)

# Métricas Runge-Kutta
mae_rk = mean_absolute_error(u_ana, u_rk)
rmse_rk = np.sqrt(mean_squared_error(u_ana, u_rk))
r2_rk = r2_score(u_ana, u_rk)
mre_rk = mean_relative_error(u_ana, u_rk)

print("\n" + "="*50)
print(f"{'EQUAÇÃO DE CALOR 1D':^50}")
print("="*50)
print(f"{'Métrica':<15} | {'Crank-Nicolson':<15} | {'Runge-Kutta':<15}")
print("-"*50)
print(f"{'R² (Score)':<15} | {r2_cn:<15.6f} | {r2_rk:<15.6f}")
print(f"{'Erro Relat. %':<15} | {mre_cn:<15.2f}% | {mre_rk:<15.2f}%")
print(f"{'MAE':<15} | {mae_cn:<15.6f} | {mae_rk:<15.6f}")
print(f"{'RMSE':<15} | {rmse_cn:<15.6f} | {rmse_rk:<15.6f}")
print("="*50)

# --- 4. Plotagem do Gráfico 2D ---

# Criando a figura com tamanho maior e DPI alto
# figsize=(8, 10) garante que o gráfico seja "alto" e preencha a página
plt.rcParams.update({'font.size': 12, 'font.family': 'serif'})
plt.figure(figsize=(10, 6), dpi=210)
plt.plot(x, u_ana, 'r-', linewidth=2.5, label=f'Analítico (Referência)')
plt.plot(x, u_cn, 'b--', linewidth=2, label=f'CN (R²: {r2_cn:.4f})')
plt.plot(x, u_rk, 'g:', linewidth=2, label=f'RK (R²: {r2_rk:.4f})')

plt.xlabel('Posição (x)')
plt.ylabel('u')
plt.legend()
plt.grid(True, alpha=0.3)
plt.tight_layout()
plt.show()

# exec_time_df = time.time() - start_time_df

# print("oi")

# start_time_rk_cuda = time.time()

# testar = RK.SERKF45_cuda(resultado[0], ['t'], resultado[1], PDES1.ic, 0, 1, 100, 1, len(PDES1.sp_vars))
# print("Resultado Numérico RK CUDA:")
# print(testar[1][0][-1])
# exec_time_rk_cuda = time.time() - start_time_rk_cuda

# print("oi")

# exec_time_total = time.time() - start_time_df

# start_time_rk = time.time()
# # testar = SERKF45.SERKF45(resultado[0], ['t'], resultado[1], PDES1.ic, 0, 0.1, 400, 1, len(PDES1.sp_vars))
# exec_time_rk = time.time() - start_time_rk


# print(f"Tempo de execução df: {exec_time_df:.2f} segundos")
# print(f"Tempo de execução RK: {exec_time_rk:.2f} segundos")
# print(f"Tempo de execução RK CUDA: {exec_time_rk_cuda:.2f} segundos")
# print(f"Tempo de execução total: {exec_time_total:.2f} segundos")


# print("Resultado Numérico:")
# print(testar[1][0][-1])

# print("Erro Absoluto:")
# erro_absoluto = np.abs(np.array(testar[1][0][-1]) - np.array(resultado_analitico))
# print(erro_absoluto)

# print("Erro Relativo (%):")
# erro_relativo = (erro_absoluto / np.array(resultado_analitico)) * 100
# print(erro_relativo)

# print("Erro Médio Absoluto:")
# erro_medio_absoluto = np.mean(erro_absoluto)
# print(erro_medio_absoluto)

# print("Erro Médio Relativo:")
# erro_medio_relativo = np.mean(erro_relativo)
# print(erro_medio_relativo)

# df = pd.DataFrame(testar[1][0][-1], columns=["Valores"])

# x = np.linspace(0, 1, disc_n)
# y = np.linspace(0, 1, disc_n)
# X, Y = np.meshgrid(x, y, indexing='xy')  # X[y,x], Y[y,x]

# F_anal_mat = 2*np.sin(X) + 3*np.cos(Y)  # t=1
# # Se o seu empacotamento linear usa ordem 'F':
# F_anal_vec = np.array(F_anal_mat, dtype=float).reshape(-1, order='F')

# F_num_vec = np.array(testar[1][0][-1], dtype=float)
# # Erros coerentes
# erro_abs = np.abs(F_num_vec - F_anal_vec)
# erro_rel = erro_abs / np.maximum(np.abs(F_anal_vec), 1e-15) * 100
# print("||erro||_inf:", erro_abs.max())
# print("MAE:", erro_abs.mean())
# print("R^2", 1 - np.sum(erro_abs**2) / np.sum((F_anal_vec - np.mean(F_anal_vec))**2))

# df = pd.DataFrame(testar[1][0][-1], columns=["Valores"])

# df.to_excel("saida.xlsx", index=False)

# if len(PDES1.sp_vars) == 2:
#     vetor = np.array(testar[1][0][-1], dtype=float)  # Convertendo para um array NumPy
#     vetor = vetor.reshape((disc_n, disc_n), order='F')  # Reshape para matriz 2D (coluna maior que linha)
#     plt.imshow(vetor, cmap='RdYlBu_r', interpolation='bilinear', extent=(0, 1, 0, 1), origin='lower')
#     plt.show()

#     lista_vetores = testar[1][0]
#     n_frames = len(lista_vetores)

#     # 1) Empilha e reshape para (n_frames, disc_n, disc_n)
#     data3d = np.array(lista_vetores, dtype=float) \
#         .reshape((n_frames, disc_n, disc_n), order='F')

#     # 2) Plot 3D
#     # cria grelha de coordenadas X, Y no domínio [0,1]×[0,1]
#     x = np.linspace(0, 1, disc_n)
#     y = np.linspace(0, 1, disc_n)
#     X, Y = np.meshgrid(x, y)

#     # 2) Setup inicial da figura 3D
#     fig = plt.figure()
#     ax = fig.add_subplot(111, projection='3d')
#     ax.set_xlim(0, 1)
#     ax.set_ylim(0, 1)
#     ax.set_zlim(np.min(data3d), np.max(data3d))
#     ax.set_xlabel('X')
#     ax.set_ylabel('Y')
#     ax.set_zlabel('Valor')

#     # plota a primeira fatia
#     surf = ax.plot_surface(X, Y, data3d[0], cmap='RdYlBu_r', edgecolor='none')


#     # 3) Função de update
#     def update(frame):
#         global surf
#         # remove superfície antiga
#         surf.remove()
#         # desenha nova
#         surf = ax.plot_surface(
#             X, Y, data3d[frame],
#             cmap='RdYlBu_r', edgecolor='none'
#         )
#         ax.set_title(f"Frame {frame + 1}/{n_frames}")


#     # 4) Cria animação
#     ani = FuncAnimation(
#         fig, update,
#         frames=n_frames,
#         interval=20,  # ms entre frames
#         repeat=False,
#         blit=False  # blit não funciona bem em 3D
#     )

#     plt.show()
#     from matplotlib.animation import FFMpegWriter

#     # configura o writer:
#     writer = FFMpegWriter(
#         fps=10,  # quadros por segundo
#         metadata=dict(artist='Você'),
#         bitrate=1800  # taxa de bits (quanto maior, melhor qualidade/tamanho)
#     )


from PDES import PDES
import PDE
from Disc_tokenfix import df
import SERKF45_corrigido as RK
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.animation import FuncAnimation
import time

# ============================================================
# Funções auxiliares de erro (robustas em torno de zero)
# ============================================================

def compute_scale(y_true, kind="rms", eps=1e-12):
    """
    Define uma escala S do problema para normalização de erros.
    kind="rms": S = ||y_true||2 / sqrt(N)  (amplitude RMS)
    kind="max": S = max |y_true|
    kind="range": S = max(y_true) - min(y_true)
    Retorna 1.0 se a escala for muito pequena, garantindo robustez perto de 0.
    """
    y_true = np.asarray(y_true, dtype=float)
    N = y_true.size if y_true.size > 0 else 1
    if kind == "rms":
        S = np.linalg.norm(y_true) / np.sqrt(N)
    elif kind == "max":
        S = np.max(np.abs(y_true)) if y_true.size else 0.0
    elif kind == "range":
        S = (np.max(y_true) - np.min(y_true)) if y_true.size else 0.0
    else:
        raise ValueError("kind deve ser 'rms', 'max' ou 'range'")
    return S if S > eps else 1.0


def symmetric_relative_error(y_hat, y_true, eps=1e-12):
    """
    Erro relativo simétrico estável em 0:
    |y_hat - y_true| / (|y_hat| + |y_true| + eps)
    """
    y_hat = np.asarray(y_hat, dtype=float)
    y_true = np.asarray(y_true, dtype=float)
    return np.abs(y_hat - y_true) / (np.abs(y_hat) + np.abs(y_true) + eps)


def norm_errors(err_vec, ref_vec=None, eps=1e-12):
    """
    Normas globais do erro com normalização robusta:
    - L2_abs  = ||e||2
    - L2_rel  = ||e||2 / (||ref||2 + eps)
    - Linf_abs= ||e||∞
    - Linf_rel= ||e||∞ / (||ref||∞ + eps)
    """
    e = np.asarray(err_vec, dtype=float)
    L2_abs = np.linalg.norm(e)
    Linf_abs = np.max(np.abs(e)) if e.size else 0.0
    if ref_vec is not None:
        r = np.asarray(ref_vec, dtype=float)
        L2_ref = np.linalg.norm(r)
        Linf_ref = np.max(np.abs(r)) if r.size else 0.0
        L2_rel = L2_abs / (L2_ref + eps)
        Linf_rel = Linf_abs / (Linf_ref + eps)
    else:
        L2_rel, Linf_rel = np.nan, np.nan
    return dict(L2_abs=L2_abs, L2_rel=L2_rel, Linf_abs=Linf_abs, Linf_rel=Linf_rel)


# ============================================================
# Valores Iniciais e configuração do problema
# ============================================================

disc_n = 5
nu = 0.1

PDE1 = PDE.PDE(
    f'dT/dt = -dT/dx - dT/dy + {nu}*(d2T/dx2 + d2T/dy2) + '
    'sech(t)**2*(sin(pi*x)+cos(pi*y)) + '
    'tanh(t)*(pi*cos(pi*x) - pi*sin(pi*y) + {nu}*pi**2*(sin(pi*x)+cos(pi*y)))'.format(nu=nu),
    ['T'], ['x', 'y'], [disc_n, disc_n], [(0, 1), (0, 1)], '0'
)

tfinal = 1.0
resultado_analitico = []
for i in range(disc_n):
    for j in range(disc_n):
        x_ = i / (disc_n - 1)
        y_ = j / (disc_n - 1)
        F_analitico = float(np.tanh(tfinal) * (np.sin(np.pi * x_) + np.cos(np.pi * y_)))
        resultado_analitico.append(F_analitico)

print("Resultado Analítico:")
# print(resultado_analitico)

PDES1 = PDES([PDE1], ['x', 'y'], ['T'])
time_start_df = time.time()
resultado = df(
    PDES1, [disc_n, disc_n],
    west_bd="Dirichlet", east_bd="Dirichlet",
    south_bd="Dirichlet", north_bd="Dirichlet",
    west_func_bd='tanh(t)*(cos(pi*y))',  # x=0
    east_func_bd='tanh(t)*(cos(pi*y))',  # x=1
    south_func_bd='tanh(t)*(sin(pi*x)+1)',  # y=0
    north_func_bd='tanh(t)*(sin(pi*x)-1)',  # y=1
    method="backward"  # upwind (u=v=+1)
)
print(resultado[0])
time_exec_df = time.time() - time_start_df
time_start_rk = time.time()
testar = RK.SERKF45_cuda(
    resultado[0], ['t'], resultado[1], PDES1.ic,
    0.0, tfinal, 1000, 1, len(PDES1.sp_vars)
)
time_exec_rk = time.time() - time_start_rk
# ============================================================
# Cálculo das métricas de erro (robustas quando x* ~ 0)
# ============================================================

print("Resultado Numérico:")
y_num = np.array(testar[1][0][-1], dtype=float)
print(y_num)

y_true = np.array(resultado_analitico, dtype=float)

# 1) Erro absoluto ponto a ponto
erro_abs = np.abs(y_num - y_true)

# 2) Erro relativo clássico (com proteção para zero)
eps_rel = 1e-12
erro_rel = erro_abs / (np.abs(y_true) + eps_rel)

# 3) Erro relativo simétrico (estável em 0)
erro_rel_sim = symmetric_relative_error(y_num, y_true, eps=eps_rel)

# 4) Erro escalado por uma escala S do problema
S = compute_scale(y_true, kind="rms", eps=1e-12)
erro_escalado = erro_abs / S

# 5) Normas globais
normas = norm_errors(erro_abs, ref_vec=y_true, eps=1e-12)

# 6) R² (coeficiente de determinação)
ss_res = np.sum((y_true - y_num)**2)
ss_tot = np.sum((y_true - np.mean(y_true))**2)
R2 = 1 - ss_res/(ss_tot + 1e-15)

# 7) Sumários
MAE = np.mean(erro_abs)
MedAE = np.median(erro_abs)
mean_rel = np.mean(erro_rel)
mean_scaled = np.mean(erro_escalado)
mean_sym_rel = np.mean(erro_rel_sim)

print("\n--- Métricas de erro robustas ---")
print(f"Escala S (kind='rms'): {S:.6e}")
print(f"L2_abs  = {normas['L2_abs']:.6e} ; L2_rel  = {normas['L2_rel']:.6e}")
print(f"Linf_abs= {normas['Linf_abs']:.6e} ; Linf_rel= {normas['Linf_rel']:.6e}")
print(f"R² = {R2:.6f}")

print("\nResumo ponto a ponto:")
print(f"MAE (média |erro|)       = {MAE:.6e}")
print(f"Mediana |erro|            = {MedAE:.6e}")
print(f"Média erro relativo       = {mean_rel:.6e}")
print(f"Média erro escalado       = {mean_scaled:.6e}")
print(f"Média erro rel. simétrico = {mean_sym_rel:.6e}")
print(f"Tempo de execução df: {time_exec_df:.2f} segundos")
print(f"Tempo de execução RK: {time_exec_rk:.2f} segundos")
print(f"Tempo de execução total: {time_exec_df + time_exec_rk:.2f} segundos")

# (Opcional) DataFrame para inspeção
df_erros = pd.DataFrame({
    "true": y_true,
    "num": y_num,
    "abs_err": erro_abs,
    "rel_err": erro_rel,
    "scaled_err": erro_escalado,
    "sym_rel_err": erro_rel_sim
})
print("\nPrimeiras linhas do DataFrame de erros:")
print(df_erros.head())


# ============================================================
# Visualização
# ============================================================

x = np.linspace(0, 1, disc_n)
y = np.linspace(0, 1, disc_n)
X, Y = np.meshgrid(x, y, indexing='xy')

if len(PDES1.sp_vars) == 2:
    vetor = np.array(testar[1][0][-1], dtype=float)  # vetor solução final
    vetor = vetor.reshape((disc_n, disc_n), order='F')
    plt.imshow(vetor, cmap='RdYlBu_r', interpolation='bilinear',
               extent=(0, 1, 0, 1), origin='lower')
    plt.title("Solução numérica (t = tfinal)")
    plt.colorbar()
    plt.show()

    lista_vetores = testar[1][0]
    n_frames = len(lista_vetores)

    # Empilha e reshape para (n_frames, disc_n, disc_n)
    data3d = np.array(lista_vetores, dtype=float) \
        .reshape((n_frames, disc_n, disc_n), order='F')

    # Plot 3D
    x = np.linspace(0, 1, disc_n)
    y = np.linspace(0, 1, disc_n)
    X, Y = np.meshgrid(x, y)

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.set_zlim(np.min(data3d), np.max(data3d))
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Valor')

    surf_holder = [ax.plot_surface(X, Y, data3d[0], cmap='RdYlBu_r', edgecolor='none')]


    def update(frame):
        # remove o surface anterior (a referência está em surf_holder[0])
        surf_holder[0].remove()
        # cria o novo surface e atualiza a referência
        surf_holder[0] = ax.plot_surface(X, Y, data3d[frame], cmap='RdYlBu_r', edgecolor='none')
        # opcional: atualizar limites se a escala variar muito
        # ax.set_zlim(np.min(data3d[frame]), np.max(data3d[frame]))
        return (surf_holder[0],)


    ani = FuncAnimation(fig, update, frames=n_frames, interval=20, repeat=False, blit=False)
    plt.show()

    from matplotlib.animation import FFMpegWriter

    writer = FFMpegWriter(fps=10, metadata=dict(artist='Você'), bitrate=1800)

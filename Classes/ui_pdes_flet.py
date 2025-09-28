# ui_pdes_flet.py
# ---------------------------------------------------------------
# Flet interface to configure, run and visualize your PDE setup.
# Place this file in the SAME folder where your modules live:
#   PDE.py, Disc.py, SERKF45.py, RKF45_novo.py, etc.
# Then run:  flet run ui_pdes_flet.py
# ---------------------------------------------------------------

import base64
import os
import threading
import time
import traceback

import flet as ft

# Matplotlib for file-based plots (headless)
import matplotlib
import io
matplotlib.use("Agg")  # use non-GUI backend
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, FFMpegWriter

import numpy as np
import pandas as pd

# --- Your project modules ---
# NOTE: We intentionally DO NOT import PDES.py because the user's file
#       appears to execute code at import-time. To avoid side effects,
#       we re-define the minimal PDES class here, mirroring the user's version.
#       (You can later refactor PDES.py to guard with `if __name__ == "__main__":`.)
import PDE
from Disc import df
import SERKF45
import RKF45_novo as RK


class PDES:
    """
    Minimal mirror of the user's PDES class (to avoid importing PDES.py which
    runs a script section on import). It aggregates equations (eqs), independent
    vars (ivars), functions (funcs), initial conditions (ic), and stores the
    spatial variables list (sp_vars) for downstream use.
    """
    eqs = []
    ivars = []
    sp_vars = []
    funcs = []
    ic = []
    n_sp = 1
    n_temp = 1

    def __init__(self, pdes, sp_vars, n_sp=1, n_temp=1):
        self.eqs = [pde.eq for pde in pdes]
        self.ivars = [pde.ivar[i] for pde in pdes for i in range(len(pde.ivar))]
        self.funcs = [pde.func[i] for pde in pdes for i in range(len(pde.func))]
        self.ic = [pde.ic[i] for pde in pdes for i in range(len(pde.ic))]
        self.sp_vars = sp_vars
        self.n_sp = n_sp
        self.n_temp = n_temp

    def xs(self, vars):
        nvars = vars.copy()
        for i in range(len(nvars)):
            nvars[i] = f"XX{i}"
        return nvars


# -------------------- Helper Parsers --------------------

def parse_list_csv(s: str):
    """Parse a comma/space separated list into a Python list of stripped strings."""
    if not s:
        return []
    items = [x.strip() for x in s.replace(";", ",").split(",")]
    return [x for x in items if x]


def parse_domain_pair(s: str):
    """Parse a pair like "0,1" -> (0.0, 1.0)."""
    parts = [p.strip() for p in s.split(",")]
    if len(parts) != 2:
        raise ValueError("Domínio deve ser no formato a,b (ex.: 0,1)")
    return float(parts[0]), float(parts[1])


# -------------------- Core Runner --------------------

def run_simulation(
    eq_text: str,
    funcs_list: list[str],
    sp_vars_list: list[str],
    disc_n: int,
    dom_x: tuple[float, float],
    dom_y: tuple[float, float],
    ic_text: str,
    method: str,
    t0: float,
    t1: float,
    nsteps: int,
    use_cuda: bool,
    save_excel: bool,
    save_heatmap: bool,
    save_mp4: bool,
    log: callable,
    image_path: str,
    north_bd: str, 
    south_bd: str, 
    east_bd: str, 
    north_func_bd: str, 
    south_func_bd: str, 
    east_func_bd: str, 
    north_alpha_bd: str, 
    south_alpha_bd: str, 
    east_alpha_bd: str, 
    north_beta_bd: str, 
    south_beta_bd: str, 
    east_beta_bd: str
):
    log("→ Iniciando simulação…")

    # Construct PDE and PDES aggregators (mirroring the user's code)
    disc = [disc_n, disc_n]
    domain = [dom_x, dom_y]

    log("• Construindo PDE…")
    eq_text = eq_text.split(",")  # single line
    eq_text1 = eq_text[0]  # single line
    eq_text2 = eq_text[1]  # single line
    pde1 = PDE.PDE(
        eq_text1,
        funcs_list,
        sp_vars_list,
        disc,
        domain,
        ic_text,
    )
    pde2 = PDE.PDE(
        eq_text2,
        funcs_list,
        sp_vars_list,
        disc,
        domain,
        ic_text,
    )

    log("• Agregando em PDES…")
    pdes = PDES([pde1,pde2], sp_vars_list)
    print(pdes.eqs)
    # Build discrete system via Disc.df
    log(f"• Discretizando com df(…): método={method} | disc={disc}")
    resultado = df(pdes, disc, inlet=[[0 for _ in range(disc_n)]], method=method, 
                    north_bd=north_bd, south_bd=south_bd, east_bd=east_bd, north_func_bd=north_func_bd, 
                    south_func_bd=south_func_bd, east_func_bd=east_func_bd, north_alpha_bd=north_alpha_bd, 
                    south_alpha_bd=south_alpha_bd, east_alpha_bd=east_alpha_bd, north_beta_bd=north_beta_bd, 
                    south_beta_bd=south_beta_bd, east_beta_bd=east_beta_bd)

    # Solve in time using provided kernels
    log("• Integrando no tempo…")
    solver_used = "RKF45_novo.SERKF45_cuda" if use_cuda else "SERKF45.SERKF45"
    try:
        if use_cuda:
            sol = RK.SERKF45_cuda(
                resultado[0], ['t'], resultado[1], pdes.ic,
                t0, t1, nsteps, 1, len(pdes.sp_vars)
            )
        else:
            sol = SERKF45.SERKF45(
                resultado[0], ['t'], resultado[1], pdes.ic,
                t0, t1, nsteps, 1, len(pdes.sp_vars)
            )
    except Exception as e:
        if use_cuda:
            log("! Falha no kernel CUDA; tentando fallback CPU (SERKF45)…")
            solver_used = "SERKF45.SERKF45 (fallback)"
            sol = SERKF45.SERKF45(
                resultado[0], ['t'], resultado[1], pdes.ic,
                t0, t1, nsteps, 1, len(pdes.sp_vars)
            )
        else:
            raise

    last_vec = sol[1][0][-1]
    log(f"• Último vetor (len={len(last_vec)}): amostra = {last_vec[:min(6, len(last_vec))]}")

    # Save Excel of final vector
    if save_excel:
        try:
            df_out = pd.DataFrame(last_vec, columns=["Valores"])  # 1D vector
            df_out.to_excel("saida.xlsx", index=False)
            log("✔ Excel salvo: ./saida.xlsx")
        except Exception:
            log("✖ Falha ao salvar Excel:\n" + traceback.format_exc())

    # Visuals for 2D fields
    if len(sp_vars_list) == 2 and disc_n > 1:
        # Heatmap from final state
        if save_heatmap:
            try:
                vec = np.array(last_vec, dtype=float).reshape((disc_n, disc_n), order='F')
                plt.clf()
                print("Gerando heatmap…")
                plt.figure()
                plt.imshow(
                    vec,
                    cmap='RdYlBu_r',
                    interpolation='bilinear',
                    extent=(dom_x[0], dom_x[1], dom_y[0], dom_y[1]),
                    origin='lower'
                )
                plt.colorbar(label='Valor')
                plt.title('Campo final (heatmap)')
                plt.xlabel(sp_vars_list[0])
                plt.ylabel(sp_vars_list[1])
                plt.tight_layout()
                plt.savefig(image_path, dpi=144)
                plt.close('all')
                log(f"✔ Heatmap salvo: ./{os.path.basename(image_path)}")
            except Exception:
                log("✖ Falha ao gerar heatmap:\n" + traceback.format_exc())

        # 3D animation across time
        if save_mp4:
            try:
                frames_list = sol[1][0]
                n_frames = len(frames_list)
                data3d = np.array(frames_list, dtype=float).reshape((n_frames, disc_n, disc_n), order='F')
                x = np.linspace(dom_x[0], dom_x[1], disc_n)
                y = np.linspace(dom_y[0], dom_y[1], disc_n)
                X, Y = np.meshgrid(x, y)

                fig = plt.figure()
                ax = fig.add_subplot(111, projection='3d')
                ax.set_xlim(dom_x[0], dom_x[1])
                ax.set_ylim(dom_y[0], dom_y[1])
                ax.set_zlim(np.min(data3d), np.max(data3d))
                ax.set_xlabel(sp_vars_list[0])
                ax.set_ylabel(sp_vars_list[1])
                ax.set_zlabel('Valor')

                surf = ax.plot_surface(X, Y, data3d[0], cmap='RdYlBu_r', edgecolor='none')

                def update(frame):
                    nonlocal surf
                    surf.remove()
                    surf = ax.plot_surface(X, Y, data3d[frame], cmap='RdYlBu_r', edgecolor='none')
                    ax.set_title(f"Frame {frame+1}/{n_frames}")

                ani = FuncAnimation(fig, update, frames=n_frames, interval=200, blit=False)

                writer = FFMpegWriter(fps=10, metadata=dict(artist='Você'), bitrate=1800)
                ani.save('animacao3d.mp4', writer=writer)
                plt.close(fig)
                log("✔ Animação 3D salva: ./animacao3d.mp4")
            except Exception:
                log(
                    "✖ Falha ao gerar MP4. Verifique se o FFmpeg está instalado e no PATH.\n" +
                    traceback.format_exc()
                )

    log(f"✔ Concluído com sucesso. Solver: {solver_used}")


# -------------------- Flet UI --------------------

def main(page: ft.Page):
    page.title = "PDE Solver & Visualizer"
    page.scroll = ft.ScrollMode.AUTO
    page.padding = 16
    page.theme_mode = ft.ThemeMode.LIGHT

    # Defaults (from your snippet)
    default_disc_n = 5
    default_eq = (
        "dT/dt = -dT/dx - dT/dy + 10*sech(t)**2 * x *(sin(x)+cos(y))"
        "+10*tanh(t)*(sin(x)+cos(y)+x*cos(x)-x *sin(y))"
    )
    default_funcs = "T"
    default_sp_vars = "x,y"
    default_dom_x = "0,1"
    default_dom_y = "0,1"
    default_ic = "0"

    def on_change_text(e):
        funcs_list = parse_list_csv(e.control.value)
        if len(funcs_list) == 0:
            e.control.error_text = "Esse campo não pode estar vazio."
        else:
            print(len(funcs_list))
            e.control.error_text = None
        e.control.update()

    def onchange_east_bc(e):
        hint_str = ""
        selected = e.control.value
        if selected == "Dirichlet":
            hint_str = "u(X,t) = g(X,t)"
            east_alpha_text.visible = False
            east_beta_text.visible = False
        if selected == "Neumann":
            hint_str = "∂u/∂n(X,t) = h(X,t)"    
            east_alpha_text.visible = False
            east_beta_text.visible = False
        if selected == "Robin":
            hint_str = "α(X)u + β(X)∂u/∂n = f(X,t) "
            east_alpha_text.visible = True
            east_beta_text.visible = True
            
        east_text.hint_text = hint_str
        page.update()
            
    def onchange_north_bc(e):
        hint_str = ""
        selected = e.control.value
        if selected == "Dirichlet":
            hint_str = "u(X,t) = g(X,t)"
        if selected == "Neumann":
            hint_str = "∂u/∂n(X,t) = h(X,t)"    
        if selected == "Robin":
            hint_str = "α(X)u + β(X)∂u/∂n = f(X,t) "
        north_text.hint_text = hint_str
        page.update()
            
    def onchange_south_bc(e):
        hint_str = ""
        selected = e.control.value
        if selected == "Dirichlet":
            hint_str = "u(X,t) = g(X,t)"
        if selected == "Neumann":
            hint_str = "∂u/∂n(X,t) = h(X,t)"    
        if selected == "Robin":
            hint_str = "α(X)u + β(X)∂u/∂n = f(X,t) "
            
        south_text.hint_text = hint_str
        page.update()
        
    east_text = ft.TextField(label="", hint_text="∂u/∂n(X,t) = h(X,t)", width=250)
    north_text = ft.TextField(label="", hint_text="∂u/∂n(X,t) = h(X,t)", width=250)
    south_text = ft.TextField(label="", hint_text="∂u/∂n(X,t) = h(X,t)", width=250)
    north_alpha_text = ft.TextField(label="", hint_text="α(X)", width=250, visible=False)
    north_beta_text = ft.TextField(label="", hint_text="β(X)", width=250, visible=False)
    south_alpha_text = ft.TextField(label="", hint_text="α(X)", width=250, visible=False)
    south_beta_text = ft.TextField(label="", hint_text="β(X)", width=250, visible=False)
    east_alpha_text = ft.TextField(label="", hint_text="α(X)", width=250, visible=False)
    east_beta_text = ft.TextField(label="", hint_text="β(X)", width=250, visible=False)
    
    combo_box_list_north = [ft.Dropdown(
                                                label="Condição de contorno norte",
                                                options=[ft.dropdown.Option("Dirichlet"), ft.dropdown.Option("Neumann"), ft.dropdown.Option("Robin")],
                                                value="Neumann",
                                                width=250,
                                                on_change=onchange_north_bc,
                                            ), north_text, north_alpha_text, north_beta_text]
    combo_box_list_south = [ft.Dropdown(
                                                label="Condição de contorno sul",
                                                options=[ft.dropdown.Option("Dirichlet"), ft.dropdown.Option("Neumann"), ft.dropdown.Option("Robin")],
                                                value="Neumann",
                                                width=250,
                                                on_change=onchange_south_bc,
                                            ), south_text, south_alpha_text, south_beta_text]
    combo_box_list_east = [ft.Dropdown(
                                                label="Condição de contorno leste",
                                                options=[ft.dropdown.Option("Dirichlet"), ft.dropdown.Option("Neumann"), ft.dropdown.Option("Robin")],
                                                value="Neumann",
                                                width=250,
                                                on_change=onchange_east_bc,
                                            ), east_text, east_alpha_text, east_beta_text]
    
    dropdown_column_north = ft.Column()
    dropdown_column_south = ft.Column()
    dropdown_column_east = ft.Column()

    dropdown_column_north.controls = combo_box_list_north
    dropdown_column_south.controls = combo_box_list_south
    dropdown_column_east.controls = combo_box_list_east

    def on_change_eqs(e):
        eqs_list = parse_list_csv(e.control.value)
        if len(eqs_list) == 0:
            e.control.error_text = "Esse campo não pode estar vazio."
        else:
            while len(eqs_list) > len(combo_box_list_north)-3:
                combo_box_list_north.insert(-5, ft.Dropdown(
                                                label=f"Condição de contorno norte {eqs_list[len(combo_box_list_north)-3]}",
                                                options=[ft.dropdown.Option("Dirichlet"), ft.dropdown.Option("Neumann"), ft.dropdown.Option("Robin")],
                                                value="Neumann",
                                                width=250,
                                                on_change=onchange_north_bc,
                                            ))
                print(len(combo_box_list_north))
            while len(eqs_list) > len(combo_box_list_south)-3:
                combo_box_list_south.insert(-5, ft.Dropdown(
                                                label=f"Condição de contorno sul {eqs_list[len(combo_box_list_south)-3]}",
                                                options=[ft.dropdown.Option("Dirichlet"), ft.dropdown.Option("Neumann"), ft.dropdown.Option("Robin")],
                                                value="Neumann",
                                                width=250,
                                                on_change=onchange_south_bc,
                                            ))
            while len(eqs_list) > len(combo_box_list_east)-3:
                combo_box_list_east.insert(-5, ft.Dropdown(
                                                label=f"Condição de contorno leste {eqs_list[len(combo_box_list_east)-3]}",
                                                options=[ft.dropdown.Option("Dirichlet"), ft.dropdown.Option("Neumann"), ft.dropdown.Option("Robin")],
                                                value="Neumann",
                                                width=250,
                                                on_change=onchange_east_bc,
                                            ))
            while len(eqs_list) < len(combo_box_list_north)-3:
                combo_box_list_north.pop(-5)
                print(len(combo_box_list_north))
            while len(eqs_list) < len(combo_box_list_south)-3:
                combo_box_list_south.pop(-5)
                print(len(combo_box_list_south))    
            while len(eqs_list) < len(combo_box_list_east)-3:
                combo_box_list_east.pop(-5)

            dropdown_column_north.controls = combo_box_list_north
            dropdown_column_south.controls = combo_box_list_south
            dropdown_column_east.controls = combo_box_list_east

            # print(len(eqs_list))
            e.control.error_text = None
        page.update()
        e.control.update()


    
    
    # Controls
    disc_n = ft.TextField(label="disc_n (malha por eixo)", value=str(default_disc_n), width=200)
    eq = ft.TextField(label="Equação PDE (texto)", value=default_eq, multiline=True, min_lines=3, max_lines=6, expand=True, on_change=on_change_eqs)
    funcs = ft.TextField(label="Funções (lista)", value=default_funcs, hint_text="ex.: T ou C,D", width=260, on_change=on_change_text)
    sp_vars = ft.TextField(label="Variáveis espaciais (lista)", value=default_sp_vars, hint_text="ex.: x,y", width=260, on_change=on_change_text)
    domx = ft.TextField(label="Domínio x: a,b", value=default_dom_x, width=200)
    domy = ft.TextField(label="Domínio y: a,b", value=default_dom_y, width=200)
    ic = ft.TextField(label="Condição inicial (texto)", value=default_ic, width=200, on_change=on_change_text)

    method = ft.Dropdown(
        label="Esquema de discretização (Disc.df)",
        options=[ft.dropdown.Option("backward"), ft.dropdown.Option("forward"), ft.dropdown.Option("central")],
        value="backward",
        width=200,
    )
    
    
    east_bc = ft.Dropdown(
        label="Condição de contorno leste",
        options=[ft.dropdown.Option("Dirichlet"), ft.dropdown.Option("Neumann"), ft.dropdown.Option("Robin")],
        value="Neumann",
        width=250,
        on_change=onchange_east_bc,
    )
    north_bc = ft.Dropdown(
        label="Condição de contorno norte",
        options=[ft.dropdown.Option("Dirichlet"), ft.dropdown.Option("Neumann"), ft.dropdown.Option("Robin")],
        value="Neumann",
        width=250,
        on_change=onchange_north_bc,
    )
    south_bc = ft.Dropdown(
        label="Condição de contorno sul",
        options=[ft.dropdown.Option("Dirichlet"), ft.dropdown.Option("Neumann"), ft.dropdown.Option("Robin")],
        value="Neumann",
        width=250,
        on_change=onchange_south_bc,
    )
    

    
    t0 = ft.TextField(label="t0", value="0", width=120)
    t1 = ft.TextField(label="t1", value="1", width=120)
    nsteps = ft.TextField(label="n passos (tempo)", value="10", width=160)

    use_cuda = ft.Switch(label="Usar kernel CUDA (RKF45_novo)", value=True)
    save_excel = ft.Switch(label="Salvar Excel do estado final (saida.xlsx)", value=True)
    save_heatmap = ft.Switch(label="Gerar Heatmap do estado final (heatmap.png)", value=True)
    save_mp4 = ft.Switch(label="Exportar animação 3D (animacao3d.mp4)", value=False)

    run_btn = ft.ElevatedButton(text="Rodar")
    busy = ft.ProgressRing(visible=False)

    log_view = ft.TextField(
        label="Logs",
        value="",
        multiline=True,
        min_lines=10,
        max_lines=20,
        read_only=True,
        expand=True,
    )

    img = ft.Image(src=None, width=560, height=420, fit=ft.ImageFit.CONTAIN, visible=False)
    
    def fig_to_base64(fig) -> str:
        buf = io.BytesIO()
        fig.savefig(buf, format="png", bbox_inches="tight", dpi=200)
        buf.seek(0)
        return base64.b64encode(buf.read()).decode("utf-8")
    
    def log(msg: str):
        log_view.value += (msg + "\n")
        log_view.update()

    def set_running(running: bool):
        run_btn.disabled = running
        busy.visible = running
        run_btn.update(); busy.update()

    def on_run_click(e: ft.ControlEvent):
        # Spawn the simulation in a background thread to keep UI responsive
        def worker():
            try:
                set_running(True)
                log_view.value = ""  # clear
                log("Configuração recebida. Validando…")

                # Parse inputs
                try:
                    discn = int(disc_n.value)
                    funcs_list = parse_list_csv(funcs.value)
                    spv_list = parse_list_csv(sp_vars.value)
                    dom_x = parse_domain_pair(domx.value)
                    dom_y = parse_domain_pair(domy.value)
                    t0f = float(t0.value)
                    t1f = float(t1.value)
                    nst = int(nsteps.value)
                except Exception as pe:
                    log("Erro de parsing de entradas: " + str(pe))
                    return

                image_path = os.path.abspath("heatmap.png")

                # Run
                run_simulation(
                    eq_text=eq.value,
                    funcs_list=funcs_list,
                    sp_vars_list=spv_list,
                    disc_n=discn,
                    dom_x=dom_x,
                    dom_y=dom_y,
                    ic_text=ic.value,
                    method=method.value,
                    t0=t0f,
                    t1=t1f,
                    nsteps=nst,
                    use_cuda=use_cuda.value,
                    save_excel=save_excel.value,
                    save_heatmap=save_heatmap.value,
                    save_mp4=save_mp4.value,
                    log=log,
                    image_path=image_path,
                    north_bd=north_bc.value,
                    south_bd=south_bc.value,
                    east_bd=east_bc.value,
                    north_func_bd=north_text.value,
                    south_func_bd=south_text.value,
                    east_func_bd=east_text.value,
                    north_alpha_bd=north_alpha_text.value,
                    south_alpha_bd=south_alpha_text.value,
                    east_alpha_bd=east_alpha_text.value,
                    north_beta_bd=north_beta_text.value,
                    south_beta_bd=south_beta_text.value,
                    east_beta_bd=east_beta_text.value
                )

                # Update image preview if we created it
                if save_heatmap.value and os.path.exists(image_path):
                    try:
                        with open(image_path, "rb") as f:
                            b64 = base64.b64encode(f.read()).decode("utf-8")
                        img.clean()
                        img.src = None                 # garante que não fica preso no src antigo
                        img.src_base64 = b64           # carrega direto no controle
                        img.visible = True
                        img.update()
                        page.update()
                    except Exception:
                        log("✖ Falha ao carregar heatmap no preview:\n" + traceback.format_exc())
                else:
                    img.visible = False
                    img.update()

            except Exception:
                log("Erro inesperado:\n" + traceback.format_exc())
            finally:
                set_running(False)

        threading.Thread(target=worker, daemon=True).start()

    run_btn.on_click = on_run_click

    # Layout
    page.add(
        ft.Row([
            ft.Text("PDE Solver & Visualizer", size=22, weight=ft.FontWeight.BOLD),
            ft.Container(expand=True),
            busy,
        ]),
        ft.Divider(),
        ft.ResponsiveRow([
            ft.Container(eq, col={'xs': 12, 'md': 12, 'lg': 12}),
        ]),
        ft.ResponsiveRow([
            ft.Container(disc_n, col={'xs': 6, 'md': 3, 'lg': 2}),
            ft.Container(funcs, col={'xs': 12, 'md': 4, 'lg': 3}),
            ft.Container(sp_vars, col={'xs': 12, 'md': 4, 'lg': 3}),
        ]),
        ft.ResponsiveRow([
            ft.Container(domx, col={'xs': 6, 'md': 2, 'lg': 2}),
            ft.Container(domy, col={'xs': 6, 'md': 2, 'lg': 2}),
            ft.Container(ic, col={'xs': 6, 'md': 2, 'lg': 2}),
            ft.Container(method, col={'xs': 6, 'md': 3, 'lg': 3}),
        ]),
        ft.ResponsiveRow([
            ft.Container(t0, col={'xs': 4, 'md': 2, 'lg': 2}),
            ft.Container(t1, col={'xs': 4, 'md': 2, 'lg': 2}),
            ft.Container(nsteps, col={'xs': 4, 'md': 2, 'lg': 2}),
        ]),
        ft.Divider(),
        ft.Row([
            dropdown_column_north,
            dropdown_column_south,
            dropdown_column_east,
            ],  alignment = ft.MainAxisAlignment.START),
        ft.Divider(),
        ft.ResponsiveRow([
            ft.Container(use_cuda, col={'xs': 12, 'md': 4, 'lg': 3}),
            ft.Container(save_excel, col={'xs': 12, 'md': 4, 'lg': 3}),
            ft.Container(save_heatmap, col={'xs': 12, 'md': 4, 'lg': 3}),
            ft.Container(save_mp4, col={'xs': 12, 'md': 4, 'lg': 3}),
        ]),
        ft.Row([run_btn]),
        ft.Divider(),
        ft.Row([
            ft.Container(img, expand=False),
            ft.Container(log_view, expand=True),
        ]),
    )


if __name__ == "__main__":
    ft.app(target=main)

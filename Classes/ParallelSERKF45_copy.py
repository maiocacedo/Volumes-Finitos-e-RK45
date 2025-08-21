import sympy as sp
from sympy.parsing.sympy_parser import parse_expr
import cupy as cp
from FuncAux import symbol_references

def SERKF45_cuda(oldexpr, ivar, funcs, yn, x0, xn, n, n_funcs, sp_vars):
    """
    RKF45 com passo adaptativo por componente (h[j]), replicando a lógica do SERKF45 (CPU):
      - t_j = x0 + i*h[j]
      - estágios k* calculados a partir do mesmo y0 (Jacobi)
      - y <- y4 só ao final do passo
      - se algum diffy>tol, repete o passo (i não avança)
    """
    # ----------------------------
    # 1) Símbolos consistentes
    # ----------------------------
    olddvar = symbol_references(funcs)   # ['y1','y2',...], strings (mantido)
    oldivar = symbol_references(ivar)    # ['t'], string (mantido)
    # crio Symbols e uso os MESMOS no parse e no lambdify
    sym_map = {name: sp.Symbol(name) for name in (oldivar + olddvar)}
    t_sym   = sym_map[oldivar[0]]
    y_syms  = [sym_map[name] for name in olddvar]

    # parse com local_dict para garantir que 't' e 'y' são os mesmos Symbols
    exprs = [parse_expr(e, local_dict=sym_map, evaluate=False) for e in oldexpr]
    m = len(exprs)

    # ----------------------------
    # 2) lambdify (uma função por EDO)
    #    mapeando p/ CuPy quando possível
    # ----------------------------
    cupy_map = {
        'sin': cp.sin, 'cos': cp.cos, 'tan': cp.tan,
        'asin': cp.arcsin, 'acos': cp.arccos, 'atan': cp.arctan,
        'sinh': cp.sinh, 'cosh': cp.cosh, 'tanh': cp.tanh,
        'exp': cp.exp, 'log': cp.log, 'sqrt': cp.sqrt,
        'Abs': cp.abs, 'sign': cp.sign,
        'Max': cp.maximum, 'Min': cp.minimum, 'mod': cp.mod,
        'floor': cp.floor, 'ceil': cp.ceil,
    }
    f_list = [sp.lambdify((t_sym, *y_syms), exprs[j], modules=[cupy_map, cp])
              for j in range(m)]

    # helper: avalia f_j(t, y_stage) -> cp.float64 escalar
    def f_j_val(j, t_val, y_stage):
        # (para evitar conversões implícitas estranhas)
        t_py = float(cp.asnumpy(t_val))
        y_py = [float(cp.asnumpy(y_stage[k])) for k in range(m)]
        v = f_list[j](t_py, *y_py)  # pode vir como float Python
        return cp.float64(v)

    # ----------------------------
    # 3) Estados, h e histórico
    # ----------------------------
    y  = cp.asarray(yn, dtype=cp.float64).reshape(m,)
    if y.size != m:
        raise ValueError(f"len(yn) ({y.size}) != número de EDOs ({m}).")

    # h e s "vetoriais" (mesma ideia do seu sp.zeros + loop)
    h = cp.full((m,), (xn - x0)/max(int(n),1), dtype=cp.float64)
    s = cp.zeros((m,), dtype=cp.float64)
    tol = 1e-3  # mesmo valor do seu código

    # estrutura de histórico no seu formato
    final_list = [[] for _ in range(n_funcs)]
    if n_funcs > 0 and (m % n_funcs == 0):
        n_elements = m // n_funcs
        if sp_vars == 2:
            # snapshot inicial (igual ao seu if sp_vars==2 no início)
            y_host = y.get().reshape((n_funcs, n_elements))
            for jgrp in range(n_funcs):
                final_list[jgrp].append(y_host[jgrp].tolist())
    else:
        n_elements = m  # fallback

    # ----------------------------
    # 4) Laço de passos (com repetição)
    # ----------------------------
    i = 0
    max_retries = 10_000_000  # guarda
    retries = 0

    while i < n:
        checkh = False

        # y0 é o "yn" do começo do passo (Jacobi)
        y0 = y.copy()

        # buffers por passo
        k1 = cp.zeros((m,), dtype=cp.float64)
        k2 = cp.zeros((m,), dtype=cp.float64)
        k3 = cp.zeros((m,), dtype=cp.float64)
        k4 = cp.zeros((m,), dtype=cp.float64)
        k5 = cp.zeros((m,), dtype=cp.float64)
        k6 = cp.zeros((m,), dtype=cp.float64)
        y4 = cp.zeros((m,), dtype=cp.float64)
        y5 = cp.zeros((m,), dtype=cp.float64)

        # --- loop em j: exatamente como no seu código ---
        for j in range(m):
            hj = h[j]
            tj = cp.float64(x0) + cp.float64(i) * hj  # t_j = x0 + i*h[j]

            # k1
            k1[j] = hj * f_j_val(j, tj, y0)

            # k2 : y_stage = y0 + (1/4)*k1[j] (esc. somado a TODAS as variáveis)
            y2 = y0 + (0.25 * k1[j])
            k2[j] = hj * f_j_val(j, tj + 0.25*hj, y2)

            # k3 : y0 + (3/32)*k1[j] + (9/32)*k2[j]
            y3 = y0 + ((3.0/32.0)*k1[j] + (9.0/32.0)*k2[j])
            k3[j] = hj * f_j_val(j, tj + (3.0/8.0)*hj, y3)

            # k4 : y0 + (1932/2197)*k1[j] - (7200/2197)*k2[j] + (7296/2197)*k3[j]
            y4s = y0 + ((1932.0/2197.0)*k1[j] - (7200.0/2197.0)*k2[j] + (7296.0/2197.0)*k3[j])
            k4[j] = hj * f_j_val(j, tj + (12.0/13.0)*hj, y4s)

            # k5 : y0 + (439/216)*k1[j] - 8*k2[j] + (3680/513)*k3[j] - (845/4104)*k4[j]
            y5s = y0 + ((439.0/216.0)*k1[j] - 8.0*k2[j] + (3680.0/513.0)*k3[j] - (845.0/4104.0)*k4[j])
            k5[j] = hj * f_j_val(j, tj + hj, y5s)

            # k6 : y0 - (8/27)*k1[j] + 2*k2[j] - (3544/2565)*k3[j] + (1859/4104)*k4[j] - (11/40)*k5[j]
            y6s = y0 + (-(8.0/27.0)*k1[j] + 2.0*k2[j] - (3544.0/2565.0)*k3[j] + (1859.0/4104.0)*k4[j] - (11.0/40.0)*k5[j])
            k6[j] = hj * f_j_val(j, tj + 0.5*hj, y6s)

            # y4[j], y5[j]
            y4[j] = y0[j] + (25.0/216.0)*k1[j] + (1408.0/2565.0)*k3[j] + (2197.0/4104.0)*k4[j] - (1.0/5.0)*k5[j]
            y5[j] = y0[j] + (16.0/135.0)*k1[j] + (6656.0/12825.0)*k3[j] + (28561.0/56430.0)*k4[j] - (9.0/50.0)*k5[j] + (2.0/55.0)*k6[j]

        # igual ao seu: salva "yns.append(yn)" ANTES de y <- y4
        # (aqui não mantenho uma lista 'yns', mas replico o bloco sp_vars==2)
        if sp_vars == 2 and n_funcs > 0 and (m % n_funcs == 0):
            y0_host = y0.get().reshape((n_funcs, m // n_funcs))
            for jgrp in range(n_funcs):
                final_list[jgrp].append(y0_host[jgrp].tolist())

        # y <- y4 (Jacobi)
        y = y4.copy()

        # erro por componente e atualização de h[j] SÓ se diff>tol
        diff = cp.abs(y5 - y4)
        # s[j] = 0.84 * (tol*h[j]/diff[j])**0.25 (onde diff>0)
        # para diff==0, podemos deixar s[j]=0 (não usado) ou 2.0 (não usado se diff<=tol)
        mask_pos = diff > 0
        s = cp.where(mask_pos, 0.84 * (tol * h / diff) ** 0.25, s)
        # clipes usuais (só interessam quando diff>tol)
        s = cp.clip(s, 0.1, 4.0)

        # onde diff>tol, ajusta h[j] e marca checkh
        mask_bad = diff > tol
        if bool(cp.any(mask_bad).item()):
            h = cp.where(mask_bad, s * h, h)
            checkh = True
        if not bool(cp.any(mask_bad).item()):
            print(i)

        # histórico do passo (após y <- y4)
        if sp_vars == 2 and n_funcs > 0 and (m % n_funcs == 0):
            y_host = y.get().reshape((n_funcs, m // n_funcs))
            for jgrp in range(n_funcs):
                final_list[jgrp].append(y_host[jgrp].tolist())

        # repetição do passo se necessário (equivalente ao seu "i=i-1")
        if checkh:
            retries += 1
            if retries > max_retries:
                raise RuntimeError("Muitas repetições de passo; verifique tolerância/dinâmica.")
            continue  # não avança i
        else:
            retries = 0
            i += 1     # avança para o próximo passo

    # retorno
    return y.get(), final_list

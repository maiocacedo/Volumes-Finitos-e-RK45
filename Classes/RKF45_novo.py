import sympy as sp
from sympy.parsing.sympy_parser import parse_expr
import cupy as cp
from FuncAux import symbol_references

def SERKF45_cuda(oldexpr, ivar, funcs, yn, x0, xn, n, n_funcs, sp_vars):
    """
    RKF45 com passo adaptativo por componente (h[j]), lógica Jacobi.
    Versão CuPy otimizada (sem cópias host↔device no miolo do passo).
    """

    # =========================
    # 1) Símbolos consistentes
    # =========================
    olddvar = symbol_references(funcs)   # ['y1','y2',...]
    oldivar = symbol_references(ivar)    # ['t']
    sym_map = {name: sp.Symbol(name) for name in (oldivar + olddvar)}
    t_sym   = sym_map[oldivar[0]]
    y_syms  = [sym_map[name] for name in olddvar]

    exprs = [parse_expr(e, local_dict=sym_map, evaluate=False) for e in oldexpr]
    m = len(exprs)
    if m == 0:
        raise ValueError("Lista de EDOs vazia.")

    # =========================
    # 2) lambdify sem Matrix  ⟶  tuple
    # =========================
    cupy_map = {
        'sin': cp.sin, 'cos': cp.cos, 'tan': cp.tan,
        'asin': cp.arcsin, 'acos': cp.arctan, 'atan': cp.arctan, 'atan2': cp.arctan2,
        'sinh': cp.sinh, 'cosh': cp.cosh, 'tanh': cp.tanh,
        'exp': cp.exp, 'log': cp.log, 'sqrt': cp.sqrt,
        'Abs': cp.abs, 'sign': cp.sign,
        'Max': cp.maximum, 'Min': cp.minimum, 'mod': cp.mod,
        'floor': cp.floor, 'ceil': cp.ceil, 'sech': lambda x: 1.0/cp.cosh(x),
        # se usar Piecewise, garanta forma avaliável no device
    }

    # Em vez de sp.Matrix(exprs), use tuple(exprs) para evitar cp.array(...) interno
    F = sp.lambdify((t_sym, *y_syms), tuple(exprs), modules=[cupy_map, cp])

    def _to_vec(o, template):
        """Garante vetor CuPy compatível com t_vec (broadcast)."""
        arr = cp.asarray(o)
        if arr.ndim == 0:
            return cp.full_like(template, arr)
        return arr

    def F_diag(t_vec, y_vec):
        """
        Avalia (f1,...,fm) com t_vec (m,) e y_vec (m,) fixo,
        empilha e extrai a diagonal M[i,i] = f_i(t_i, y_vec).
        """
        # y_vec[k] são 0-D cupy scalars; broadcasta com t_vec (m,)
        out = F(t_vec, *[y_vec[k] for k in range(m)])
        if not isinstance(out, (tuple, list)) or len(out) != m:
            raise RuntimeError("F retornou forma inesperada; espere tuple/list de tamanho m.")
        M = cp.stack([_to_vec(o, t_vec) for o in out], axis=0)  # (m, m)
        idx = cp.arange(m)
        return M[idx, idx]  # (m,)

    # =========================
    # 3) Estados e histórico
    # =========================
    dtype = cp.float64
    y  = cp.asarray(yn, dtype=dtype).reshape(m,)
    if y.size != m:
        raise ValueError(f"len(yn) ({y.size}) != número de EDOs ({m}).")

    h = cp.full((m,), (xn - x0)/max(int(n), 1), dtype=dtype)
    s = cp.zeros((m,), dtype=dtype)
    tol = dtype(1e-3)

    final_list = [[] for _ in range(n_funcs)]
    if n_funcs > 0 and (m % n_funcs == 0):
        n_elements = m // n_funcs
        if sp_vars == 2:
            y_host = y.get().reshape((n_funcs, n_elements))
            for jgrp in range(n_funcs):
                final_list[jgrp].append(y_host[jgrp].tolist())
    else:
        n_elements = m

    # =========================
    # 4) Buffers reutilizados
    # =========================
    y0   = cp.empty_like(y)
    k1   = cp.empty_like(y); k2 = cp.empty_like(y); k3 = cp.empty_like(y)
    k4   = cp.empty_like(y); k5 = cp.empty_like(y); k6 = cp.empty_like(y)
    y4   = cp.empty_like(y); y5 = cp.empty_like(y)
    ytmp = cp.empty_like(y)

    t1 = cp.empty_like(y); t2 = cp.empty_like(y); t3 = cp.empty_like(y)
    t4 = cp.empty_like(y); t5 = cp.empty_like(y); t6 = cp.empty_like(y)

    @cp.fuse()
    def _s_and_h(diff, h, tol):
        s = cp.where(diff > 0, 0.84 * (tol * h / diff) ** 0.25, 0.0)
        s = cp.clip(s, 0.1, 4.0)
        hnew = cp.where(diff > tol, s * h, h)
        return s, hnew

    # =========================
    # 5) Laço de passos
    # =========================
    i = 0
    max_retries = 10_000_000
    retries = 0
    x0_scalar = float(x0)

    while i < n:
        checkh = False

        y0[...] = y
        t1[...] = x0_scalar + (i * h)

        # k1
        f1 = F_diag(t1, y0);                  k1[...] = h * f1

        # k2
        cp.multiply(k1, 0.25, out=ytmp);      cp.add(ytmp, y0, out=ytmp)
        t2[...] = t1 + 0.25 * h
        f2 = F_diag(t2, ytmp);                k2[...] = h * f2

        # k3
        cp.multiply(k1, 3.0/32.0, out=y4)
        cp.multiply(k2, 9.0/32.0, out=y5)
        cp.add(y4, y5, out=ytmp);             cp.add(ytmp, y0, out=ytmp)
        t3[...] = t1 + (3.0/8.0) * h
        f3 = F_diag(t3, ytmp);                k3[...] = h * f3

        # k4
        cp.multiply(k1, 1932.0/2197.0, out=y4)
        cp.multiply(k2, -7200.0/2197.0, out=y5)
        cp.add(y4, y5, out=ytmp)
        cp.multiply(k3, 7296.0/2197.0, out=y4)
        cp.add(ytmp, y4, out=ytmp);           cp.add(ytmp, y0, out=ytmp)
        t4[...] = t1 + (12.0/13.0) * h
        f4 = F_diag(t4, ytmp);                k4[...] = h * f4

        # k5
        cp.multiply(k1, 439.0/216.0, out=y4)
        cp.multiply(k2, -8.0, out=y5)
        cp.add(y4, y5, out=ytmp)
        cp.multiply(k3, 3680.0/513.0, out=y4)
        cp.add(ytmp, y4, out=ytmp)
        cp.multiply(k4, -845.0/4104.0, out=y4)
        cp.add(ytmp, y4, out=ytmp);           cp.add(ytmp, y0, out=ytmp)
        t5[...] = t1 + h
        f5 = F_diag(t5, ytmp);                k5[...] = h * f5

        # k6
        cp.multiply(k1, -(8.0/27.0), out=y4)
        cp.multiply(k2, 2.0, out=y5)
        cp.add(y4, y5, out=ytmp)
        cp.multiply(k3, -3544.0/2565.0, out=y4)
        cp.add(ytmp, y4, out=ytmp)
        cp.multiply(k4, 1859.0/4104.0, out=y4)
        cp.add(ytmp, y4, out=ytmp)
        cp.multiply(k5, -11.0/40.0, out=y4)
        cp.add(ytmp, y4, out=ytmp);           cp.add(ytmp, y0, out=ytmp)
        t6[...] = t1 + 0.5 * h
        f6 = F_diag(t6, ytmp);                k6[...] = h * f6

        # y4 (ordem 4)
        cp.multiply(k1, 25.0/216.0, out=y4)
        cp.multiply(k3, 1408.0/2565.0, out=y5);       cp.add(y4, y5, out=y4)
        cp.multiply(k4, 2197.0/4104.0, out=ytmp);     cp.add(y4, ytmp, out=y4)
        cp.multiply(k5, -1.0/5.0, out=ytmp);          cp.add(y4, ytmp, out=y4)
        cp.add(y4, y0, out=y4)

        # y5 (ordem 5)
        cp.multiply(k1, 16.0/135.0, out=y5)
        cp.multiply(k3, 6656.0/12825.0, out=ytmp);    cp.add(y5, ytmp, out=y5)
        cp.multiply(k4, 28561.0/56430.0, out=ytmp);   cp.add(y5, ytmp, out=y5)
        cp.multiply(k5, -9.0/50.0, out=ytmp);         cp.add(y5, ytmp, out=y5)
        cp.multiply(k6, 2.0/55.0, out=ytmp);          cp.add(y5, ytmp, out=y5)
        cp.add(y5, y0, out=y5)

        # histórico ANTES de y <- y4
        if sp_vars == 2 and n_funcs > 0 and (m % n_funcs == 0):
            y0_host = y0.get().reshape((n_funcs, m // n_funcs))
            for jgrp in range(n_funcs):
                final_list[jgrp].append(y0_host[jgrp].tolist())

        # aceitar y4 provisoriamente
        y[...] = y4

        # erro e ajuste de h
        diff = cp.abs(y5 - y4)
        s, h_candidate = _s_and_h(diff, h, tol)
        bad_mask = diff > tol
        bad_count = int(cp.count_nonzero(bad_mask).get())

        # histórico APÓS y <- y4
        if sp_vars == 2 and n_funcs > 0 and (m % n_funcs == 0):
            y_host = y.get().reshape((n_funcs, m // n_funcs))
            for jgrp in range(n_funcs):
                final_list[jgrp].append(y_host[jgrp].tolist())

        if bad_count > 0:
            h = cp.where(bad_mask, h_candidate, h)  # diminui passo onde necessário
            retries += 1
            if retries > max_retries:
                raise RuntimeError("Muitas repetições de passo; verifique tolerância/dinâmica.")
            continue  # repete o passo (i não avança)
        else:
            print(i)
            retries = 0
            i += 1

    return y.get(), final_list

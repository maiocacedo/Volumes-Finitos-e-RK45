import numpy as np
import sympy as sp
from sympy.parsing.sympy_parser import parse_expr

def thomas_solver(a, b, c, d):
    """Resolve o sistema tridiagonal Ax = d de forma vetorizada."""
    n = len(d)
    cp = np.zeros(n-1)
    dp = np.zeros(n)
    x = np.zeros(n)
    
    cp[0] = c[0] / b[0]
    dp[0] = d[0] / b[0]
    for i in range(1, n):
        denom = b[i] - a[i] * cp[i-1]
        if i < n-1:
            cp[i] = c[i] / denom
        dp[i] = (d[i] - a[i] * dp[i-1]) / denom
    
    x[-1] = dp[-1]
    for i in range(n-2, -1, -1):
        x[i] = dp[i] - cp[i] * x[i+1]
    return x

def cn_1d_optimized(flat_list, d_vars, nt, dt, u_init_array):
    """
    Solver Crank-Nicolson 1D com vetorização completa.
    Resolve sistemas lineares/não-lineares extraindo coeficientes via SymPy.
    """
    n = len(d_vars)
    u = u_init_array.copy().astype(np.float64)
    dx = 1.0 / (n - 1)
    
    t_sym = sp.Symbol('t')
    sym_list = [sp.Symbol(v) for v in d_vars]
    parsed_eqs = [parse_expr(eq_str) for eq_str in flat_list]
    
    # --- COMPILAÇÃO VETORIZADA ---
    # Geramos funções que aceitam o vetor 'u' completo de uma vez
    func_map = []
    for i in range(n):
        expr = parsed_eqs[i]
        # Extrai Coeficientes para a, b, c (sub, diag, super)
        # i-1: a, i: b, i+1: c
        c_i_m1 = expr.coeff(sym_list[i-1]) if i > 0 else 0
        c_i    = expr.coeff(sym_list[i])
        c_i_p1 = expr.coeff(sym_list[i+1]) if i < n-1 else 0
        fonte  = expr.as_coeff_Add()[0]
        
        # Lambdify com 'numpy' permite passar o array u_args
        f_a = sp.lambdify((t_sym, *sym_list), c_i_m1, 'numpy')
        f_b = sp.lambdify((t_sym, *sym_list), c_i, 'numpy')
        f_c = sp.lambdify((t_sym, *sym_list), c_i_p1, 'numpy')
        f_s = sp.lambdify((t_sym, *sym_list), fonte, 'numpy')
        
        func_map.append((f_a, f_b, f_c, f_s))

    # --- LOOP DE TEMPO ---
    for passo in range(nt):
        t_atual = passo * dt
        u_args = tuple(u) # Desempacota o vetor para as funções lambdify
        
        a_vec, b_vec, c_vec, s_vec = np.zeros(n), np.zeros(n), np.zeros(n), np.zeros(n)
        
        # Avaliação dos coeficientes (Linearização de Picard no tempo n)
        for i in range(n):
            fa, fb, fc, fs = func_map[i]
            a_val = fa(t_atual, *u_args)
            b_val = fb(t_atual, *u_args)
            c_val = fc(t_atual, *u_args)
            s_val = fs(t_atual, *u_args)
            
            if i == 0 or i == n-1: # Borda
                b_vec[i] = 1.0
                s_vec[i] = fs(t_atual, *u_args)
            else:
                # Termos implícitos (LHS): [I - (dt/2) * Matriz_Espacial]
                a_vec[i] = -(dt / 2.0) * a_val
                b_vec[i] = 1.0 - (dt / 2.0) * b_val
                c_vec[i] = -(dt / 2.0) * c_val
                s_vec[i] = s_val

        # --- MONTAGEM DO LADO DIREITO (RHS) ---
        rhs = np.zeros(n)
        rhs[0] = s_vec[0]
        rhs[-1] = s_vec[-1]
        
        # RHS = [I + (dt/2) * Matriz_Espacial] * u_n + dt * fonte
        # Usamos os mesmos coeficientes a, b, c calculados acima
        for i in range(1, n-1):
            # Lembre-se: a_vec, b_vec, c_vec já estão multiplicados por -(dt/2)
            # Para o RHS precisamos de +(dt/2), então invertemos o sinal
            term_espacial_n = (
                -a_vec[i] * u[i-1] + 
                (1.0 - b_vec[i]) * u[i] + # (1 - b_implícito) = (dt/2)*coeff
                -c_vec[i] * u[i+1]
            )
            rhs[i] = u[i] + term_espacial_n + dt * s_vec[i]

        # Resolve
        u = thomas_solver(a_vec, b_vec, c_vec, rhs)
        
    return u
import numpy as np
import sympy as sp
from sympy.parsing.sympy_parser import parse_expr

def thomas_solver(a, b, c, d):
    """Resolve o sistema tridiagonal Ax = d."""
    n = len(d)
    cp, dp, x = np.zeros(n-1), np.zeros(n), np.zeros(n)
    
    # Eliminação progressiva
    cp[0] = c[0] / b[0]
    dp[0] = d[0] / b[0]
    for i in range(1, n):
        denom = b[i] - a[i] * cp[i-1]
        if i < n-1:
            cp[i] = c[i] / denom
        dp[i] = (d[i] - a[i] * dp[i-1]) / denom
    
    # Substituição regressiva
    x[-1] = dp[-1]
    for i in range(n-2, -1, -1):
        x[i] = dp[i] - cp[i] * x[i+1]
    return x

def cn_1d(flat_list, d_vars, nt, dt, u_init_val):
    """
    Solver Crank-Nicolson 1D otimizado com Lambdify.
    Resolve EDPs Lineares e Não-Lineares (via linearização de Picard).
    """
    n = len(d_vars)
    u = np.full(n, u_init_val, dtype=np.float64)
    
    # Preparação de Símbolos
    t_sym = sp.Symbol('t')
    sym_list = [sp.Symbol(v) for v in d_vars]
    parsed_eqs = [parse_expr(eq_str) for eq_str in flat_list]
    
    # --- ETAPA DE COMPILAÇÃO (FORA DO LOOP) ---
    # Guardamos funções que calculam coeficientes e termos de fonte
    mapa_coeficientes = {} 
    
    
    for i, expr in enumerate(parsed_eqs):
        mapa_coeficientes[i] = {'coeffs': [], 'fonte': None}
        
        # Borda (Dirichlet)
        if i == 0 or i == n - 1:
            mapa_coeficientes[i]['fonte'] = sp.lambdify((t_sym, *sym_list), expr)
        else:
            # Pontos internos: extrair dF/dt = Coeff*F + Fonte
            # Note: Para não-lineares, o coeff_sym ainda terá símbolos
            for j, sym in enumerate(sym_list):
                coeff_sym = expr.coeff(sym)
                if coeff_sym != 0:
                    func_coeff = sp.lambdify((t_sym, *sym_list), coeff_sym)
                    mapa_coeficientes[i]['coeffs'].append((j, func_coeff))
            
            # Extrair termo independente (fonte)
            fonte_sym = expr.as_coeff_Add()[0]
            mapa_coeficientes[i]['fonte'] = sp.lambdify((t_sym, *sym_list), fonte_sym)
    
    # --- LOOP DE TEMPO (ALTA PERFORMANCE) ---
    tempo_total = nt * dt
    for passo in range(nt):
        tempo_atual = passo * dt
        a_diag, b_diag, c_diag = np.zeros(n), np.zeros(n), np.zeros(n)
        fontes_s = np.zeros(n)
        rhs = np.zeros(n)
        
        # u_args ajuda a passar o array para o lambdify de uma vez
        u_args = tuple(u)

        for i in range(n):
            if i == 0 or i == n - 1:
                b_diag[i] = 1.0
                fontes_s[i] = mapa_coeficientes[i]['fonte'](tempo_atual, *u_args)
            else:
                # Extrai coeficientes numéricos (linearização no tempo n)
                for j, func_c in mapa_coeficientes[i]['coeffs']:
                    c_val = func_c(tempo_atual, *u_args)
                    val_implicito = -(dt / 2.0) * c_val
                    
                    if i == j: b_diag[i] = 1.0 + val_implicito
                    elif j == i - 1: a_diag[i] = val_implicito
                    elif j == i + 1: c_diag[i] = val_implicito
                
                # Termo de fonte numérico
                fontes_s[i] = mapa_coeficientes[i]['fonte'](tempo_atual, *u_args)

        # Montagem do Lado Direito (RHS) do Crank-Nicolson
        for i in range(n):
            if i == 0 or i == n - 1:
                rhs[i] = fontes_s[i]
            else:
                # RHS = [I - (dt/2)A] * u^n + dt * fonte
                rhs[i] = u[i] * (2.0 - b_diag[i])
                if i > 0:   rhs[i] -= a_diag[i] * u[i-1]
                if i < n-1: rhs[i] -= c_diag[i] * u[i+1]
                rhs[i] += dt * fontes_s[i] 

        # Resolve o sistema linear tridiagonal
        u = thomas_solver(a_diag, b_diag, c_diag, rhs)
        
    return u
# cn_2d_adi: Solver Crank-Nicolson 2D ADI
def cn_2d_adi(flat_list, d_vars, u_init, nt, dt):
    """
    flat_list: Lista de strings discretizadas gerada pela sua função df()
    d_vars: Lista de nomes das variáveis simbólicas (XX0_i_j)
    """
    nx, ny = u_init.shape
    u = u_init.copy()
    symbols = {v: sp.Symbol(v) for v in d_vars}
    
    # --- EXTRAÇÃO DOS COEFICIENTES p, q, r, s DO QUADRO ---
    # Pegamos uma equação central da flat_list para extrair a lógica espacial
    # Nota: Assumimos que a física é uniforme na malha
    exemplo_eq = parse_expr(flat_list[len(flat_list)//2])
    
    # Coeficientes espaciais puros (L) da sua engine
    # i+1, i, i-1
    c_p1 = float(exemplo_eq.coeff(symbols[d_vars[len(d_vars)//2 + 1]]))
    c_i  = float(exemplo_eq.coeff(symbols[d_vars[len(d_vars)//2]]))
    c_m1 = float(exemplo_eq.coeff(symbols[d_vars[len(d_vars)//2 - 1]]))
    # Termo de fonte s (termo constante na expressão)
    c_s  = float(exemplo_eq.as_coeff_Add()[0])

    # De acordo com seu quadro:
    p = -(dt / 2.0) * c_p1
    q = 1.0 - (dt / 2.0) * c_i
    r = -(dt / 2.0) * c_m1
    s = -(dt / 2.0) * c_s
    
    # Coeficiente q_chapéu do lado direito (RHS)
    q_hat = 1.0 + (dt / 2.0) * c_i

    # Loop Temporal
    for t in range(nt):
        # --- MEIO-PASSO 1: X Implícito ---
        u_meio = u.copy()
        for j in range(1, ny - 1):
            rhs = np.zeros(nx - 2)
            for i in range(1, nx - 1):
                # Aplicando a lógica do lado direito do quadro:
                # RHS = -p*F_i+1 + q_hat*F_i - r*F_i-1 - s
                rhs[i-1] = -p*u[i, j+1] + q_hat*u[i, j] - r*u[i, j-1] - 2*s
            
            # Resolve a linha com os coeficientes p, q, r do LHS
            u_meio[1:-1, j] = thomas_solver(np.full(nx-2, r), 
                                            np.full(nx-2, q), 
                                            np.full(nx-2, p), rhs)

        # --- MEIO-PASSO 2: Y Implícito ---
        u_novo = u_meio.copy()
        for i in range(1, nx - 1):
            rhs = np.zeros(ny - 2)
            for j in range(1, ny - 1):
                # RHS usando u_meio
                rhs[j-1] = -p*u_meio[i+1, j] + q_hat*u_meio[i, j] - r*u_meio[i-1, j] - 2*s
            
            u_novo[i, 1:-1] = thomas_solver(np.full(ny-2, r), 
                                            np.full(ny-2, q), 
                                            np.full(ny-2, p), rhs)
        u = u_novo.copy()
        
    return u
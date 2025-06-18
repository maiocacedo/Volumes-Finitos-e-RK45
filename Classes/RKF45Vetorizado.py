import sympy as sp
from sympy.parsing.sympy_parser import parse_expr
import cupy as cp


def symbol_references(in_list):
    slist = []
    for e in in_list:
        globals()[e] = sp.Symbol(e)
        slist.append(e)
    return slist

def SERKF45_cuda(oldexpr, ivar, funcs, yn, x0, xn, n, nfuncs, sp_vars):
    # 1) símbolos e lambdify para CuPy
    olddvar = symbol_references(funcs)
    oldivar = symbol_references(ivar)
    exprs = [parse_expr(e) for e in oldexpr]
    f_gpu = sp.lambdify(
        (oldivar[0], *olddvar),
        exprs,
        modules=[cp]
    )
    from inspect import signature
    
    print("len(funcs)       =", len(funcs))
    print("len(initial_values) =", len(yn))
    print("f_gpu signature  =", signature(f_gpu))
    
    n_elements = len(yn) // nfuncs
    final_list = [ [] for i in range(nfuncs)] 
    
    
    

    # 2) inicialização
    m = len(exprs)
    y  = cp.array(yn, dtype=cp.float64)
    h  = cp.full((m,), (xn - x0)/n, dtype=cp.float64)
    tol = 0.001
    
    y_np = cp.asnumpy(y).reshape((nfuncs, n_elements))
    for j in range(nfuncs):
        final_list[j].append(y_np[j].tolist())
        
    # buffers RKF
    k1 = cp.zeros(m, dtype=cp.float64)
    k2 = cp.zeros(m, dtype=cp.float64)
    k3 = cp.zeros(m, dtype=cp.float64)
    k4 = cp.zeros(m, dtype=cp.float64)
    k5 = cp.zeros(m, dtype=cp.float64)
    k6 = cp.zeros(m, dtype=cp.float64)

    # helper corrigido: preenche manualmente um cupy.ndarray
    def get_deriv(x_arr, y_arr):
        raw = f_gpu(x_arr, *[y_arr[i] for i in range(m)])
        arr = cp.empty(m, dtype=cp.float64)
        for i, val in enumerate(raw):
            arr[i] = val   # aceita Python float ou cupy scalar
        return arr
        

    # 3) loop adaptativo
    for step in range(n):
        x_base = x0 + step * h

        k1 = h * get_deriv(x_base,                    y)
        k2 = h * get_deriv(x_base + 0.25*h,           y + 0.25*k1)
        k3 = h * get_deriv(x_base + (3/8)*h,         y + (3/32)*k1 + (9/32)*k2)
        k4 = h * get_deriv(x_base + (12/13)*h,
                          y + (1932/2197)*k1 - (7200/2197)*k2 + (7296/2197)*k3)
        k5 = h * get_deriv(x_base +      h,
                          y + (439/216)*k1 - 8*k2 + (3680/513)*k3 - (845/4104)*k4)
        k6 = h * get_deriv(x_base + 0.5*h,
                          y - (8/27)*k1 + 2*k2 - (3544/2565)*k3 + (1859/4104)*k4 - (11/40)*k5)

        yn4 = y + (25/216)*k1 + (1408/2565)*k3 + (2197/4104)*k4 - (1/5)*k5
        yn5 = y + (16/135)*k1 + (6656/12825)*k3 + (28561/56430)*k4 - (9/50)*k5 + (2/55)*k6

        diff = cp.abs(yn5 - yn4)
        s    = 0.84 * (tol * h / diff)**0.25
        h    = cp.where(diff > tol, s * h, h)

        y = yn4
        
        y_np = cp.asnumpy(y).reshape((nfuncs, n_elements))
        
        for j in range(nfuncs):
            final_list[j].append(y_np[j].tolist())
            
    return cp.asnumpy(y), final_list

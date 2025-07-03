import cupy as cp
import sympy as sp
from sympy.parsing.sympy_parser import parse_expr

def SERKF45_cuda(oldexpr, ivar, funcs, yn, x0, xn, n, nfuncs, ):
    exprs = [parse_expr(e) for e in oldexpr]

    def format_expr(expr, y_prefix="y"):
        ccode = sp.ccode(expr)
        for j, f in enumerate(funcs):
            ccode = ccode.replace(f, f'{y_prefix}[{j}]')
        ccode = ccode.replace(ivar[0], 'x')
        return ccode

    def deriv_code(y_label, out_label):
        return '\n'.join([
            f'        {out_label}[{i}] = h * ({format_expr(expr, y_label)});'
            for i, expr in enumerate(exprs)
        ])

    # Geração do kernel completo com todas as etapas RKF45
    kernel = f"""
    extern "C" _global_
    void rkf45_kernel(double* y_in, double x0, double h, int nsteps, int nfuncs, double tol, double* y_out) {{
        int idx = threadIdx.x;
        if (idx != 0) return;

        double y[{nfuncs}];
    """
    for i in range(nfuncs):
        kernel += f"        y[{i}] = y_in[{i}];\n"

    kernel += f"""
        double x = x0;
        double k1[{nfuncs}], k2[{nfuncs}], k3[{nfuncs}], k4[{nfuncs}], k5[{nfuncs}], k6[{nfuncs}];
        double yn4[{nfuncs}], yn5[{nfuncs}], diff[{nfuncs}], s[{nfuncs}], yt[{nfuncs}];

        for (int step = 0; step < nsteps; ++step) {{
{deriv_code('y', 'k1')}

            for (int i = 0; i < nfuncs; ++i) yt[i] = y[i] + 0.25 * k1[i];
            x = x0 + step * h + 0.25 * h;
{deriv_code('yt', 'k2')}

            for (int i = 0; i < nfuncs; ++i)
                yt[i] = y[i] + (3.0/32.0)*k1[i] + (9.0/32.0)*k2[i];
            x = x0 + step * h + (3.0/8.0) * h;
{deriv_code('yt', 'k3')}

            for (int i = 0; i < nfuncs; ++i)
                yt[i] = y[i] + (1932.0/2197.0)*k1[i] - (7200.0/2197.0)*k2[i] + (7296.0/2197.0)*k3[i];
            x = x0 + step * h + (12.0/13.0) * h;
{deriv_code('yt', 'k4')}

            for (int i = 0; i < nfuncs; ++i)
                yt[i] = y[i] + (439.0/216.0)*k1[i] - 8.0*k2[i] + (3680.0/513.0)*k3[i] - (845.0/4104.0)*k4[i];
            x = x0 + step * h + h;
{deriv_code('yt', 'k5')}

            for (int i = 0; i < nfuncs; ++i)
                yt[i] = y[i] - (8.0/27.0)*k1[i] + 2.0*k2[i] - (3544.0/2565.0)*k3[i] +
                        (1859.0/4104.0)*k4[i] - (11.0/40.0)*k5[i];
            x = x0 + step * h + 0.5 * h;
{deriv_code('yt', 'k6')}

            for (int i = 0; i < nfuncs; ++i) {{
                yn4[i] = y[i] + (25.0/216.0)*k1[i] + (1408.0/2565.0)*k3[i] +
                         (2197.0/4104.0)*k4[i] - (1.0/5.0)*k5[i];
                yn5[i] = y[i] + (16.0/135.0)*k1[i] + (6656.0/12825.0)*k3[i] +
                         (28561.0/56430.0)*k4[i] - (9.0/50.0)*k5[i] + (2.0/55.0)*k6[i];
                diff[i] = fabs(yn5[i] - yn4[i]);
                s[i] = 0.84 * pow(tol * h / diff[i], 0.25);
            }}
            for (int i = 0; i < nfuncs; ++i) {{
                h = diff[i] > tol ? s[i] * h : h;
                y[i] = yn4[i];
            }}
        }}
    """
    for i in range(nfuncs):
        kernel += f"        y_out[{i}] = y[{i}];\n"
    kernel += "    }\n"

    # Compilação do kernel
    mod = cp.RawModule(code=kernel)
    rkf_kernel = mod.get_function("rkf45_kernel")

    # Execução do kernel
    y_in = cp.array(yn, dtype=cp.float64)
    y_out = cp.zeros_like(y_in)
    h = cp.float64((xn - x0) / n)
    x0 = cp.float64(x0)
    nsteps = cp.int32(n)
    nfuncs = cp.int32(nfuncs)
    tol = cp.float64(0.001)

    rkf_kernel((1,), (1,), (y_in, x0, h, nsteps, nfuncs, tol, y_out))

    return cp.asnumpy(y_out)
import PDES 


disc_n = 10
PDE1 = PDES.PDE.PDE('dC/dt = -0.01*(exp(-15400/(8.34*T)) * C**0.524)',  ['C'], ['x','y'],[disc_n,disc_n], [(0,1),(0,1)], 'x*y')
PDE2 = PDES.PDE.PDE('dT/dt = -0.1*d2T/dx2 - 0.2*d2T/dy2',  ['T'], ['x','y'],[disc_n,disc_n], [(0,1),(0,1)], 'sin(x)*cos(y)')

PDES1 = PDES.PDES([PDE1, PDE2], ['x','y'], ['T'])

print(PDES1.ic)
resultado = PDES1.df([disc_n,disc_n], inlet=[[0 for i in range(disc_n)]], method="foward")

# print(SERKF45.SERKF45(resultado[0], ['t'], resultado[1], initial_values, 0, 0.1, 10,2,len(PDES1.sp_vars)))
testar = PDES.SERKF45_cuda(resultado[0], ['t'], resultado[1], PDES1.ic, 0, 1, 100, 1, len(PDES1.sp_vars))
print(testar[1][0][-1])
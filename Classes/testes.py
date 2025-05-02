import PDE
import PDES

PDE1 = PDE.PDE('dT/dt = 0.1*d2T/dx2 ', ['t', 'x'], ['T'])
# PDE2 = PDE.PDE('u_t = u_xx + v_yy', ['t', 'x', 'y'], ['u','v'])
PDES1 = PDES.PDES([PDE1], ['x'], ['T'])

# print(PDES1.eqs)
# print(PDES1.ivars)
# print(PDES1.funcs)
print(PDES1.df([4]))



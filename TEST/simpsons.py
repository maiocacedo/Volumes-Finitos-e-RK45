# Implementação do metodo de 3/8 de simpson

"""

1. Ler equações
2. Implementar a integral
3. Implementar a lógica do metodo

"""
#import numpy as np 
import sympy as sym 

x = sym.Symbol("x")   
y = sym.Symbol("y")

def simp(y0, ini, fim, expr):
    h = ((fim-ini)/3)
    x0 = ini
    x1 = x0+h
    x2 = x1+h
    x3 = x2+h
    
    nsubi = 3 
    subi = 1
    
    
    f0 = expr.subs([(x,x0), (y,y0)]) 
    y1 = y0 + h*f0
    print("y1: ", y1)
    
    f1 = expr.subs([(x,x1), (y,y1)])
    y2 = y1 + h*f1
    #print("y2: ", y2)
    
    f2 = expr.subs([(x,x2), (y,y2)])
    y3 = y2 + h*f2
    #print("y3 por taylor: ",y3)
    
    f3 = expr.subs([(x,x3), (y,y3)])
    
    y3 = sym.solve(-y+ y0 + (3/8)*h*(f0 + 3*(f1+f2) + expr.subs(x,x3)))[0]
    #print("y3 por simpson: ", y3)
    return y3

#simp(2,1,1.3, (y*y)+1)    
    
def EDO(x0,xn,y0,n,expr):
    print(y0)
    h0=((xn-x0)/n)
    for i in range(n):
        y0 = simp(y0,x0,x0+h0,expr)
        x0 = x0+h0
        print(y0)
    
EDO(1,1.3,2,600,(y*y)+1)

    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
import numpy as np    
    
def int_simp(nsubi):
    x = 0
    j=1
    nsubi 
    subi = 1
    soma = 0
    
    
    for i in np.linspace(1.2/nsubi, 1.2, nsubi):
        x=i
        k=3
        if (j == 3):
            k=2
            j=0
        expr = k*(x**2 + x**3)
        soma +=expr
        print("sub-intervalo",subi,":", expr)
        j+=1
        subi+=1

    print("Soma:", soma)
    print(3/8*((1.2-1.2/90)/90)*soma) 
    
import sympy as sp
from sympy.parsing.sympy_parser import parse_expr

# Função para receber todos os simbolos e retornar uma lista com eles
def symbol_references(in_list):
    slist = []

    for e in in_list:
        globals()[e] = sp.Symbol(e)
        slist.append(e)
    return slist

# Função para resolver um sistema de EDOs utilizando o método de Runge-Kutta-Fehlberg de ordem 4 e 5
def SERKF45(oldexpr, ivar, funcs, yn, x0, xn, n):
    olddvar = symbol_references(funcs)
    oldivar = symbol_references(ivar)
    expr = [parse_expr(oldexpr[i]) for i in range(len(oldexpr))]
        
 
    allvar = list()
    allvar.append(oldivar)
    print(f'expr:{expr}')
    for i in range(len(olddvar)):
        allvar.append(olddvar[i])
    print(f'allvar{allvar}')
        
    # Inicialização do h e s como vetores de zeros
    h = sp.zeros(len(expr))
    s = sp.zeros(len(expr))
    
    # Atribuição de valores para h
    for j in range(len(expr)):
        h[j] = (xn - x0)/n
    
    # Cópias de yn para yn4 e yn5
    yn4 = yn.copy()
    yn5 = yn.copy()
    tol = 0.001 # Definindo a tolerância
    
    # Inicialização dos vetores k1, k2, k3, k4, k5 e k6 como vetores de zeros
    k1 = sp.zeros(len(expr))
    k2 = sp.zeros(len(expr))
    k3 = sp.zeros(len(expr))
    k4 = sp.zeros(len(expr))
    k5 = sp.zeros(len(expr))
    k6 = sp.zeros(len(expr))

    # print(expr[0])
    
    # Percorrendo as expressões e calculando os valores de k1, k2, k3, k4, k5 e k6 no intervalo [x0, xn]
    checkh = False
    for i in range(n):
        if (checkh == True):
            i = i - 1
            
        checkh = False
        
        for j in range(len(expr)):
            
            # Calculando k1
            set1 = {}
            set1.update({ivar[0]:x0 + i*h[j]}) # Atualizando o valor de x
            for k in range(1, len(allvar)):
                set1.update({allvar[k]:yn[k-1]}) # Atualizando os valores de y
             
            k1[j] = h[j]*expr[j].subs(set1).evalf() # Calculando k1 subtituindo os valores de x e y 

            
            # Calculando k2
            set2 = {}
            set2.update({ivar[0]:x0 + i*h[j] +(1/4)*h[j]}) # Atualizando o valor de x
            for k in range(1, len(allvar)):
                set2.update({allvar[k]:yn[k-1]+ (1/4)*k1[j]}) # Atualizando os valores de y

            k2[j] = h[j]*expr[j].subs(set2).evalf() # Calculando k2 subtituindo os valores de x e y
            
            
            # Calculando k3
            set3 = {}
            set3.update({ivar[0]:x0 + i*h[j] + (3/8)*h[j]}) # Atualizando o valor de x
            for k in range(1, len(allvar)):
                set3.update({allvar[k]:yn[k-1]+(3/32)*k1[j] + (9/32)*k2[j]}) # Atualizando os valores de y

            k3[j] = h[j]*expr[j].subs(set3).evalf() # Calculando k3 subtituindo os valores de x e y
            
            
            # Calculando k4
            set4 = {}
            set4.update({ivar[0]:x0 + i*h[j] + (12/13)*h[j]}) # Atualizando o valor de x
            for k in range(1, len(allvar)):
                set4.update({allvar[k]:yn[k-1] + (1932/2197)*k1[j] - (7200/2197)*k2[j] + (7296/2197)*k3[j]}) # Atualizando os valores de y 
    
            k4[j] = h[j]*expr[j].subs(set4).evalf() # Calculando k4 subtituindo os valores de x e y
            
            
            # Calculando k5
            set5 = {}
            set5.update({ivar[0]:x0 + i*h[j] + h[j]}) # Atualizando o valor de x
            for k in range(1, len(allvar)):
                set5.update({allvar[k]:yn[k-1]+(439/216)*k1[j] - 8*k2[j] + (3680/513)*k3[j] - (845/4104)*k4[j]}) # Atualizando os valores de y
 
            k5[j] = h[j]*expr[j].subs(set5).evalf() # Calculando k5 subtituindo os valores de x e y
            
            
            # Calculando k6
            set6 = {}
            set6.update({ivar[0]:x0 + i*h[j] + (1/2)*h[j]}) # Atualizando o valor de x
            for k in range(1, len(allvar)):
                set6.update({allvar[k]:yn[k-1]- (8/27)*k1[j] + 2*k2[j]-(3544/2565)*k3[j] + (1859/4104)*k4[j]-(11/40)*k5[j]}) # Atualizando os valores de y
  
            k6[j] = h[j]*expr[j].subs(set6).evalf() # Calculando k6 subtituindo os valores de x e y
            
            # Calculando yn4 e yn5 utilizando os valores de k1, k2, k3, k4, k5 e k6
            yn5[j] = yn[j] + (16/135)*k1[j] + (6656/12825)*k3[j] + (28561/56430)*k4[j] - (9/50)*k5[j] + (2/55)*k6[j]
            yn4[j] = yn[j] + (25/216)*k1[j] +  (1408/2565)*k3[j] + (2197/4104)*k4[j] - (1/5)*k5[j]
        

            # Verficando se a necessidade de atualizar o valor de h de acordo com a tolerância definida
            diffy = abs(yn5[j]-yn4[j])
            s[j]=0.84*(tol*h[j]/diffy)**0.25
            if diffy > tol:    
                h[j] = s[j]*h[j]
                checkh = True
        
        yn = yn4.copy()
        print(yn)
    return yn

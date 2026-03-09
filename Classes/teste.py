import numpy as np

disc_n = 5
resultado_analitico = []
for i in range(disc_n):
    t = i / (disc_n - 1)
    lista_analitico = []
    for j in range(disc_n):
        x_ = j / (disc_n - 1)
        
        F_analitico = float(x_ + np.exp(t * x_))  # t=1
        # G_analitico = t + np.cosh(x_ - y_)
        resultado_analitico.append(F_analitico)
        # resultado_analitico.append(G_analitico)
        lista_analitico.append(F_analitico)
    print(f"t={t}: {lista_analitico}")

print("Resultado Analítico:")
print(resultado_analitico)

r = (1/(disc_n-1))/(2*(1/(disc_n - 1))**2)
dt = 1/(disc_n-1)
a1 = [[1 + 2*r - dt/2, - r, 0], [-r , 1 + 2*r - dt/2, -r], [0, -r, 1 + 2*r - dt/2]] 
b1 = [[1 - 2*r + dt/2, r, 0], [r , 1 - 2*r + dt/2, r], [0, r, 1 - 2*r + dt/2]]
print(f"r = {r}")

print("Matriz A:")
print(np.array(a1))

print("Matriz B:")
print(np.array(b1))

print("termo-fonte:")
fonte = [[4.0561669], [0.124452],[8.765397]]
inicial = [[1.25],[1.5], [1.75]]
direita = np.array(b1)@np.array(inicial) + np.array(fonte)

print("Direita:")
print(np.array(direita))
inversa_a1 = np.linalg.inv(a1)

resultado = inversa_a1 @ direita

print("Resultado:")
print(resultado)

#! entender o termo fonte para dirichlet = F no CN.
from RKF45 import SERKF45

# Função para gerar variáveis dependentes da equação diferencial parcial (EDP)
def xs(vars):
    nvars = vars.copy()
    for i in range(len(nvars)): nvars[i] = f'XX{i}'
    return nvars

def df(eq, dvar, discvar, npart):
    # Gerando lista das variáveis dependentes
    xdvar = xs(dvar)
    eqd = eq.split('=')[1]
    allvar = ''
    
    # Gerando lista com todas as variáveis de discretização
    for i in range(len(discvar)):
        allvar = allvar + discvar[i]
    
    # Adicionando as variáveis de discretização como indices das variáveis dependentes, para que possam ser identificadas
    for i in range(len(xdvar)):
        eqd = eqd.replace(f'{dvar[i]}', f'{xdvar[i]}{allvar}')
    
    # Substituindo as derivadas parciais por diferenças finitas adiantadas
    for k in range(len(allvar)):
        for i in range(len(xdvar)): 
            eqd = eqd.replace(f'd2{xdvar[i]}{allvar}/d{allvar[k]}2', f'({xdvar[i]}i+1j - 2*{xdvar[i]}ij + {xdvar[i]}i-1j)/ h{allvar[k]} ** 2')
            eqd = eqd.replace(f'd{xdvar[i]}{allvar}/d{allvar[k]}', f'({xdvar[i]}i+1j - {xdvar[i]}ij)/ h{allvar[k]}')
            eqd = eqd.replace(f'{xdvar[i]}{allvar}', f'{xdvar[i]}ij')
    
    print(eqd)
    print(allvar)
    # Substituindo o indice i por seu valor correspondente
    lista = [eqd.replace('i+1', str(i+1)).replace('i-1', str(i-1)).replace('i',str(i)) for i in range(1,npart[0])]
    lista2 = []
    for i in range(5):
        lista2.append('0')

    # Substituindo o indice j por seu valor correspondente
    for i in range(len(lista)):
        for j in range(0,npart[1]+1):
            lista2.append(lista[i].replace('j+1', str(j+1)).replace('j-1', str(j-1)).replace('j',str(j)))
            
    dvars = []
    
    # Gerando lista de variáveis dependentes
    for l in range(len(lista2)):
        for k in range(len(dvar)):
            for i in range(0,npart[0]+1):
                for j in range(0,npart[1]+1):
                    if not(f'XX{k}{i}{j}' in dvars):
                        dvars.append(f'XX{k}{i}{j}')
    print('dvars')
    print(dvars)
    print(f'tamanho dvars: {len(dvars)}')
    
    # Adicionando os valores de h na lista de equações
    for i in range(len(lista2)):
        for j in range(len(npart)):
            lista2[i] = lista2[i].replace(f'h{allvar[j]}', str(1/npart[j]))
    
    # Adicionando as equações da fronteira
    for i in range(0,npart[1]+1):
        lista2.append(f'1.25*XX0{npart[0]}{i} - 1.25*XX0{npart[0]-1}{i}')

    
    
    print(f'Tamanho lista 2:{len(lista2)}') 
    
    return lista2, dvars

# Gerando os valores iniciais
initial_values = []
for i in range(0, 5):
    initial_values.append(100)
    
for i in range (0, 20):
    initial_values.append(80)
    
print(initial_values)

resultadodisc = df('dT/dt = 0.1*d2T/dx2 ', ['T'], ['x','y'], [4,4])     

print(resultadodisc[0], resultadodisc[1])
SERKF45(resultadodisc[0], ['t'], resultadodisc[1], initial_values, 0, 0.1, 10)   

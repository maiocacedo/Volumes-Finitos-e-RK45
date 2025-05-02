def xs(vars):
    nvars = vars.copy()
    for i in range(len(nvars)): nvars[i] = f'XX{i}'
    return nvars

class Disc:
    def df(eq, dvar, discvar, npart):
    
        # Gerando lista das variáveis dependentes
        xdvar = xs(dvar)
        eqd = eq.split('=')[1]
        allvar = ''
        print(f'xdvar: {xdvar}')
        print(f'dvar: {dvar}')
        # Gerando lista com todas as variáveis de discretização
        for i in range(len(discvar)):
            allvar = allvar + discvar[i]
        
        # Adicionando as variáveis de discretização como indices das variáveis dependentes, para que possam ser identificadas
        for i in range(len(xdvar)):
            eqd = eqd.replace(f'{dvar[i]}', f'{xdvar[i]}{allvar}')
        
        # Substituindo as derivadas parciais por diferenças finitas adiantadas
        for k in range(len(allvar)):
            if k == 0:
                for i in range(len(xdvar)): 
                    eqd = eqd.replace(f'd2{xdvar[i]}{allvar}/d{allvar[k]}2', f'({xdvar[i]}i+1j - 2*{xdvar[i]}ij + {xdvar[i]}i-1j)/ h{allvar[k]} ** 2')
                    eqd = eqd.replace(f'd{xdvar[i]}{allvar}/d{allvar[k]}', f'({xdvar[i]}i+1j - {xdvar[i]}ij)/ h{allvar[k]}')
                    eqd = eqd.replace(f'{xdvar[i]}{allvar}', f'{xdvar[i]}ij')
            elif k == 1:
                for i in range(len(xdvar)): 
                    eqd = eqd.replace(f'd2{xdvar[i]}{allvar}/d{allvar[k]}2', f'({xdvar[i]}ij+1 - 2*{xdvar[i]}ij + {xdvar[i]}ij-1)/ h{allvar[k]} ** 2')
                    eqd = eqd.replace(f'd{xdvar[i]}{allvar}/d{allvar[k]}', f'({xdvar[i]}ij+1 - {xdvar[i]}ij)/ h{allvar[k]}')
                    eqd = eqd.replace(f'{xdvar[i]}{allvar}', f'{xdvar[i]}ij')
            
        print(eqd)
        print(allvar)
        
        # Substituindo o indice i por seu valor correspondente
        lista = [eqd.replace('i+1', str(i+1)).replace('i-1', str(i-1)).replace('i',str(i)) for i in range(1,npart[0])]
        lista2 = []
        if len(allvar) == 2:
            for i in range(npart[1]):
                lista2.append('0')
                
        lista2.append('0')
        
        # Substituindo o indice j por seu valor correspondente
        if len(allvar) == 2:
            for i in range(len(lista)):
                for j in range(0,npart[1]+1):
                    lista2.append(lista[i].replace('j+1', str(j+1)).replace('j-1', str(j-1)).replace('j',str(j)))
        
        elif len(allvar) == 1:
            for i in range(len(lista)):
                lista2.append(lista[i].replace('j', str(0)))
                
        dvars = []
        
        # Gerando lista de variáveis dependentes
        for l in range(len(lista2)):
            for k in range(len(dvar)):
                for i in range(0,npart[0]+1):
                    if len(allvar) == 2:
                        for j in range(0,npart[1]+1):
                            if not(f'XX{k}{i}{j}' in dvars):
                                dvars.append(f'XX{k}{i}{j}')
                    elif len(allvar) == 1:
                        if not(f'XX{k}{i}0' in dvars):
                            dvars.append(f'XX{k}{i}0')
        print('dvars')
        print(dvars)
        print(f'tamanho dvars: {len(dvars)}')
        
        
        # Adicionando os valores de h na lista de equações
        if len(allvar) == 2:
            # Adicionando as equações da fronteira
            for i in range(0,npart[1]+1):
                for j in range(len(xdvar)):
                    lista2.append(f'({xdvar[j]}{npart[0]}{i} - {xdvar[j]}{npart[0]-1}{i})/h{allvar[1]}**2')
                
        elif len(allvar) == 1:
            # Adicionando as equações da fronteira
            for j in range(len(xdvar)):
                lista2.append(f'({xdvar[j]}{npart[0]}0 - {xdvar[j]}{npart[0]-1}0)/h{allvar[0]}**2')
            
            
        for i in range(len(lista2)):
            for j in range(len(allvar)):
                lista2[i] = lista2[i].replace(f'h{allvar[j]}', str(1/npart[j]))
        
        print('Equações discretizadas:')
        for i in range(len(lista2)):
            print(lista2[i])
        print(f'Tamanho lista 2:{len(lista2)}') 
        
        return lista2, dvars


    def vf():
        return
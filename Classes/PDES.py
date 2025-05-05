import PDE


#*                 IDEIAS:
        
#*                 variavel como XX0i0000....? Dessa forma, substitui de acordo com o número de variáveis de discretização.
#*                 Então, se temos ij: XX0i0. Se temos ijk: XX0i00. Se temos ijkl: XX0i000.
#*                 Ou seja, o número de zeros é igual ao número de variáveis de discretização - 1.
#*                 Generalizando o número de variáveis de discretização.
        
#*                Trate lista de eq como uma matriz.
        
#*                


class PDES:
    eqs = []
    ivars = []
    sp_vars = []
    funcs = []
    n_sp = 1
    n_temp = 1
    def __init__(self, pdes, sp_vars, n_sp=1, n_temp=1):
        self.eqs = [pde.eq for pde in pdes]
        self.ivars = list([pde.ivar[i] for pde in pdes for i in range(len(pde.ivar))]) # Tirei o set.
        self.funcs = list([pde.func[i] for pde in pdes for i in range(len(pde.func))])
        self.sp_vars = sp_vars
        self.n_sp = n_sp
        self.n_temp = n_temp
        
        pass
    
    
    def xs(self,vars):
        nvars = vars.copy()
        for i in range(len(nvars)): nvars[i] = f'XX{i}'
        return nvars


    def df(self, npart):
        
        # Gerando lista das variáveis dependentes
        xdvar = self.xs(self.funcs)

        eqrs = [eq.split('=')[1] for eq in self.eqs]
        str_sp_vars = ''

        
        # Gerando lista com todas as variáveis de discretização
        for i in range(len(self.sp_vars)):
            str_sp_vars = str_sp_vars + self.sp_vars[i]
        
       
        # Adicionando as variáveis de discretização como indices das variáveis dependentes, para que possam ser identificadas
        for j in range(len(eqrs)):
            for i in range(len(xdvar)):
                eqrs[j] = eqrs[j].replace(f'{str(self.funcs[i])}', f'{xdvar[i]}{str_sp_vars}')

        
      
        
        for j in range(len(eqrs)):
            # Substituindo as derivadas parciais por diferenças finitas adiantadas
            for k in range(len(str_sp_vars)):
                if k == 0:
                    for i in range(len(xdvar)): 
                        eqrs[j] = eqrs[j].replace(f'd2{xdvar[i]}{str_sp_vars}/d{str_sp_vars[k]}2', f'({xdvar[i]}i+1j - 2*{xdvar[i]}ij + {xdvar[i]}i-1j)/ h{str_sp_vars[k]} ** 2')
                        eqrs[j] = eqrs[j].replace(f'd{xdvar[i]}{str_sp_vars}/d{str_sp_vars[k]}', f'({xdvar[i]}i+1j - {xdvar[i]}ij)/ h{str_sp_vars[k]}')
                        eqrs[j] = eqrs[j].replace(f'{xdvar[i]}{str_sp_vars}', f'{xdvar[i]}ij')
                elif k == 1:
                    for i in range(len(xdvar)): 
                        eqrs[j] = eqrs[j].replace(f'd2{xdvar[i]}{str_sp_vars}/d{str_sp_vars[k]}2', f'({xdvar[i]}ij+1 - 2*{xdvar[i]}ij + {xdvar[i]}ij-1)/ h{str_sp_vars[k]} ** 2')
                        eqrs[j] = eqrs[j].replace(f'd{xdvar[i]}{str_sp_vars}/d{str_sp_vars[k]}', f'({xdvar[i]}ij+1 - {xdvar[i]}ij)/ h{str_sp_vars[k]}')
                        eqrs[j] = eqrs[j].replace(f'{xdvar[i]}{str_sp_vars}', f'{xdvar[i]}ij')
                
        print(f'\nEquações: {eqrs}\n')
        
        
        #* A partir daqui: Substituição dos indices i e j por seus valores correspondentes em eqrs:
        
        #*                 [[centros eq1][centros de eq2]]
        #*                 [W0, W1, W2, W3, S1, ce0, ce1, N1, S2, ce2, ce3, N2, E0, E1, E2, E3]
        
        #*                 Generalizando para 2 variáveis de discretização.
        
        
        
        lista = []
        

        
        #* South-Center-North
        
        for j in range(len(eqrs)):
            lista.append([eqrs[j].replace('i+1', str(i+1)).replace('i-1', str(i-1)).replace('i',str(i)) for i in range(1,npart[0]-1)])
            
        print(f'Lista: {lista}')
 
        lista2 = []
        lista2 = [[] for i in range(len(lista))]
        
        #* West
        
        if len(str_sp_vars) == 1:
            for j in range(len(lista2)):
                lista2[j].append('W0')
            
            for j in range(len(lista2)):
                for k in range(1,npart[0]):
                    lista2[j].append("S" + str(k))
                    for i in range(0, npart[1]-2):
                        lista2[j].append('Ce' + str(i))
                    lista2[j].append('N' + str(k))
                    
            for j in range(len(lista2)):
                lista2[j].append('E0')
                    
        
        if len(str_sp_vars) == 2:
            for j in range(len(lista2)):
                l = 0
                
                for i in range(npart[1]):
                    lista2[j].append('W' + str(i))
                    
                for k in range(1,npart[0]-1):
                    lista2[j].append("S" + str(k))
                    for i in range(l, l + npart[1]-2):
                        lista2[j].append('Ce' + str(i))
                    lista2[j].append('N' + str(k))
                    l += npart[1]-2
                
                for i in range(npart[1]):    
                    lista2[j].append('E' + str(i))
                
        # if len(str_sp_vars) == 2:
        #     for j in range(len(lista)):
        #         for i in range(len(lista[j])):
        #             for k in range(npart[1]):
        #                 lista2[j].append(lista[j][i].replace('j+1', str(k+1)).replace('j-1', str(k-1)).replace('j',str(k)))
        #     for i in range(len(lista2)):
        #         print(f'Lista2[{i}]: {lista2[i]}')
        
        
        # elif len(str_sp_vars) == 1:
        #     for j in range(len(lista)):
        #         for i in range(len(lista2[j])):
        #             lista2[j][i] = lista2[j][i].replace('j', str(0))
                
        dvars = []
        
        # # Gerando lista de variáveis dependentes
        # for l in range(len(lista2)):
        #     for k in range(len(self.funcs)):
        #         for i in range(0,npart[0]+1):
        #             if len(str_sp_vars) == 2:
        #                 for j in range(0,npart[1]+1):
        #                     if not(f'XX{k}{i}{j}' in dvars):
        #                         dvars.append(f'XX{k}{i}{j}')
        #             elif len(str_sp_vars) == 1:
        #                 if not(f'XX{k}{i}0' in dvars):
        #                     dvars.append(f'XX{k}{i}0')
        # print('dvars')
        # print(dvars)
        # print(f'tamanho dvars: {len(dvars)}')
        
        
        # * East and h
        listeast = [[] for i in range(len(lista2))]
        if len(str_sp_vars) == 2:
            for i in range(len(lista2)):
                for j in range(npart[i], 0, -1):
                    listeast[i].append(lista2[i][-j])
                    print(j)
            
            for i in range(len(listeast)):
                for j in range(len(listeast[i])):
                    lista2[i].append(listeast[i][j])    
                    
                
        
        elif len(str_sp_vars) == 1:
            for j in range(len(lista2)):
                lista2[j].append(lista2[j][-1])
        
        # if len(str_sp_vars) == 2:
        #     # Adicionando as equações da fronteira
        #     for i in range(0,npart[1]+1):
        #         for j in range(len(xdvar)):
        #             lista2.append(f'({xdvar[j]}{npart[0]}{i} - {xdvar[j]}{npart[0]-1}{i})/h{str_sp_vars[1]}**2')
                
        # elif len(str_sp_vars) == 1:
        #     # Adicionando as equações da fronteira
        #     for j in range(len(xdvar)):
        #         lista2.append(f'({xdvar[j]}{npart[0]}0 - {xdvar[j]}{npart[0]-1}0)/h{str_sp_vars[0]}**2')
            
            
        # for i in range(len(lista2)):
        #     for j in range(len(str_sp_vars)):
        #         lista2[i] = lista2[i].replace(f'h{str_sp_vars[j]}', str(1/npart[j]))
        
        # print('Equações discretizadas:')
        # for i in range(len(lista2)):
        #     print(lista2[i])
        # print(f'Tamanho lista 2:{len(lista2)}') 
        
        for i in range(len(lista2)):
            for j in range(len(lista2[i])):
                print(f'Lista2[{i}]: {lista2[i][j]}')
        print(f'Tamanho lista 2: {len(lista2[0])}')
        print(f'Tamanho lista 2: {len(lista2[1])}')
        return lista2, dvars

PDE1 = PDE.PDE('dC/dt = 0.1*d2C/dx2 + T', ['t', 'x'], ['C'])
PDE2 = PDE.PDE('dT/dt = 0.1*d2T/dx2 ', ['t', 'x'], ['T'])


# PDE2 = PDE.PDE('u_t = u_xx + v_yy', ['t', 'x', 'y'], ['u','v'])
PDES1 = PDES([PDE2, PDE1], ['x', 'y'], ['T', 'C'])

# print(PDES1.eqs)
# print(PDES1.ivars)
# print(PDES1.funcs)
print(PDES1.df([4,4])) 

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
        self.ivars = list([pde.ivar[i] for pde in pdes for i in range(len(pde.ivar))])
        self.funcs = list([pde.func[i] for pde in pdes for i in range(len(pde.func))])
        self.sp_vars = sp_vars
        self.n_sp = n_sp
        self.n_temp = n_temp
        
        pass
    
    
    def xs(self,vars):
        nvars = vars.copy()
        for i in range(len(nvars)): nvars[i] = f'XX{i}'
        return nvars


    def df(self, n_part, list_east = []):
       
        # Gerando lista das variáveis dependentes
        xd_var = self.xs(self.funcs)

        eqrs = [eq.split('=')[1] for eq in self.eqs]
        str_sp_vars = ''

        
        # Gerando lista com todas as variáveis de discretização
        for i in range(len(self.sp_vars)):
            str_sp_vars = str_sp_vars + self.sp_vars[i]
        
       
        # Adicionando as variáveis de discretização como indices das variáveis dependentes, para que possam ser identificadas
        for j in range(len(eqrs)):
            for i in range(len(xd_var)):
                eqrs[j] = eqrs[j].replace(f'{str(self.funcs[i])}', f'{xd_var[i]}{str_sp_vars}')

        
      
        
        for j in range(len(eqrs)):
            # Substituindo as derivadas parciais por diferenças finitas adiantadas
            for k in range(len(str_sp_vars)):
                if k == 0: # se a variável de discretização for a primeira
                    for i in range(len(xd_var)): 
                        eqrs[j] = eqrs[j].replace(f'd2{xd_var[i]}{str_sp_vars}/d{str_sp_vars[k]}2', f'({xd_var[i]}i+1j - 2*{xd_var[i]}ij + {xd_var[i]}i-1j)/ h{str_sp_vars[k]} ** 2') 
                        eqrs[j] = eqrs[j].replace(f'd{xd_var[i]}{str_sp_vars}/d{str_sp_vars[k]}', f'({xd_var[i]}i+1j - {xd_var[i]}ij)/ h{str_sp_vars[k]}')
                        eqrs[j] = eqrs[j].replace(f'{xd_var[i]}{str_sp_vars}', f'{xd_var[i]}ij')
                elif k == 1: # se a variável de discretização for a segunda
                    for i in range(len(xd_var)): 
                        eqrs[j] = eqrs[j].replace(f'd2{xd_var[i]}{str_sp_vars}/d{str_sp_vars[k]}2', f'({xd_var[i]}ij+1 - 2*{xd_var[i]}ij + {xd_var[i]}ij-1)/ h{str_sp_vars[k]} ** 2')
                        eqrs[j] = eqrs[j].replace(f'd{xd_var[i]}{str_sp_vars}/d{str_sp_vars[k]}', f'({xd_var[i]}ij+1 - {xd_var[i]}ij)/ h{str_sp_vars[k]}')
                        eqrs[j] = eqrs[j].replace(f'{xd_var[i]}{str_sp_vars}', f'{xd_var[i]}ij')
                

        
        
        #* A partir daqui: Substituição dos indices i e j por seus valores correspondentes em eqrs:
        
        #*                 [[centros eq1][centros de eq2]]
        #*                 [W0, W1, W2, W3, S1, ce0, ce1, N1, S2, ce2, ce3, N2, E0, E1, E2, E3]
        
        #*                 Generalizando para 2 variáveis de discretização.
        

        partial_list_eq = []

        # * Gerando lista de equações com indices parcialmente substituídos
        for j in range(len(eqrs)):
            partial_list_eq.append([eqrs[j].replace('i+1', str(i+1)).replace('i-1', str(i-1)).replace('i',str(i)) for i in range(1,n_part[0]-1)])
            

        list_eq = []
        list_eq = [[] for i in range(len(partial_list_eq))]
        # * Gerando lista de equações com indices totalmente substituídos para duas variáveis de discretização
        if len(str_sp_vars) == 2:
            for j in range(len(partial_list_eq)): # para cada lista de equações discretizadas
                for i in range(len(partial_list_eq[j])): # para cada equação da lista
                    for k in range(1,n_part[1]-1): # indice j da segunda variável de discretização
                        list_eq[j].append(partial_list_eq[j][i].replace('j+1', str(k+1)).replace('j-1', str(k-1)).replace('j',str(k)))
          
        
        # * Gerando lista de equações com indices totalmente substituídos para uma única variável de discretização
        elif len(str_sp_vars) == 1:
            for j in range(len(partial_list_eq)): # para cada lista de equações discretizadas
                for i in range(len(list_eq[j])): # para cada equação da lista
                    list_eq[j][i] = list_eq[j][i].replace('j', str(0))
      
        # * Gerando Lista de Posições para duas variáveis de discretização
        list_positions = []
        
        if len(str_sp_vars) == 2:
            for func in range(len(self.funcs)):
                list_aux = []
                for i in range(n_part[0]):
                    for j in range(n_part[1]):
                        # Se a posição for da primeira parede do domínio (West)
                        if i == 0:
                            list_aux.append(f'W{func}{i}{j}') 
                        # Se a posição for da última parede do domínio (East)
                        elif i == n_part[0]-1:
                            list_aux.append(f'E{func}{i}{j}')
                        # Se a posição for da parede de baixo do domínio (South)
                        elif j == 0:
                            list_aux.append(f'S{func}{i}{j}')
                        # Se a posição for da parede de cima do domínio (North)
                        elif j == n_part[1]-1:
                            list_aux.append(f'N{func}{i}{j}')
                        # Se a posição for do centro do domínio (Center)
                        else:
                            list_aux.append(f'Ce{func}{i}{j}')
                list_positions.append(list_aux)    
                
        # * Gerando lista de posições para uma única variável de discretização
        elif len(str_sp_vars) == 1:    
            for func in range(len(self.funcs)):
                list_aux = []
                # Se a posição for da primeira parede do domínio (West)
                list_aux.append(f'W{func}00')
                
                # Se a posição for do centro do domínio (Center)
                for i in range(1,n_part[0]-1):
                    list_aux.append(f'Ce{func}{i}0')
                    
                # Se a posição for da última parede do domínio (East)
                list_aux.append(f'E{func}{n_part[0]-1}0')
                
                list_positions.append(list_aux)    
                

        
        # * Substituindo centros na lista de posições
        if len(str_sp_vars) == 2:
            for ii in range(len(list_positions)): # para cada lista de posições
                for jj in range(len(list_positions[ii])): # para cada posição da lista
                    for i in range(1,n_part[0]-1): # para cada índice i da primeira variável de discretização
                        for j in range(1,n_part[1]-1): # para cada índice j da segunda variável de discretização
                            # Substituindo centros na lista de posições por suas respectivas equações
                            list_positions[ii][jj] = list_positions[ii][jj].replace(f'Ce{ii}{i}{j}', list_eq[ii][(n_part[0]-2)*(i-1)+j-1])
            
            # * Gerando lista norte e sul
            list_north = [[] for i in range(len(list_eq))]
            list_south = [[] for i in range(len(list_eq))]
            for i in range(len(list_eq)): # para cada lista de equações discretizadas
                for j in range(0, len(list_eq[i])): # para cada equação da lista
                    if j % (n_part[1]-2) == 0: # se o índice j for múltiplo de (n_part[1]-2)
                        list_south[i].append(list_eq[i][j])
                        list_north[i].append(list_eq[i][j+n_part[1]-3])
                    

            
            # * Gerando lista leste
            if list_east == []:
                list_east = [[] for i in range(len(list_eq))]
                for ii in range(len(list_eq)):
                    list_east[ii].append(list_south[ii][-1])
                    for i in range((n_part[1]-2)*(n_part[0]-2)-(n_part[1]-2), (n_part[1]-2)*(n_part[0]-2)): 
                        list_east[ii].append(list_eq[ii][i])
                    list_east[ii].append(list_north[ii][-1])        
            
        
        if len(str_sp_vars) == 1:
            # * Substituindo centros na lista de posições
            for ii in range(len(list_positions)):
                for jj in range(len(list_positions[ii])):
                    for i in range(1,n_part[0]-1):
                        list_positions[ii][jj] = list_positions[ii][jj].replace(f'Ce{ii}{i}0', list_eq[ii][(n_part[0]-2)*(i-1)])

            # * Gerando lista leste
            if list_east == []: # Verifica se a lista leste já foi criada pelo usuáriow
                list_east = [[] for i in range(len(list_eq))]
                for ii in range(len(list_eq)):
                    list_east[ii].append(list_eq[ii][n_part[0]])
        
        # * Substituindo as posições por equações
        final_list_eq = []
        if len(str_sp_vars) == 2:
            for func in range(len(list_positions)):
                S = 0
                N = 0
                E = 0
                for len_list in range(len(list_positions[func])):
                    if "S" in list_positions[func][len_list]:
                        list_positions[func][len_list] = list_south[func][S]
                        S += 1
                    elif "N" in list_positions[func][len_list]:
                        list_positions[func][len_list] = list_north[func][N]
                        N += 1
                    elif "E" in list_positions[func][len_list]:
                        list_positions[func][len_list] = list_east[func][E]
                        E += 1
                                                                                    
        d_vars = []
    
        # Gerando lista de variáveis dependentes
        for l in range(len(list_positions)):
            for k in range(len(self.sp_vars)):
                for i in range(0,n_part[0]+1):
                    if len(self.sp_vars) == 2:
                        for j in range(0,n_part[1]+1):
                            if not(f'XX{k}{i}{j}' in d_vars):
                                d_vars.append(f'XX{k}{i}{j}')
                    elif len(str_sp_vars) == 1:
                        if not(f'XX{k}{i}0' in d_vars):
                            d_vars.append(f'XX{k}{i}0')
            

        flat_list_positions = []
        
        for list in list_positions:
            flat_list_positions.extend(list)
            
        print(f'flat_list_positions: {flat_list_positions}')
        
        return list_positions, d_vars

PDE1 = PDE.PDE('dC/dt = 0.1*d2C/dx2 + T', ['t', 'x'], ['C'])
PDE2 = PDE.PDE('dT/dt = 0.1*d2T/dx2 ', ['t', 'x'], ['T'])

PDES1 = PDES([PDE2, PDE1], ['x', 'y'], ['T', 'C'])

print(PDES1.df([4,4]))

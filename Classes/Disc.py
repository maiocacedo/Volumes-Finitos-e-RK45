def df(pdes, n_part, inlet = [], method="forward", north_bd = "neumann", south_bd = "neumann", east_bd = "neumann",
       north_func_bd = "0", south_func_bd = "0", east_func_bd = "0", north_alpha_bd = "0", south_alpha_bd = "0", east_alpha_bd = "0", north_beta_bd = "1", south_beta_bd = "1", east_beta_bd = "1" ):

        # Gerando lista das variáveis dependentes
        xd_var = pdes.xs(pdes.funcs)

        eqrs = [eq.split('=')[1] for eq in pdes.eqs]
        str_sp_vars = ''

        
        # Gerando lista com todas as variáveis de discretização
        for i in range(len(pdes.sp_vars)):
            str_sp_vars = str_sp_vars + pdes.sp_vars[i]
        
       
       
        # Adicionando as variáveis de discretização como indices das variáveis dependentes, para que possam ser identificadas
        for j in range(len(eqrs)):
            for i in range(len(xd_var)):
                eqrs[j] = eqrs[j].replace(f'{str(pdes.funcs[i])}', f'{xd_var[i]}{str_sp_vars}')

        if (method == "forward"):
            for j in range(len(eqrs)):
                # Substituindo as derivadas parciais por diferenças finitas adiantadas
                for k in range(len(str_sp_vars)):
                    if k == 0:  # primeira variável de discretização
                        for i in range(len(xd_var)):
                            eqrs[j] = eqrs[j].replace(
                                f'd2{xd_var[i]}{str_sp_vars}/d{str_sp_vars[k]}2',
                                f'({xd_var[i]}_i+1_j - 2*{xd_var[i]}_ii_j + {xd_var[i]}_i-1_j)/ h{xd_var[i]}_ ** 2'
                            )
                            eqrs[j] = eqrs[j].replace(
                                f'd{xd_var[i]}{str_sp_vars}/d{str_sp_vars[k]}',
                                f'({xd_var[i]}_i+1_j - {xd_var[i]}_ii_j)/ h{xd_var[i]}_'
                            )
                    elif k == 1:  # segunda variável de discretização
                        for i in range(len(xd_var)):
                            eqrs[j] = eqrs[j].replace(
                                f'd2{xd_var[i]}{str_sp_vars}/d{str_sp_vars[k]}2',
                                f'({xd_var[i]}_ii_j+1 - 2*{xd_var[i]}_ii_j + {xd_var[i]}_ii_j-1)/ h{xd_var[i]}_ ** 2'
                            )
                            eqrs[j] = eqrs[j].replace(
                                f'd{xd_var[i]}{str_sp_vars}/d{str_sp_vars[k]}',
                                f'({xd_var[i]}_ii_j+1 - {xd_var[i]}_ii_j)/ h{xd_var[i]}_'
                            )
                            eqrs[j] = eqrs[j].replace(
                                f'{xd_var[i]}{str_sp_vars}', f'{xd_var[i]}_ii_j'
                            )
            for j in range(len(eqrs)):
                eqrs[j] = eqrs[j].replace(f'{str_sp_vars[0]}', f'ii * h{str_sp_vars[0]}_')
                eqrs[j] = eqrs[j].replace(f'{str_sp_vars[1]}', f'j * h{str_sp_vars[1]}_')




        elif (method == "central"):

            for j in range(len(eqrs)):

                # Substituindo as derivadas parciais por diferenças finitas centrais

                for k in range(len(str_sp_vars)):

                    if k == 0:  # primeira variável de discretização

                        for i in range(len(xd_var)):
                            eqrs[j] = eqrs[j].replace(

                                f'd2{xd_var[i]}{str_sp_vars}/d{str_sp_vars[k]}2',

                                f'({xd_var[i]}_i+1_j - 2*{xd_var[i]}_ii_j + {xd_var[i]}_i-1_j)/ h{xd_var[i]}_ ** 2'

                            )

                            # 1ª derivada central correta

                            eqrs[j] = eqrs[j].replace(

                                f'd{xd_var[i]}{str_sp_vars}/d{str_sp_vars[k]}',

                                f'({xd_var[i]}_i+1_j - {xd_var[i]}_i-1_j)/(2* h{xd_var[i]}_)'

                            )

                    elif k == 1:  # segunda variável de discretização

                        for i in range(len(xd_var)):
                            eqrs[j] = eqrs[j].replace(

                                f'd2{xd_var[i]}{str_sp_vars}/d{str_sp_vars[k]}2',

                                f'({xd_var[i]}_ii_j+1 - 2*{xd_var[i]}_ii_j + {xd_var[i]}_ii_j-1)/ h{xd_var[i]}_ ** 2'

                            )

                            # 1ª derivada central correta

                            eqrs[j] = eqrs[j].replace(

                                f'd{xd_var[i]}{str_sp_vars}/d{str_sp_vars[k]}',

                                f'({xd_var[i]}_ii_j+1 - {xd_var[i]}_ii_j-1)/(2* h{xd_var[i]}_)'

                            )

                            eqrs[j] = eqrs[j].replace(

                                f'{xd_var[i]}{str_sp_vars}', f'{xd_var[i]}_ii_j'

                            )

            for j in range(len(eqrs)):
                eqrs[j] = eqrs[j].replace(f'{str_sp_vars[0]}', f'ii * h{xd_var[0]}_')

                eqrs[j] = eqrs[j].replace(f'{str_sp_vars[1]}', f'j * h{xd_var[0]}_')


        elif (method == "backward"):

            for j in range(len(eqrs)):

                # Substituindo as derivadas parciais por diferenças finitas atrasadas

                for k in range(len(str_sp_vars)):

                    if k == 0:  # primeira variável de discretização

                        for i in range(len(xd_var)):
                            eqrs[j] = eqrs[j].replace(

                                f'd2{xd_var[i]}{str_sp_vars}/d{str_sp_vars[k]}2',

                                f'({xd_var[i]}_i+1_j - 2*{xd_var[i]}_ii_j + {xd_var[i]}_i-1_j)/ h{xd_var[i]}_ ** 2'

                            )

                            eqrs[j] = eqrs[j].replace(

                                f'd{xd_var[i]}{str_sp_vars}/d{str_sp_vars[k]}',

                                f'({xd_var[i]}_ii_j - {xd_var[i]}_i-1_j)/ h{xd_var[i]}_'

                            )


                    elif k == 1:  # segunda variável de discretização

                        for i in range(len(xd_var)):
                            eqrs[j] = eqrs[j].replace(

                                f'd2{xd_var[i]}{str_sp_vars}/d{str_sp_vars[k]}2',

                                f'({xd_var[i]}_ii_j+1 - 2*{xd_var[i]}_ii_j + {xd_var[i]}_ii_j-1)/ h{xd_var[i]}_ ** 2'

                            )

                            eqrs[j] = eqrs[j].replace(

                                f'd{xd_var[i]}{str_sp_vars}/d{str_sp_vars[k]}',

                                f'({xd_var[i]}_ii_j - {xd_var[i]}_ii_j-1)/ h{xd_var[i]}_'

                            )

                            eqrs[j] = eqrs[j].replace(

                                f'{xd_var[i]}{str_sp_vars}', f'{xd_var[i]}_ii_j'

                            )

            for j in range(len(eqrs)):
                eqrs[j] = eqrs[j].replace(f'{str_sp_vars[0]}', f'ii * h{str_sp_vars[0]}_')

                eqrs[j] = eqrs[j].replace(f'{str_sp_vars[1]}', f'j * h{str_sp_vars[1]}_')

        else:
            raise ValueError("Método de discretização inválido. Use 'forward', 'central' ou 'backward'.")
        
        #* A partir daqui: Substituição dos indices i e j por seus valores correspondentes em eqrs:
        
        #*                 [[centros eq1][centros de eq2]]
        #*                 [W0, W1, W2, W3, S1, ce0, ce1, N1, S2, ce2, ce3, N2, E0, E1, E2, E3]
        
        #*                 Generalizando para 2 variáveis de discretização.
        

        partial_list_eq = []

        # * Gerando lista de equações com indices parcialmente substituídos
        for j in range(len(eqrs)):
            partial_list_eq.append([eqrs[j].replace('i+1', str(i+1)).replace('i-1', str(i-1)).replace('i-2', str(i-2)).replace('i+2', str(i+2)).replace('ii',str(i)) for i in range(1,n_part[0]-1)])
        

            
        list_eq = []
        list_eq = [[] for i in range(len(partial_list_eq))]
        
        if len(str_sp_vars) == 3:
            pass
        
        # * Gerando lista de equações com indices totalmente substituídos para duas variáveis de discretização
        elif len(str_sp_vars) == 2:
            for j in range(len(partial_list_eq)): # para cada lista de equações discretizadas
                for i in range(len(partial_list_eq[j])): # para cada equação da lista
                    for k in range(1,n_part[1]-1): # indice j da segunda variável de discretização
                        list_eq[j].append(partial_list_eq[j][i].replace('j+1', str(k+1)).replace('j-1', str(k-1)).replace('j-2', str(k-2)).replace('j+2', str(k+2)).replace('j',str(k)))
          
        
        # * Gerando lista de equações com indices totalmente substituídos para uma única variável de discretização
        elif len(str_sp_vars) == 1:
            for j in range(len(partial_list_eq)): # para cada lista de equações discretizadas
                for i in range(len(partial_list_eq[j])): # para cada equação da lista
                    list_eq[j].append(partial_list_eq[j][i].replace('j', str(0)))
        
        # * Gerando Lista de Posições para duas variáveis de discretização
        list_positions = []
        
        if len(str_sp_vars) == 3:
            pass
        
        elif len(str_sp_vars) == 2:
            for func in range(len(pdes.funcs)):
                list_aux = []
                for i in range(n_part[0]):
                    for j in range(n_part[1]):
                        # Se a posição for da primeira parede do domínio (West)
                        if i == 0:
                            list_aux.append(f'W{func}_{i}_{j}') 
                        # Se a posição for da última parede do domínio (East)
                        elif i == n_part[0]-1:
                            list_aux.append(f'E{func}_{i}_{j}')
                        # Se a posição for da parede de baixo do domínio (South)
                        elif j == 0:
                            list_aux.append(f'S{func}_{i}_{j}')
                        # Se a posição for da parede de cima do domínio (North)
                        elif j == n_part[1]-1:
                            list_aux.append(f'N{func}_{i}_{j}')
                        # Se a posição for do centro do domínio (Center)
                        else:
                            list_aux.append(f'Ce{func}_{i}_{j}')
                list_positions.append(list_aux)    
        
            
        # * Gerando lista de posições para uma única variável de discretização
        elif len(str_sp_vars) == 1:    
            for func in range(len(pdes.funcs)):
                list_aux = []
                for i in range(0,n_part[0]):
                    # Se a posição for da primeira parede do domínio (West)
                    if i == 0:
                        list_aux.append(f'W{func}{i}0')
                
                
                    elif i == n_part[0]-1:
                        # Se a posição for da última parede do domínio (East)
                        list_aux.append(f'E{func}{i}0')
                    else:
                        # Se a posição for do centro do domínio (Center)
                        list_aux.append(f'Ce{func}{i}0')
                list_positions.append(list_aux)    
                
        if len(str_sp_vars) == 3:
            print("não está implementado.")
            pass
        
        # * Substituindo centros na lista de posições
        elif len(str_sp_vars) == 2:
        
            for func in range(len(list_positions)):
                C = 0
                for i in range(len(list_positions[func])):
                    if 'C' in list_positions[func][i]:
                        list_positions[func][i] = list_eq[func][C]
                        C+=1
            
            # * Gerando lista norte e sul
            # list_north = [[] for i in range(len(list_eq))]
            # list_south = [[] for i in range(len(list_eq))]

            if south_bd.lower() == 'dirichlet':
                list_south = dirichlet(south_func_bd, list_eq, 'south', n_part)
            elif south_bd.lower() == 'neumann':
                list_south = neumann(south_func_bd, list_eq, 'south', n_part, xd_var, str_sp_vars)
            elif south_bd.lower() == 'robin':
                list_south = robin(south_func_bd, list_eq, 'south', south_alpha_bd, south_beta_bd, n_part, xd_var)

            if north_bd.lower() == 'dirichlet':
                list_north = dirichlet(north_func_bd, list_eq, 'north', n_part)
            elif north_bd.lower() == 'neumann':
                list_north = neumann(north_func_bd, list_eq, 'north', n_part, xd_var, str_sp_vars)
            elif north_bd.lower() == 'robin':
                list_north = robin(north_func_bd, list_eq, 'north', north_alpha_bd, north_beta_bd, n_part, xd_var)

            # * Gerando lista leste
            list_east = [[] for i in range(len(list_eq))]
            for func in range(len(list_eq)):
                list_east[func].append(list_south[func][-1])
                # for i in range((n_part[1]-2)*(n_part[0]-2)-(n_part[1]-2), (n_part[1]-2)*(n_part[0]-2)): 
                #     list_east[func].append(list_eq[func][i])
                if east_bd.lower() == 'dirichlet':
                    centro = dirichlet(east_func_bd, list_eq, 'east', n_part)[func]
                elif east_bd.lower() == 'neumann':
                    centro = neumann(east_func_bd,list_eq, 'east', n_part, xd_var, str_sp_vars)[func]
                elif east_bd.lower() == 'robin':
                    centro = robin(east_func_bd, list_eq, 'east', east_alpha_bd, east_beta_bd, n_part, xd_var)[func]
                        
                for i in range(len(centro)):
                    list_east[func].append(centro[i])       
                list_east[func].append(list_north[func][-1])        

                
        elif len(str_sp_vars) == 1:
           
            for func in range(len(list_positions)):
                C = 0
                for i in range(len(list_positions[func])):
                    if 'C' in list_positions[func][i]:
                        list_positions[func][i] = list_eq[func][C]
                        C+=1                        
            # * Gerando lista leste

        # * Substituindo as posições por equações
        final_list_eq = []
        
        if len(str_sp_vars) == 3:
            pass
        
        elif len(str_sp_vars) == 2:
            for func in range(len(list_positions)):
                # Contadores
                S = 0
                N = 0
                E = 0
                W = 0
                for len_list in range(len(list_positions[func])):
                    # Substitui South's
                    if "S" in list_positions[func][len_list]:
                        list_positions[func][len_list] = list_south[func][S]
                        S += 1
                    # Substitui North's
                    elif "N" in list_positions[func][len_list]:
                        list_positions[func][len_list] = list_north[func][N] # -T_ij + Tw
                        N += 1
                    # Substitui East's
                    elif "E" in list_positions[func][len_list]:
                        list_positions[func][len_list] = list_east[func][E]
                        E += 1
                    elif "W" in list_positions[func][len_list]:
                        list_positions[func][len_list] = str(inlet[func][W])
                        W += 1
                              
        elif len(str_sp_vars) == 1:
            for func in range(len(list_positions)):
                # Contadores
                E = 0
                W = 0
                
                list_positions[func][-1] = list_east[func][0]
                
                list_positions[func][0] = str(inlet[func][0])
        
        d_vars = []
    
        if len(str_sp_vars) == 2:
            for func in range(len(list_positions)):
                for i in range(0,n_part[0]):
                    for j in range(0,n_part[1]):
                            if not(f'XX{func}_{i}_{j}' in d_vars):
                                d_vars.append(f'XX{func}_{i}_{j}')
        elif len(str_sp_vars) == 1:
            if not(f'XX{func}_{i}_0' in d_vars):
                d_vars.append(f'XX{func}_{i}_0')
            

        flat_list_positions = []
        
        for list in list_positions:
            flat_list_positions.extend(list)
        
        
        for _ in range(len(str_sp_vars)):
            for i in range(len(flat_list_positions)):
                flat_list_positions[i] = flat_list_positions[i].replace(f'h{str_sp_vars[_]}_', str(1/n_part[_]))
                for j in range(len(xd_var)):
                    flat_list_positions[i] = flat_list_positions[i].replace(f'h{xd_var[j]}_', str(1/n_part[_]))
                    
        # print(flat_list_positions)
        return flat_list_positions, d_vars


def dirichlet(bd_func, list_eq, bd, n_part, str_sp_vars = ''):
    if bd.lower() == 'north':
        list_north = [[] for i in range(len(list_eq))]
        for func in range(len(list_eq)): # para cada lista de equações discretizadas
            for j in range(n_part[func]): # para cada equação da lista
                list_north[func].append(bd_func)
                list_north[func][-1] = list_north[func][-1].replace(f'{str_sp_vars[0]}', f'{int(j / (n_part[1]-2)+n_part[1]-3)} * h{str_sp_vars[0]}_')
                list_north[func][-1] = list_north[func][-1].replace(f'{str_sp_vars[1]}', f'{int(j / (n_part[1]-2)+n_part[1]-3)} * h{str_sp_vars[1]}_')
        return list_north        
            
    elif bd.lower() == 'south':
        list_south = [[] for i in range(len(list_eq))]
        for func in range(len(list_eq)): # para cada lista de equações discretizadas
            for j in range(n_part[func]): # para cada equação da lista
                list_south[func].append(bd_func)
                list_south[func][-1] = list_south[func][-1].replace(f'{str_sp_vars[0]}', f'{int(j / (n_part[1]-2))} * h{str_sp_vars[0]}_')
                list_south[func][-1] = list_south[func][-1].replace(f'{str_sp_vars[1]}', f'{int(j / (n_part[1]-2))} * h{str_sp_vars[1]}_')
        return list_south
    
    
    elif bd.lower() == 'east':
        list_east = [[] for i in range(len(list_eq))]
        for func in range(len(list_eq)): # para cada lista de equações discretizadas
            for j in range(n_part[func]): # para cada equação da lista
                list_east[func].append(bd_func)
                list_east[func][-1] = list_east[func][-1].replace(f'{str_sp_vars[0]}', f'{i} * h{str_sp_vars[0]}_')
                list_east[func][-1] = list_east[func][-1].replace(f'{str_sp_vars[1]}', f'{i} * h{str_sp_vars[1]}_')
        return list_east
    
    else:
        print("Invalid boundary. Try 'east', 'north' or 'south' ")
        return []


def robin(bd_func, list_eq, bd, alpha, beta, n_part, xd_var, str_sp_vars = ''):
    if bd.lower() == 'north':
        list_north = [[] for i in range(len(list_eq))]
        for func in range(len(list_eq)): # para cada lista de equações discretizadas
            for j in range(0, len(list_eq[func])): # para cada equação da lista
                if j % (n_part[1]-2) == 0: # se o índice j for múltiplo de (n_part[1]-2)          
                    list_north[func].append(f'(h{xd_var[func]}_*{bd_func}-{list_eq[func][j+n_part[1]-3]}*({alpha}*h{xd_var[func]}_-{beta}))/{beta}')
                    list_north[func][-1] = list_north[func][-1].replace(f'{str_sp_vars[0]}', f'{int(j / (n_part[1]-2)+n_part[1]-3)} * h{str_sp_vars[0]}_')
                    list_north[func][-1] = list_north[func][-1].replace(f'{str_sp_vars[1]}', f'{int(j / (n_part[1]-2)+n_part[1]-3)} * h{str_sp_vars[1]}_')
        return list_north        
            
    elif bd.lower() == 'south':
        list_south = [[] for i in range(len(list_eq))]
        for func in range(len(list_eq)): # para cada lista de equações discretizadas
            for j in range(0, len(list_eq[func])): # para cada equação da lista
                if j % (n_part[1]-2) == 0: # se o índice j for múltiplo de (n_part[1]-2)  
                    list_south[func].append(f'(h{xd_var[func]}_*{bd_func}-{list_eq[func][j]}*({alpha}*h{xd_var[func]}_-{beta}))/{beta}')
                    list_south[func][-1] = list_south[func][-1].replace(f'{str_sp_vars[0]}', f'{int(j / (n_part[1]-2))} * h{str_sp_vars[0]}_')
                    list_south[func][-1] = list_south[func][-1].replace(f'{str_sp_vars[1]}', f'{int(j / (n_part[1]-2))} * h{str_sp_vars[1]}_')
        return list_south
    
    
    elif bd.lower() == 'east':
        list_east = [[] for i in range(len(list_eq))]
        for func in range(len(list_eq)):
            
            for i in range((n_part[1]-2)*(n_part[0]-2)-(n_part[1]-2), (n_part[1]-2)*(n_part[0]-2)): 
                list_east[func].append(f'(h{xd_var[func]}_*{bd_func}-{list_eq[func][i]}*({alpha}*h{xd_var[func]}_-{beta}))/{beta}')
                list_east[func][-1] = list_east[func][-1].replace(f'{str_sp_vars[0]}', f'{i} * h{str_sp_vars[0]}_')
                list_east[func][-1] = list_east[func][-1].replace(f'{str_sp_vars[1]}', f'{i} * h{str_sp_vars[1]}_')
            
        return list_east
    
    else:
        print("Invalid boundary. Try 'east', 'north' or 'south' ")
        return []


def neumann(bd_func, list_eq, bd, n_part, xd_var, str_sp_vars = ''):
    
    if bd.lower() == 'north':
        list_north = [[] for i in range(len(list_eq))]
        for func in range(len(list_eq)): # para cada lista de equações discretizadas
            for j in range(0, len(list_eq[func])): # para cada equação da lista
                if j % (n_part[1]-2) == 0: # se o índice j for múltiplo de (n_part[1]-2)          
                    list_north[func].append(f'h{xd_var[func]}_*{bd_func}+{list_eq[func][j+n_part[1]-3]}')

                    list_north[func][-1] = list_north[func][-1].replace(f'{str_sp_vars[0]}', f'{int(j / (n_part[0]-2)+n_part[0]-3)} * h{xd_var[0]}_')
                    list_north[func][-1] = list_north[func][-1].replace(f'{str_sp_vars[1]}', f'{n_part[1]-1} * h{xd_var[0]}_')

        return list_north

    elif bd.lower() == 'south':
        list_south = [[] for i in range(len(list_eq))]
        for func in range(len(list_eq)): # para cada lista de equações discretizadas
            for j in range(0, len(list_eq[func])): # para cada equação da lista
                if j % (n_part[1]-2) == 0: # se o índice j for múltiplo de (n_part[1]-2)  
                    list_south[func].append(f'h{xd_var[func]}_*{bd_func}+{list_eq[func][j]}')
                    list_south[func][-1] = list_south[func][-1].replace(f'{str_sp_vars[0]}', f'{int(j / (n_part[0]-2))} * h{xd_var[0]}_')
                    list_south[func][-1] = list_south[func][-1].replace(f'{str_sp_vars[1]}', f'0 * h{xd_var[0]}_')

        return list_south
    
    
    elif bd.lower() == 'east':
        list_east = [[] for i in range(len(list_eq))]
        for func in range(len(list_eq)):
            
            for i in range((n_part[1]-2)*(n_part[0]-2)-(n_part[1]-2), (n_part[1]-2)*(n_part[0]-2)): 
                list_east[func].append(f'h{xd_var[func]}_*{bd_func}+{list_eq[func][i]}')
                list_east[func][-1] = list_east[func][-1].replace(f'{str_sp_vars[0]}', f'{n_part[0]-1} * h{xd_var[0]}_')
                list_east[func][-1] = list_east[func][-1].replace(f'{str_sp_vars[1]}', f'{i} * h{xd_var[0]}_')

            
            
        return list_east
    
    else:
        print("Invalid boundary. Try 'east', 'north' or 'south' ")
        return []
    


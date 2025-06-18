import PDE
import SERKF45
import RKF45Vetorizado as RK
import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.mplot3d import Axes3D  # habilita o 3D em matplotlib
from matplotlib.animation import FuncAnimation

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


    def df(self, n_part, list_east = [], inlet = [], method="foward"):
       
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

        
      
        if (method == "foward"):
            for j in range(len(eqrs)):
                # Substituindo as derivadas parciais por diferenças finitas adiantadas
                for k in range(len(str_sp_vars)):
                    if k == 0: # se a variável de discretização for a primeira
                        for i in range(len(xd_var)): 
                            eqrs[j] = eqrs[j].replace(f'd2{xd_var[i]}{str_sp_vars}/d{str_sp_vars[k]}2', f'({xd_var[i]}i+1j - 2*{xd_var[i]}ij + {xd_var[i]}i-1j)/ h{str_sp_vars[k]} ** 2') 
                            eqrs[j] = eqrs[j].replace(f'd{xd_var[i]}{str_sp_vars}/d{str_sp_vars[k]}', f'({xd_var[i]}i+1j - {xd_var[i]}ij)/ h{str_sp_vars[k]}')
                            
                    elif k == 1: # se a variável de discretização for a segunda
                        for i in range(len(xd_var)): 
                            
                            eqrs[j] = eqrs[j].replace(f'd2{xd_var[i]}{str_sp_vars}/d{str_sp_vars[k]}2', f'({xd_var[i]}ij+1 - 2*{xd_var[i]}ij + {xd_var[i]}ij-1)/ h{str_sp_vars[k]} ** 2')
                            eqrs[j] = eqrs[j].replace(f'd{xd_var[i]}{str_sp_vars}/d{str_sp_vars[k]}', f'({xd_var[i]}ij+1 - {xd_var[i]}ij)/ h{str_sp_vars[k]}')
                            eqrs[j] = eqrs[j].replace(f'{xd_var[i]}{str_sp_vars}', f'{xd_var[i]}ij')
        elif (method == "central"):            
            for j in range(len(eqrs)):
                # Substituindo as derivadas parciais por diferenças finitas adiantadas
                for k in range(len(str_sp_vars)):
                    if k == 0: # se a variável de discretização for a primeira
                        for i in range(len(xd_var)): 
                            eqrs[j] = eqrs[j].replace(f'd2{xd_var[i]}{str_sp_vars}/d{str_sp_vars[k]}2', f'({xd_var[i]}i+1j - 2*{xd_var[i]}ij + {xd_var[i]}i-1j)/ h{str_sp_vars[k]} ** 2') 
                            eqrs[j] = eqrs[j].replace(f'd{xd_var[i]}{str_sp_vars}/d{str_sp_vars[k]}', f'(({xd_var[i]}i+1j - {xd_var[i]}ij)-({xd_var[i]}ij - {xd_var[i]}i-1j))/ h{str_sp_vars[k]}')
                            
                    elif k == 1: # se a variável de discretização for a segunda
                        for i in range(len(xd_var)): 
                            
                            eqrs[j] = eqrs[j].replace(f'd2{xd_var[i]}{str_sp_vars}/d{str_sp_vars[k]}2', f'({xd_var[i]}ij+1 - 2*{xd_var[i]}ij + {xd_var[i]}ij-1)/ h{str_sp_vars[k]} ** 2')
                            eqrs[j] = eqrs[j].replace(f'd{xd_var[i]}{str_sp_vars}/d{str_sp_vars[k]}', f'(({xd_var[i]}ij+1 - {xd_var[i]}ij)-({xd_var[i]}ij - {xd_var[i]}ij-1))/ h{str_sp_vars[k]}')
                            eqrs[j] = eqrs[j].replace(f'{xd_var[i]}{str_sp_vars}', f'{xd_var[i]}ij')
        elif (method == "backward"):
            for j in range(len(eqrs)):
                # Substituindo as derivadas parciais por diferenças finitas adiantadas
                for k in range(len(str_sp_vars)):
                    if k == 0: # se a variável de discretização for a primeira
                        for i in range(len(xd_var)): 
                            eqrs[j] = eqrs[j].replace(f'd2{xd_var[i]}{str_sp_vars}/d{str_sp_vars[k]}2', f'({xd_var[i]}i+1j - 2*{xd_var[i]}ij + {xd_var[i]}i-1j)/ h{str_sp_vars[k]} ** 2') 
                            eqrs[j] = eqrs[j].replace(f'd{xd_var[i]}{str_sp_vars}/d{str_sp_vars[k]}', f'({xd_var[i]}ij - {xd_var[i]}i-1j)/ h{str_sp_vars[k]}')
                            
                    elif k == 1: # se a variável de discretização for a segunda
                        for i in range(len(xd_var)): 
                            
                            eqrs[j] = eqrs[j].replace(f'd2{xd_var[i]}{str_sp_vars}/d{str_sp_vars[k]}2', f'({xd_var[i]}ij+1 - 2*{xd_var[i]}ij + {xd_var[i]}ij-1)/ h{str_sp_vars[k]} ** 2')
                            eqrs[j] = eqrs[j].replace(f'd{xd_var[i]}{str_sp_vars}/d{str_sp_vars[k]}', f'({xd_var[i]}ij - {xd_var[i]}ij-1)/ h{str_sp_vars[k]}')
                            eqrs[j] = eqrs[j].replace(f'{xd_var[i]}{str_sp_vars}', f'{xd_var[i]}ij')
        else:
            raise ValueError("Método de discretização inválido. Use 'foward', 'central' ou 'backward'.")
        
        #* A partir daqui: Substituição dos indices i e j por seus valores correspondentes em eqrs:
        
        #*                 [[centros eq1][centros de eq2]]
        #*                 [W0, W1, W2, W3, S1, ce0, ce1, N1, S2, ce2, ce3, N2, E0, E1, E2, E3]
        
        #*                 Generalizando para 2 variáveis de discretização.
        

        partial_list_eq = []

        # * Gerando lista de equações com indices parcialmente substituídos
        for j in range(len(eqrs)):
            partial_list_eq.append([eqrs[j].replace('i+1', str(i+1)).replace('i-1', str(i-1)).replace('i-2', str(i-2)).replace('i+2', str(i+2)).replace('i',str(i)) for i in range(1,n_part[0]-1)])
        

        list_eq = []
        list_eq = [[] for i in range(len(partial_list_eq))]
        
        if len(str_sp_vars) == 3:
            pass
        
        # * Gerando lista de equações com indices totalmente substituídos para duas variáveis de discretização
        elif len(str_sp_vars) == 2:
            for j in range(len(partial_list_eq)): # para cada lista de equações discretizadas
                for i in range(len(partial_list_eq[j])): # para cada equação da lista
                    for k in range(1,n_part[1]-1): # indice j da segunda variável de discretização
                        list_eq[j].append(partial_list_eq[j][i].replace('j+1', str(k+1)).replace('j-1', str(k-1)).replace('j-2', str(k-2)).replace('j+2', str(k-2)).replace('j',str(k)))
          
        
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
            pass
        
        # * Substituindo centros na lista de posições
        elif len(str_sp_vars) == 2:
        #     for func in range(len(list_positions)): # para cada lista de posições
        #         for jj in range(len(list_positions[func])): # para cada posição da lista
        #             for i in range(1,n_part[0]-1): # para cada índice i da primeira variável de discretização
        #                 for j in range(1,n_part[1]-1): # para cada índice j da segunda variável de discretização
        #                     # Substituindo centros na lista de posições por suas respectivas equações
        #                     print((n_part[0]-2)*(i-1)+j-1)
        #                     list_positions[func][jj] = list_positions[func][jj].replace(f'Ce{func}{i}{j}', list_eq[func][(n_part[0]-2)*(i-1)+j-1])
            
            for func in range(len(list_positions)):
                C = 0
                for i in range(len(list_positions[func])):
                    if 'C' in list_positions[func][i]:
                        list_positions[func][i] = list_eq[func][C]
                        C+=1
            
            
            
            # * Gerando lista norte e sul
            list_north = [[] for i in range(len(list_eq))]
            list_south = [[] for i in range(len(list_eq))]
            for func in range(len(list_eq)): # para cada lista de equações discretizadas
                i = 0
                k = 1
                for j in range(0, len(list_eq[func])): # para cada equação da lista
                    if j % (n_part[1]-2) == 0: # se o índice j for múltiplo de (n_part[1]-2)
                        list_south[func].append(list_eq[func][j])
                        if func == 1:
                            list_north[func].append(f'0') #! a
                            # list_north[func].append(f'2.78*(100 - XX1{n_part[1]-1}{k})') #! a
                            print(list_north[func][-1])
                            k+=1
                        else:
                            list_north[func].append("(("+list_eq[func][j+n_part[1]-3]+")-(" +"0"+"))") #! a
                            print("(("+list_eq[func][j+n_part[1]-3]+")-(" +"0"+"))")
                        i+=1
                    

            
            # * Gerando lista leste
            if list_east == []:
                list_east = [[] for i in range(len(list_eq))]
                for func in range(len(list_eq)):
                    list_east[func].append(list_south[func][-1])
                    for i in range((n_part[1]-2)*(n_part[0]-2)-(n_part[1]-2), (n_part[1]-2)*(n_part[0]-2)): 
                        list_east[func].append(list_eq[func][i])
                    list_east[func].append(list_north[func][-1])        
            
        
        elif len(str_sp_vars) == 1:
            # * Substituindo centros na lista de posições
            # for func in range(len(list_positions)):
            #     for jj in range(len(list_positions[func])):
            #         for i in range(1,n_part[0]-1):
            #             list_positions[func][jj] = list_positions[func][jj].replace(f'Ce{func}{i}0', list_eq[func][(n_part[0]-2)*(i-1)])
            for func in range(len(list_positions)):
                C = 0
                for i in range(len(list_positions[func])):
                    if 'C' in list_positions[func][i]:
                        list_positions[func][i] = list_eq[func][C]
                        C+=1                        
            # * Gerando lista leste
            if list_east == []: # Verifica se a lista leste já foi criada pelo usuário            
                list_east = [[] for i in range(len(list_eq))]
                for func in range(len(list_positions)):              
                    
                    list_east[func].append(list_eq[func][-1])
        
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
                    
                    # # Substitui East's
                    # if "E" in list_positions[func][len_list]:
                    #     list_positions[func][len_list] = list_east[func][E]
                    #     E += 1
                    # elif "W" in list_positions[func][len_list]:
                    #     list_positions[func][len_list] = str(inlet[func][W])
                    #     W += 1
        
        
        d_vars = []
    

        if len(str_sp_vars) == 2:
            for func in range(len(list_positions)):
                for i in range(0,n_part[0]):
                    for j in range(0,n_part[1]):
                            if not(f'XX{func}{i}{j}' in d_vars):
                                d_vars.append(f'XX{func}{i}{j}')
        # elif len(str_sp_vars) == 1:
        #     if not(f'XX{func}{i}0' in d_vars):
        #                 d_vars.append(f'XX{func}{i}0')
            

        flat_list_positions = []
        
        for list in list_positions:
            flat_list_positions.extend(list)
        
        
        for _ in range(len(str_sp_vars)):
            for i in range(len(flat_list_positions)):
                flat_list_positions[i] = flat_list_positions[i].replace(f'h{str_sp_vars[_]}', str(1/n_part[_]))
        
        
        print(flat_list_positions)
        return flat_list_positions, d_vars

#! Valores Iniciais

PDE2 = PDE.PDE('dT/dt = 0.1*d2T/dx2 + 0.9*d2T/dy2', ['t', 'x', 'y'], ['T'])
PDE1 = PDE.PDE('dC/dt = 0.1*d2C/dx2', ['t', 'x', 'y'], ['C','T'])
PDE3 = PDE.PDE('dD/dt = -15400/(8.34*T) * C**0.524', ['t'], ['D','T'])

PDES1 = PDES([PDE1, PDE2, PDE3], ['x','y'], ['C','T','D'])


disc_n=9
resultado = PDES1.df([disc_n,disc_n], inlet=[[0 for i in range(disc_n)],[1 for i in range(disc_n)], [0 for i in range(disc_n)]], method="foward")
print(resultado[1])
initial_values = []
for i in range(0,disc_n):
    initial_values.append(0)
for i in range(disc_n,disc_n**2):
    initial_values.append(0)

for i in range(0,disc_n):
    initial_values.append(100)
for i in range(disc_n,disc_n**2):
    if (i+1) % disc_n == 0:
        initial_values.append(100)
    else:
        initial_values.append(20)

for i in range(0,disc_n):
    initial_values.append(0)
for i in range(disc_n,disc_n**2):
    initial_values.append(0)
# for i in range(0, disc_n[0]+1):
#     initial_values.append(1)
# for i in range(disc_n[0]+1,(disc_n[1]+1)*(disc_n[0]+1)):
#     initial_values.append(1)

# for i in range(0, disc_n[0]+1):
#     initial_values.append(100)
# for i in range(disc_n[0]+1,(disc_n[1]+1)*(disc_n[0]+1)):
#     initial_values.append(20)    

# initial_values.append(100)        
# for i in range (1, disc_n):
#     initial_values.append(80)

# print(initial_values)
# print(resultado[1])
# print(SERKF45.SERKF45(resultado[0], ['t'], resultado[1], initial_values, 0, 0.1, 10,2,len(PDES1.sp_vars)))
testar = RK.SERKF45_cuda(resultado[0], ['t'], resultado[1], initial_values, 0, 0.5, 100, 3, len(PDES1.sp_vars))

if len(PDES1.sp_vars) == 2:
    
    vetor = np.array(testar[1][0][-1], dtype=float)  # Convertendo para um array NumPy
    vetor = vetor.reshape((disc_n, disc_n), order='F')  # Reshape para matriz 2D (coluna maior que linha)
    plt.imshow(vetor, cmap='RdYlBu_r', interpolation='bilinear', extent=(0, 1, 0, 1), origin='lower')
    plt.show()
    
    lista_vetores = testar[1][0]  
    n_frames = len(lista_vetores)

    # 1) Empilha e reshape para (n_frames, disc_n, disc_n)
    data3d = np.array(lista_vetores, dtype=float) \
             .reshape((n_frames, disc_n, disc_n), order='F')
    
    # 2) Plot 3D
    # cria grelha de coordenadas X, Y no domínio [0,1]×[0,1]
    x = np.linspace(0, 1, disc_n)
    y = np.linspace(0, 1, disc_n)
    X, Y = np.meshgrid(x, y)

    # 2) Setup inicial da figura 3D
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.set_xlim(0,1)
    ax.set_ylim(0,1)
    ax.set_zlim(np.min(data3d), np.max(data3d))
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Valor')

    # plota a primeira fatia
    surf = ax.plot_surface(X, Y, data3d[0], cmap='RdYlBu_r', edgecolor='none')

    # 3) Função de update
    def update(frame):
        global surf
        # remove superfície antiga
        surf.remove()
        # desenha nova
        surf = ax.plot_surface(
            X, Y, data3d[frame],
            cmap='RdYlBu_r', edgecolor='none'
        )
        ax.set_title(f"Frame {frame+1}/{n_frames}")

    # 4) Cria animação
    ani = FuncAnimation(
        fig, update,
        frames=n_frames,
        interval=200,      # ms entre frames
        blit=False         # blit não funciona bem em 3D
    )

    plt.show()
    from matplotlib.animation import FFMpegWriter

    # configura o writer:
    writer = FFMpegWriter(
        fps=10,                # quadros por segundo
        metadata=dict(artist='Você'),
        bitrate=1800           # taxa de bits (quanto maior, melhor qualidade/tamanho)
    )

    # salva a animação:
    ani.save('animacao3d.mp4', writer=writer)

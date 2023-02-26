import os
import time
import random

class jogo_da_velha:
    def __init__(self):
        #o fluxo do jogo todo vai ser feito aqui, as demais funções não são necessárias de alterar.
        self._status = [[0 for i in range(3)] for j in range(3) ]
        self.rodada=0
        print('Bem-Vindo ao jogo da velha, escreva que espaço quer preencher!')
        print(self._desenho_partida())
        while self._jogo_nao_terminou():
            self._status = self._jogada_player()
            self._limpar_tela()
            self.rodada+=1
            print(self._desenho_partida())
            if(self.rodada<9):
                print('O Bot está pensando...')
                self.rodada+=1
                time.sleep(2) 
                self._limpar_tela()
                self._status = self._jogada_bot()
                print(self._desenho_partida())
                

    def _limpar_tela(self):
        os.system('clear')

    def _jogo_nao_terminou(self):
        if self._alguem_ganhou():
            print('Jogo terminado!')
            return False
        elif self.rodada==9:
            print('Jogo Empatado!')
        else:
            return True

    def _alguem_ganhou(self):
        return self._check_colunas() or self._check_diagonais() or self._check_linhas()


    def _check_diagonais(self):
        #checando diagonais
        diagonal_1 = self._status[0][0] + self._status[1][1] + self._status[2][2]
        diagonal_2 = self._status[0][2] + self._status[1][1] + self._status[2][0]
        if diagonal_1 == 3 or diagonal_2 == 3:
            print('PLAYER GANHOU!')
            return True
        elif diagonal_1 ==-3 or diagonal_2 == -3:
            print('BOT GANHOU!')
            return True
        else:
            return False
        
    def _check_colunas(self):
        #checando colunas
        for linha in range(3):
            if self._status[linha][0] + self._status[linha][1] + self._status[linha][2] == 3:
                print('PLAYER GANHOU!')
                return True
            elif self._status[linha][0] + self._status[linha][1] + self._status[linha][2] == -3:
                print('BOT GANHOU!')
                return True
        return False
 
 
    def _check_linhas(self):
        #checando linhas
        for coluna in range(3):
            if self._status[0][coluna] + self._status[1][coluna] + self._status[2][coluna] == 3:
                print('PLAYER GANHOU!')
                return True
            elif self._status[0][coluna] + self._status[1][coluna] + self._status[2][coluna] == -3:
                print('BOT GANHOU!')
                return True
        return False


    def _converter_status_para_desenho(self):
        preenche__ = [[' ' if x==0 else x for x in self._status[y][:]] for y in range(3)]
        preenche_X = [['X' if x==1 else x for x in preenche__[y][:]] for y in range(3)]
        preenche_O = [['0' if x==-1 else x for x in preenche_X[y][:]] for y in range(3)]
        return preenche_O
    
    def _desenho_partida(self):
        print('Rodada {}'.format(self.rodada))
        desenho_status = self._converter_status_para_desenho()
        desenho = '   1   2   3\na  {a1} | {a2} | {a3} \n  ---|---|---\nb  {b1} | {b2} | {b3} \n  ---|---|---\nc  {c1} | {c2} | {c3} '.format(
            a1=desenho_status[0][0],
            a2=desenho_status[0][1],
            a3=desenho_status[0][2],
            b1=desenho_status[1][0],
            b2=desenho_status[1][1],
            b3=desenho_status[1][2],
            c1=desenho_status[2][0],
            c2=desenho_status[2][1],
            c3=desenho_status[2][2],
        )
        return desenho
    
    def _jogada_player(self):
        x=input('Que espaço quer preencher?\n')
        while not self._jogada_eh_valida(x,1):
           x = input('Entre com um espaço válido!\n')    
        return self._status
         
    def _jogada_eh_valida(self,input,tag):
        list_espaco_valido = ['a1','a2','a3','b1','b2','b3','c1','c2','c3']
        if input in list_espaco_valido:
            if self._espaco_ja_preenchido_se_nao_preenche(input,tag):
                print('Esse espaço já esta ocupado!')
                tag_eh_valido = False 
            else:
                tag_eh_valido = True
        else:
            print('Entre com algum espaço válido!\nExemplo: a1,b2,c3...')
            tag_eh_valido = False

        return tag_eh_valido        

    def _espaco_ja_preenchido_se_nao_preenche(self,input,tag):
        tag_espaco_preenchido = False
        if input == 'a1' and self._status[0][0]==0:
            self._status[0][0] = tag
        elif input == 'a2' and self._status[0][1]==0:
            self._status[0][1] = tag
        elif input == 'a3' and self._status[0][2]==0:
            self._status[0][2] = tag
        elif input == 'b1' and self._status[1][0]==0:
            self._status[1][0] = tag
        elif input == 'b2' and self._status[1][1]==0:
            self._status[1][1] = tag
        elif input == 'b3' and self._status[1][2]==0:
            self._status[1][2] = tag
        elif input == 'c1' and self._status[2][0]==0:
            self._status[2][0] = tag
        elif input == 'c2' and self._status[2][1]==0:
            self._status[2][1] = tag  
        elif input == 'c3' and self._status[2][2]==0:
            self._status[2][2] = tag
        else:
            tag_espaco_preenchido = True
        return tag_espaco_preenchido

    
    def _jogada_bot(self):
        #Simulando aleatoriedade
        list_espaco_valido = ['a1','a2','a3','b1','b2','b3','c1','c2','c3'] 
        input = random.sample(list_espaco_valido,k=1)[0]
        while self._espaco_ja_preenchido_se_nao_preenche(input=input,tag=-1):
            list_espaco_valido.remove(input)
            input = random.sample(list_espaco_valido,k=1)[0]
        return self._status
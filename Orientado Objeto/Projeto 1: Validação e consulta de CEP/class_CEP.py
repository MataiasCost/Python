'''
Nota do Autor
Nesse código será declarada a classe CEP, onde será feita a validação e obtenção de informações com base no número de um CEP
'''
## Importanto os pacotes necessários
# Pacote para acessar a API do site: https://viacep.com.br/
import requests
# Pacote para usar regex, e com isso criar padrões de busca/consulta
import re

class CEP:
    def __init__(self,num_cep):          
        if self.format_cep(str(num_cep)):
            self.num_cep = str(num_cep)
            if 'erro' in self.request_API().keys():
                print('O CEP informado não existe no banco de dados do site utilizado.\nPara mais informações https://viacep.com.br/')
            else:
                r = self.request_API()
                # print(r)
                print("Informações obtidas sobre o CEP:\nEstado:{}\nBairro:{}\nRua:{}".format(r['uf'],
                                                                                              r['bairro'],
                                                                                              r['logradouro']),
                '\nPara mais informações: https://viacep.com.br/'
                )


    def __str__(self):
         return self.mask_cep()

    def format_cep(self,num_cep):
        padrao_cep = '^[0-9]{5}[-]?[0-9]{3}$'
        #Testa se o CEP tem um dos dois formatos, de 8 ou 9 termos
        # Ex: 12345678 ou 12345-678
        if re.findall(padrao_cep,num_cep):
            return num_cep
        #Caso não, o programa para e retorna que o formato é inválido
        else:
            print('O formato do CEP é inválido!\nInforme CEP no formato 12345-678 ou 12345678')
            return False

    def mask_cep(self):
        padrao = '[0-9]{5}-[0-9]{2}'
        if re.findall(padrao,self.num_cep):
            return self.num_cep
        else:
            return '{}-{}'.format(self.num_cep[:5],self.num_cep[5:])
    def request_API(self):
        res = re.sub('[-]', '', self.num_cep)
        #Criando o link de acesso
        url = 'https://viacep.com.br/ws/{}/json/'.format(res)
        # Pedindo request ao site, através do pacote request
        r = requests.get(url)
        #Retorna o valor com base na resposta da comunicação:
        #Informational responses (100 – 199)
        #Successful responses (200 – 299)
        #Redirection messages (300 – 399)
        #Client error responses (400 – 499)
        #Server error responses (500 – 599)
        #Mais informações https://developer.mozilla.org/en-US/docs/Web/HTTP/Status
        if 200 <= r.status_code <= 299:
            return r.json()
        else:
            raise Exception('Desculpe, mas não conseguimos acessar a API do site. Código: {}'.format(r.status_code)) 
#O que ele faz:
# #-Le:
#     -1.Tabela com as observacoes do xmm 
# que estao no sdss
# Para cada objeto, procura o redshift no NED
# (primeiro por nome, depois por ra,dec)
# -Atualiza 1. com esses valores e marca com uma flag
# os que nao tem redshift

#---------------------------------------------------------------
#FUNCOES
#---------------------------------------------------------------

def Compara(obj, nome_url, z_url):
  # Vamos ver o comprimento de nome_url. A 1st casa eh pro nome da 
  #variavel. as demais, pros valores. entao se so tiver um obj, len = 2
  # Caso tenha mais de um dado, temos que escolher entre eles. Os criterios sao:
  # 1)maior semelhanca no nome
  # 2) menor z

  comp = len(nome_url)
  if comp > 2:

    #compara nome:faz uma lista com os caracteres do nome
    #a comparacao eh feita de frente pra tras e no reverso
    #a cada caracter igual, 'igual +=1'.
    igual = zeros(comp)
    i = 1
    l_obj = list(obj)
    l_obj_rev = list(reversed(l_obj))


    while i< comp:
      l_NED = list(nome_url[i])
      l_NED_rev = list(reversed(l_NED))
      nome_url[i] = nome_url[i].replace(' ','') #tira so espacos
      letras = min(len(l_obj),len(nome_url[i])) 
      for c in range(letras):
        if l_obj[c] == l_NED[c]:
          igual[i] +=1
        if l_obj_rev[c] == l_NED_rev[c]:
          igual[i] +=1
      i+=1

    if any(item!=0 for item in igual):
      dado_certo = argmax(igual)
    else:
      #se nao achar o mais proximo pelo nome, vai pelo menor z
      dado_certo = argmin(z_url)

  else:
    dado_certo = 1
      
  return dado_certo

#---------------------------------------------------------------
#MAIN
#---------------------------------------------------------------
#from scypy import *
from numpy import *
import urllib3 as ur
pool = ur.PoolManager()
import astropy.io.ascii as ap
from astropy.table import *


#pergunta qual o arquivo que sera lido
#arq = input('nome do arquivo: ')
#path = input('path do arquivo: ')
arq = 'xmm_clean.csv'
path = '/home/natalia/Dados/XMM'

#Le o arquivo
le = path+'/'+arq
xmm = ap.read('%s'%le,delimiter=',',guess=False)

OBJ = xmm["TARGET"]
RA = xmm['RA(deci)']
DEC = xmm['DEC(deci)']
RAH = xmm['RA(h)']
DECH = xmm['DEC(h)']

#esses dois arrays serao usados pra guardar as infos
#buscadas nos urls do NED

nome_NED=[]
ra_NED = []
dec_NED = []
z_NED=[]

sem_nome = 0
i = 0

for i in range(len(OBJ)):
#for i in range(1):
  obj = OBJ[i].replace(' ','')
  ra = RA[i]
  dec = DEC[i]
  rah = RAH[i]
  dech = DECH[i]
  raio = 2.0
  urlnome = "https://ned.ipac.caltech.edu/cgi-bin/objsearch?objname="+obj+"&radius="+'%s'%raio+"&ot_include=ANY&in_objtypes1=GGroups&in_objtypes1=GClusters&ex_objtypes1=Galaxies&nmp_op=ANY&out_csys=Equatorial&out_equinox=J2000.0&obj_sort=Distance+to+search+center&of=ascii_bar&zv_breaker=30000.0&list_limit=5&img_stamp=YES&search_type=Near+Name+Search"
  #urlpos = "https://ned.ipac.caltech.edu/cgi-bin/objsearch?in_csys=Equatorial&in_equinox=J2000.0&lon="+'%f'%ra+"d&lat="+'%f'%dec+"d&radius=2.0&ot_include=ANY&in_objtypes1=GGroups&in_objtypes1=GClusters&ex_objtypes1=Galaxies&nmp_op=ANY&out_csys=Equatorial&out_equinox=J2000.0&obj_sort=Distance+to+search+center&of=ascii_bar&zv_breaker=30000.0&list_limit=5&img_stamp=YES&search_type=Near+Position+Search"
  urlpos = "https://ned.ipac.caltech.edu/cgi-bin/objsearch?in_csys=Equatorial&in_equinox=J2000.0&lon="+'%s'%rah+"&lat="+'%s'%dech+"&radius="+'%s'%raio+"&ot_include=ANY&in_objtypes1=GGroups&in_objtypes1=GClusters&ex_objtypes1=Galaxies&nmp_op=ANY&out_csys=Equatorial&out_equinox=J2000.0&obj_sort=Distance+to+search+center&of=ascii_bar&zv_breaker=30000.0&list_limit=5&img_stamp=YES&search_type=Near+Position+Search"

  #procurando o url ( faremos ate 5 tentativas)
  # BUSCA POR NOME
  chance = 0
  while chance < 5:
    # try: se um erro acontece no bloco do 'try', pula
    #pro bloco do 'exception'. 
    #vou colocar a informacao do url num txt pra
    #poder ler como uma tabela
    try:
      info=pool.request('GET',urlnome).data
      o = open('out.txt','wb')
      o.write(info)
      o.close()
      break
    #Agora, se der um ERRO
    except ur.URLError as err:
      chance +=1
      print (type(err))
      print ('\n')

  #TESTE DE ERRO => NUM DE LINHAS
  #NAO ENTENDI AAAAAAAAAAAA
  #pass: quando vc nao quer uma saida explicita, 
  #mas precia escrever alguma acao pro comando
  #No 'j' vao ficar guardado o tamanho da ultima linha
  #do 'out', e no 'l', seu conteudo 
  with open('out.txt') as o:
    for j,l in enumerate(o):
      pass
  #achou pelo nome:
  #pega no 'out.txt' os dados
  if not j+1<=22:
    nome_url = genfromtxt('out.txt',delimiter = '|',skip_header=23, usecols=[1],dtype="str")
    ra_url = genfromtxt('out.txt',delimiter = '|',skip_header=23, usecols=[2])
    dec_url = genfromtxt('out.txt',delimiter = '|',skip_header=23, usecols=[3])
    z_url = genfromtxt('out.txt',delimiter = '|',skip_header=23, usecols=[6])
        
    obj_certo = Compara(obj,nome_url,z_url)
    nome_data = nome_url[obj_certo]
    ra_data = ra_url[obj_certo]
    dec_data = dec_url[obj_certo]
    z_data = z_url[obj_certo]

    nome_NED.append(nome_data)
    ra_NED.append(ra_data)
    dec_NED.append(dec_data)
    z_NED.append(z_data)



  #Nao encontrou pelo nome (obj)
  else:
    print ('objeto', i, 'nao foi encontrado pelo nome\n')
    sem_nome +=1

    #BUSCA POR POSICAO
    chance = 0
    while chance < 5:
        try:
          info=pool.request('GET',urlpos).data
          o = open('out.txt','wb')
          o.write(info)
          o.close()
          break
        except url.URLError as err:
          chance +=1
          print (type(err))

    with open('out.txt') as o:
      for j,l in enumerate(o):
        pass
    #achou pela pos:
    if not j+1<=22:
      nome_url = genfromtxt('out.txt',delimiter = '|',skip_header=23, usecols=[1],dtype="str")
      ra_url = genfromtxt('out.txt',delimiter = '|',skip_header=23, usecols=[2])
      dec_url = genfromtxt('out.txt',delimiter = '|',skip_header=23, usecols=[3])
      z_url = genfromtxt('out.txt',delimiter = '|',skip_header=23, usecols=[6])

      obj_certo = Compara(obj,nome_url,z_url)
      nome_data = nome_url[obj_certo]
      ra_data = ra_url[obj_certo]
      dec_data = dec_url[obj_certo]
      z_data = z_url[obj_certo]

      nome_NED.append(nome_data)
      ra_NED.append(ra_data)
      dec_NED.append(dec_data)
      z_NED.append(z_data)
    #Nao encontrou pela pos
    else:
      print ('objeto', i, 'nao foi encontrado pela posicao\n')
      nome_data = 'nao encontrado'
      ra_data = 999999
      dec_data = 999999
      z_data = 99999

      nome_NED.append(nome_data)
      ra_NED.append(ra_data)
      dec_NED.append(dec_data)
      z_NED.append(z_data)


#Colocando esse resultado no mesmo array de entrada

xmm['nome_NED'] = nome_NED
xmm['ra_NED'] = ra_NED
xmm['dec_NED'] = dec_NED
xmm['z_NED'] = z_NED

#Com isso pronto, vamos separar os dados em duas tabelas:
#uma contendo os objetos que tem z, e outra dos que nao tem


headers = xmm.colnames
tipo = xmm.dtype
h_tipo = []

i = 0
for i in range(len(headers)):
  aux = tipo[i]
  h_tipo.append(aux)


xmm_sz = Table(names = headers, dtype = h_tipo)
remove = []

i = 0
for i in range(len(z_NED)):
  if z_NED[i]>100 or isnan(z_NED[i])==True:
    remove.append(i)
    passa = xmm[:][i]
    xmm_sz.add_row(passa)

xmm.remove_rows(remove)

#Por fim, vamos ordenar a tabela xmm por redshift (crescente)
xmm.sort('z_NED')


#fala onde o arquivo ser asalvo e com qual nome
arq_sai = 'xmm_ned.csv'
arq_sai_sz = 'xmm_ned_sem_z.csv'
le_sai = path+'/'+arq_sai
le_sai_sz = path+'/'+arq_sai_sz

xmm.write('%s'%le_sai, delimiter = ',', overwrite = True)
xmm_sz.write('%s'%le_sai_sz, delimiter = ',', overwrite = True)





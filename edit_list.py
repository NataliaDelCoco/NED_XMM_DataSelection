#---------------------------------------------------------------
#MAIN
#---------------------------------------------------------------
#from scypy import *
from numpy import *
import urllib3 as ur
import astropy.io.ascii as ap
from astropy.table import * 
import ast

#pergunta qual o arquivo que sera lido
arq_xmm = input('nome do arquivo xmm: ')
#path_xmm = input('path do arquivo xmm: ')
arq_tof = input('nome do arquivo True or False: ')
#path_tof = input('path do arquivo True or False: ')
#arq_xmm = '1_NXSA-Results.csv'
path_xmm = '/home/natalia/Dados/XMM'
#arq_tof = '1_true_false.csv'
path_tof = '/home/natalia/Dados/XMM'

le_xmm = path_xmm+'/'+arq_xmm
le_tof = path_tof+'/'+arq_tof

xmm = ap.read('%s'%le_xmm,delimiter=',',guess=False)
tof = ap.read('%s'%le_tof,delimiter=',',guess=False)

#mudando os nomes das colunas, e retirar aquelas que nao queremos
aux = xmm.colnames
xmm.rename_column(aux[0] ,'ID')
xmm.rename_column(aux[1], 'TARGET')
xmm.rename_column(aux[2], 'RA(deci)')
xmm.rename_column(aux[3], 'DEC(deci)')
xmm.rename_column(aux[4], 'DURATION')
xmm.rename_column(aux[5], 'TYPE')

# xmm.remove_columns(['PROPOSAL.PI_SURNAME','PROPOSAL.TYPE',
#   'OBSERVATION.PROPRIETARY_END_DATE','OBSERVATION.COORD_OBS'])

#Primeiro, vamos retirar aqueles com FALSE no footprint do SDSS
count = len(tof)
i = 0

#o axis = 1 no delete indica que ta apagando a linha i inteira de todas as colunas
#quando so tem 1 coluna, usa o axis = 0
while i<count:
  if 'True' in tof[i]:
    i += 1
  else:
    tof.remove_row(i)
    xmm.remove_row(i)
    count -= 1


#Vamos retirar todos os objetos que nao sao CLUSTERS
#comparacao pelo nome
count = len(xmm['TYPE'])
i = 0
while i<count:
  if 'CLUSTER' in xmm['TYPE'][i]:
    i += 1
  else:
    xmm.remove_row(i)
    count -= 1

#retirando objetos com menos de 20000s de duracao
count = len(xmm['ID'])
i = 0
while i<count:
  if xmm['DURATION'][i] >= 20000:
    i += 1
  else:
    xmm.remove_row(i)
    count -= 1

#Agora, vamos retirar as observacoes repetidas
#olharemos num raio de  0.3d ao redor de cada objeto. os Targets 
#encontrados dentro desse raio serao considerados repetidos.
#Dentre os repetidos, aquele com maior duracao sera mantido
count = len(xmm['ID'])
i = 0
marker = 0

while i<count:
  j = i + 1 
  while j< count: 
    dif_ra = abs(xmm['RA(deci)'][i] - xmm['RA(deci)'][j])
    dif_dec = abs(xmm['DEC(deci)'][i] - xmm['DEC(deci)'][j])
    dist = sqrt(dif_ra**2 + dif_dec**2)

    if dist <= 0.3:
      if xmm['DURATION'][i] < xmm['DURATION'][j]:
        xmm.remove_row(i)
        count -= 1
        marker = 2 # esse marker eh pra nao somar i=i+1, e continuar com o mesmo i de antes, pois a informacao contida na linha i agora eh outra, ja que apaguei a antiga linha i.
        break # o break tira desse loop. ou seja, vai pro "if marker".
      else:
        xmm.remove_row(j)
        count -= 1
        marker = 1 # pra quando sair do loop, fazer marker==1, pra poder mudar o i
        continue #o continue leva pro inicio desse loop "j<count"
    else:
      j += 1
      marker = 1
  if marker==1:
    i += 1


#Convertendo ra, dec pra h,min,sec e deg,min,sec:
#RA
ra_h = xmm['RA(deci)']//15
ra_h = ra_h.astype(int)
restoh =  xmm['RA(deci)']/15 - ra_h
ra_min = restoh*60
ra_min = ra_min.astype(int)
restom = restoh*60 - ra_min
ra_sec = restom*60


#DEC
dec_d = xmm['DEC(deci)'].astype(int)
restod = xmm['DEC(deci)'] - dec_d
dec_min = restod*60
dec_min = dec_min.astype(int)
restomin = restod*60 - dec_min
dec_sec = restomin*60 

RAh, DECh = [],[]
count = len(xmm['ID'])
i = 0

while i<count:
  aux_ra = '%d'%(ra_h[i])+'h'+ '%d'%(ra_min[i])+'m'+ '%.2f'%(ra_sec[i])+'s'
  aux_dec = '%d'%(dec_d[i])+'d'+ '%d'%(dec_min[i])+'m'+ '%.2f'%(dec_sec[i])+'s'
  RAh.append("%s"%(aux_ra))
  DECh.append("%s"%(aux_dec))
  i +=1

col_ra = Column(name = 'RA(h)', data = [RAh])
col_ra = col_ra.reshape(count,)
col_dec = Column(name = 'DEC(h)', data = [DECh])
col_dec = col_dec.reshape(count,)

xmm.add_columns([col_ra, col_dec])

#Uniformizando os nomes
count = len(xmm['TARGET'])
i = 0
for i in count:
  word='Abell'
  if word in xmm['TARGET'][i]:
    xmm['TARGET'][i].replace(word,'A')


#salva o arquivo
arq_sai = 'xmm_clean.csv'
le2 = path_xmm+'/'+arq_sai

xmm.write('%s'%le2, delimiter = ',', overwrite = True)










#esse programa le uma tabela em csv, filtra os dados pela 
#incerteza no redshift (zErr<=2sigmas), e cria uma nova tabela
# (data_mapa.txt), 

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import csv
import sys
from astropy.cosmology import FlatLambdaCDM
cosmo = FlatLambdaCDM(H0=70, Om0=0.3)
import astropy.units as un
from scipy import  integrate
import pandas as pd
from cor_k import calc_kcor

def Luminosity(Mag):
  L=(10**(-Mag/2.5))*4*np.pi*100
  return L

def abs_mag(ap_mag, z, K_corr): #Michael R. Blanton and Sam Roweis, 2006
  #valor de M
  print('Cálculo da Mag_absoluta')
  M=[]
  for d in range(len(z)):
    DL=(cosmo.luminosity_distance(z[d]).value)*1e6 #1e6 para passar de Mpc para pc
    DM=5*np.log10(DL/10)
    M_aux = ap_mag[d] - DM - K_corr[d] 

    M.append(M_aux)

  return M

#=======================================================
#entradas

def filtra_z(file_in,z_ag):

  # file_in = 'A2631.csv'
  # z_ag = 0.27300
  fil=pd.read_csv(file_in, delimiter=',')

  #============================================
  #LIMITES DE Z
  sigma=0.08

  #corta fatia do aglomerado e joga fora valores com erro superior ao limite
  x=fil.query('Column3 < 0.03*(1+Column2)')
  xx=x.query('(@z_ag - @sigma) <= Column2 <= (@z_ag + @sigma)')

  #calcula correção-k, Mags absoultas e cores
  K_g_gr_corr=[]
  for q in range(len(xx['dered_g'].values)):
    aux_cor=xx['dered_g'].values[q]-xx['dered_r'].values[q]
    corr=calc_kcor('g', xx['Column2'].values[q], 'g - r', aux_cor)
    K_g_gr_corr.append(corr)

  # cor r
  K_r_gr_corr=[]
  for q in range(len(xx['dered_g'].values)):
    aux_cor=xx['dered_g'].values[q]-xx['dered_r'].values[q]
    corr=calc_kcor('r', xx['Column2'].values[q], 'g - r', aux_cor)
    K_r_gr_corr.append(corr)

  # cor i
  K_i_gi_corr=[]
  for q in range(len(xx['dered_g'].values)):
    aux_cor=xx['dered_g'].values[q]-xx['dered_i'].values[q]
    corr=calc_kcor('i', xx['Column2'].values[q], 'g - i', aux_cor)
    K_i_gi_corr.append(corr)

  # Magnitude absoluta corrigida das galaxias
  xx['Kcorr_g']=K_g_gr_corr
  xx['Kcorr_i']=K_i_gi_corr
  xx['Kcorr_i']=K_r_gr_corr
  xx['Mg']= abs_mag(xx['dered_g'].values,  xx['Column2'].values, K_g_gr_corr)
  xx['Mi'] = abs_mag(xx['dered_i'].values,  xx['Column2'].values, K_i_gi_corr)
  xx['Mr'] = abs_mag(xx['dered_r'].values,  xx['Column2'].values, K_r_gr_corr)
  xx['gr']=xx['Mg'].values-xx['Mr'].values
  xx['gi']=xx['Mg'].values-xx['Mi'].values
  xx['ri']=xx['Mr'].values-xx['Mi'].values

  # Luminosidades
  xx['Lumr']=Luminosity(xx.Mr.values)
  xx['Lumi']=Luminosity(xx.Mi.values)
  xx['Lumg']=Luminosity(xx.Mg.values)

  #===============================================
  #SALVA

  out_1="data_med_mapa_n_v2.csv"
  xx.to_csv(out_1)

  out_2="data_med_mapa_n_2_v2.txt"
  cols=[xx['ra'].values,xx['dec'].values]
  list_df = pd.DataFrame([ pd.Series(value) for value in cols ]).transpose()
  list_df.to_csv(out_2, sep=' ', header=None, index=False)

  return








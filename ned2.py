#!/usr/bin/python
#coding: utf-8
#Ultima atualização 20/09/2017

# ...
# Esse script faz: Encontra Redshift no NED(database) para grupos/aglomerado
# Entrada: um arquivo com Target_Name, RA, DEC
# Como funciona: Procura primeiro o nome de um grupo/aglomerado no NED se não encontrar
# Procura por RA,DEC.
# ...

import os
import astropy.io.ascii as at
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit
import math
import urllib2



## ================================================
# Carregando o arquivo
#entrada = 'chandra_sdss.csv'
entrada = 'chaser.csv'
saida = "chaser_z.csv"

chandra = at.read('%s'%entrada,delimiter=',',guess=False)
n = len(chandra["RA"])

z = []
name = []

for i in range(n):

    #Pega 1 nome , RA e DEC do arquivo

    alvo = chandra["Target_Name"][i]
    RA = str(chandra["RA"][i])
    DEC = str(chandra["Dec"][i])

    #Define o URL de busca por nome e por RA,DEC
    url1 = "https://ned.ipac.caltech.edu/cgi-bin/objsearch?objname="+alvo+"&radius=1.0&ot_include=ANY&in_objtypes1=GGroups&in_objtypes1=GClusters&ex_objtypes1=Galaxies&nmp_op=ANY&out_csys=Equatorial&out_equinox=J2000.0&obj_sort=Distance+to+search+center&of=ascii_bar&zv_breaker=30000.0&list_limit=5&img_stamp=YES&search_type=Near+Name+Search"
    url2 = "https://ned.ipac.caltech.edu/cgi-bin/objsearch?in_csys=Equatorial&in_equinox=J2000.0&lon="+ RA +"d&lat="+ DEC +"d&radius=15.0&ot_include=ANY&in_objtypes1=GGroups&in_objtypes1=GClusters&ex_objtypes1=Galaxies&nmp_op=ANY&out_csys=Equatorial&out_equinox=J2000.0&obj_sort=Distance+to+search+center&of=ascii_bar&zv_breaker=30000.0&list_limit=5&img_stamp=YES&search_type=Near+Position+Search"

    #Algoritimo de url, parecido com wget para bash shell
    attempts = 0
    while attempts < 5:
        try:
            response = urllib2.urlopen(url1, timeout = 5)
            content = response.read()
            f = open( "saida_ned_name.txt", 'w' )
            f.write( content )
            f.close()
            break
        except urllib2.URLError as e:
            attempts += 1
            print type(e)

    #Teste de erro : Numero de linhas
    with open("saida_ned_name.txt") as f:
        for j, l in enumerate(f):
            pass
    #Se o objeto foi encontrado por nome
    if not j+1 <= 22:
        data1 = np.genfromtxt("saida_ned_name.txt",delimiter='|',skip_header=23, usecols=[1],dtype="str")
        data2 = np.genfromtxt("saida_ned_name.txt",delimiter='|',skip_header=23, usecols=[6])
        out1 = data1[1]
        out2 = data2[1]
        name.append(out1)
        z.append(out2)
        print( out2 )
    #Se o objeto não foi encontrado por nome
    else:
        print( "name not found" )
        # Busca por RA, DEC
        attempts = 0
        while attempts < 5:
            try:
                response = urllib2.urlopen(url2, timeout = 5)
                content = response.read()
                f = open( "saida_ned_RADEC.txt", 'w' )
                f.write( content )
                f.close()
                break
            except urllib2.URLError as e:
                attempts += 1
                print type(e)
        #Teste de erro : Numero de linhas
        with open("saida_ned_RADEC.txt") as f:
            for j, l in enumerate(f):
                pass
        # Se objeto não foi encontrado por nome nem RA,DEC
        if j+1 <= 22:
            print("object not found")
            out3 = "no_name"
            out4 = 100
            name.append(out3)
            z.append(out4)
        #Se o objeto foi encontrado por RA,DEC
        else:
            data1 = np.genfromtxt("saida_ned_RADEC.txt",delimiter='|',skip_header=23, usecols=[1],dtype="str")
            data2 = np.genfromtxt("saida_ned_RADEC.txt",delimiter='|',skip_header=23, usecols=[6])
            res2 = data2[1]
            res1 = data1[1]
            z.append(res2)
            name.append(res1)
            print( res2 )

#Guardando o resultado em um array
z = np.array(z)
name = np.array(name)

print(len(z))
print(len(name))

#Colocando os resultados na tabela de entrada
chandra["Redshift"] = z
chandra["NED_Name"] = name

# Criando o novo catalogo
with open(saida, "w") as file:
	file.write("Obs_ID,RA,Dec,Exposure_,Redshift,Target_Name,NED_Name \n")
	n = len(chandra["Obs_ID"])
	for j in range(n):
		file.write("%i,%.3f,%.3f,%.2f,%.5f,%s,%s \n" %(chandra["Obs_ID"][j],chandra["RA"][j],chandra["Dec"][j],chandra["Exposure_"][j],chandra["Redshift"][j],chandra["Target_Name"][j],chandra["NED_Name"][j] ))


# entrada = 'chandra_sdss_z.csv'
# chandra = at.read('%s'%entrada,delimiter=',',guess=False)
# n = len(chandra["RA"])
# a = chandra["Redshift"]
#
# with open("chandra_sdss_z.csv", "w") as file:
#     file.write("Obs_ID,RA,Dec,Target_Name,Exposure_,Redshift \n")
#     for i in range(len(a)):
#         if not a[i] == 100.0 :
#             file.write("%i,%.3f,%.3f,%s,%.2f,%.5f \n" %(chandra["Obs_ID"][i],chandra["RA"][i],chandra["Dec"][i],chandra["Target_Name"][i],chandra["Exposure_"][i],a[i]))
#             print(i)

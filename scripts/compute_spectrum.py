#! /usr/bin/python3
import numpy as np
from scipy import linalg
import matplotlib.pyplot as plt
import argparse
import csv

#############################
#        Arguments          #
#############################

parser = argparse.ArgumentParser()
parser.add_argument("--ni",type=int)
parser.add_argument("--geo",choices=['emboite','non_emboite'],type=str)
parser.add_argument("--type",choices=[0,1],type=int)
args = parser.parse_args()
ni = args.ni
geo = args.geo
type = args.type

###########################
#  Chargement matrices    #
###########################
print("Chargement matrices...")
A = np.loadtxt('../output/matrices/A_'+geo+'_'+str(ni)+"_"+str(type)+'.txt').view(complex)
nc = int(np.sqrt(A.shape[0]))
print(nc)
A = np.reshape(A,(nc,nc))

####
M = np.loadtxt('../output/matrices/M_'+geo+'_'+str(ni)+"_"+str(type)+'.txt')
nc = int(np.sqrt(M.shape[0]))
M = np.reshape(M,(nc,nc))

####
P = np.loadtxt('../output/matrices/P_'+geo+'_'+str(ni)+"_"+str(type)+'.txt')
nc = int(np.sqrt(P.shape[0]))
P = np.reshape(P,(nc,nc))

###########################
#  Calcul valeurs propres #
###########################
print("Calcul valeurs propres...\n")
alphas = [-1,-0.5,0,0.5]



with open("../output/csv/spectrum_"+geo+"_"+str(ni)+"_"+str(type)+".csv", 'w') as csv_file:
    writer = csv.writer(csv_file)

    # Columns head
    writer.writerow([str(alpha)+"_"+xy for alpha in alphas for xy in ["Real","Imag"]])
    
    data = []
    for alpha in alphas:
        print(alpha)
        MTF = -0.5*(A-alpha*M-(1-alpha)*P)

        val = linalg.eig(MTF,P,right=False)
        x = val.real
        y = val.imag

        data.append(x)
        data.append(y)


    writer.writerows([list(i) for i in zip(*data)])


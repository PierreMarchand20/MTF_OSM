#! /usr/bin/python3
import numpy as np
from scipy import linalg
from scipy import sparse
from scipy import io
import matplotlib.pyplot as plt
from scipy.sparse.linalg import gmres
from scipy.sparse import coo_matrix
import argparse

#############################
#        Arguments          #
#############################

parser = argparse.ArgumentParser()
parser.add_argument("--ni",type=int)
parser.add_argument("--geo",choices=['emboite','non_emboite'],type=str)
parser.add_argument("--type",choices=[0,1,2],type=int)
args = parser.parse_args()
ni = args.ni
geo = args.geo
type = args.type

###########################
#  Chargement matrices    #
###########################
print("Chargement matrices...")
h=0.1
A = np.loadtxt('../output/matrices/A_'+geo+'_'+str(ni)+"_"+str(type)+'.txt').view(complex)
nc = int(np.sqrt(A.shape[0]))
A = np.reshape(A,(nc,nc))
####
M = np.loadtxt('../output/matrices/M_'+geo+'_'+str(ni)+"_"+str(type)+'.txt')
nc = int(np.sqrt(M.shape[0]))
M = np.reshape(M,(nc,nc))
####
P = np.loadtxt('../output/matrices/P_'+geo+'_'+str(ni)+"_"+str(type)+'.txt')
nc = int(np.sqrt(P.shape[0]))
P = np.reshape(P,(nc,nc))
####
io.savemat("../output/matrices/P_"+geo+"_"+str(ni)+"_"+str(type)+".mat",mdict={"P_"+geo+"_"+str(ni)+"_"+str(type):P})
###########################
#  Calcul valeurs propres #
###########################
print("Calcul valeurs propres...\n")
plt.rc('text', usetex=True)
plt.rc('font', family='serif')
fig = plt.figure()
ax  = fig.gca()
count=0
for alpha in [-1,-0.5,0,0.5]:
    print(alpha)
    MTF = -0.5*(A-alpha*M-(1-alpha)*P)
    io.savemat("../output/matrices/MTF_"+geo+"_"+str(ni)+"_"+str(type)+"_"+str(count)+".mat",mdict={"MTF_"+geo+"_"+str(ni)+"_"+str(type)+"_"+str(count):MTF})
    count+=1

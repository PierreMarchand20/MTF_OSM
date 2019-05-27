#! /usr/bin/python3
import numpy as np
from scipy import linalg
from scipy.sparse.linalg import gmres
import scipy.sparse.linalg as spla
from scipy.sparse import coo_matrix
import scipy.sparse.linalg as spla
from numpy import linalg as la
import argparse
import copy as cp
import pandas

#############################
#  Compteur d'iterations    #
#############################
class gmres_counter(object):
    res = []
    def __init__(self, disp=True):
        self._disp = disp
        self.niter = 0
    def __call__(self, rk=None):
        self.niter += 1
        if self._disp:
            #print('iter %3i\trk = %s' % (self.niter, str(rk)))
            self.res.append(rk)
    def reset(self):
        self.niter = 0
        self.res.clear()



#############################
#        Arguments          #
#############################

parser = argparse.ArgumentParser()
parser.add_argument("--nbr",type=int)
parser.add_argument("--geo",choices=['emboite','non_emboite'],type=str)
parser.add_argument("--type",choices=[0,1,2,3,4,5],type=int)
args = parser.parse_args()
nbr = args.nbr
geo = args.geo
type = args.type


#############################
#         Calculs           #
#############################
dico = {}
for alpha in [-1,-0.5,0,0.5]:
    # iterations_no_precond = np.zeros((nbr,1))
    dico[str(alpha)] = np.zeros(nbr)


    for ni in range(1,len(dico[str(alpha)])+1):
        print(alpha,ni)
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
        b   = np.loadtxt('../output/matrices/f_'+geo+'_'+str(ni)+"_"+str(type)+'.txt').view(complex)
        uex = np.loadtxt('../output/matrices/uex_'+geo+'_'+str(ni)+"_"+str(type)+'.txt').view(complex)


        ###
        row = np.loadtxt('../output/matrices/Prec_'+geo+'_'+str(ni)+"_"+str(type)+'.txt',usecols=(0))
        col = np.loadtxt('../output/matrices/Prec_'+geo+'_'+str(ni)+"_"+str(type)+'.txt',usecols=(1))
        val = np.loadtxt('../output/matrices/Prec_'+geo+'_'+str(ni)+"_"+str(type)+'.txt',usecols=(2,3)).view(complex)
        val = val[:,0]
        Prec = coo_matrix((val,(row,col)), shape=(nc,nc)).tocsr()

        ###
        row = np.loadtxt('../output/matrices/Mass_'+geo+'_'+str(ni)+"_"+str(type)+'.txt',usecols=(0))
        col = np.loadtxt('../output/matrices/Mass_'+geo+'_'+str(ni)+"_"+str(type)+'.txt',usecols=(1))
        val = np.loadtxt('../output/matrices/Mass_'+geo+'_'+str(ni)+"_"+str(type)+'.txt',usecols=(2,3)).view(complex)
        val = val[:,0]
        Mass = coo_matrix((val,(row,col)), shape=(nc,nc)).tocsr()

        ###
        temp = Mass.dot(uex)
        norm=np.vdot(temp,uex)

        ###########################
        #   Lancement solveur     #
        ###########################
        print("Lancement solveur...")
        R_x = lambda x: spla.spsolve(Prec, x)
        R = spla.LinearOperator((nc, nc), R_x)
        counter = gmres_counter()
        counter.reset()
        MTF = -0.5*(A-alpha*M-(1-alpha)*P)
        x,info = gmres(MTF, b, M=R, restart=200, callback=counter, maxiter=1e4)
        dico[str(alpha)][ni-1]=cp.deepcopy(counter.niter)
        err = Mass.dot(x-uex)
        print(np.vdot(err,x-uex)/norm,counter.niter)
        # counter.reset()
        # x,info = gmres(MTF, b, restart=200, callback=counter, maxiter=1e4)
        # iterations_no_precond[ni-1]=cp.deepcopy(counter.niter)
        # err = Mass.dot(x-uex)
        # print(np.vdot(err,x-uex)/norm)


df = pandas.DataFrame(dico)
df.to_csv('../output/csv/gmres_iterations_'+geo+'_'+str(nbr)+"_"+str(type)+'.csv',index=False)
#! /usr/bin/python3
import argparse
import itertools
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
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
parser.add_argument("--hsize",type=float,default=7)
parser.add_argument("--wsize",type=float,default=10)
parser.add_argument("--show",type=int,choices=[0,1],default=1)
parser.add_argument("--save",type=str,default="")
args = parser.parse_args()
nbr = args.nbr
geo = args.geo
type = args.type


#############################
#          Plot             #
#############################
plt.style.use('articles')
marker = itertools.cycle(('s','P','^','v'))
fig, axes = plt.subplots(1,1,figsize=(args.wsize,args.hsize))
axes.xaxis.set_major_locator(MaxNLocator(integer=True))
axes.yaxis.set_major_locator(MaxNLocator(integer=True))

df = pandas.read_csv("../output/csv/gmres_iterations_"+geo+"_"+str(nbr)+"_"+str(type)+".csv")

#  Get alpha
alphas = set()
for header in list(df):
    alphas.add(header)

# Plot
for alpha in sorted(alphas,key=float):
    print(alpha)

    axes.plot(range(1,len(df.loc[:,str(alpha)])+1),df.loc[:,str(alpha)],label=r"$\alpha = $"+alpha,marker=next(marker))

fig.legend()
axes.grid(True)
axes.set_xlabel("Number of interfaces")
axes.set_ylabel("Number of iterations")
bottom, top = axes.get_ylim()
left, right = axes.get_xlim()
if top>4000:
    axes.set_ylim(bottom=0,top=4000)
if right>1000:
    axes.set_xlim(left=0,right=1000)

if args.show:
    plt.show()

if args.save!="":
    fig.savefig(args.save+'/iterations_'+geo+'_'+str(nbr)+"_"+str(type)+'.png',bbox_inches = "tight")

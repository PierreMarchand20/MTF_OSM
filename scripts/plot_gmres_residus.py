#! /usr/bin/python3

import argparse
import itertools
import matplotlib.pyplot as plt
import pandas





#############################
#        Arguments          #
#############################

parser = argparse.ArgumentParser()
parser.add_argument("--ni",type=int)
parser.add_argument("--geo",choices=['emboite','non_emboite'],type=str)
parser.add_argument("--type",choices=[0,1,2,3,4,5],type=int)
parser.add_argument("--hsize",type=float,default=7)
parser.add_argument("--wsize",type=float,default=10)
parser.add_argument("--show",type=int,choices=[0,1],default=1)
parser.add_argument("--save",type=str,default="")
args = parser.parse_args()
ni = args.ni
geo = args.geo
type = args.type

#############################
#           plot            #
#############################
plt.style.use('articles')
marker = itertools.cycle(('s','P','^','v'))
fig, axes = plt.subplots(1,1,figsize=(args.wsize,args.hsize))

df = pandas.read_csv("../output/csv/gmres_residus_"+geo+"_"+str(ni)+"_"+str(type)+".csv")

#  Get alpha
alphas = set()
for header in list(df):
    alphas.add(header.split("_")[0])

# Plot
for alpha in sorted(alphas,key=float):
    print(alpha)
    current_marker = next(marker)

    axes.plot(range(1,len(df.loc[:,str(alpha)+"_preconditioned"])+1),df.loc[:,str(alpha)+"_preconditioned"],label="Preconditioned "+r"$\alpha = $"+alpha,marker=current_marker)

    axes.plot(range(1,len(df.loc[:,str(alpha)+"_no_preconditioned"])+1),df.loc[:,str(alpha)+"_no_preconditioned"],label="No preconditioner "+r"$\alpha = $"+alpha,marker=current_marker)


fig.legend()
axes.set_yscale('log', basey = 10)
axes.grid(True)
axes.set_xlabel("Number of iterations")
axes.set_ylabel("Residual")
bottom, top = axes.get_ylim()
left, right = axes.get_xlim()

if right>1000:
    axes.set_xlim(left=0,right=1000)


if args.show:
    plt.show()

if args.save!="":
    fig.savefig(args.save+'/residus_'+geo+'_'+str(ni)+"_"+str(type)+'.png',bbox_inches = "tight")

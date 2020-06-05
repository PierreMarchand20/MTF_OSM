#! /usr/bin/python3

import argparse
import pandas
import matplotlib.pyplot as plt
from matplotlib.markers import MarkerStyle





#############################
#        Arguments          #
#############################

parser = argparse.ArgumentParser()
parser.add_argument("--ni",type=int)
parser.add_argument("--geo",choices=['emboite','non_emboite'],type=str)
parser.add_argument("--type",choices=[0,1,2],type=int)
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
fig, axes = plt.subplots(1,1,figsize=(args.wsize,args.hsize))


df = pandas.read_csv("../output/csv/spectrum_"+geo+"_"+str(ni)+"_"+str(type)+".csv")
print(list(df))

#  Get alpha
alphas = set()
for header in list(df):
    alphas.add(header.split("_")[0])



# Plot
for alpha in sorted(alphas,key=float):
    print(alpha)
    axes.scatter(df.loc[:,str(alpha)+"_Real"],df.loc[:,str(alpha)+"_Imag"],label=r"$\alpha = $"+alpha,marker='x')

fig.legend()
axes.grid(True)
axes.set_xlabel(r"Real($\lambda$)")
axes.set_ylabel(r"Imag($\lambda$)")

if args.show:
    plt.show()

if args.save!="":
    fig.savefig(args.save+'/spectrum_'+geo+'_'+str(ni)+"_"+str(type)+'.png',bbox_inches = "tight")

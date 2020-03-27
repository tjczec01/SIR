# -*- coding: utf-8 -*-
"""
Created on Thu Mar 26 22:55:37 2020

@author: Travis Czechorski tjcze01@gmail.com
"""

import matplotlib.pyplot as plt 
from scipy.integrate import solve_ivp
import os

clear = lambda: os.system('cls')
cwd = os.getcwd()
dir_path = os.path.dirname(os.path.realpath(__file__))
path_fol = "{}\SEIR Model for Spread of Disease".format(dir_path)
try:
       os.mkdir(path_fol)
except:
       pass

def R0(α, β, μ, λ):
       R_0 = (α/(μ + α))*(β/(μ + λ))
       return R_0

def SEIR(t, y, *args):
       a, b, u, lamd, gamm,  np = args
       BN = float(b/np)
       dsdt = gamm - u*y[0] - BN*y[2]*y[0]
       dedt = BN*y[2]*y[0] - (u + a)*y[1]
       didt = a*y[1] - (u + lamd)*y[2]
       drdt = lamd*y[2] - u*y[3]
       return [dsdt, dedt, didt, drdt]

b = 0.45/100.0 # Chance of getting infected per contact
r = 0.5 # Contacts per day
rb = r*b
Λ = 0.0 # Birth rate
μ = 0.0 # Death rate
Tc = 1.0/r # Typical time between contacts
β = 1.0/Tc
Tr = 11.1 # Typical time until recovery
λ = 1.0/Tr
Ti = 5.2 # Average incubation period
α = Tr/Tc
Infi = 100 # Initial infected
Daysnn = 100
NP = 329436928
S0 = NP - Infi
its = 100000
itern = Daysnn/its
Days = [0.0,Daysnn]
Time = [i for i in range(0, int(Daysnn + 1), 1)]
tt = list(range(0,its,1))
Time_f = [i*itern for i in tt]
Y0 = [NP, Infi, 0.0, 0.0]

answer = solve_ivp(SEIR, Days, Y0, t_eval=Time_f, method = 'RK45', args=(α, β, μ, λ, Λ, NP))

Sn = answer.y[0]
En = answer.y[1]
In = answer.y[2]
Rn = answer.y[3]
Spb = answer.y[0]
Epb = answer.y[1]
Ipb = answer.y[2]
Rpb = answer.y[3]

Sp = [(i/NP)*100.0 for i in answer.y[0]]
Ep = [(i/NP)*100.0 for i in answer.y[1]]
Ip = [(i/NP)*100.0 for i in answer.y[2]]
Rp = [(i/NP)*100.0 for i in answer.y[3]]

m = max(In)
mi = (In.tolist()).index(max(In))
mip = mi/its
peakn = round(Daysnn*mip)
my = max(Ip)
myi = (Ip).index(max(Ip))
myp = myi/its
peakyn = round(Daysnn*myp)

fig = plt.figure()
plt.plot(Time_f, Sn, 'b-', label=r'$\it{Susceptible}$')
plt.plot(Time_f, En, 'c-', label=r'$\it{Exposed}$')
plt.plot(Time_f, In, 'r-', label=r'$\it{Infected}$')
plt.plot(Time_f, Rn, 'g-', label=r'$\it{Recovered}$')
plt.axvline(x=int(Time_f[mi]), color='k', linestyle='--')
plt.legend([r'$\it{Susceptible}$', r'$\it{Exposed}$', r'$\it{Infected}$', r'$\it{Recovered}$', r'$\it{}$'.format("{}{}{}".format('{',"Peak \ = {}{}{} \ Days".format('{',r'\ {}'.format(peakn),'}'),'}'))],loc="best", fontsize=15)
plt.xlim((0,Daysnn))
plt.ylim((0,NP))
plt.gca().get_yaxis().get_major_formatter().set_scientific(False)
plt.gca().set_yticklabels([r'{:,}'.format(int(x)) for x in plt.gca().get_yticks()])
plt.xlabel(r'$\bf{Time \ [Days]}$', fontsize=15)
plt.ylabel(r'$\bf{Number \ of \ people}$', fontsize=15)
plt.title(r'$\bf{SEIR \ Method  \ for \ Spread \ of \ Disease}$', fontsize=18)
plt.grid()
fig.savefig(r"{}\SEIR Population.pdf".format(path_fol), bbox_inches='tight')
fig.savefig(r"{}\SEIR Population.svg".format(path_fol), bbox_inches='tight')

fig = plt.figure()
plt.plot(Time_f, Sp, 'b-', label=r'$\it{Susceptible}$')
plt.plot(Time_f, Ep, 'c-', label=r'$\it{Exposed}$')
plt.plot(Time_f, Ip, 'r-', label=r'$\it{Infected}$')
plt.plot(Time_f, Rp, 'g-', label=r'$\it{Recovered}$')
plt.axvline(x=int(Time_f[myi]), color='k', linestyle='--')
plt.legend([r'$\it{Susceptible}$', r'$\it{Exposed}$', r'$\it{Infected}$', r'$\it{Recovered}$', r'$\it{}$'.format("{}{}{}".format('{',"Peak \ = {}{}{} \ Days".format('{',r'\ {}'.format(peakn),'}'),'}'))],loc="best", fontsize=15)
plt.xlim((0,Daysnn))
plt.ylim((0,100))
plt.gca().set_yticklabels([f'{int(y)}%' for y in plt.gca().get_yticks()])
plt.gca().set_xticklabels([f'{int(x)}' for x in plt.gca().get_xticks()])
plt.xlabel(r'$\bf{Time \ [Days]}$', fontsize=15)
plt.ylabel(r'$\bf{Percentage \ of \ population}$', fontsize=15)
plt.title(r'$\bf{SEIR \ Method  \ for \ Spread \ of \ Disease}$', fontsize=18)
plt.grid()
fig.savefig(r"{}\SEIR Percent.pdf".format(path_fol), bbox_inches='tight')
fig.savefig(r"{}\SEIR Percent.svg".format(path_fol), bbox_inches='tight')

# -*- coding: utf-8 -*-
"""
Created on Fri Feb 21 07:45:26 2020

@author: Travis Czechorski tjcze01@gmail.com
"""

import matplotlib.pyplot as plt 
import sympy as sp
from scipy.integrate import solve_ivp
from tkinter import Tk, StringVar, ttk, N, W, E, S
import os

clear = lambda: os.system('cls')
cwd = os.getcwd()
dir_path = os.path.dirname(os.path.realpath(__file__))
path_fol = "{}\SIR Model for Spread of Disease".format(dir_path)
try:
       os.mkdir(path_fol)
except:
       pass

inits = []

def close_window():
      global entry
      entry = 1.0/float(bs.get())
      inits.append(entry)
      entry2 = 1.0/float(ks.get())
      inits.append(entry2)
      entry3 = float(Ns.get())
      inits.append(entry3)
      entry4 = float(Is.get())
      inits.append(entry4)
      entry5 = float(Ds.get())
      inits.append(entry5)
      root.destroy()

root = Tk()
root.title("SIR Method Initial Values")
mainframe = ttk.Frame(root, padding="3 3 12 12")
mainframe.grid(column=0, row=0, sticky=(N, W, E, S))
root.columnconfigure(0, weight=1)
root.rowconfigure(0, weight=1)	
Bb = StringVar()
bs = ttk.Entry(mainframe, width=7, textvariable=Bb)
bs.grid(column=2, row=1, sticky=(W, E))
ttk.Label(mainframe, text="Average number of days an infected person makes infecting contact").grid(column=1, row=1, sticky=W)
Kk = StringVar()
ks = ttk.Entry(mainframe, width=7, textvariable=Kk)
ks.grid(column=2, row=2, sticky=(W, E))
ttk.Label(mainframe, text="Average duration of infection (Days)").grid(column=1, row=2, sticky=W)
Nn = StringVar()
Ns = ttk.Entry(mainframe, width=7, textvariable=Nn)
Ns.grid(column=2, row=3, sticky=(W, E))
ttk.Label(mainframe, text="Population number").grid(column=1, row=3, sticky=W)
In = StringVar()
Is = ttk.Entry(mainframe, width=7, textvariable=In)
Is.grid(column=2, row=4, sticky=(W, E))
ttk.Label(mainframe, text="Initial Infected").grid(column=1, row=4, sticky=W)
Dn = StringVar()
Ds = ttk.Entry(mainframe, width=7, textvariable=Dn)
Ds.grid(column=2, row=5, sticky=(W, E))
ttk.Label(mainframe, text="Number of Days").grid(column=1, row=5, sticky=W)
B = ttk.Button(root, text = "OK", command = close_window).grid(column=3, row=1)
root.mainloop()

def SIR(t, y, *args):
       b, k = args
       dsdt = -b*y[0]*y[1]
       didt = b*y[0]*y[1] - k*y[1]
       drdt = k*y[1]
       
       return [dsdt, didt, drdt]

bb = inits[0]
kk = inits[1]
Nu = inits[2]
In = inits[3]
Daysn = inits[4]
# bb = 1.0/2.0
# kk = 1.0/14.0
# Nu = 329433166
# In = 10.0
# Daysn = 120
argsl = [bb, kk]
Y0 = [1.0, In/float(Nu), 0.0]
Days = [0.0,Daysn]
t = sp.symbols('t')
Y = [sp.symbols('S'), sp.symbols('I')]
Time = [i for i in range(0, int(Daysn + 1), 1)]
answer = solve_ivp(lambda t, Y: SIR(t, Y, *argsl), Days, Y0, t_eval=Time, method = 'RK45')
Ss = answer.y[0]*Nu
I = answer.y[1]*Nu
R = answer.y[2]*Nu
Sp = answer.y[0]*100
Ip = answer.y[1]*100
Rp = answer.y[2]*100

m = max(I)
mi = (I.tolist()).index(max(I.tolist()))
my = max(Ip)
myi = (Ip.tolist()).index(max(Ip.tolist()))

fig = plt.figure()
plt.plot(Time, Ss, 'b-', label=r'$\it{Susceptible}$')
plt.plot(Time, I, 'r-', label=r'$\it{Infected}$')
plt.plot(Time, R, 'g-', label=r'$\it{Recovered}$')
plt.axvline(x=int(mi), color='k', linestyle='--')
plt.legend([r'$\it{Susceptible}$', r'$\it{Infected}$', r'$\it{Recovered}$', r'$\it{}$'.format("{}{}{}".format('{',"Peak \ = {}{}{} \ Days".format('{',r'\ {}'.format(mi),'}'),'}'))],loc="best", fontsize=15)
plt.xlim((0,Daysn))
plt.ylim((0,Nu))
plt.gca().get_yaxis().get_major_formatter().set_scientific(False)
plt.gca().set_yticklabels([r'{:,}'.format(int(x)) for x in plt.gca().get_yticks()])
plt.xlabel(r'$\bf{Time \ [Days]}$', fontsize=15)
plt.ylabel(r'$\bf{Number \ of \ people}$', fontsize=15)
plt.title(r'$\bf{SIR \ Method  \ for \ Spread \ of \ Disease}$', fontsize=18)
plt.grid()
fig.savefig(r"{}\SIR Population.pdf".format(path_fol), bbox_inches='tight')
fig.savefig(r"{}\SIR Population.svg".format(path_fol), bbox_inches='tight')

fig = plt.figure()
plt.plot(Time, Sp, 'b-', label=r'$\it{Susceptible}$')
plt.plot(Time, Ip, 'r-', label=r'$\it{Infected}$')
plt.plot(Time, Rp, 'g-', label=r'$\it{Recovered}$')
plt.axvline(x=int(mi), color='k', linestyle='--')
plt.xlim((0,Daysn))
plt.ylim((0,100))
plt.xlabel(r'$\bf{Time \ [Days]}$', fontsize=15)
plt.ylabel(r'$\bf{Percentage \ of \ population}$', fontsize=15)
plt.legend([r'$\it{Susceptible}$', r'$\it{Infected}$', r'$\it{Recovered}$', r'$\it{}$'.format("{}{}{}".format('{',"Peak \ = {}{}{} \ Days".format('{',r'\ {}'.format(mi),'}'),'}'))],loc="best", fontsize=15)
plt.title(r'$\bf{SIR \ Method  \ for \ Spread \ of \ Disease}$', fontsize=18)
plt.gca().set_yticklabels([f'{int(y)}%' for y in plt.gca().get_yticks()])
plt.gca().set_xticklabels([f'{int(x)}' for x in plt.gca().get_xticks()])
plt.grid()
fig.savefig(r"{}\SIR Percent.pdf".format(path_fol), bbox_inches='tight')
fig.savefig(r"{}\SIR Percent.svg".format(path_fol), bbox_inches='tight')

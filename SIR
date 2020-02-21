# -*- coding: utf-8 -*-
"""
Created on Fri Feb 21 07:45:26 2020

@author: tjcze
"""

import matplotlib.pyplot as plt 
import sympy as sp
from scipy.integrate import solve_ivp
from tkinter import Tk, StringVar, ttk, N, W, E, S

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

def SIR(t, y, *args):
       b, k = args
       dsdt = -b*y[0]*y[1]
       didt = b*y[0]*y[1] - k*y[1]
       drdt = k*y[1]
       
       return [dsdt, didt, drdt]

root = Tk()
root.title("SIR Method Initial Values")
mainframe = ttk.Frame(root, padding="3 3 12 12")
mainframe.grid(column=0, row=0, sticky=(N, W, E, S))
root.columnconfigure(0, weight=1)
root.rowconfigure(0, weight=1)	
Bb = StringVar()
bs = ttk.Entry(mainframe, width=7, textvariable=Bb)
bs.grid(column=2, row=1, sticky=(W, E))
ttk.Label(mainframe, text="Average number of conctacts a day with susceptible people").grid(column=1, row=1, sticky=W)
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

bb = inits[0]
kk = inits[1]
N = inits[2]
In = inits[3]
Daysn = inits[4]
argsl = [bb, kk]
Y0 = [1.0, 10/N, 0.0]
Days = [0.0,Daysn]
t = sp.symbols('t')
Y = [sp.symbols('S'), sp.symbols('I')]
Time = [i for i in range(0, int(Daysn), 1)]
answer = solve_ivp(lambda t, Y: SIR(t, Y, *argsl), Days, Y0, t_eval=Time, method = 'RK45')
S = answer.y[0]*N
I = answer.y[1]*N
R = answer.y[2]*N

fig = plt.figure()
plt.plot(Time, S, 'b--', label='Susceptible')
plt.plot(Time, I, 'g--', label='Infected')
plt.plot(Time, R, 'r--', label='Recovered')
plt.xlabel(r'Time [Days]')
plt.ylabel(r'Number of people')
plt.legend(loc='best')
plt.title('SIR Method')
plt.grid()

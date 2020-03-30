# -*- coding: utf-8 -*-
"""
Created on Thu Mar 26 22:55:37 2020

@author: Travis Czechorski tjcze01@gmail.com
"""

import matplotlib.pyplot as plt 
from scipy.integrate import solve_ivp
import numpy as np
import math
import pandas as pd
import os
import time

start = time.time() #Real time when the program starts to run

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

def aplha(day):
       if int(day) <= 23:
              return 0
       elif int(day) > 23 and int(day) <= 29:
              return 0.4239
       else:
              return 0.8478

def Beta(σ, b, D, N, k):
       B = b*(1 - σ)*((1 - D/N)**k)
       return B

def SEIR(t, y, *args):
       σ, β, γ, μ, Λ, F, α, d, κ, λ = args
       β_t = Beta(α, β, y[5], y[4], κ)
       dsdt = Λ - μ*y[0] - ((β*F)/y[4])*y[2]*y[0] - (β_t/y[4])*y[2]*y[0] - μ*y[0]
       dedt = ((β*F)/y[4])*y[2]*y[0] + (β_t/y[4])*y[2]*y[0] - (μ + σ)*y[1]
       didt = σ*y[1] - (μ + γ)*y[2]
       drdt = γ*y[2] - μ*y[3]
       dndt = -μ*y[4]
       dDdt = d*γ*y[2] - λ*y[5]
       dcdt = σ*y[1]
       return [dsdt, dedt, didt, drdt, dndt, dDdt, dcdt]

def jacobian(t, y, *args):
        σ, β, γ, μ, Λ, F, α, d, κ, λ = args
        β_t = Beta(α, β, y[5], y[4], κ)
        return [[-F*y[2]*β/y[4] - y[2]*β_t/y[4] - 2*μ,      0, -F*y[0]*β/y[4] - y[0]*β_t/y[4],  0,  F*y[2]*y[0]*β/y[4]**2 + y[2]*y[0]*β_t/y[4]**2,  0, 0],
                [       F*y[2]*β/y[4] + y[2]*β_t/y[4], -μ - σ,  F*y[0]*β/y[4] + y[0]*β_t/y[4],  0, -F*y[2]*y[0]*β/y[4]**2 - y[2]*y[0]*β_t/y[4]**2,  0, 0],
                [                       0,      σ,             -γ - μ,  0,                            0,  0, 0],
                [                       0,      0,                  γ, -μ,                            0,  0, 0],
                [                       0,      0,                  0,  0,                           -μ,  0, 0],
                [                       0,      0,                d*γ,  0,                            0, -λ, 0],
                [                       0,      σ,                  0,  0,                            0,  0, 0]]       

def roundup(x, places):
    return int(math.ceil(x / int(places))) * int(places)

b_0 = 0.45/100.0 # Chance of getting infected per contact
r = 10.0 # Contacts per day
Λ = 0.0 # Birth rate
μ = 0.0 # Death rate
# Λ = 0.01 # Birth rate
# μ = 0.0205 # Death rate
Tc = 2.0 # Typical time between contacts
# β = 0.5944 #1.0/Tc 
β = 1.68
# Tr = 11.1 # Typical time until recovery
Tr = 14.0
γ = 1.0/Tr
iis = [5.2, 5.2, 6.1, 5.5, 4.8, 5.0, 6.5, 4.8]
Ti = sum(iis)/len(iis)
σ = Ti**-1
rb = (Tr/Tc)*b_0
F = 0
α = 0.0 #
# α = 0.4239
# α = 0.8478
d = 0.2
# k = 1117.3
# k = 100
k = 0
λb = 11.2
λ = λb**-1
Infi = 10 # Initial infected
Daysnn = 100
NP = 329436928 # 1437904257
S0 = NP - Infi
its = 10000
itern = Daysnn/its
Days = [0.0, Daysnn]
Time = [i for i in range(0, int(Daysnn + 1), 1)]
tt = list(range(0,its,1))
Time_f = [i*itern for i in tt]
Y0 = [NP, 0.0, Infi, 0.0, NP, 0.0, Infi]
Ro = R0(σ, β, μ, γ)
# print(Ro)
# print('Λ') 
# print('μ') 
# print(α)
# print(β)
# print(RB)
# print(λ)

answer = solve_ivp(SEIR, Days, Y0, t_eval=Time_f, method = 'Radau', args=(σ, β, γ, μ, Λ, F, α, d, k, λ), jac=jacobian, rtol=1E-10, atol=1E-10)
ts = answer.t

Bs = [Beta(σ, β, i, j, k) for i,j in zip(answer.y[5],answer.y[4])]

Sn = answer.y[0]
En = answer.y[1]
In = answer.y[2]
Rn = answer.y[3]
Nn = answer.y[4]
Dn = answer.y[5]
Cn = answer.y[6]
Spb = answer.y[0]/NP
Epb = answer.y[1]/NP
Ipb = answer.y[2]/NP
Rpb = answer.y[3]/NP
Npb = answer.y[4]/NP
Dpb = answer.y[5]/NP
Cpb = answer.y[6]/NP
Sp = [i*100.0 for i in Spb]
Ep = [i*100.0 for i in Epb]
Ip = [i*100.0 for i in Ipb]
Rp = [i*100.0 for i in Rpb]
Np = [i*100.0 for i in Npb]
Dp = [i*100.0 for i in Dpb]
Cp = [i*100.0 for i in Cpb]

m = max(In)
mi = (In.tolist()).index(max(In))
mip = mi/its
peakn = round(Daysnn*mip)
my = max(Ip)
myi = (Ip).index(max(Ip))
myp = myi/its
peakyn = round(Daysnn*myp)

PEAK = [int(round(Daysnn*(mi/its)))]
nPEAK = np.array(PEAK, ndmin=2)
Tdata = np.array((Time_f, Sn, En, In, Rn))
TTdata = np.array((Time_f, Spb, Epb, Ipb, Rpb, Sp, Ep, Ip, Rp))
Tdatal = Tdata.tolist()

if its <= 16384:
       writer = pd.ExcelWriter(r'{}\SIR Population.xlsx', engine='xlsxwriter')
       writerp = pd.ExcelWriter(r'{}\SIR Percent.xlsx', engine='xlsxwriter')
       indexes =  ['Time [Days]', 'Susceptible', 'Exposed', 'Infected', 'Recovered', 'Peak [Day]']
       indexest = ['Time [Days]', 'Susceptible', 'Exposed', 'Infected', 'Recovered' , 'Susceptible [%]', 'Exposed [%]', 'Infected [%]', 'Recovered [%]', 'Peak [Day]']
       df = pd.DataFrame([Time_f, Sn, En, In, Rn, PEAK], index=[*indexes])
       dft = pd.DataFrame([Time_f, Spb, Epb, Ipb, Rpb, Sp, Ep, Ip, Rp, PEAK], index=[*indexest])
       df.to_excel(r"{}\SIR Population.xlsx".format(path_fol), sheet_name="SIR Population.xlsx", header=True, startrow=1)
       dft.to_excel(r"{}\SIR Percent.xlsx".format(path_fol), sheet_name="SIR Percent.xlsx", header=True, startrow=1)

elif its > 16384 and its <= 1048576: 
       writer = pd.ExcelWriter(r'{}\SIR Population.xlsx', engine='xlsxwriter')
       writerp = pd.ExcelWriter(r'{}\SIR Percent.xlsx', engine='xlsxwriter')
       indexesb =  ['Time [Days]', 'Susceptible', 'Exposed', 'Infected', 'Recovered']
       indexestb = ['Time [Days]', 'Susceptible', 'Exposed', 'Infected', 'Recovered' , 'Susceptible [%]', 'Exposed [%]', 'Infected [%]', 'Recovered [%]']
       df = pd.DataFrame(Tdata.T, columns=[*indexesb])
       df.T
       df.loc[:,'Peak [Day]'] = pd.Series([PEAK])
       dft = pd.DataFrame(TTdata.T, columns=[*indexestb])
       dft.T
       dft.loc[:,'Peak [Day]'] = pd.Series([PEAK])
       df.to_excel(r"{}\SIR Population.xlsx".format(path_fol), sheet_name="SIR Population.xlsx", header=True, startrow=1)
       dft.to_excel(r"{}\SIR Percent.xlsx".format(path_fol), sheet_name="SIR Percent.xlsx", header=True, startrow=1)
else:
       print("To many data points, couldn't save to excel file")
       print('Maximum excel dimensions are 1,048,576 rows by 16,384 columns')

fig = plt.figure()
plt.plot(Time_f, Sn, 'b-', label=r'$\it{Susceptible}$')
plt.plot(Time_f, En, 'c-', label=r'$\it{Exposed}$')
plt.plot(Time_f, In, 'r-', label=r'$\it{Infected}$')
plt.plot(Time_f, Rn, 'g-', label=r'$\it{Recovered}$')
plt.axvline(x=int(Time_f[mi]), color='k', linestyle='--')
plt.legend([r'$\it{Susceptible}$', r'$\it{Exposed}$', r'$\it{Infected}$', r'$\it{Recovered}$', r'$\it{}$'.format("{}{}{}".format('{',"Peak \ = {}{}{} \ Days".format('{',r'\ {}'.format(peakn),'}'),'}'))],loc="best", fontsize=15)
plt.xlim((0,Daysnn))
plt.ylim((0,NP))
plt.yticks([roundup(i*(NP/10.0), 1000000) for i in range(11)])
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
plt.yticks([i*10 for i in range(11)])
plt.gca().set_yticklabels([f'{int(y)}%' for y in plt.gca().get_yticks()])
plt.gca().set_xticklabels([f'{int(x)}' for x in plt.gca().get_xticks()])
plt.xlabel(r'$\bf{Time \ [Days]}$', fontsize=15)
plt.ylabel(r'$\bf{Percentage \ of \ population}$', fontsize=15)
plt.title(r'$\bf{SEIR \ Method  \ for \ Spread \ of \ Disease}$', fontsize=18)
plt.grid()
fig.savefig(r"{}\SEIR Percent.pdf".format(path_fol), bbox_inches='tight')
fig.savefig(r"{}\SEIR Percent.svg".format(path_fol), bbox_inches='tight')
plt.show()

end = time.time() #Time when it finishes, this is real time

def timer(start,end):
    hours, rem = divmod(end-start, 3600)
    minutes, seconds = divmod(rem, 60)
    print("Completion Time: {} Hours {} Minutes {} Seconds".format(int(hours),int(minutes),int(seconds)))

timer(start,end) # Prints the amount of time passed since starting the program

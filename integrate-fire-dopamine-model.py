#Integrate-and-Fire Model of Dopamine to inform Drug Addiction

#This code is based on the paper: "A mathematical model of reward-mediated learning in drug adduction" (https://doi.org/10.1063/5.0082997)
#The model was adapted for a final project in BME 4409 Quantitative Physiology
#This program contains functions for plotting a graph similar to that found in Fig 2 of the paper and a pair of w(t|theta) and RPE plots over time.

import numpy as np
import matplotlib.pyplot as plt
from scipy import integrate
from array import array

#This function plots a graph similar to that found in Fig 2 of the referenced paper
def variedBeta(alpha, Ta, beta_values, Tb, delta):
    # Time grid
    t = np.linspace(0, 1000, 10000)
    
    # Solve for w(t|theta) and plot
    for beta in beta_values:
        solution = (Ta * delta / (alpha - 1)) * (
            (1 - Tb / (beta - 1)) * np.exp(-t) - (1 - Tb / (beta - alpha)) * np.exp(-alpha * t) - (Tb / (beta - alpha) - Tb / (beta - 1)) * np.exp(-beta * t)
        )
        plt.plot(t, solution, label=f"Tolerant beta={beta}")

    plt.xlabel("Time")
    plt.ylabel("w(t|$\Theta$)")
    plt.legend()
    plt.xlim(0, 20)
    plt.axhline(0, color='grey', linestyle='--')
    plt.show()

#This function plots the w(t|theta) and RPE responses
def plotwRPE(alpha, Ta, beta, Tb, delta, T, freqs, days, G, B):
    k = freq*days + 1  # Number of dosages the patient will take +1
    for i in range(k-1):
        b2 = beta[i] * (1 - B * delta)
        Tb2 = Tb[i] * (1 + G * delta)
        beta.append(b2)
        Tb.append(Tb2)
        T.append(freq*(i+1))
    #print(beta)
    # Time grid
    steps = 1000
    t = np.linspace(0, freq*k, steps*k)  # Covering a longer time span for better visualization
    tlarge = np.linspace(0, freqs[0]*k, steps*k)
    wt = [0]
    n=1
    save = wt[len(wt)-1]
    last = 0
    RPE = [0]
    for i in range(steps*k):
        if n==k:
            break
        elif t[i]>=T[n-1] and t[i]<T[n]:
            wt.append(save + (Ta * delta / (alpha - 1)) * (
                (1 - Tb[n-1] / (beta[n-1] - 1)) *
                np.exp(-(t[i]-T[n-1])) -
                (1 - Tb[n-1] / (beta[n-1] - alpha)) *
                np.exp(-alpha * (t[i]-T[n-1])) -
                (Tb[n-1] / (beta[n-1] - alpha) - Tb[n-1] / (beta[n-1] - 1)) *
                np.exp(-beta[n-1] * (t[i]-T[n-1]))
            ))
            x = t[len(wt)-2:len(wt)]
            RPE.append(RPE[len(RPE)-1] + integrate.trapezoid(wt[len(wt)-2:len(wt)], x))
        else:
            save = wt[len(wt)-1]
            n += 1
            last = i 
    if len(t) < len(tlarge):
        t = np.interp(tlarge, t, t)
        RPE = np.interp(tlarge, t, RPE)
    plt.subplot(211)
    plt.plot(t[range(0,len(wt))], wt, label=f"B={B}") #change label as necessary
    plt.xlabel("Time")
    plt.ylabel("w(t|$\Theta$)")
    plt.axhline(0, color='grey', linestyle='--')
    plt.tick_params(axis='both', which='major', direction='in', length=4, width=1)
    plt.tick_params(axis='both', which='minor', direction='in', length=2, width=1)
    plt.minorticks_on()
    plt.legend()
    plt.subplot(212)
    plt.plot(t[range(0,len(RPE))], RPE, label=f"B={B}") #change label as necessary
    plt.xlabel("Time")
    plt.ylabel("RPE")
    plt.axhline(0, color='grey', linestyle='--')
    plt.tick_params(axis='both', which='major', direction='in', length=4, width=1)
    plt.tick_params(axis='both', which='minor', direction='in', length=2, width=1)
    plt.minorticks_on()
    plt.legend()

# Initialize variables
alpha = 0.3
Ta = 1
beta_values = [0.5, 0.17149999999999996, 0.058824499999999974]
Tb = 0.1
delta = 3

variedBeta(alpha, Ta, beta_values, Tb, delta)

# Initialize variables
# naive values
alpha = 0.3
Ta = 1
beta = [0.5]
Tb = [0.1]
delta = 1  # Dosage in mg
T = [0]  # Frequency of intake
freqs = [6]
days = 3
k = 19  # Number of dosages the patient will take +1
G = [0,0.05,0.1]  # Genetic predisposition
B = 0.05

# tolerant values
alpha = 0.3
Ta = 1
beta = [0.5]
Tb = [0.1]
k = 19  # Number of dosages the patient will take +1
delta = 3  # Dosage in mg
T = [0]  # Frequency of intake
freqs = [2]
days = 3
G = 0.1  # Genetic predisposition
B = [0,0.05,0.1]  


# dynamic dopamine response and RPE vs time

for B in B: #change loop to plot other values
    print(B)
    # Updates values of beta and Tb for every new intake
    freq = freqs[0]
    beta = [0.5]
    Tb = [0.1]
    T = [0]
    k = freq*days + 1  # Number of dosages the patient will take +1
    plotwRPE(alpha, Ta, beta, Tb, delta, T, freqs, days, G, B)
plt.show()





# -*- coding: utf-8 -*-
"""
Created on Tue Jun  5 17:15:39 2018

@author: atkachen
"""

import numpy as np
import matplotlib.pyplot as plt
from math import log10, log, exp

file_dens = open("density_model1001.txt", "r").readlines() # we are going to read the file as a list of strings
file_el_pre = open("electronpressure_model1001.txt", "r").readlines()
file_pr_pre = open("protonpressure_model1001.txt", "r").readlines()
file_tot_pre = open("gaspressure_model1001.txt", "r").readlines()



height = []                    #these are our lists of relevant parameters in file 
electron_dens = []
proton_dens = []
hydr_dens = []
tot_dens = []

el_pres = []
ion_pres = []
tot_pres = []


for line_numb in range(2, len(file_dens)):    #we start reading the data and writing it into our lists,
    file_dens[line_numb] = file_dens[line_numb].split()#we split the strings from 2nd to last into a list of values  
    height.append(float(file_dens[line_numb][0]))
    electron_dens.append(float(file_dens[line_numb][1]))
    proton_dens.append(float(file_dens[line_numb][2]))
    hydr_dens.append(float(file_dens[line_numb][3]))
    tot_dens.append(float(file_dens[line_numb][1]))
 
for line_numb in range(2, len(file_el_pre)):    #we start reading the data and writing it into our lists,
    file_el_pre[line_numb] = file_el_pre[line_numb].split()#we split the strings from 2nd to last into a list of values  
    el_pres.append(float(file_el_pre[line_numb][5]))
    
for line_numb in range(2, len(file_pr_pre)):    #we start reading the data and writing it into our lists,
    file_pr_pre[line_numb] = file_pr_pre[line_numb].split()#we split the strings from 2nd to last into a list of values  
    ion_pres.append(float(file_pr_pre[line_numb][5]))
    
for line_numb in range(2, len(file_tot_pre)):    #we start reading the data and writing it into our lists,
    file_tot_pre[line_numb] = file_tot_pre[line_numb].split()#we split the strings from 2nd to last into a list of values  
    tot_pres.append(float(file_tot_pre[line_numb][5]))
        
"""for ind in range(len(height)):
    if(height[ind] == height[ind-1]):
        height[ind] += 1.5*10**5 """

###############interpolation of the values dens(height)
x = np.linspace(2.0*10**8, -1.0*10**7, 101)
n_e_interp = []
n_i_interp = []
n_n_interp = []

for val in x:
    for i in range(len(height)):
#        if(val>=0):
        if(val < height[i] and val > height[i+1]):
            k = i
            n_e_interp.append(exp(log(electron_dens[k])+(log(electron_dens[k+1])-log(electron_dens[k]))*(val-height[k])/(height[k+1]-height[k])))  #the interpolation formula       
            n_i_interp.append(exp(log(proton_dens[k])+(log(proton_dens[k+1])-log(proton_dens[k]))*(val-height[k])/(height[k+1]-height[k])))
            n_n_interp.append(exp(log(hydr_dens[k]) + (log(hydr_dens[k+1])-log(hydr_dens[k]))*(val-height[k])/(height[k+1]-height[k])))    
        elif(val == height[i]):
            n_e_interp.append(electron_dens[k])
            n_i_interp.append(proton_dens[k])
            n_n_interp.append(hydr_dens[k])
################################################
            
            
def pressure_derivative(h, component_pressure):
    deriv = []
    for i in range(len(h)-1):
        deriv.append(float((component_pressure[i+1]-component_pressure[i]))/(abs(h[i+1])-abs(h[i])))
    deriv.append(deriv[len(h)-2])    
    return deriv
    
###################reverse order of values
h = []
n_e = []
n_i = []
n_n = []
p_e = []
p_i = []
p_n = []

for d in range(1, len(height)+1):
    h.append(height[len(height)-d])
    n_e.append(electron_dens[len(height)-d])
    n_i.append(proton_dens[len(height)-d])
    n_n.append(hydr_dens[len(height)-d])
#    p_e.append(el_pres[len(height)-i])
#    p_i.append(ion_pres[len(height)-i])
#    p_n.append(tot_pres[len(height)-i]-ion_pres[len(height)-i]-el_pres[len(height)-i])
###########################################

############constants and gravit. acceleration (z)

m_e = 9.10938356*10**(-28)
m_p = 1.6726219*10**(-24)
m_n = m_p + m_e
G = 6.67259*10**(-8)
M_sun = 2*10**33
g = 27400
#g = []
#for i in range(len(x)):
#    g.append(G*M_sun/(7*10**10+x[i])**2)
#################################################
 
#############integrating the pressure with interpolated values of density
def trapezoidal(arr, a, b, m):        
    I = (m*g*arr[a]+m*g*arr[b])*(abs(x[b]-x[a]))/2            
    return I

el_pres_integr = []
ion_pres_integr = []
hydr_pres_integr = []

el_pres_integr.append(el_pres[0])
ion_pres_integr.append(ion_pres[0])
hydr_pres_integr.append(tot_pres[0] - el_pres[0] - ion_pres[0])

for l in range(len(x)-1):
    el_pres_integr.append(el_pres_integr[l] + trapezoidal(n_e_interp, l, l+1, m_e)) 
    ion_pres_integr.append(ion_pres_integr[l] + trapezoidal(n_i_interp, l, l+1, m_p))
    hydr_pres_integr.append(hydr_pres_integr[l] + trapezoidal(n_n_interp, l, l+1, m_n))

######################################################################    
rho_e = []
rho_i = []
rho_n = []

for j in range(len(n_e_interp)):
    rho_e.append(m_e*n_e_interp[j])
    rho_i.append(m_p*n_i_interp[j])
    rho_n.append(m_n*n_n_interp[j])

####################################################writing to a file
file_wr = open("Atmospheric_model_log_lin_interp_const_g.dat", "w")
file_wr.write("100")
file_wr.write('\n')
file_wr.write(str(-7.9*10**6)+"   ")
file_wr.write(str(2.0*10**8)+"   ")
file_wr.write(str(2.1*10**6)+"   ")
file_wr.write('\n')
for z in range(1, len(x)):
    file_wr.write("%+1.5e" %x[len(x)-1-z] + "   ")
    file_wr.write("%-1.5e" %rho_e[len(x)-1-z] + "   ")
    file_wr.write("%-1.5e" %rho_i[len(x)-1-z] + "   ")
    file_wr.write("%-1.5e" %rho_n[len(x)-1-z] + "   ")
    file_wr.write("%-1.5e" %el_pres_integr[len(x)-1-z] + "   ")
    file_wr.write("%-1.5e" %ion_pres_integr[len(x)-1-z] + "   ")
    file_wr.write("%-1.5e" %hydr_pres_integr[len(x)-1-z] + "   ")
    file_wr.write('\n')

file_wr.close()    
    


#####################################################################
    
#el_pres_grad = pressure_derivative(height, pres_integr)
#deviat = []

#for i in range(len(height)):
#    deviat.append(log10(abs(el_pres_grad[i] - m_e*electron_dens[i]*g[i])))
    
#deviat_log = []

#for i in range(len(height)):
#    deviat_log.append(log10(abs(deviat[i])))

#diff = []
#diff.append(0)
#for i in range(1, len(height)):
#    diff.append(log10(abs(el_pres[i]-pres_integr[i])))
#plt.plot(height, deviat)
#plt.plot(height, diff)
#plt.xlim([0, 2*10**8])
#plt.ylim([0, 0.3])

#height_cm = np.array(height) #creating numpy arrays of data
#tot_dens_cm_3 = np.array(tot_dens)
#electr_dens_cm_3 = np.array(electron_dens)
#proton_dens_cm_3 = np.array(proton_dens)
#hydr_dens_cm_3 = np.array(hydr_dens)



#height_kkm = height_cm / 10 ** 8



#plt.scatter(temperat_K, height_kkm, s = electr_, c = "g", alpha = 0.8)
#plt.xscale("log")
#plt.yticks([0.0, 0.5, 1.0, 1.5, 2.0], ["0", "50M", "100M", "150M", "200M"])
#plt.xlabel("Height, cm")
#plt.ylabel("dp/dz - rho*g, log scale")
#plt.title("Deviation from equilibrium")
#plt.grid(True)
#plt.savefig("Deviation_from_equilibrium_for_electron_pressure_in_Fontenla_model_for_integrated.png")
#plt.show()

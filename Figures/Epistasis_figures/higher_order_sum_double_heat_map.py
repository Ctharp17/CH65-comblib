import numpy as np
import matplotlib.pyplot as plt
import csv
import sys
import seaborn as sns
from scipy.special import comb
from sklearn.preprocessing import PolynomialFeatures
from matplotlib.colors import LogNorm
import matplotlib as mpl
import pandas as pd
import itertools as it
import matplotlib.lines as lines
from matplotlib.patches import Patch
from matplotlib.tri import Triangulation

def powerset(iterable):
    "powerset([1,2,3]) --> () (1,) (2,) (3,) (1,2) (1,3) (2,3) (1,2,3)"
    s = list(iterable)
    return it.chain.from_iterable(it.combinations(s, r) for r in range(len(s)+1))

mpl.rc_file_defaults()
plt.rcParams.update({'font.size': 7})
plt.rcParams['font.sans-serif'] = "Arial"
plt.rcParams['font.family'] = "sans-serif"
plt.rcParams['axes.linewidth'] = 0.5
plt.rcParams["xtick.major.size"] = 2
plt.rcParams["ytick.major.size"] = 2



num_mutations_G189E = 16
order_G189E = 5
num_term_list_G189E = np.array([int(comb(num_mutations_G189E,i)) for i in range(1,order_G189E+1)])
total_params_G189E = sum(num_term_list_G189E)
order_start_indices_G189E = list(np.cumsum(num_term_list_G189E)+1)
order_start_indices_G189E.insert(0,1)

num_mutations_MA90 = 16
order_MA90 = 4
num_term_list_MA90 = np.array([int(comb(num_mutations_MA90,i)) for i in range(1,order_MA90+1)])
total_params_MA90 = sum(num_term_list_MA90)
order_start_indices_MA90 = list(np.cumsum(num_term_list_MA90)+1)
order_start_indices_MA90.insert(0,1)

num_mutations_SI06 = 15
order_SI06 = 5
num_term_list_SI06 = np.array([int(comb(num_mutations_SI06,i)) for i in range(1,order_SI06+1)])
total_params_SI06 = sum(num_term_list_SI06)
order_start_indices_SI06 = list(np.cumsum(num_term_list_SI06)+1)
order_start_indices_SI06.insert(0,1)


full_mut_names = ['N26D','S29R','Y35N','Y48C','D49Y','V98I','G31D','Y33H','M34I','H35N','N52H','G57D','L83V','S84N','R85G','R87K']

ep_type = 'biochem'

coefs_G189E = np.zeros(total_params_G189E+1)
names_G189E = []
sig_G189E = np.full((total_params_G189E+1),0)
stderr_G189E = np.zeros(total_params_G189E+1)
cis_G189E = np.zeros((total_params_G189E+1,2))

with open('../../Epistasis_Inference/G189E/CH65_G189E_102022_'+str(order_G189E)+'order_'+ep_type+'.txt','r') as readfile:
    coef_reader = csv.reader(readfile,delimiter='\t')
    num_params = int(next(coef_reader)[-1])
    r2_train = float(next(coef_reader)[-1])
    header = next(coef_reader)
    print(header)
    for i in range(total_params_G189E+1):
        row = next(coef_reader)
        names_G189E.append(row[0])
        coefs_G189E[i] = float(row[1])  
        if i >= 1:
            stderr_G189E[i] = float(row[2])
            cis_G189E[i,0] = float(row[4])
            cis_G189E[i,1] = float(row[5])
            if float(row[4])*float(row[5]) > 0:
                sig_G189E[i] = 1
    readfile.close()
            
print(len(coefs_G189E))    


coefs_MA90 = np.zeros(total_params_MA90+1)
names_MA90 = []
sig_MA90 = np.full((total_params_MA90+1),0)
stderr_MA90 = np.zeros(total_params_MA90+1)
cis_MA90 = np.zeros((total_params_MA90+1,2))

with open('../../Epistasis_Inference/MA90/CH65_MA90_102022_'+str(order_MA90)+'order_'+ep_type+'.txt','r') as readfile:
    coef_reader = csv.reader(readfile,delimiter='\t')
    num_params = int(next(coef_reader)[-1])
    r2_train = float(next(coef_reader)[-1])
    header = next(coef_reader)
    print(header)
    for i in range(total_params_MA90+1):
        row = next(coef_reader)
        names_MA90.append(row[0])
        coefs_MA90[i] = float(row[1])  
        if i >= 1:
            stderr_MA90[i] = float(row[2])
            cis_MA90[i,0] = float(row[4])
            cis_MA90[i,1] = float(row[5])
            if float(row[4])*float(row[5]) > 0:
                sig_MA90[i] = 1
    readfile.close()
            
print(len(coefs_MA90))   


coefs_SI06 = np.zeros(total_params_SI06+1)
names_SI06 = []
sig_SI06 = np.full((total_params_SI06+1),0)
stderr_SI06 = np.zeros(total_params_SI06+1)
cis_SI06 = np.zeros((total_params_SI06+1,2))

with open('../../Epistasis_Inference/SI06/CH65_SI06_102022_'+str(order_MA90)+'order_'+ep_type+'.txt','r') as readfile:
    coef_reader = csv.reader(readfile,delimiter='\t')
    num_params = int(next(coef_reader)[-1])
    r2_train = float(next(coef_reader)[-1])
    header = next(coef_reader)
    print(header)
    for i in range(total_params_SI06+1):
        row = next(coef_reader)
        names_SI06.append(row[0])
        coefs_SI06[i] = float(row[1])  
        if i >= 1:
            stderr_SI06[i] = float(row[2])
            cis_SI06[i,0] = float(row[4])
            cis_SI06[i,1] = float(row[5])
            if float(row[4])*float(row[5]) > 0:
                sig_SI06[i] = 1
    readfile.close()
            
print(len(coefs_SI06))   
# initialize matrices to store values

# total (lower diagonal)
total_epistasis = np.zeros((16,16),dtype=float)
epistasis_sum_MA90 = np.zeros((16,1),dtype=float)
epistasis_sum_G189E = np.zeros((16,1),dtype=float)
epistasis_sum_SI06 = np.zeros((16,1),dtype=float)

epistasis_first_MA90 = np.zeros((16,1),dtype=float)
epistasis_first_G189E = np.zeros((16,1),dtype=float)
epistasis_first_SI06 = np.zeros((16,1),dtype=float)
for i in range(16):  
    total_epistasis[i,i] = np.nan


# add up all coefficients
for i in range(1,len(coefs_G189E)):

    muts_involved = [int(x)-1 for x in names_G189E[i].split(',')]
    
    # only consider 3rd order and higher
    if len(muts_involved) >= 2:
        # only consider significant terms
        if sig_G189E[i]:
            for j in range(len(muts_involved)):
                epistasis_sum_G189E[muts_involved[j]] += coefs_G189E[i]  
                for k in range(j+1,len(muts_involved)):
                    total_epistasis[muts_involved[k],muts_involved[j]] += coefs_G189E[i]
                    #total_epistasis[muts_involved[j],muts_involved[k]] += np.abs(coefs_G189E[i])

# add up all SI06 coefficients
for i in range(1,len(coefs_SI06)):

    muts_involved = [int(x)-1 for x in names_SI06[i].split(',')]
    
    # only consider 3rd order and higher
    if len(muts_involved) >= 2:
        # only consider significant terms
        if sig_SI06[i]:
            for j in range(len(muts_involved)):
                epistasis_sum_SI06[muts_involved[j]] += coefs_SI06[i]  

                         
# add up all MA90 coefficients
for i in range(1,len(coefs_MA90)):

    muts_involved = [int(x)-1 for x in names_MA90[i].split(',')]
    
    # only consider 3rd order and higher
    if len(muts_involved) >= 2:
        # only consider significant terms
        if sig_MA90[i]:
            for j in range(len(muts_involved)):
                epistasis_sum_MA90[muts_involved[j]] += coefs_MA90[i]   

for i in range(len(epistasis_first_MA90)):
    epistasis_first_MA90[i] = coefs_MA90[i+1]
    epistasis_first_G189E[i] = coefs_G189E[i+1]

for i in range(10):
    epistasis_first_SI06[i] = coefs_SI06[i+1]

for i in range(11,16):
    epistasis_first_SI06[i] = coefs_SI06[i+1]

epistasis_first_SI06[10] = 3
epistasis_sum_SI06[10] = 3

ep_first_order = np.concatenate((epistasis_first_MA90,epistasis_first_SI06,epistasis_first_G189E),axis=1) 
ep_higher_order = np.concatenate((epistasis_sum_MA90,epistasis_sum_SI06,epistasis_sum_G189E),axis=1)  
    
ep_first_order =  np.flipud(ep_first_order) 
ep_higher_order =  np.flipud(ep_higher_order)
full_mut_names =np.flipud(full_mut_names)

M = 3
N = 16
x = np.arange(M + 1)
y = np.arange(N + 1)
xs, ys = np.meshgrid(x, y)
zs = (xs * ys) % 10
zs = zs[:-1, :-1].ravel()

fig, ax = plt.subplots(figsize=(2*.9,4*.9))
triangles1 = [(i + j*(M+1), i+1 + j*(M+1), i + (j+1)*(M+1)) for j in range(N) for i in range(M)]
triangles2 = [(i+1 + j*(M+1), i+1 + (j+1)*(M+1), i + (j+1)*(M+1)) for j in range(N) for i in range(M)]
triang1 = Triangulation(xs.ravel(), ys.ravel(), triangles1)
triang2 = Triangulation(xs.ravel(), ys.ravel(), triangles2)
img1 = plt.tripcolor(triang1, ep_higher_order.ravel(), cmap='bwr', vmin=-2.5,vmax=2.5,linewidths=0)
img2 = plt.tripcolor(triang2, ep_first_order.ravel(), cmap='bwr', vmin=-2.5,vmax=2.5,linewidths=0)


ax.text(0, 16.4, 'MA90', size=6, color='black',fontweight='bold')
ax.text(1.2, 16.4, 'SI06', size=6, color='black',fontweight='bold')
ax.text(2, 16.4, 'G189E', size=6, color='black',fontweight='bold')
for _, spine in ax.spines.items():
    spine.set_visible(False)

plt.colorbar(img2, pad=0.05)
plt.yticks(np.arange(0.6,num_mutations_G189E+0.4,1),full_mut_names,rotation='0')
plt.xticks(color='w')
plt.tick_params(length=0,pad=2)

plt.tight_layout()
plt.savefig('CH65_biochem_combo.eps')
plt.show()

HA_cont = []
CD_cont = []
with open('structural_info.csv','r') as readfile:
    kd_reader = csv.reader(readfile)
    header = next(kd_reader)
    for row in kd_reader:

        HA_contact = row[1]
        CD_contact = row[2]
        
        if len(HA_contact) != 0:  
            HA_cont.append(float(HA_contact))
            CD_cont.append(float(CD_contact))
    readfile.close()
    
HA_cont = np.array(HA_cont)  
CD_cont = np.array(CD_cont) 

HA_cont = HA_cont.reshape(16,1)
CD_cont = CD_cont.reshape(16,1)


contact = np.concatenate((HA_cont,CD_cont),axis=1) 
full_mut_names =np.flipud(full_mut_names)

sns.set_style({"axes.facecolor": "k"})       
fig, ax = plt.subplots(figsize=(1.5*0.85,2.5*.9))
sns.heatmap(contact,cmap='inferno',cbar_kws={"pad": 0.1}, linewidths=0.0) #,vmin=0.0,vmax=np.nanmax(total_epistasis)) 
plt.yticks(np.arange(0.6,num_mutations_G189E+0.4,1),full_mut_names,rotation='0',size=5)
plt.tick_params(length=0,pad=2)
plt.tick_params(length=0,pad=2)

ax.text(0, -0.6, '     HA ', size=4.5, color='black',fontweight='bold')
ax.text(0, -0.15, ' contact  contact', size=4.5, color='black',fontweight='bold')
ax.text(1.2, -0.6, 'CDR3  ', size=4.5, color='black',fontweight='bold')
for _, spine in ax.spines.items():
    spine.set_visible(False)

plt.xticks(color='w')

plt.tight_layout()
plt.savefig('mut_contact.eps')
plt.show()

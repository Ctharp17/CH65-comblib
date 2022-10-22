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
ep_type_long = 'biochemical'
antigen = 'G189E'



if antigen == 'MA90':
    coefs = np.zeros(total_params_MA90+1)
    names = []
    sig = np.full((total_params_MA90+1),0)
    stderr = np.zeros(total_params_MA90+1)
    cis = np.zeros((total_params_MA90+1,2))
    
    with open('../../Epistasis_Inference/MA90/'+ep_type_long+'/CH65_MA90_102022_'+str(order_MA90)+'order_'+ep_type+'.txt','r') as readfile:
        coef_reader = csv.reader(readfile,delimiter='\t')
        num_params = int(next(coef_reader)[-1])
        r2_train = float(next(coef_reader)[-1])
        header = next(coef_reader)
        print(header)
        for i in range(total_params_MA90+1):
            row = next(coef_reader)
            names.append(row[0])
            coefs[i] = float(row[1])  
            if i >= 1:
                stderr[i] = float(row[2])
                cis[i,0] = float(row[4])
                cis[i,1] = float(row[5])
                if float(row[4])*float(row[5]) > 0:
                    sig[i] = 1
        readfile.close()


if antigen == 'SI06': 
    coefs = np.zeros(total_params_SI06+1)
    names = []
    sig = np.full((total_params_SI06+1),0)
    stderr = np.zeros(total_params_SI06+1)
    cis = np.zeros((total_params_SI06+1,2))
    
    with open('../../Epistasis_Inference/SI06/'+ep_type_long+'/CH65_SI06_102022_'+str(order_SI06)+'order_'+ep_type+'.txt','r') as readfile:
        coef_reader = csv.reader(readfile,delimiter='\t')
        num_params = int(next(coef_reader)[-1])
        r2_train = float(next(coef_reader)[-1])
        header = next(coef_reader)
        print(header)
        for i in range(total_params_SI06+1):
            row = next(coef_reader)
            names.append(row[0])
            coefs[i] = float(row[1])  
            if i >= 1:
                stderr[i] = float(row[2])
                cis[i,0] = float(row[4])
                cis[i,1] = float(row[5])
                if float(row[4])*float(row[5]) > 0:
                    sig[i] = 1
        readfile.close()


if antigen == 'G189E':    
    # read model coefficients
    coefs = np.zeros(total_params_G189E+1)
    names = []
    sig = np.full((total_params_G189E+1),0)
    stderr = np.zeros(total_params_G189E+1)
    cis = np.zeros((total_params_G189E+1,2))
    
    with open('../../Epistasis_Inference/G189E/'+ep_type_long+'/CH65_G189E_102022_'+str(order_G189E)+'order_'+ep_type+'.txt','r') as readfile:
        coef_reader = csv.reader(readfile,delimiter='\t')
        num_params = int(next(coef_reader)[-1])
        r2_train = float(next(coef_reader)[-1])
        header = next(coef_reader)
        print(header)
        for i in range(total_params_G189E+1):
            row = next(coef_reader)
            names.append(row[0])
            coefs[i] = float(row[1])  
            if i >= 1:
                stderr[i] = float(row[2])
                cis[i,0] = float(row[4])
                cis[i,1] = float(row[5])
                if float(row[4])*float(row[5]) > 0:
                    sig[i] = 1
        readfile.close()
 
        
# initialize matrices to store values

total_epistasis = np.zeros((16,16),dtype=float)
epistasis_sum = np.zeros((16,1),dtype=float)
epistasis_higher = np.zeros((16,1),dtype=float)
epistasis_1 = np.zeros((16,1),dtype=float)

for i in range(16):  
    total_epistasis[i,i] = np.nan

    
 
# add up all MA90 coefficients
for i in range(1,len(coefs)):

    muts_involved = [int(x)-1 for x in names[i].split(',')]
       
    if len(muts_involved) == 1:
        if sig[i]:
            for j in range(len(muts_involved)):
                epistasis_1[muts_involved[j]] = coefs[i]  
                 
    # only consider 2nd 
    if len(muts_involved) == 2:
        # only consider significant terms
        if sig[i]:
            for j in range(len(muts_involved)):
                epistasis_sum[muts_involved[j]] += coefs[i]  
                for k in range(j+1,len(muts_involved)):
                     total_epistasis[muts_involved[k],muts_involved[j]] += coefs[i]
                     total_epistasis[muts_involved[j],muts_involved[k]] += coefs[i]   
                                        
    if len(muts_involved) > 2:
        if sig[i]:
            for j in range(len(muts_involved)):
                epistasis_higher[muts_involved[j]] += coefs[i]  

    

#finding min and max value to set scale bars
if np.abs(np.nanmin(epistasis_1)) >= np.nanmax(epistasis_1):
    val_1 = np.abs(np.nanmin(epistasis_1)) + 0.2
else:
    val_1 = np.nanmax(epistasis_1) + 0.2

if np.abs(np.nanmin(total_epistasis)) >= np.nanmax(total_epistasis):
    val_2 = np.abs(np.nanmin(total_epistasis)) + 0.2
else:
    val_2 = np.nanmax(total_epistasis) + 0.2

if np.abs(np.nanmin(epistasis_higher)) >= np.nanmax(epistasis_higher):
    val_h = np.abs(np.nanmin(epistasis_higher)) + 0.2
else:
    val_h = np.nanmax(epistasis_higher) + 0.2

if antigen == 'SI06':
    total_epistasis[10,0:10] = np.nanmax(total_epistasis) + 1  
    total_epistasis[11:16,10] = np.nanmax(total_epistasis) + 1 
    epistasis_1[10] = np.nanmax(epistasis_1) + 1  
    epistasis_higher[10] = np.nanmax(epistasis_higher) + 1  

    
sns.set_style({"axes.facecolor": "k"})       
fig, ax = plt.subplots(figsize=(3.0*.9,2.4*.9))

sns.heatmap(total_epistasis,cmap='bwr',xticklabels=full_mut_names,yticklabels=full_mut_names,vmin=-val_2,vmax=val_2,linewidths=0,linecolor='w')
plt.tick_params(length=0,pad=3)
for _, spine in ax.spines.items():
    spine.set_visible(True)
line1 = lines.Line2D([0,num_mutations_MA90], [0,num_mutations_MA90],lw=0.5, color='white', axes=ax)
line2 = lines.Line2D([0,6],[6],lw=1, color='gray', axes=ax)
line3 = lines.Line2D([6],[0,6],lw=1, color='gray', axes=ax)
line4 = lines.Line2D([6],[6,16],lw=1, color='tan', axes=ax)
line5 = lines.Line2D([6,16],[6],lw=1, color='tan', axes=ax)

ax.xaxis.labelpad = 20
ax.add_line(line1)
ax.add_line(line2)
ax.add_line(line3)
ax.add_line(line4)
ax.add_line(line5)
ax.get_xticklabels()[-6].set_weight("bold")
ax.get_xticklabels()[-8].set_weight("bold")
ax.get_xticklabels()[-10].set_weight("bold")
ax.get_yticklabels()[-6].set_weight("bold")
ax.get_yticklabels()[-8].set_weight("bold")
ax.get_yticklabels()[-10].set_weight("bold")
ax.text(0, -.5, '%s' %antigen, size=8, color='black',fontweight='bold')
ax.tick_params(pad=2)


plt.tight_layout()
plt.savefig('CH65_'+antigen+'_'+ep_type+'_2_sum.eps')
plt.show()


sns.set_style({"axes.facecolor": "k"})       
fig, ax = plt.subplots(figsize=(10.0*.1,2.4*.9))
sns.heatmap(epistasis_1,cmap='bwr',cbar_kws={"pad": 0.2},yticklabels=full_mut_names,vmin=-val_1,vmax=val_1,linewidths=0,linecolor='w')
for _, spine in ax.spines.items():
    spine.set_visible(True)
plt.tick_params(length=0,pad=2)
plt.xticks(color='w')
plt.tight_layout()
plt.yticks(size=5.5)
plt.savefig('CH65_'+antigen+'_'+ep_type+'_1.eps')
plt.show()


sns.set_style({"axes.facecolor": "k"})       
fig, ax = plt.subplots(figsize=(9.0*.1,2.4*.9))
sns.heatmap(epistasis_higher,cmap='PuOr_r',cbar_kws={"pad": 0.2},yticklabels=full_mut_names,vmin=-val_h,vmax=val_h,linewidths=0,linecolor='w')
for _, spine in ax.spines.items():
    spine.set_visible(True)
plt.tick_params(length=0,pad=2)
plt.tight_layout()
plt.xticks(color='w')
plt.yticks(size=5.5)
plt.savefig('CH65_'+antigen+'_'+ep_type+'_higher.eps')
plt.show()

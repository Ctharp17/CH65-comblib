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
import matplotlib.colors as colors

from mutation_info import *

def powerset(iterable):
    "powerset([1,2,3]) --> () (1,) (2,) (3,) (1,2) (1,3) (2,3) (1,2,3)"
    s = list(iterable)
    return it.chain.from_iterable(it.combinations(s, r) for r in range(len(s)+1))

plt.rcParams.update({'font.size': 7})
plt.rcParams['font.sans-serif'] = "Arial"
plt.rcParams['font.family'] = "sans-serif"
plt.rcParams['axes.linewidth'] = 0.5
plt.rcParams["xtick.major.size"] = 2
plt.rcParams["ytick.major.size"] = 2

def truncate_colormap(cmap, minval=0.0, maxval=1.0, n=100):
    new_cmap = colors.LinearSegmentedColormap.from_list(
        'trunc({n},{a:.2f},{b:.2f})'.format(n=cmap.name, a=minval, b=maxval),
        cmap(np.linspace(minval, maxval, n)))
    return new_cmap

# set some things
num_mutations_MA90 = 16
order_MA90 = 5
num_term_list_MA90 = np.array([int(comb(num_mutations_MA90,i)) for i in range(1,order_MA90+1)])
total_params_MA90 = sum(num_term_list_MA90)
order_start_indices_MA90 = list(np.cumsum(num_term_list_MA90)+1)
order_start_indices_MA90.insert(0,1)
print(num_term_list_MA90,total_params_MA90)
print(order_start_indices_MA90)

num_mutations_G189E = 16
order_G189E = 6
num_term_list_G189E = np.array([int(comb(num_mutations_G189E,i)) for i in range(1,order_G189E+1)])
total_params_G189E = sum(num_term_list_G189E)
order_start_indices_G189E = list(np.cumsum(num_term_list_G189E)+1)
order_start_indices_G189E.insert(0,1)
print(num_term_list_G189E,total_params_G189E)
print(order_start_indices_G189E)

num_mutations_SI06 = 15
order_SI06 = 5
num_term_list_SI06 = np.array([int(comb(num_mutations_SI06,i)) for i in range(1,order_SI06+1)])
total_params_SI06 = sum(num_term_list_SI06)
order_start_indices_SI06 = list(np.cumsum(num_term_list_SI06)+1)
order_start_indices_SI06.insert(0,1)
print(num_term_list_SI06,total_params_SI06)
print(order_start_indices_SI06)



mut_names = ['30','35','36','57','64','65','66','79','82','83','84','85','92','95','103','113']

MA90_color = '#e8735c'
SI06_color = '#72c2a6'
G189E_color = '#5482a7'



ep_type = 'stat'




coefs_MA90 = np.zeros(total_params_MA90+1)
stderr_MA90 = np.zeros(total_params_MA90+1)
names_MA90 = []
sig_MA90 = np.full((total_params_MA90+1),0)
cis_MA90 = np.zeros((total_params_MA90+1,2))

with open('CH65_MA90_102022_'+str(order_MA90)+'order_'+ep_type+'.txt','r') as readfile:
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
    


coefs_SI06 = np.zeros(total_params_SI06+1)
stderr_SI06 = np.zeros(total_params_SI06+1)
names_SI06 = []
sig_SI06 = np.full((total_params_SI06+1),0)
cis_SI06 = np.zeros((total_params_SI06+1,2))

with open('CH65_SI06_102022_'+str(order_SI06)+'order_'+ep_type+'.txt','r') as readfile:
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



coefs_G189E = np.zeros(total_params_G189E+1)
stderr_G189E = np.zeros(total_params_G189E+1)
names_G189E = []
sig_G189E = np.full((total_params_G189E+1),0)
cis_G189E = np.zeros((total_params_G189E+1,2))

with open('CH65_G189E_102022_'+str(order_G189E)+'order_'+ep_type+'.txt','r') as readfile:
    coef_reader = csv.reader(readfile,delimiter='\t')
    num_params = int(next(coef_reader)[-1])
    r2_train = float(next(coef_reader)[-1])
    print(r2_train)
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


# read in genos and phenos
df = pd.read_csv('20221008_CH65_QCfilt_REPfilt.csv',dtype={"variant": str})
#df = pd.read_csv('20210427_HA_unadj_fil_merg.csv',dtype={"variant": str})

# G189E
df_G189E = df.dropna(subset=['G189E_mean'])
genos_G189E = np.array(df_G189E[['pos'+x for x in G189E_mutations]].copy(),dtype=np.float64)
phenos_G189E = df_G189E[['G189E_mean']].values.flatten()
print(genos_G189E.shape,phenos_G189E.shape)
genos_G189E = 2*genos_G189E-1


poly_current_G189E = PolynomialFeatures(order_G189E,interaction_only=True)
genos_current_G189E = poly_current_G189E.fit_transform(genos_G189E)

phenos_pred = np.tensordot(genos_current_G189E,coefs_G189E,axes=1)
pearsonr = np.corrcoef(phenos_pred,phenos_G189E)[1,0]
print(pearsonr,pearsonr**2)


order_r2_G189E = []
for i in range(1,order_G189E+1):
    phenos_pred = np.tensordot(genos_current_G189E[:,:order_start_indices_G189E[i]],coefs_G189E[:order_start_indices_G189E[i]],axes=1)
    pearsonr = np.corrcoef(phenos_pred,phenos_G189E)[1,0]
    order_r2_G189E.append(pearsonr**2)
print('order r2:', order_r2_G189E)

var_expl = np.array(order_r2_G189E)/order_r2_G189E[-1]

print('order -1:', order_r2_G189E[-1])

print('var expl:', var_expl)

total_ep = 1-var_expl[0]

print('total ep:', total_ep)


delta_ep = list(var_expl)
delta_ep.insert(0,0)
print('delta ep:', delta_ep)


delta_ep = [delta_ep[i]-delta_ep[i-1] for i in range(1,6)]
print(delta_ep)

plt.figure(figsize=(1.6,2))
plt.bar(range(1,6),delta_ep,width=0.8,color=G189E_color,label='G189E')
ax = plt.gca()
ax.set_facecolor('white')
plt.ylim([0,1])
plt.ylabel('Fractional variance explained')
plt.xlabel('Order of interaction')
plt.xticks([1,2,3,4,5],['1','2','3','4','5'])
plt.legend(facecolor="white")
plt.tight_layout()
plt.savefig('CH65_G189E_var_partitioning_102022.eps')
plt.show()


# read in genos and phenos
df = pd.read_csv('20221008_CH65_QCfilt_REPfilt.csv',dtype={"variant": str})
# H3
# # for H3, filter for the three required mutations and remove them
df_SI06 = df.dropna(subset=['SI06_mean'])
for mut in SI_required_mutations:
    df_SI06 = df_SI06.loc[df_SI06['pos'+mut] == 1]

genos_SI06 = np.array(df_SI06[['pos'+x for x in SI_mutations]].copy(),dtype=np.float64)
phenos_SI06 = df_SI06[['SI06_mean']].values.flatten()
print(genos_SI06.shape,phenos_SI06.shape)
genos_SI06 = 2*genos_SI06-1


poly_current_SI06 = PolynomialFeatures(order_SI06,interaction_only=True)
genos_current_SI06 = poly_current_SI06.fit_transform(genos_SI06)

phenos_pred = np.tensordot(genos_current_SI06,coefs_SI06,axes=1)
pearsonr = np.corrcoef(phenos_pred,phenos_SI06)[1,0]
print(pearsonr,pearsonr**2)


order_r2_SI06 = []
for i in range(1,order_SI06+1):
    phenos_pred = np.tensordot(genos_current_SI06[:,:order_start_indices_SI06[i]],coefs_SI06[:order_start_indices_SI06[i]],axes=1)
    pearsonr = np.corrcoef(phenos_pred,phenos_SI06)[1,0]
    order_r2_SI06.append(pearsonr**2)
print(order_r2_SI06)
var_expl = np.array(order_r2_SI06)/order_r2_SI06[-1]

total_ep = 1-var_expl[0]
print(total_ep)

delta_ep = list(var_expl)
delta_ep.insert(0,0)
print(delta_ep)

delta_ep = [delta_ep[i]-delta_ep[i-1] for i in range(1,order_SI06+1)]
print(delta_ep)



plt.figure(figsize=(1.4,2))
plt.bar(range(1,order_SI06+1),delta_ep,width=0.8,color=SI06_color,label='SI06')
ax = plt.gca()
ax.set_facecolor('white')
plt.ylim([0,1])
plt.ylabel('Fractional variance explained')
plt.xlabel('Order of interaction')
plt.xticks([1,2,3,4,5],['1','2','3','4','5'])
plt.legend(facecolor="white")
plt.tight_layout()
plt.savefig('CH65_SI06_var_partitioning_102022.eps')
plt.show()



# read in genos and phenos
df = pd.read_csv('20221008_CH65_QCfilt_REPfilt.csv',dtype={"variant": str})
#df = pd.read_csv('20210427_HA_unadj_fil_merg.csv',dtype={"variant": str})

# MA90
df_MA90 = df.dropna(subset=['MA90_mean'])
genos_MA90 = np.array(df_MA90[['pos'+x for x in MA_mutations]].copy(),dtype=np.float64)
phenos_MA90 = df_MA90[['MA90_mean']].values.flatten()
print(genos_MA90.shape,phenos_MA90.shape)
genos_MA90 = 2*genos_MA90-1


poly_current_MA90 = PolynomialFeatures(order_MA90,interaction_only=True)
genos_current_MA90 = poly_current_MA90.fit_transform(genos_MA90)

phenos_pred = np.tensordot(genos_current_MA90,coefs_MA90,axes=1)
pearsonr = np.corrcoef(phenos_pred,phenos_MA90)[1,0]
print(pearsonr,pearsonr**2)


order_r2_MA90 = []
for i in range(1,order_MA90+1):
    phenos_pred = np.tensordot(genos_current_MA90[:,:order_start_indices_MA90[i]],coefs_MA90[:order_start_indices_MA90[i]],axes=1)
    pearsonr = np.corrcoef(phenos_pred,phenos_MA90)[1,0]
    order_r2_MA90.append(pearsonr**2)
print('order r2:', order_r2_MA90)

var_expl = np.array(order_r2_MA90)/order_r2_MA90[-1]

print('order -1:', order_r2_MA90[-1])

print('var expl:', var_expl)

total_ep = 1-var_expl[0]

print('total ep:', total_ep)


delta_ep = list(var_expl)
delta_ep.insert(0,0)
print('delta ep:', delta_ep)


delta_ep = [delta_ep[i]-delta_ep[i-1] for i in range(1,6)]
print(delta_ep)

plt.figure(figsize=(1.6,2))
plt.bar(range(1,6),delta_ep,width=0.8,color=MA90_color,label='MA90')
ax = plt.gca()
ax.set_facecolor('white')
plt.ylim([0,1])
plt.ylabel('Fractional variance explained')
plt.xlabel('Order of interaction')
plt.xticks([1,2,3,4,5],['1','2','3','4','5'])
plt.legend(facecolor="white")
plt.tight_layout()
plt.savefig('CH65_MA90_var_partitioning_102022.eps')
plt.show()
import numpy as np
import csv
# import sys
from sklearn.preprocessing import PolynomialFeatures
import statsmodels.api as sm

# choose statistical or biochemical epistasis
ep_type = 'biochem' 
# ep_type = 'stat'


# read in rsq data from CV folds        
num_folds = 10
if ep_type == 'stat':
 	num_folds = 9
rsq_train_list_G189E = []
rsq_test_list_G189E = []
for f in range(num_folds):
 	rsq_train_list_G189E.append([])
 	rsq_test_list_G189E.append([])
 	if ep_type == 'stat':
         fold = f #f+1?
 	else:
         fold = f
 	with open('/n/desai_lab/users/caelanb/G189E/temp_CH65_G189E_CV_rsquared_newdata_'+str(fold)+'_'+ep_type+'.csv','r') as readfile:
         rsq_reader = csv.reader(readfile)
         header = next(rsq_reader)
         for row in rsq_reader:
             rsq_train_list_G189E[f].append(float(row[1]))
             rsq_test_list_G189E[f].append(float(row[2]))
        
rsq_train_list_G189E = np.array(rsq_train_list_G189E)
rsq_test_list_G189E = np.array(rsq_test_list_G189E)

# average over folds to get optimal order
mean_rsq_train_G189E = np.mean(rsq_train_list_G189E,axis=0)
stdev_rsq_train_G189E = np.std(rsq_train_list_G189E,axis=0)
mean_rsq_test_G189E = np.mean(rsq_test_list_G189E,axis=0)
stdev_rsq_test_G189E = np.std(rsq_test_list_G189E,axis=0)

optimal_G189E_order = np.argmax(mean_rsq_test_G189E)
# print('Optimal order, CH65: ',optimal_G189E_order,file=sys.stdout,flush=True)

# print(mean_rsq_test_G189E)

# write averages to new file
with open('/n/desai_lab/users/caelanb/G189E/CH65_CV_rsquared_G189E_newdata_'+ep_type+'.csv','w') as writefile:
# with open('CH65_CV_rsquared_G189E_'+ep_type+'.csv','w') as writefile:
 	rsq_writer = csv.writer(writefile)
 	rsq_writer.writerow(['Optimal order: '+str(optimal_G189E_order)])
 	rsq_writer.writerow(['Type','Order','Mean','Std'])
 	for i in range(len(mean_rsq_train_G189E)):
         rsq_writer.writerow(['Train',str(i),mean_rsq_train_G189E[i],stdev_rsq_train_G189E[i]])
 	for i in range(len(mean_rsq_test_G189E)):
         rsq_writer.writerow(['Test',str(i),mean_rsq_test_G189E[i],stdev_rsq_test_G189E[i]])
 	writefile.close()
 

# read in data

geno_vectors_G189E = []
phenos_G189E = []

mutations_G189E = [str(x) for x in range(1,17)]

with open('20220603_CH65_QCfilt_REPfilt.csv','r') as readfile:
    kd_reader = csv.reader(readfile)
    header = next(kd_reader)
    for row in kd_reader:
        geno = row[0]
        
        geno_vec = np.array([float(x) for x in geno])

        pheno_G189E = row[11]
        
            
        if len(pheno_G189E) != 0:  
            geno_vectors_G189E.append(geno_vec)
            phenos_G189E.append(float(pheno_G189E))
    readfile.close()

phenos_G189E = np.array(phenos_G189E)

genos_G189E = np.empty((len(phenos_G189E),len(geno_vectors_G189E[0])))
for i in range(len(phenos_G189E)):
    genos_G189E[i] = geno_vectors_G189E[i][:]
    


if ep_type == 'stat':
    genos_G189E = 2*(genos_G189E-0.5)

# print(genos_G189E.shape,phenos_G189E.shape)


# # Fit final models

np.random.seed(2112)
indices_permuted_G189E = np.random.permutation(np.arange(len(genos_G189E)))

# fit models of increasing order
for order in range(1,optimal_G189E_order+1):
# for order in range(1,2):
    # print(order)
    genos_G189E_permuted = genos_G189E[indices_permuted_G189E]
    phenos_G189E_permuted = phenos_G189E[indices_permuted_G189E]
    # print('Order: ',str(order),file=sys.stdout,flush=True)
    poly_G189E_current = PolynomialFeatures(order,interaction_only=True)
    genos_G189E_current = poly_G189E_current.fit_transform(genos_G189E_permuted)

    # fit
    reg_G189E_current = sm.OLS(phenos_G189E_permuted,genos_G189E_current).fit()
    reg_G189E_coefs_current = reg_G189E_current.params
    reg_G189E_CIs_current = reg_G189E_current.conf_int(alpha=0.05/float(len(reg_G189E_coefs_current)), cols=None)
    reg_G189E_stderr = reg_G189E_current.bse
    reg_G189E_pvalues = reg_G189E_current.pvalues
    
    num_sig = len(np.where(reg_G189E_pvalues < 0.05/float(len(reg_G189E_coefs_current)))[0])

    predicted_phenos_permuted_G189E = reg_G189E_current.predict(genos_G189E_current)
    rsquared_G189E_current = reg_G189E_current.rsquared
    # print('Params: ',len(reg_G189E_coefs_current),file=sys.stdout,flush=True)
    # print('Performance: ',rsquared_G189E_current,file=sys.stdout,flush=True)
    # print(num_sig,file=sys.stdout,flush=True)
	 

    # write model to file
    coef_names = poly_G189E_current.get_feature_names(input_features = mutations_G189E)
    with open('/n/desai_lab/users/caelanb/G189E/CH65_G189E_newdata_'+str(order)+'order_'+ep_type+'.txt','w') as writefile:
        coef_writer = csv.writer(writefile,delimiter='\t')
        coef_writer.writerow(['Params: ',len(reg_G189E_coefs_current)])
        coef_writer.writerow(['Performance: ',rsquared_G189E_current])
        coef_writer.writerow(['Term','Coefficient','Standard Error','p-value','95% CI lower','95% CI upper'])
        coef_writer.writerow(['Intercept',reg_G189E_coefs_current[0]])
        for i in range(1,len(reg_G189E_coefs_current)):
            coef_writer.writerow([','.join(coef_names[i].split(' ')),reg_G189E_coefs_current[i],reg_G189E_stderr[i],
                                  reg_G189E_pvalues[i],reg_G189E_CIs_current[i][0],reg_G189E_CIs_current[i][1]])
        writefile.close()



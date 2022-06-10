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
rsq_train_list_SI06 = []
rsq_test_list_SI06 = []
for f in range(num_folds):
 	rsq_train_list_SI06.append([])
 	rsq_test_list_SI06.append([])
 	if ep_type == 'stat':
         fold = f #f+1?
 	else:
         fold = f
 	with open('/n/desai_lab/users/caelanb/SI06/temp_CH65_SI06_fil_CV_rsquared_newdata_'+str(fold)+'_'+ep_type+'.csv','r') as readfile:
         rsq_reader = csv.reader(readfile)
         header = next(rsq_reader)
         for row in rsq_reader:
             rsq_train_list_SI06[f].append(float(row[1]))
             rsq_test_list_SI06[f].append(float(row[2]))
        
rsq_train_list_SI06 = np.array(rsq_train_list_SI06)
rsq_test_list_SI06 = np.array(rsq_test_list_SI06)

# average over folds to get optimal order
mean_rsq_train_SI06 = np.mean(rsq_train_list_SI06,axis=0)
stdev_rsq_train_SI06 = np.std(rsq_train_list_SI06,axis=0)
mean_rsq_test_SI06 = np.mean(rsq_test_list_SI06,axis=0)
stdev_rsq_test_SI06 = np.std(rsq_test_list_SI06,axis=0)

optimal_SI06_order = np.argmax(mean_rsq_test_SI06)
# print('Optimal order, CH65: ',optimal_SI06_order,file=sys.stdout,flush=True)

# print(mean_rsq_test_SI06)

# write averages to new file
with open('/n/desai_lab/users/caelanb/SI06/CH65_CV_rsquared_SI06_newdata_'+ep_type+'.csv','w') as writefile:
#with open('CH65_CV_rsquared_SI06_fil_'+ep_type+'.csv','w') as writefile:
 	rsq_writer = csv.writer(writefile)
 	rsq_writer.writerow(['Optimal order: '+str(optimal_SI06_order)])
 	rsq_writer.writerow(['Type','Order','Mean','Std'])
 	for i in range(len(mean_rsq_train_SI06)):
         rsq_writer.writerow(['Train',str(i),mean_rsq_train_SI06[i],stdev_rsq_train_SI06[i]])
 	for i in range(len(mean_rsq_test_SI06)):
         rsq_writer.writerow(['Test',str(i),mean_rsq_test_SI06[i],stdev_rsq_test_SI06[i]])
 	writefile.close()
 

# read in data

geno_vectors_SI06 = []
phenos_SI06 = []

mutations_SI06 = ['1','2','3','4','5','6','7','8','9','10','12','13','14','15','16']

with open('20220603_CH65_QCfilt_REPfilt.csv','r') as readfile:
    kd_reader = csv.reader(readfile)
    header = next(kd_reader)
    for row in kd_reader:
        geno = row[0]
        
        geno_vec = np.array([float(x) for x in geno])

        pheno_SI06 = row[7]
        
            
        if len(pheno_SI06) != 0 and row[23] == '1':  
            geno_vectors_SI06.append(geno_vec)
            phenos_SI06.append(float(pheno_SI06))
    readfile.close()

phenos_SI06 = np.array(phenos_SI06)

genos_SI06 = np.empty((len(phenos_SI06),len(geno_vectors_SI06[0])))
for i in range(len(phenos_SI06)):
    genos_SI06[i] = geno_vectors_SI06[i][:]
    
genos_SI06 = np.delete(genos_SI06, 10, 1)


if ep_type == 'stat':
    genos_SI06 = 2*(genos_SI06-0.5)

# print(genos_SI06.shape,phenos_SI06.shape)


# # Fit final models

np.random.seed(2112)
indices_permuted_SI06 = np.random.permutation(np.arange(len(genos_SI06)))

# fit models of increasing order
for order in range(1,optimal_SI06_order+1):
# for order in range(1,2):
    # print(order)
    genos_SI06_permuted = genos_SI06[indices_permuted_SI06]
    phenos_SI06_permuted = phenos_SI06[indices_permuted_SI06]
    # print('Order: ',str(order),file=sys.stdout,flush=True)
    poly_SI06_current = PolynomialFeatures(order,interaction_only=True)
    genos_SI06_current = poly_SI06_current.fit_transform(genos_SI06_permuted)

    # fit
    reg_SI06_current = sm.OLS(phenos_SI06_permuted,genos_SI06_current).fit()
    reg_SI06_coefs_current = reg_SI06_current.params
    reg_SI06_CIs_current = reg_SI06_current.conf_int(alpha=0.05/float(len(reg_SI06_coefs_current)), cols=None)
    reg_SI06_stderr = reg_SI06_current.bse
    reg_SI06_pvalues = reg_SI06_current.pvalues
    
    num_sig = len(np.where(reg_SI06_pvalues < 0.05/float(len(reg_SI06_coefs_current)))[0])

    predicted_phenos_permuted_SI06 = reg_SI06_current.predict(genos_SI06_current)
    rsquared_SI06_current = reg_SI06_current.rsquared
    # print('Params: ',len(reg_SI06_coefs_current),file=sys.stdout,flush=True)
    # print('Performance: ',rsquared_SI06_current,file=sys.stdout,flush=True)
    # print(num_sig,file=sys.stdout,flush=True)
	 

    # write model to file
    coef_names = poly_SI06_current.get_feature_names(input_features = mutations_SI06)
    with open('/n/desai_lab/users/caelanb/SI06/CH65_SI06_newdata_'+str(order)+'order_'+ep_type+'.txt','w') as writefile:
        coef_writer = csv.writer(writefile,delimiter='\t')
        coef_writer.writerow(['Params: ',len(reg_SI06_coefs_current)])
        coef_writer.writerow(['Performance: ',rsquared_SI06_current])
        coef_writer.writerow(['Term','Coefficient','Standard Error','p-value','95% CI lower','95% CI upper'])
        coef_writer.writerow(['Intercept',reg_SI06_coefs_current[0]])
        for i in range(1,len(reg_SI06_coefs_current)):
            coef_writer.writerow([','.join(coef_names[i].split(' ')),reg_SI06_coefs_current[i],reg_SI06_stderr[i],
                                  reg_SI06_pvalues[i],reg_SI06_CIs_current[i][0],reg_SI06_CIs_current[i][1]])
        writefile.close()



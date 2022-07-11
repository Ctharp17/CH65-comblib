import numpy as np

num_muts_total = 16
mutations = np.array([str(x) for x in range(1,num_muts_total+1)])
mut_names = np.array(['N26D','S29R','Y35N','Y48C','D49Y','V98I','G31D','Y33H','M34I','H35N','N52H','G57D','L83V','S84N','R85G','R87K'])

num_muts_MA = 16
MA_indices = np.arange(num_muts_total,dtype=int)
MA_mutations = mutations
MA_mut_names = mut_names

num_muts_SI = 15
SI_required_indices = [10]
SI_required_mutations = [str(x+1) for x in SI_required_indices]
SI_remaining_indices = [0,1,2,3,4,5,6,7,8,9,11,12,13,14,15]
SI_mutations = mutations[SI_remaining_indices]
SI_mut_names = mut_names[SI_remaining_indices]

num_muts_G189E = 16
G189E_indices = np.arange(num_muts_total,dtype=int)
G189E_mutations = mutations
G189E_mut_names = mut_names

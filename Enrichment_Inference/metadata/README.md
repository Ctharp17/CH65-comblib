# Settings and sequences

TODO Write something more descriptive here


`sample_info.tsv`
-----------------
- `concentration`: Concentration of the solute, string used to label the tubes
- `concentration_float`: Measured concentration of the solute
- `bin`: bin name used to label the tubes

**notes:

change CH65a_uns fastq file to be called '*_reseq'

before running, copy directory to holyscratch:
cd ../..
scp -r desai_lab/users/aphillips/babyYODA/titeseq/CH65_full_titeseq_analysis/ holyscratch01/desai_lab/Users/aphillips/babyYODA/
then run via sbatch
cd holyscratch01/desai_lab/Users/aphillips/babyYODA/CH65_full_titeseq_analysis
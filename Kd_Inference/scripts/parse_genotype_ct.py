import regex
import pandas as pd
from Bio.Seq import Seq
import numpy as np
import gzip

geno_for_regex = regex.compile(open(snakemake.input.geno_for_regex).read())

variants_barcodes = pd.read_csv(config["variant_barcode_table"], dtype = 'str')
geno_binary = dict(zip(variants_barcodes['barcode'], variants_barcodes['binary_genotypes']))

class BarcodeDoesntMatchError(ValueError):
    pass


def mutations_to_binary_geno(match):
    barcode = match.group('barcode').upper()
    s = ''
    if variants_barcodes['barcode'].str.contains(barcode).any():
        s = geno_binary[barcode]

    else:
        raise BarcodeDoesntMatchError

    return s

unmatched_reads = 0
reads_with_err_in_barcode = 0
reads_with_no_err = 0
reads_with_indels = 0
reads_with_one_sub = 0
reads_with_two_or_more_subs = 0
total_accepted = 0
nb_reads = 0

with \
     open(snakemake.input.tsv, 'r') as f, \
     open(snakemake.output.tsv, 'w') as fw, \
     gzip.open(snakemake.output.bad_reads, 'wt') as fbad:

    f.readline() # remove header line
    fw.write("\t".join(['UMI', 'geno']) + "\n")
    fbad.write("\t".join(['read_1', 'read_2']) + "\n")
    for line in f:
        UMI, read_1 = line.strip().split("\t")

        match_1 = geno_for_regex.match(read_1)
        nb_reads += 1

        if match_1 is None:
            unmatched_reads += 1
            fbad.write(f"{read_1}\n")
            continue

        try:
            geno = mutations_to_binary_geno(match_1)
        except BarcodeDoesntMatchError:
            reads_with_err_in_barcode += 1
            # TODO we can fix the error if it is closer to one codon than the other
            continue

        sub1, in1, del1 = match_1.fuzzy_counts
        num_indels = in1 + del1
        num_subs = sub1

        if num_indels != 0:
            reads_with_indels += 1
        elif num_subs == 1:
            reads_with_one_sub += 1
        elif num_subs >= 2:
            reads_with_two_or_more_subs += 1
        else:
            reads_with_no_err += 1

        total_accepted += 1
        fw.write("\t".join((UMI, geno)) + "\n")



# Write stats
(pd.DataFrame({
    'sample' : [snakemake.wildcards.sample],
    'unmatched_reads' : [unmatched_reads],
    'reads_with_err_in_barcode' : [reads_with_err_in_barcode],
    'reads_with_indels' : [reads_with_indels],
    'reads_with_one_sub' : [reads_with_one_sub],
    'reads_with_two_or_more_sub' : [reads_with_two_or_more_subs],
    'reads_with_no_err' : [reads_with_no_err],
    'total_accepted' : [total_accepted],
    'fraction_of_geno_accepted' : [(total_accepted) / nb_reads] if nb_reads != 0 else np.nan,
})
 .set_index('sample')
 .to_csv(snakemake.output.stats, sep='\t')
)

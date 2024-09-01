import math
from Bio.Seq import Seq
import pandas as pd
import regex
import numpy as np


def rc(x):
    return str(Seq(x).reverse_complement())

UMI_length = snakemake.config['umi_length']
col_inline_idx = snakemake.config['col_inline_indices']
row_inline_idx = snakemake.config['row_inline_indices']

num_idx_subs = snakemake.config['num_substitutions_tolerated_in_inline_idx']
num_substitutions_tolerated_per_ten_bp = snakemake.config['num_substitutions_tolerated_per_ten_bp']

def make_or_regex(seqs, name, dist=0):
    '''Returns a regex string that is an OR of all the seqs in seqs.
    The capture is given a name `name`.
    Substitutions only are accepted within edit distance `dist`.
    '''
    ret = f"(?P<{name}>{'|'.join(seqs)})"
    if dist != 0:
        ret += '{s<=%d}' % dist
    return ret


UMI_regex = "(?P<UMI>[ACGT]{%d})" % UMI_length

cols_regex = make_or_regex(col_inline_idx, name='index', dist=num_idx_subs)
rows_regex = make_or_regex(row_inline_idx, name='index', dist=num_idx_subs)
read_seq_regex = "(?P<read_seq>[ACGTN]+)"

read_1_regex = cols_regex + read_seq_regex
read_2_regex = UMI_regex + rows_regex + read_seq_regex

with open(snakemake.output.read_1_regex, 'w') as f:
    f.write(read_1_regex)

with open(snakemake.output.read_2_regex, 'w') as f:
    f.write(read_2_regex)


for_read_constant_region_1 = config['const_reg_1']

for_read_constant_region_2 = config['const_reg_2']
#following line is not needed for barcode-variant table approach and is thus commented out
#rev_read_constant_regions = list(map(rc, reversed(constant_regions)))

barcode_regex = "(?P<barcode>[ACGT]{3}AT[ACGT]{5}AT[ACGT]{5}AT[ACGT]{3})" 

print(for_read_constant_region_1)
print(for_read_constant_region_2)

def make_fuzzy_regex(search, num_err, name='', bestmatch=True, allow_indels=False):
    if search == '':
        return ''
    r = ''
    if bestmatch and num_err != 0:
        r += '(?e)'
    r += '('
    if name != '':
        r += f'?P<{name}>'
    r += search
    r += ')'
    if num_err != 0:
        r += '{%c<=%d}' % ('e' if allow_indels else 's', num_err)
    return r


def assemble_regex(const_region, num):
    regex_parts = []
    num_mut_tolerated = math.ceil(
        (len(const_region) / 10) * num_substitutions_tolerated_per_ten_bp)
    is_indel_tolerated = are_indels_tolerated_in_primers

    regex_parts.append(make_fuzzy_regex(const_region,
                                        num_err=num_mut_tolerated,
                                        allow_indels=is_indel_tolerated,
                                        name=f'const_{num}'))
    print('\n'.join(regex_parts))
    return ''.join(regex_parts)


geno_for_regex = assemble_regex(for_read_constant_region_1, 1) + barcode_regex + assemble_regex(for_read_constant_region_2, 2)

#following line is also not needed for barcode-variant table approach
#geno_rev_regex = assemble_regex(rev_read_constant_regions)

with open(snakemake.output.geno_for_regex, 'w') as f:
    f.write(geno_for_regex)

#following line is also not needed for barcode-variant table approach
#with open(snakemake.output.geno_rev_regex, 'w') as f:
    #f.write(geno_rev_regex)

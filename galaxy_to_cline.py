if False:"""
projfile --altprojfile NONE --genbank NONE --dump_alignments NONE --dump_records NONE 
--threads NONE --fasta*** NONE --fasta_chartype NONE --feature*** NONE --featureis_fasname FALSE --qualis_fasname NONE
 --qualis_seqid NONE --qualis_seqdesc NONE --expng_seqandmeta FALSE --write_loci_csv NONE --read_loci_csv NONE
    --add_locus*** NONE --banish_loci FALSE --write_meta_csv NONE --read_meta_csv NONE --start FALSE --overwrite FALSE
    --verbose_report NONE --version_ctrl FALSE
    
*** get a string "a%%B%%c", turn to a B c
"""
import sys

galaxy_cline = sys.argv[2:]
cline = "python base_reprophylo.py %s"%sys.argv[1]
pairs = zip(galaxy_cline[::2], galaxy_cline[1::2])
for pair in pairs:
    if pair[0] in ['--featureis_fasname', '--expng_seqandmeta', '--banish_loci',
                   '--start', '--overwrite', '--version_ctrl']:
        if pair[1] == 'TRUE':
            cline += ' %s'%pair[0]
    elif pair[0] in ['--feature','--add_locus', '--fasta']:
        if not pair[1] == 'NONE':
            cline += ' %s'%pair[0]
            for i in pair[1].split('%%'):
                cline += ' %s'%i
    elif not pair[1] == 'NONE':
        cline += ' %s %s'%pair

from subprocess import call
print 'ReproPhylo has been callsed as follows:'
print cline
output = call(cline, shell=True) 
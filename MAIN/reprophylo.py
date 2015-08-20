reprophyloversion=1.0
############################################################################################
if False:
    """
    ReproPhylo version 1 
    
    General purpose phylogenetics package for reproducible and experimental analysis
    
    Amir Szitenebrg
    A.Szitenberg@Hull.ac.uk
    Szitenberg@gmail.com
    
    David H Lunt
    D.H.Lunt@Hull.ac.uk
    
    EvoHull.org
    University of Hull
    
    
    Developed with:
    CPython 2.7.6
    IPython 1.2.1
    ete2 2.2rev1056
    biopython 1.64
    dendropy 3.12.0
    cloud 2.8.5
    numpy 1.8.2
    matplotlib 1.3.1
    pandas
    
    RAxML 8
    Phylobayes 3
    Trimal 1
    Muscle 
    Mafft 7
    Pal2nal 14
    """
##############################################################################################



from Bio import SeqIO
import os, csv, sys, dendropy, re, time, random, glob, platform, warnings, rpgit, ast, gb_syn,css
import HTML, inspect, shutil
import subprocess as sub
from Bio.Seq import Seq
import numpy as np
import matplotlib.pyplot as plt
from Bio.Alphabet import IUPAC
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation, CompoundLocation
from Bio.Align.Applications import MafftCommandline, MuscleCommandline
from StringIO import StringIO 
from Bio import AlignIO 
from Bio.Phylo.Applications import RaxmlCommandline
from Bio.Align import MultipleSeqAlignment
from Bio.SeqUtils import GC
from ete2 import *
from collections import Counter
#import pandas as pd
import math
import __builtin__




##############################################################################################
class Locus:
##############################################################################################

    """ Configure the loci stored in the ReproPhylo Project.
        
    >>> locus = Locus('dna', 'CDS', 'coi', ['cox1','COX1','coi','COI','CoI'])
    >>> print(locus)
    Locus(char_type=dna, feature_type=CDS, name=coi, aliases=cox1; COX1; coi; COI; CoI)
    """

    char_type = 'NotSet'
    feature_type = 'NotSet'
    name = 'NotSet'
    aliases = []

    def __init__(self, char_type=char_type, feature_type=feature_type,
                 name=name, aliases=aliases):

        
        self.char_type = char_type
        self.feature_type = feature_type
        self.name = name
        self.aliases = aliases
        
        valid = ['dna','prot']
        if not self.char_type in valid:
            raise ValueError('self.char_type should be \'dna\' or \'prot\'')
        if not type(self.feature_type) is str:
            raise ValueError('self.feature_type should be a string')
        if not type(self.name) is str:
            raise ValueError('self.name should be a string')
        if not type(self.aliases) is list:
            raise ValueError('self.aliases should be a list')
        else:
            for a in self.aliases:
                if not type(a) is str:
                    raise ValueError('aliases in self.aliases have to be strings')
            


    def __str__(self):
        aliases_str = ('; ').join(self.aliases)
        return ('Locus(char_type='+self.char_type+', feature_type='+self.feature_type+
                ', name='+self.name+', aliases='+aliases_str+')')

    
    
##############################################################################################
class Concatenation:
##############################################################################################

    """This class is used to configure concatenations given loci and rules.
    
    >>> coi = Locus('dna', 'CDS', 'coi', ['cox1','COX1','coi','COI','CoI'])
    >>> ssu = Locus('dna', 'rRNA', '18S', ['18S rRNA','SSU rRNA'])
    >>> bssu = Locus('dna', 'rRNA', '16S', ['16S rRNA'])
    >>> lsu = Locus('dna', 'rRNA', '28S', ['28S rRNA', 'LSU rRNA'])
    >>> alg11 = Locus('dna', 'CDS', 'ALG11', ['ALG11'])
    >>> loci = [coi, ssu, bssu, lsu, alg11]
    >>> concatenation = Concatenation(name='combined', loci=loci,
    ...                                otu_meta='OTU_name',
    ...                                otu_must_have_all_of=['coi'],
    ...                                otu_must_have_one_of =[['16S','28S'],['ALG11','18S']],
    ...                                define_trimmed_alns=["MuscleDefaults@dummyTrimMethod"])
    >>> print(str(concatenation))
    Concatenation named combined, with loci coi,18S,16S,28S,ALG11,
    of which coi must exist for all species
    and at least one of each group of [ 16S 28S ][ ALG11 18S ] is represented.
    Alignments with the following names: MuscleDefaults@dummyTrimMethod are prefered
    """
    
    otu_must_have_all_of = []
    otu_must_have_one_of = 'any'
    define_trimmed_alns = [] #should be Locus_name@Alignment_method_name@Trimming_mathod_name
    
    feature_id_dict = {}
    
    def __init__(self,
                 name,
                 loci,
                 otu_meta,
                 otu_must_have_all_of = otu_must_have_all_of,
                 otu_must_have_one_of = otu_must_have_one_of,
                 define_trimmed_alns = define_trimmed_alns):
        self.name = name
        self.loci = loci
        self.otu_meta = otu_meta
        self.otu_must_have_all_of = otu_must_have_all_of
        self.otu_must_have_one_of = otu_must_have_one_of
        if isinstance(otu_must_have_all_of,str):
            raise IOError('The keyword \'otu_must_have_all_of\' has to be a list')
        if isinstance(otu_must_have_one_of[0],str) and not otu_must_have_one_of == 'any':
            raise IOError('The keyword \'otu_must_have_one_of\' has to be a list of lists')
        if self.otu_must_have_one_of == 'any':
            self.otu_must_have_one_of = [[l.name for l in self.loci]]
        self.feature_id_dict = {} # Will hold the feature_id list for each otu
        self.define_trimmed_alns = define_trimmed_alns # To choose between alternative
        # alignments of the same locus
        self.used_trimmed_alns = {} #To hold the alignment chosen for each locus
        
        # Validate loci list
        seen = []
        for locus in loci:
            if not isinstance(locus, Locus):
                raise TypeError("Expecting Locus object in loci list")
            if locus.name in seen:
                raise NameError('Locus ' + locus.name + ' appears more than once in self.loci')
            else:
                seen.append(locus.name)
      
                
                
    def __str__(self):
        loci_names = [i.name for i in self.loci]
        loci_string = ''
        for l in loci_names:
            loci_string += l+','
        loci_string = loci_string[:-1]
        must_have = ''
        for i in self.otu_must_have_all_of:
            must_have += i+','
        must_have = must_have[:-1]
        trimmed_alignmnets_spec = ''
        one_of = ''
        for i in self.otu_must_have_one_of:
            one_of += '[ '
            for j in i:
                one_of += j+' '
            one_of += ']'
        if (self.define_trimmed_alns) > 0:
            for i in self.define_trimmed_alns:
                trimmed_alignmnets_spec += i
        return ("Concatenation named %s, with loci %s,\n"
                "of which %s must exist for all species\n"
                "and at least one of each group of %s is represented.\n"
                "Alignments with the following names: %s are prefered"
                % (self.name, loci_string, must_have, one_of, trimmed_alignmnets_spec))
        
        
        
##############################################################################################
if False:
    """
    Reprophylo Project Utilities
    
    Used in the Project class but are not in the classe's methods
    """
##############################################################################################


## Git management

__builtin__.git = False


# git log template
gitline = "<<<<\n%s\nSTDOUT:\n%s\nSTDERR:%s\n>>>>\n"

def undate_git_log(pj, out, err):
    if not err:
        err = 'None'
    if not out:
        out = 'None'
    pj.git_log += gitline%(str(time.asctime()),str(out), str(err))

def start_git(pj):
    __builtin__.git = True # flag it on
    cwd = os.getcwd() 
    if os.path.isdir(cwd + '/.git'):
        # a repo exists, check it belongs to this project by checking the description
        try:
            assert open(cwd + '/.git/description').read().strip().rstrip() == pj.pickle_name.strip().rstrip()
            warnings.warn('Git repository exists for this Project')
        except:
            raise RuntimeError('The Git repository in the CWD does not belong to this project. Either the pickle'+
                              ' moved, or this is a preexsisting repo. Try one of the following: Delete the local '+
                              ' .Git dir if you don\'t need it, move the pickle and the notebook to a new work dir,'+
                              ' or if possible, move them back to their original location. You may also disable Git'+
                              ' by with stop_git().')
    else:
        # start a rep
        out, err = rpgit.gitInit()
        undate_git_log(pj, out, err)
        # write the pickle name as the repo description
        hndl = open(cwd + '/.git/description', 'wt')
        hndl.write(pj.pickle_name.strip().rstrip())
        hndl.close()
        warnings.warn('The new repository is called %s.'%open(cwd + '/.git/description', 'r').read().rstrip())
    
    # list scripts and notebooks
    import fnmatch
    matches = []
    for root, dirnames, filenames in os.walk(cwd):
        for filename in fnmatch.filter(filenames, '*.py'):
            matches.append(os.path.join(root, filename))
        for filename in fnmatch.filter(filenames, '*.ipynb'):
            matches.append(os.path.join(root, filename))
    # git add scripts and notebooks
    for match in matches:
        out, err = rpgit.gitAdd(match)
        undate_git_log(pj, out, err)
    # commit scripts and notebooks
    comment = "%i script file(s) from %s" % (len(matches), time.asctime())
    out, err = rpgit.gitCommit(comment)
    undate_git_log(pj, out, err)
    
    
    
def stop_git():
    __builtin__.git = False # flad it off
    cwd = os.getcwd()
    
    # list, git add and git commit scripts and notebooks
    # list
    import fnmatch
    matches = []
    for root, dirnames, filenames in os.walk(cwd):
        for filename in fnmatch.filter(filenames, '*.py'):
            matches.append(os.path.join(root, filename))
        for filename in fnmatch.filter(filenames, '*.ipynb'):
            matches.append(os.path.join(root, filename))
    # add
    for match in matches:
        out, err = rpgit.gitAdd(match)
        undate_git_log(pj, out, err)
    # commit
    comment = "%i script file(s) from %s" % (len(matches), time.asctime())
    out, err = rpgit.gitCommit(comment)
    undate_git_log(pj, out, err)
    
## end git management

def platform_report():
    
    """ 
    Prints machine specs, os specs and dependencies at time of execution
    
    >>> isinstance(platform_report(), list)
    True
    """
    import pkg_resources
    modules = [] # and their versions
    for i in ('ete2','biopython','dendropy','cloud'):
        try:
            modules.append(i+' version: '+
                           pkg_resources.get_distribution(i).version)
        except:
            pass
    modules.append('reprophylo version %s'%str(reprophyloversion))
                   
    return(['Platform: '+platform.platform(aliased=0, terse=0),
            'Processor: '+platform.processor(),
            'Python build: '+platform.python_build()[0] + platform.python_build()[1],
            'Python compiler: '+platform.python_compiler(),
            'Python implementation: ' +platform.python_implementation(),
            'Python version: ' + platform.python_version()]+
             modules+
            ['User: ' +platform.uname()[1]])



def write_alns(pj, format = 'fasta'):
    """
    Writes untrimmed sequence alignment files that are in pj in a biopython format
    """
    
    if len(pj.alignments.keys()) == 0:
        raise IOError('Align the records first')
    else:
        for key in pj.alignments:
            AlignIO.write(pj.alignments[key], key+'_aln.'+format, format)



def keep_feature(feature, loci):
    
    """ Returns true if a feature's type is in one of the loci and if the gene
    or product qualifiers is in the aliases of one of the loci, for data collection
    from a genbank or embl file
    
    # making a dummy feature
    >>> coi = Locus('dna','CDS','coi', ['cox1','COX1','coi','COI','CoI'])
    >>> location = FeatureLocation(1,100)
    >>> feature = SeqFeature()
    >>> feature.location = location
    >>> feature.type = 'CDS'
    >>> feature.qualifiers['gene'] = ['CoI']
    
    # testing if fits any of the Project Locus objects
    >>> a = keep_feature(feature, [coi])
    >>> print(a)
    True"""
    
    keep = 0
    for g in loci:
        if not g.name in g.aliases:
            g.aliases.append(g.name)
        if feature.type == 'source':
            keep = 1
        elif feature.type == g.feature_type:
            qual = None
            if 'gene' in feature.qualifiers.keys():
                qual = 'gene'
            elif 'product' in feature.qualifiers.keys():
                qual = 'product'
            if qual and feature.qualifiers[qual][0] in g.aliases:
                keep = 1
    if keep == 1:
        return True
    else:
        return False
    return
  


def dwindle_record(record, loci):
    
    """ 
    Retains only features that are called by Locus objects and records with features that are
    called by Locus objects
    
    # Making a dummy locus    
    >>> coi = Locus('dna','CDS','coi', ['cox1','COX1','coi','COI','CoI'])
    
    # Making a dummy record with a feature that fits a Locus object (kept_feature)
    # and a feature that does not (dwindled_feature)
    >>> location = FeatureLocation(1,100)
    >>> kept_feature = SeqFeature()
    >>> kept_feature.location = location
    >>> kept_feature.type = 'CDS'
    >>> kept_feature.qualifiers['gene'] = ['CoI']
    >>> dwindled_feature = SeqFeature()
    >>> dwindled_feature.location = location
    >>> dwindled_feature.type = 'rRNA'
    >>> dwindled_feature.qualifiers['gene'] = ['LSU']
    >>> s = 'atgc'*1000
    >>> record = SeqRecord(seq=Seq(s, IUPAC.ambiguous_dna), id='1', description='spam')
    >>> record.features.append(kept_feature)
    >>> record.features.append(dwindled_feature)
    >>> print(len(record.features))
    2
    
    # Dwindling the record
    >>> a = dwindle_record(record, [coi])
    >>> print(len(record.features))
    1
    """
    
    dwindled_features = []
    feature_count = 0
    for feature in record.features:
        if keep_feature(feature, loci)== True:
            # determine feature id       
            if feature.type == 'source' and not 'feature_id' in feature.qualifiers.keys():
                feature.qualifiers['feature_id'] = [record.id + '_source']
            elif not 'feature_id' in feature.qualifiers.keys():
                feature.qualifiers['feature_id'] = [record.id + '_f' + str(feature_count)]
                feature_count += 1
            # determine prop ambiguity and GC content       
            if not feature.type == 'source':
                feature_seq = feature.extract(record.seq)
                degen = len(feature_seq)
                for i in ['A','T','G','C','U','a','t','g','c','u']:
                    degen -= feature_seq.count(i)
                feature.qualifiers['GC_content'] = [str(GC(feature_seq))]
                feature.qualifiers['nuc_degen_prop'] = [str(float(degen)/len(feature_seq))]
                if 'translation' in feature.qualifiers.keys():
                    transl = feature.qualifiers['translation'][0]
                    degen = 0
                    for i in ['B', 'X', 'Z', 'b', 'x', 'z']:
                        degen += transl.count(i)
                    feature.qualifiers['prot_degen_prop'] = [str(float(degen)/len(transl))]                    
            dwindled_features.append(feature)
    record.features = dwindled_features
    return record
            
 
    
def is_embl_or_gb(input_filename):
    suffixes = ['.gb','.embl']
    gb = False
    for s in suffixes:
        if s in input_filename:
            gb = True
    return gb



def parse_input(input_filename, fmt):
    return SeqIO.parse(input_filename, fmt)



def list_to_string(List):
    
    """
    Handles list printing as a nice string in the pj.write(format="csv") method
    
    >>> L = ['a','b','b']
    >>> print(list_to_string(L))
    a;b;b
    """
    
    string = ''
    for i in List:
        if type(i) is str and '\n' in i:
            string += lines_to_line(i).rstrip()+';'
        else:
            string += str(i)+';'
    return string[:-1]




def lines_to_line(lines):
    
    """
    Replaces newline with space in the pj.write(format="csv") method
    """
    
    lines = lines.split('\n')
    return (' ').join(lines)




def type_to_single_line_str(var):
    
    """
    Returns any type as a one line string for the pj.write(format="csv") method
    """
    
    if type(var) is str and '\n' in var:
        return lines_to_line(var)
    elif type(var) is str or type(var) is int or type(var) is float:
        return str(var)
    elif type(var) is list and len(var) == 1:
        return str(var[0])
    elif type(var) is list and len(var) > 0:
        return list_to_string(var)
    else:
        return var



def get_qualifiers_dictionary(project, feature_id):
    
    """
    Takes sequence record annotation, source qualifiers and feature qualifiers and puts them
    in a flat dictionary
    
    This is being replaced by __get_qualifiers_dictionary__ which deals with the records as a dict
    and is much faster. Eventually, records will be handled as a dict throughout, instead of as
    a list.
        
    # Making a dummy locus    
    >>> coi = Locus('dna','CDS','coi', ['cox1','COX1','coi','COI','CoI'])
    
    # Making a dummy Project
    >>> pj = Project([coi], git=False)
    
    # making a dummy record
    >>> s = 'atgc'*1000
    >>> location = FeatureLocation(1,100)
    >>> feature = SeqFeature()
    >>> feature.location = location
    >>> feature.type = 'CDS'
    >>> feature.qualifiers['gene'] = ['CoI']
    >>> feature.qualifiers['feature_id'] = ['1_f0']
    >>> source = SeqFeature()
    >>> source.location = FeatureLocation(0,3999)
    >>> source.type = 'source'
    >>> source.qualifiers['organism'] = ['Tetillda radiata']
    >>> record = SeqRecord(seq=Seq(s, IUPAC.ambiguous_dna), id='1', description='spam')
    >>> record.features.append(feature)
    >>> record.features.append(source)
    >>> record.annotations["evidence"] = 'made up'
    >>> pj.records = [record]
    
    # executing get_qualifiers_dictionary()
    >>> qual_dict = get_qualifiers_dictionary(pj, '1_f0')
    >>> qual_items = qual_dict.items()
    >>> qual_items.sort(key = lambda i: i[0])
    >>> for key, val in qual_items: print(key.ljust(20,' ') + val.ljust(20,' '))
    annotation_evidence made up             
    feature_id          1_f0                
    gene                CoI                 
    record_id           1                   
    source_organism     Tetillda radiata    
    """
    if type(feature_id) is list and len(feature_id) > 1:
        raise IOError('get_qualifiers_dictionary takes one feature_id at a time')
    if type(feature_id) is list:
        feature_id = feature_id[0]
    record_id = feature_id.rpartition('_')[0]
    qualifiers_dictionary={'record_id': record_id}
    for record in project.records:
        if record.id in feature_id:
            for annotation in record.annotations.keys():
                qualifiers_dictionary['annotation_'+annotation]=record.annotations[annotation]
            for feature in record.features:
                if feature.type == 'source':
                    for qualifier in feature.qualifiers.keys():
                        qualifiers_dictionary['source_'+qualifier]=feature.qualifiers[qualifier][0]
                elif feature.qualifiers['feature_id'][0] == feature_id:
                    for qualifier in feature.qualifiers.keys():
                        qualifiers_dictionary[qualifier]=feature.qualifiers[qualifier][0]
    return qualifiers_dictionary

def __get_qualifiers_dictionary__(project, feature_id):
    
    """
    This will replace the public version. It uses the Project._records_dict to  pull 
    the record using the record id, instead of iterating Project.records, which is very slow.
    
    It requires Project.__records_list_to_dict__() to execute beforehand.
    """
    if type(feature_id) is list and len(feature_id) > 1:
        raise IOError('get_qualifiers_dictionary takes one feature_id at a time')
    if type(feature_id) is list:
        feature_id = feature_id[0]
    record_id = feature_id.rpartition('_')[0]
    record = project._records_dict[record_id]
    qualifiers_dictionary={'record_id': record_id}
    for annotation in record.annotations.keys():
        qualifiers_dictionary['annotation_'+annotation]=record.annotations[annotation]
    for feature in record.features:
        if feature.type == 'source':
            for qualifier in feature.qualifiers.keys():
                qualifiers_dictionary['source_'+qualifier]=feature.qualifiers[qualifier][0]
        elif feature.qualifiers['feature_id'][0] == feature_id:
            for qualifier in feature.qualifiers.keys():
                qualifiers_dictionary[qualifier]=feature.qualifiers[qualifier][0]
    return qualifiers_dictionary

def seq_format_from_suffix(suffix):
    
    """
    Guesses input format from suffix
    
    >>> print(seq_format_from_suffix('gb'))
    genbank
    """
    
    suffixes = {'fasta': ['fas','fasta','fa','fna'],
                'genbank': ['gb','genbank'],
                'embl': ['embl']}
    found = False
    for key in suffixes.keys():
        if suffix in suffixes[key]:
            found = True
            return key
    if not found:
        raise RuntimeError(suffix+' is not a recognised suffix of an unaligned sequence file')



def read_feature_quals_from_tab_csv(csv_filename):
                   
    """
    This is used to update feature qualifiers from a tab delimited file
    """               
                   
    import re
    header = open(csv_filename, 'r').readlines()[0].rstrip().split('\t')
    feature_id_col = header.index('feature_id')
    taxonomy_col = header.index('taxonomy')
    seq_col = header.index('seq')
    translation_col = None
    if 'translation' in header:
        translation_col = header.index('translation')
    csv_info = {}
    for line in [l.rstrip().split('\t') for l in open(csv_filename, 'r').readlines()[1:]]:
        if not line[0] in csv_info.keys():
            csv_info[line[0]] = {'source':{},
                                 'taxonomy':[],
                                 'features':{}
                                 }
        if csv_info[line[0]]['taxonomy'] == []:
            csv_info[line[0]]['taxonomy'] = line[taxonomy_col].split(';')
        csv_info[line[0]]['features'][line[feature_id_col]] = {}
        get_source = False
        if csv_info[line[0]]['source'] == {}:
            get_source = True
        for i in range(len(header)):
            if get_source and 'source:_' in header[i]:
                qual_name = re.sub('source:_','',header[i])
                if not line[i] == 'null' and not line[i] == '' and line[i]:
                    csv_info[line[0]]['source'][qual_name] = line[i].split(';')
            elif (not 'source:_' in header[i] and not line[i] == 'null' and not line[i] == '' and line[i] and
                  not i in [seq_col, translation_col, taxonomy_col, feature_id_col]):
                csv_info[line[0]]['features'][line[feature_id_col]][header[i]] = line[i].split(';')
    return csv_info


## Alignment statistics
                   
def count_positions(aln_column):
    counts = {}
    for i in aln_column:
        if i in counts.keys():
            counts[i] += 1
        else:
            counts[i] = 1
    return counts

def global_aln_stats(aln_obj):
    total_gaps = 0
    prop_list = []
    non_uniform_count = aln_obj.get_alignment_length()
    parsimony_informative = 0
    for i in range(aln_obj.get_alignment_length()):
        total_gaps += aln_obj[:, i].count('-')
        prop_list.append(aln_obj[:, i].count('-')/float(len(aln_obj)))
        if len(count_positions(aln_obj[:, i]).keys()) == 1:
            non_uniform_count -= 1
        elif (len(count_positions(aln_obj[:, i]).keys()) == 2 and
              '-' in  count_positions(aln_obj[:, i]).keys()):
            non_uniform_count -= 1
        if len([p for p in count_positions(aln_obj[:, i]).keys() if (p != '-' and  count_positions(aln_obj[:, i])[p] > 1)]) > 1:
            parsimony_informative += 1
    mean_gap_prop = sum(prop_list)/aln_obj.get_alignment_length()
    return (mean_gap_prop, non_uniform_count, parsimony_informative)

def count_undetermined_lines(aln_obj, cutoff=0):
    count = 0
    ids = []
    
    if aln_obj.get_alignment_length() < cutoff*2:
        warnings.warn('The cutoff to exclude a sequence is more than half of the alignmnet length')
    elif aln_obj.get_alignment_length() <= cutoff:
        raise RuntimeWarning('The cutoff to exclude a sequence is as long or longer than the alignment')
    
    for seq in aln_obj:
        if str(seq.seq).count('-') >= aln_obj.get_alignment_length()-cutoff:
            count += 1
            ids.append(seq.id)
    return count, ids

def count_collapsed_aln_seqs(aln_obj):
    count = 1
    seen_seqs = [str(aln_obj[0].seq)]
    for seq in aln_obj[1:]:
        str_seq = str(seq.seq)
        if len([s for s in seen_seqs if (str_seq in s or s in str_seq or s == str_seq)]) == 0:
            count += 1
        seen_seqs.append(str_seq)
    return count

def aln_summary(aln_obj, cutoff=0):
    lines = ["Alignment length: %i" % aln_obj.get_alignment_length(),
             "Number of rows: %i" % len(aln_obj),
             "Unique sequences: %i"%count_collapsed_aln_seqs(aln_obj),
             "Average gap prop.: %f\nVariable columns: %i\nParsimony informative: %i"
             %global_aln_stats(aln_obj),
             "Undetermined sequences: %i"%(count_undetermined_lines(aln_obj, cutoff=cutoff)[0]),
             "Undetermined sequence cutoff: %i"%cutoff
             ]
    return [lines, len(aln_obj), count_undetermined_lines(aln_obj, cutoff=cutoff), count_collapsed_aln_seqs(aln_obj)]        
            
##
                   
def loci_list_from_csv(loci):
    """
    Parse the loci csv file given to Project
    """               
    
    # verify format
    if any(len(line.split(',')) >= 4 for line in open(loci, 'r').readlines()):
        pass
    else:
        raise IOError("File %s has no valid loci of format char_type,feature_type,name,aliases"%loci)
        
        
    loci_dict = {}
    loci_list = []
    for line in [line.rstrip() for line in open(loci, 'r').readlines() if len(line.rstrip()) > 0]:
        # verify format
        if len(line.split(',')) < 4:
            raise IOError("The line %s in file %s is missing arguments. Needs at least char_type,feature_type,name,aliases"%
                          (line.rstrip(), loci))
        # look for synonyms
        else:
            group = None
            try:
                group = int(line.rstrip().split(',')[-1])
            except:
                pass
            
            if group:
                locus_exists = False
                for name in loci_dict:
                    if 'group' in loci_dict[name].keys() and loci_dict[name]['group'] == group:
                        loci_dict[name]['aliases'] += line.split(',')[3:-1]
                        locus_exists = True
                if not locus_exists:
                    loci_dict[line.split(',')[2]] = {'group': int(line.rstrip().split(',')[-1]),
                                                     'char_type': line.split(',')[0],
                                                     'feature_type': line.split(',')[1],
                                                     'aliases': line.split(',')[3:-1]
                                                     }
            else:
                loci_dict[line.split(',')[2]] = {'group': None,
                                                 'char_type': line.split(',')[0],
                                                 'feature_type': line.split(',')[1],
                                                 'aliases': line.split(',')[3:]
                                                 }
                
            
            
    for name in loci_dict:
        loci_list.append(Locus(loci_dict[name]['char_type'],
                               loci_dict[name]['feature_type'],
                               name,
                               loci_dict[name]['aliases']))
    return loci_list

def parse_paup_charset(nexus_filename):
    
    """ 
    Takes a nexus file with PAUP style charset commands.
    Returns a dictionary with partition names as keys and a list of
    integers representing the start and end of the partition as a value.
    Position count starts from 0.
    
    Handles paup commands of the following format:
    CHARSET  locus_name=1-129;
    or
    charset locus_name = 1 - 129 ;
    """
    
    try:
        AlignIO.read(nexus_filename, 'nexus')
    except:
        n = len(list(AlignIO.parse(nexus_filename, 'nexus')))
        raise IOError('Cannot handle more then one matrix in %s. Got %i matrices'%
                     (nexus_filename, n))
    
    charsets = {}
    charset_lines = [l for l in open(nexus_filename,'r').readlines() if 
                     (l.startswith('CHARSET') or l.startswith('charset'))]
    if len(charset_lines) == 0:
        raise IOError("There are no CHARSET commands in %s"%nexus_filename)
    for line in charset_lines:
        try:
            info = line.split()[1].split(';')[0]
            locus_name, range = info.split('=')
            locus_name = locus_name.strip().rstrip()
            start = int(range.split('-')[0].strip().rstrip())-1
            end = int(range.split('-')[1].strip().rstrip())-1
            charsets[locus_name] = [start,end]
        except:
            raise IOError('Expects "charset set_name = start_int - end_int;"'+
                          ' (case insensitive, spaces around the "=" or "-" not mandatory). Got %s'%line)
    return charsets
        
def pj_from_nexus_w_charset(nexus_filename, output_dir, char_type,
                            feature_type, project=False, pickle=False, git=False):
    
    """ 
    Takes a nexus file with PAUP style charset commands as input.
    Creates a separate fasta file for each alignment partition
    Returns a list of fasta filenames and a list of Locus objects
    If project==True, returns a Project instance with the loci, alignments and records instead
    """
    
    from reprophylo import Locus
    from Bio import AlignIO
    
    charsets = parse_paup_charset(nexus_filename)
    
    alignment =  AlignIO.read(nexus_filename, 'nexus')
    filenames = []
    loci_list = []
    for locus_name in charsets:
        s = charsets[locus_name][0]
        e = charsets[locus_name][1]
        outname = "%s/%s.fasta"%(output_dir,locus_name)
        AlignIO.write(alignment[:, s:e], outname, 'fasta')
        filenames.append(outname)
        loci_list.append(Locus(char_type, feature_type, locus_name, [locus_name]))
    
    if project:
        from reprophylo import Project
        pj = Project(loci_list, pickle=pickle, git=git)
        i=1
        for f in filenames:
            locus_name = f.split('/')[-1].split('.')[0]
            print '%i/%i reading %s'%(i,len(filenames), locus_name)
            i += 1
            pj.read_alignment(f, char_type, feature_type, locus_name)
        return pj
            
    else:  
        return filenames, loci_list

##############################################################################################
class Project:
##############################################################################################
    
    """
    The Project class contians all the data and has methods to analyze it. It allows for
    experimental analysis by running alternative analyses and formally comparing the 
    outputs. The pickle_pj() function allows to pickle the project, including the data,
    intermediates and results, as well as a description of the methods.It allows for a rerun
    of the whole analysis as is, as well as for a reconfiguration of the analysis or addition
    of data. If git is installed, it can be called by 'import rpgit'. As a result, git can be 
    initiated using start_git(). A git repository will be created in the CWD, if it doesn't already exist. 
    Input datasets, .py, .ipynb and .pkl files in the CWD will be version controlled. 
    Version control can be paused in the middle of the script
    by calling stop_git() and restarted by calling start_git() again.
    """

    def __init__(self, loci, pickle=None, git=True):
        
        """
        # making dummy loci
        >>> coi = Locus('dna','CDS','coi',['COX1','cox1'])
        >>> ssu = Locus('dna','rRNA','18S',['18S','SSU'])
        
        # Making a Project object
        >>> pj = Project([coi,ssu], git=False)
        >>> print(str(pj))
        Project object with the loci coi,18S,
        """
        self.records = []
        self._records_dict = {}
        self.starttime = str(time.asctime())
        self.user = None
        if os.path.isfile('USER'):
            self.user = []
            for line in open('USER','r').readlines():
                key, arg = line.rstrip().split('=')
                self.user.append([key, arg])
        self.loci = loci
        self.records_by_locus = {}
        self.concatenations = []
        self.alignments = {}
        self.trimmed_alignments = {}
        self.trees = {}
        self.used_methods = {}
        self.sets = {}
        self.git_log = ''
        
        self.pickle_name=pickle
        
        if self.pickle_name and os.path.exists(self.pickle_name):
            raise IOError('Pickle %s exists. If you want to keep using it do pj=unpickle_pj(\'%s\') instead.'%(self.pickle_name,self.pickle_name))
        
        if git and not self.pickle_name:
            raise IOError('Must have pickle to run Git. Either specify pickle or git=False')
        elif git:
            start_git(self)
         
        self.defaults = {'raxmlHPC': programspath+'raxmlHPC-PTHREADS-SSE3',
                         'mafft': 'mafft',
                         'muscle': programspath+'muscle',
                         'trimal': programspath+'trimal',
                         'pb': programspath+'pb',
                         'bpcomp': programspath+'bpcomp',
                         'tracecomp': programspath+'tracecomp',
                         'fasttree': programspath+'FastTreeMP',
                         'pal2nal': programspath+'pal2nal.pl',
                         # PROGRAM PLUG
                         # 'program name': programpath+'the basic command'
                        }
    
        seen = []
        if isinstance(loci,list):
            for locus in loci:
                if not isinstance(locus, Locus):
                    raise TypeError("Expecting Locus object in loci list. "+locus+
                                    " not a Locus object")
                if locus.name in seen:
                    raise NameError('Locus ' + locus.name + ' apears more than once in self.loci')
                else:
                    seen.append(locus.name)
        elif isinstance(loci,str):
            self.loci = loci_list_from_csv(loci)
            #print 'Read the following loci from file %s:'%loci
            #for l in self.loci:
                #print str(l)
        self.aln_summaries = []
        
        if self.pickle_name:
            pickle_pj(self, self.pickle_name, track=False)
        
        if __builtin__.git and self.pickle_name:
            import rpgit
            comment = "%s from %s" % (str(self), time.asctime())
            out, err = rpgit.gitAdd(self.pickle_name)
            undate_git_log(self, out, err)
            cwd = os.getcwd()
            import fnmatch
            matches = []
            for root, dirnames, filenames in os.walk(cwd):
                for filename in fnmatch.filter(filenames, '*.py'):
                    matches.append(os.path.join(root, filename))
                for filename in fnmatch.filter(filenames, '*.ipynb'):
                    matches.append(os.path.join(root, filename))
            for match in matches:
                out, err = rpgit.gitAdd(match)
                undate_git_log(self, out, err)
            out, err = rpgit.gitCommit(comment)
            undate_git_log(self, out, err)
                
                
    def __str__(self):
        loci_string = ''
        for i in self.loci:
            loci_string += i.name+','
        return 'Project object with the loci '+loci_string


    def __records_list_to_dict__(self):
        self._records_dict = SeqIO.to_dict(self.records)
        
    def last_git_log(self):
        print self.git_log.split('<<<<')[-1]
        
    def show_commits(self):
        print rpgit.gitLog()[0]
        
    ###################################
    # Project methods for reading data
    ###################################  



    def read_embl_genbank(self, input_filenames_list):
        
        """
        Read a file from Genbank of EMBL
        
        >>> input_filenames = ['test-data/test.gb']
        >>> locus = Locus('dna', 'CDS', 'coi', ['cox1','COX1','coi','COI','CoI'])
        >>> pj = Project([locus], git=False)
        >>> print(len(pj.records))
        0
        >>> pj.read_embl_genbank(input_filenames)
        >>> print(len(pj.records))
        89
        """
        
        if __builtin__.git:
            import rpgit
        else:
            warnings.warn('Version control off')    
        generators = []
        for input_filename in input_filenames_list:
            if __builtin__.git:
                import rpgit
                out, err = rpgit.gitAdd(input_filename)
                undate_git_log(self, out, err)
            generators.append(parse_input(input_filename, 'gb'))
            for generator in generators:
                for record in generator:
                    dwindled_record = dwindle_record(record, self.loci)
                    if len(record.features) > 1:
                        self.records.append(dwindled_record)
                    elif len(record.features) == 1 and not record.features[0].type == 'source':
                        self.records.append(dwindled_record)
        if __builtin__.git:
            import rpgit
            comment = "%i genbank/embl data file(s) from %s" % (len(input_filenames_list), time.asctime()) 
            for filename in input_filenames_list:
                out, err = rpgit.gitAdd(filename)
                undate_git_log(self, out, err)
            cwd = os.getcwd()
            import fnmatch
            matches = []
            for root, dirnames, filenames in os.walk(cwd):
                for filename in fnmatch.filter(filenames, '*.py'):
                    matches.append(os.path.join(root, filename))
                for filename in fnmatch.filter(filenames, '*.ipynb'):
                    matches.append(os.path.join(root, filename))
            for match in matches:
                out, err = rpgit.gitAdd(match) 
                undate_git_log(self, out, err)
            out, err = rpgit.gitCommit(comment) 
            undate_git_log(self, out, err)
            
        
        
    def read_denovo(self, input_filenames, char_type, format = 'fasta'):
        
        """
        Include records from a fasta file. Fasta file records will be given record ids 
        of the form 'denovo1'. The record.id and record.description will be placed in a
        source feature under the 'original_id' and 'original_desc' qualifiers. Denovo sequences
        require the use of the add_feature_to_record() method in order to be included in the
        anaysis.
        

        >>> input_filenames = ['test-data/test.gb']
        >>> locus = Locus('dna', 'CDS', 'coi', ['cox1','COX1','coi','COI','CoI'])
        >>> pj = Project([locus], git=False)
        >>> print(len(pj.records))
        0
        >>> pj.read_embl_genbank(input_filenames)
        >>> print(len(pj.records))
        89
        >>> input_filenames = ['test-data/test.fasta']
        >>> pj.read_denovo(input_filenames, 'dna')
        1
        >>> print(len(pj.records))
        90
        
        # Since the denovo sequence has no feature it is not included
        >>> pj.extract_by_locus()
        >>> print(len(pj.records_by_locus['coi']))
        89
        
        # Making a feature for the denovo record.
        >>> pj.add_feature_to_record('denovo0', 'CDS',  qualifiers={'gene': 'coi'})
        'denovo0_f0'
        >>> pj.extract_by_locus()
        >>> print(len(pj.records_by_locus['coi']))
        90
        """
        
        if __builtin__.git:
            import rpgit
        else:
            warnings.warn('Version control off')
        
        count = 0
        # start the counter where it stoped the last time we read denovo things
        for record in self.records:
            if 'denovo' in record.id:
                serial = int(record.id[6:])
                if serial >= count:
                    count = serial+1
        for input_filename in input_filenames:
            if __builtin__.git:
                import rpgit
                rpgit.gitAdd(input_filename)
            denovo = SeqIO.parse(input_filename, format)
            for record in denovo:
                source = SeqFeature(FeatureLocation(0, len(record.seq)), type='source', strand=1)
                source.qualifiers['original_id'] = [record.id]
                source.qualifiers['original_desc'] = [(' ').join(record.description.split()[1:])]
                record.id = 'denovo'+str(count)
                record.name = record.id
                source.qualifiers['feature_id'] = [record.id+'_source']
                record.features = [source]
                if '-' in str(record.seq):
                    record.seq = Seq(str(record.seq).replace('-',''))
                    warnings.warn("Reseting gaps in records from %s"%input_filename)
                if char_type == 'prot':
                    record.seq.alphabet = IUPAC.protein
                    #feature = SeqFeature(FeatureLocation(0, len(record.seq)), type='Protein', strand=1)
                    #feature.qualifiers['translation'] = [str(record.seq)]
                    #feature.qualifiers['protonly']=['true']
                elif char_type == 'dna':
                    record.seq.alphabet = IUPAC.ambiguous_dna
                count += 1
                self.records.append(record)
                
        if __builtin__.git:
            import rpgit
            comment = "%i denovo data file(s) from %s" % (len(input_filenames), time.asctime())
            for filename in input_filenames:
                out, err = rpgit.gitAdd(filename)
                undate_git_log(self, out, err)
            cwd = os.getcwd()
            import fnmatch
            matches = []
            for root, dirnames, filenames in os.walk(cwd):
                for filename in fnmatch.filter(filenames, '*.py'):
                    matches.append(os.path.join(root, filename))
                for filename in fnmatch.filter(filenames, '*.ipynb'):
                    matches.append(os.path.join(root, filename))
            for match in matches:
                out, err = rpgit.gitAdd(match)
                undate_git_log(self, out, err)
            out, err = rpgit.gitCommit(comment)
            undate_git_log(self, out, err)
        return count       
    
    def read_alignment(self, filename, char_type, feature_type, locus_name, format="fasta", aln_method_name = "ReadDirectly", exclude=[]):
    
        if __builtin__.git:
            import rpgit
        else:
            warnings.warn('Version control off')
            
        if not any([locus.name == locus_name for locus in self.loci]):
            raise RuntimeError("Locus %s does not exist"%locus_name)
        elif not [locus for locus in self.loci if locus.name == locus_name][0].char_type == char_type:
            raise RuntimeError("%s is not a %s locus"%(locus_name, char_type))
        elif not [locus for locus in self.loci if locus.name == locus_name][0].feature_type == feature_type:
            raise RuntimeError("The feature_type %s is not %s"%(locus_name, feature_type))
            
        count = 0
        # start the counter where it stoped the last time we read denovo things
        for record in self.records:
            if 'denovo' in record.id:
                serial = int(record.id[6:])
                if serial >= count:
                    count = serial+1
        # Read the alignment:
        raw_aln_input = list(AlignIO.read(filename, format))
        
        
        
        # make records
        records = []
        aln_records = []
        for record in raw_aln_input:
            total_seq_len = len(record.seq)
            if (str(record.seq).count('-') + str(record.seq).count('.') + 
                str(record.seq).count('?') + str(record.seq).count('N') +
                str(record.seq).count('n') + str(record.seq).count('X') + 
                str(record.seq).count('x')) == total_seq_len:
                print 'dropping seq %s in locus %s: missing data'%(locus_name,record.id)
            elif not record.id in exclude:
                # remove gaps
                new_record = SeqRecord(seq=Seq(str(record.seq).replace('-','').replace('.','')))
                aln_record = SeqRecord(seq=Seq(str(record.seq).replace('.','-')))
                
                #set alphabet
                if char_type == 'prot':
                    new_record.seq.alphabet = IUPAC.protein
                    aln_record.seq.alphabet = IUPAC.protein
                    
                elif char_type == 'dna':
                    new_record.seq.alphabet = IUPAC.ambiguous_dna
                    aln_record.seq.alphabet = IUPAC.ambiguous_dna
                    
                # set denovo record id
                new_record.id = 'denovo%i'%count
                new_record.name = new_record.id
                
                # set source and first feature
                source = SeqFeature(FeatureLocation(0, len(new_record.seq)), type='source', strand=1)
                source.qualifiers['original_id'] = [record.id]
                source.qualifiers['original_desc'] = [(' ').join(record.description.split()[1:])]
                source.qualifiers['feature_id'] = [new_record.id+'_source']
                feature = SeqFeature(FeatureLocation(0, len(new_record.seq)), type=feature_type, strand=1)
                feature.qualifiers['feature_id'] = [new_record.id+'_f0']
                feature.qualifiers['gene'] = [locus_name]
                feature_seq = feature.extract(new_record.seq)
                degen = len(feature_seq)
                if char_type == 'dna':    
                    for i in ['A','T','G','C','U','a','t','g','c','u']:
                        degen -= feature_seq.count(i)
                    feature.qualifiers['GC_content'] = [str(GC(feature_seq))]
                    feature.qualifiers['nuc_degen_prop'] = [str(float(degen)/len(feature_seq))]
                    warnings.warn("To get translations, add a feature manually")
                elif char_type == 'prot': 
                    degen = 0
                    for i in ['B', 'X', 'Z', 'b', 'x', 'z']:
                        degen += feature_seq.count(i)
                    feature.qualifiers['prot_degen_prop'] = [str(float(degen)/len(feature_seq))]   
                new_record.features = [source, feature]
                
                aln_record.id = 'denovo%i_f0'%count
                
                count += 1
                
                records.append(new_record)
                aln_records.append(aln_record)
            
        token = "%s@%s"%(locus_name, aln_method_name)
        
        if token in self.alignments.keys():
            warnings.warn("Replacing alignment %s"%token)
        # need to add the denovo id's inside the alignment    
        self.alignments[token] = MultipleSeqAlignment(aln_records)
        self.records += records
        self.extract_by_locus()
        
        if __builtin__.git:
            import rpgit
            comment = "Alignment file %s" % (time.asctime())
            out, err = rpgit.gitAdd(filename)
            undate_git_log(self, out, err)
            cwd = os.getcwd()
            import fnmatch
            matches = []
            for root, dirnames, filenames in os.walk(cwd):
                for filename in fnmatch.filter(filenames, '*.py'):
                    matches.append(os.path.join(root, filename))
                for filename in fnmatch.filter(filenames, '*.ipynb'):
                    matches.append(os.path.join(root, filename))
            for match in matches:
                out, err = rpgit.gitAdd(match)
                undate_git_log(self, out, err)
            out, err = rpgit.gitCommit(comment)
            undate_git_log(self, out, err)

    def add_feature_to_record(self, record_id, feature_type, location='full', qualifiers={}):
    
        """
        # Making a dummy locus    
        >>> coi = Locus('dna','CDS','coi', ['cox1','COX1','coi','COI','CoI'])
    
        # Making a dummy Project
        >>> pj = Project([coi], git=False)
    
        # making a dummy record
        >>> s = 'atgc'*1000
        >>> source = SeqFeature()
        >>> source.location = FeatureLocation(0,3999)
        >>> source.type = 'source'
        >>> record = SeqRecord(seq=Seq(s, IUPAC.ambiguous_dna), id='1', description='spam')
        >>> record.features.append(source)
        >>> pj.records = [record]
        >>> print(len(pj.records[0].features))
        1
        
        # adding a feature to a record in the pj
        >>> import warnings
        >>> warnings.simplefilter("ignore")
        >>> pj.add_feature_to_record('1', 'CDS', qualifiers={'gene': 'madeuplocus'})
        '1_f0'
        >>> print(len(pj.records[0].features))
        2
        """
    
        for record in self.records:
            if record.id == record_id:
                #determine new feature id
                feature_id = None
                serials = []
                if len(record.features) > 0:
                    for feature in record.features:
                        if 'feature_id' in feature.qualifiers.keys():
                            if '_f' in feature.qualifiers['feature_id'][0]:
                                f = feature.qualifiers['feature_id'][0]
                                serials.append(int(f.split('_')[1][1:]))
                    serials.sort(reverse = True)
                if len(serials) > 0:
                    feature_id = record.id + '_f' + str(serials[0]+1)
                else:
                    feature_id = record.id + '_f0'
                feature = None
                if location == 'full':
                    feature = SeqFeature(FeatureLocation(0, len(record.seq)),
                                         type=feature_type,
                                         strand=1)
                elif isinstance(location, list):
                    for i in location:
                        if not isinstance(i, list):
                            raise RuntimeError('\'location\' takes either \'full\' or a list of lists')
                    if len(location) == 1:
                        feature = SeqFeature(FeatureLocation(int(location[0][0])-1,int(location[0][1])),
                                             type=feature_type, strand=int(location[0][2]))
                    elif len(location) > 1:
                        list_of_locations = []
                        for i in location:
                            start = int(i[0]-1)
                            end = int(i[1])
                            strand = int(i[2])
                            list_of_locations.append(FeatureLocation(start,end,strand=strand))
                        feature = SeqFeature(CompoundLocation(list_of_locations),type=feature_type)
                feature.qualifiers['feature_id'] = [feature_id]
                
                
                
                if len(qualifiers.keys()) > 0:
                    for key in qualifiers.keys():
                        feature.qualifiers[key] = [qualifiers[key]]
                if (('codon_start' in qualifiers.keys()) and
                    ('transl_table' in qualifiers.keys())):
                    cds = feature.extract(record.seq)
                    if str(qualifiers['codon_start']) == '2':
                        cds = cds[1:]
                    elif str(qualifiers['codon_start']) == '3':
                        cds = cds[2:]
                    translation = cds.translate(table=int(qualifiers['transl_table']))
                    if len(translation)*3 < float(0.9)*len(cds):
                        raise RuntimeWarning('The translation of feature '+feature_id+' uses less than 90%'+
                                             ' of the coding sequence')
                    feature.qualifiers['translation'] = [str(translation)]
                
                
                feature_seq = feature.extract(record.seq)
                degen = len(feature_seq)
                for i in ['A','T','G','C','U','a','t','g','c','u']:
                    degen -= feature_seq.count(i)
                feature.qualifiers['GC_content'] = [str(GC(feature_seq))]
                feature.qualifiers['nuc_degen_prop'] = [str(float(degen)/len(feature_seq))]
                if 'translation' in feature.qualifiers.keys():
                    transl = feature.qualifiers['translation'][0]
                    degen = 0
                    for i in ['B', 'X', 'Z', 'b', 'x', 'z']:
                        degen += transl.count(i)
                    feature.qualifiers['prot_degen_prop'] = [str(float(degen)/len(transl))]   
                
                record.features.append(feature)
                return feature_id



    ##############################################
    # Project methods for managing concatenations
    ##############################################  


 
    def add_concatenation(self, concatenation_object):
        
        """
        add a Concatenation object to the Project
        
        # making dummy loci
        >>> coi = Locus('dna','CDS','coi',['COX1','cox1'])
        >>> ssu = Locus('dna','rRNA','18S',['18S','SSU'])
        >>> lsu = Locus('dna','rRNA','28S',['28S','LSU'])
        >>> loci = [coi,ssu,lsu]
        
        # making dummy Concatenation
        >>> combined = Concatenation(name='combined', 
        ...                          loci=loci, otu_meta='OTU_dict',
        ...                          otu_must_have_all_of=['coi'],
        ...                          otu_must_have_one_of =[['18S','28S']],
        ...                          define_trimmed_alns=["MafftLinsi@Gappyout"])
        >>> print(str(combined))
        Concatenation named combined, with loci coi,18S,28S,
        of which coi must exist for all species
        and at least one of each group of [ 18S 28S ] is represented.
        Alignments with the following names: MafftLinsi@Gappyout are prefered
        
        # making a dummy Project
        >>> pj = Project(loci, git=False)
        
        # Including the Concatenation in the Project
        >>> pj.add_concatenation(combined)
        >>> print(len(pj.concatenations))
        1
        """
        
        if isinstance(concatenation_object, Concatenation):
            # correct characters offending raxml, back up original values
            meta = concatenation_object.otu_meta
            self.copy_paste_within_feature(meta, "original_%s"%meta)
            offensive = ['\t','\r','\n', "'", '"', ',', ' ',
                         ';', ':', ']','[','(',')','/']
            for r in self.records:
                for f in r.features:
                    if meta in f.qualifiers.keys():
                        for o in offensive:
                            if o in f.qualifiers[meta][0]:
                                print (("found raxml offensive char %s in OTU %s. Replacing with '_ro_'."+
                                        "Backing up original in the qualifier original_%s.")%(o, f.qualifiers[meta][0], meta))
                                f.qualifiers[meta][0] = f.qualifiers[meta][0].replace(o,'_ro_')
            
            seen = []
            for s in self.concatenations:
                seen.append(s.name)
            if concatenation_object.name in seen:
                raise NameError('Concatenation ' + concatenation_object.name +
                                ' apears more than once in self.concatenations')
            else:
                self.concatenations.append(concatenation_object)
        else:
            raise TypeError("Expecting Concatenation object")

    
    #def make_concatenation_alignments(self):
        
    #    """
    #    Concatenates a trimmed alignment based on each of the Concatenation objects and adds them
    #    to the pj.trimmed_alignments dictionary. While a trimmed alignment of an individual locus will have a key
    #    following the patten "locus_name@alignment_method_name@trimming_method_name, the key for a concatenated
    #    trimmed alignment will be the Concatenation object name attribute.
    #    """
    #    for s in self.concatenations:
    #        
    #        # get a non-redundant list of OTUs stored in 'meta', such as voucher specimen
    #        meta = s.otu_meta
    #        OTU_list = []
    #        for record in self.records:
    #            for feature in record.features:
    #                if not feature.type == 'source':
    #                    qualifiers_dictionary = get_qualifiers_dictionary(self,
    #                                                                      feature.qualifiers['feature_id'])
    #                    if (meta in qualifiers_dictionary.keys() and
    #                        not qualifiers_dictionary[meta] in OTU_list):
    #                        OTU_list.append(qualifiers_dictionary[meta])
    #                        
    #        
    #    
    #
    #        included_individuals = {} #included_individuals[otu][locus]=feautre_id
    #        
    #        #Get the correct trimmed alignment tokens
    #        keys_of_trimmed_alignments_to_use_in_concat = []
    #        for locus in s.loci:
    #            trimmed_aln = None
    #            all_locus_trimmed_alns_in_pj = []
    #            for key in self.trimmed_alignments.keys():
    #                if locus.name == key.split('@')[0]:
    #                    all_locus_trimmed_alns_in_pj.append(key)
    #            if len(all_locus_trimmed_alns_in_pj) == 1:
    #                trimmed_aln = all_locus_trimmed_alns_in_pj[0]
    #            elif len(all_locus_trimmed_alns_in_pj) == 0:
    #                raise RuntimeError('Locus '+locus.name+' have no trimmed alignments')
    #            else:
    #                s.define_trimmed_alns.sort(key = lambda i: i.count('@'), reverse=True)
    #                for definition in s.define_trimmed_alns:
    #                    if definition.count('@') == 2 and locus.name == definition.split('@')[0]:
    #                        trimmed_aln = definition
    #                    elif definition.count('@') == 1 and any([definition in i for i in all_locus_trimmed_alns_in_pj]):
    #                        trimmed_aln = locus.name+'@'+definition
    #                    else:
    #                        raise RuntimeError("Could not determine which alignment/trimming alternative to use for locus '"+
    #                                            locus.name+"' out of "+str(all_locus_trimmed_alns_in_pj))
    #            if trimmed_aln:
    #                keys_of_trimmed_alignments_to_use_in_concat.append(trimmed_aln)
    #            else:
    #                raise RuntimeError('Could not find trimmed aln for locus '+locus.name+' given the rulls '+str(s.define_trimmed_alns))
    #        
    #        #print "%i individuals will be included in the concatenations %s"%(len(included_individuals.keys()), s.name)
    #        
    #        #if len(included_individuals.keys()) < 4:
    #        #    raise RuntimeError("Concatenation %s has less than 4 OTUs and cannot be analyzed"%s.name)
    #        for otu in OTU_list:
    #            otu_features = {}
    #            use = True
    #            
    #            # Check first rule
    #            for locus in s.otu_must_have_all_of:
    #                token = [t for t in keys_of_trimmed_alignments_to_use_in_concat if "%s@"%locus in t][0]
    #                feature_ids = [r.id for r in self.trimmed_alignments[token]]
    #                feature_found = False
    #                count = 0
    #                for feature_id in feature_ids:
    #                    qualifiers = get_qualifiers_dictionary(self, feature_id)
    #                    if meta in qualifiers.keys() and otu == qualifiers[meta]:
    #                        count += 1
    #                        feature_found = True
    #                        otu_features[locus] = feature_id
    #                if count > 1:
    #                    raise RuntimeError("%s is not unique in %s"%(otu, locus))
    #                if not feature_found:
    #                    use = False
    #                    
    #            # Check second rule
    #            if use:
    #                for group in s.otu_must_have_one_of:
    #                    if isinstance(group,str):
    #                        raise IOError('The keyword \'otu_must_have_one_of\' has to be a list of lists')
    #                    feature_found = False
    #                    for locus in group:
    #                        token = [t for t in keys_of_trimmed_alignments_to_use_in_concat if "%s@"%locus in t][0]
    #                        feature_ids = [r.id for r in self.trimmed_alignments[token]]
    #                        count = 0
    #                        for feature_id in feature_ids:
    #                            qualifiers = get_qualifiers_dictionary(self, feature_id)
    #                            if meta in qualifiers.keys() and otu == qualifiers[meta]:
    #                                count += 1
    #                                feature_found = True
    #                                otu_features[locus] = feature_id
    #                        if count > 1:
    #                            raise RuntimeError("%s is not unique in %s"%(otu, locus))
    #                    if not feature_found:
    #                        use = False
    #            if use:
    #                included_individuals[otu] = otu_features
    #        
    #        # printing a table of the alignment
    #        included_indivduals_table = ''
    #        loci_names = [l.name for l in s.loci]
    #        line = 'OTU'.ljust(30,' ')
    #        for name in loci_names:
    #            line += name.ljust(20,' ')
    #        included_indivduals_table += line+'\n'
    #        for otu in included_individuals.keys():
    #            line = otu.ljust(30,' ')
    #            for locus_name in loci_names:
    #                if locus_name in included_individuals[otu].keys():
    #                    line += included_individuals[otu][locus_name].ljust(15,' ')
    #                else:
    #                    line += ''.ljust(15,' ')
    #            included_indivduals_table += line+'\n'
    #        print "Concatenation %s will have the following data"%s.name
    #        print included_indivduals_table
    #        
    #        # remove partitions with less then 4 sequences
    #        for name in loci_names:
    #            if len([otu for otu in included_individuals.keys() if name in included_individuals[otu].keys()]) < 4:
    #                print (("Locus %s has less then 4 sequences in concatenation %s and where excluded "+
    #                                     "from the concatenation")%(name,s.name))
    #                for key in keys_of_trimmed_alignments_to_use_in_concat:
    #                    if name in key:
    #                        keys_of_trimmed_alignments_to_use_in_concat.remove(key)
    #                        
    #                                    
    #        
    #        # build alignment
    #        # concat_records = []
    #        alignment = []
    #        for individual in included_individuals.keys():
    #            sequence = ''
    #            for key in keys_of_trimmed_alignments_to_use_in_concat:
    #                locus_name = key.split('@')[0]                    
    #                length = len(self.trimmed_alignments[key][0])
    #                s.used_trimmed_alns[key] = length
    #                if locus_name in included_individuals[individual].keys():
    #                    for record in self.trimmed_alignments[key]:
    #                        if record.id == included_individuals[individual][locus_name]:
    #                            sequence += str(record.seq)
    #                else:
    #                    sequence += '?'*length
    #            concat_sequence = SeqRecord(seq = Seq(sequence), id = individual, description = '')
    #            alignment.append(concat_sequence)
    #        self.trimmed_alignments[s.name] = MultipleSeqAlignment(alignment)                
    #        s.feature_id_dict = included_individuals    

    def make_concatenation_alignments(self):
        
        """
        Concatenates a trimmed alignment based on each of the Concatenation objects and adds them
        to the pj.trimmed_alignments dictionary. While a trimmed alignment of an individual locus will have a key
        following the patten "locus_name@alignment_method_name@trimming_method_name, the key for a concatenated
        trimmed alignment will be the Concatenation object name attribute.
        """
        
        self.__records_list_to_dict__()
        
        for s in self.concatenations:
            
            # get a non-redundant list of OTUs stored in 'meta', such as voucher specimen
            meta = s.otu_meta
            OTU_list = []
            for record in self.records:
                for feature in record.features:
                    if not feature.type == 'source':
                        qualifiers_dictionary = __get_qualifiers_dictionary__(self,
                                                                          feature.qualifiers['feature_id'])
                        if (meta in qualifiers_dictionary.keys() and
                            not qualifiers_dictionary[meta] in OTU_list):
                            OTU_list.append(qualifiers_dictionary[meta])
                            
            
        
    
            included_individuals = {} #included_individuals[otu][locus]=feautre_id
            
            #Get the correct trimmed alignment tokens
            keys_of_trimmed_alignments_to_use_in_concat = []
            for locus in s.loci:
                trimmed_aln = None
                all_locus_trimmed_alns_in_pj = []
                for key in self.trimmed_alignments.keys():
                    if locus.name == key.split('@')[0]:
                        all_locus_trimmed_alns_in_pj.append(key)
                if len(all_locus_trimmed_alns_in_pj) == 1:
                    trimmed_aln = all_locus_trimmed_alns_in_pj[0]
                elif len(all_locus_trimmed_alns_in_pj) == 0:
                    raise RuntimeError('Locus '+locus.name+' have no trimmed alignments')
                else:
                    s.define_trimmed_alns.sort(key = lambda i: i.count('@'), reverse=True)
                    for definition in s.define_trimmed_alns:
                        if definition.count('@') == 2 and locus.name == definition.split('@')[0]:
                            trimmed_aln = definition
                        elif definition.count('@') == 1 and any([definition in i for i in all_locus_trimmed_alns_in_pj]):
                            trimmed_aln = locus.name+'@'+definition
                        else:
                            raise RuntimeError("Could not determine which alignment/trimming alternative to use for locus '"+
                                                locus.name+"' out of "+str(all_locus_trimmed_alns_in_pj))
                if trimmed_aln:
                    keys_of_trimmed_alignments_to_use_in_concat.append(trimmed_aln)
                else:
                    raise RuntimeError('Could not find trimmed aln for locus '+locus.name+' given the rulls '+str(s.define_trimmed_alns))
            
            #print "%i individuals will be included in the concatenations %s"%(len(included_individuals.keys()), s.name)
            
            #if len(included_individuals.keys()) < 4:
            #    raise RuntimeError("Concatenation %s has less than 4 OTUs and cannot be analyzed"%s.name)
            for otu in OTU_list:
                otu_features = {}
                use = True
                
                # Check first rule
                for locus in s.otu_must_have_all_of:
                    token = [t for t in keys_of_trimmed_alignments_to_use_in_concat if "%s@"%locus in t][0]
                    feature_ids = [r.id for r in self.trimmed_alignments[token]]
                    feature_found = False
                    count = 0
                    for feature_id in feature_ids:
                        qualifiers = __get_qualifiers_dictionary__(self, feature_id)
                        if meta in qualifiers.keys() and otu == qualifiers[meta]:
                            count += 1
                            feature_found = True
                            otu_features[locus] = feature_id
                    if count > 1:
                        raise RuntimeError("%s is not unique in %s"%(otu, locus))
                    if not feature_found:
                        use = False
                        
                # Check second rule
                if use:
                    for group in s.otu_must_have_one_of:
                        if isinstance(group,str):
                            raise IOError('The keyword \'otu_must_have_one_of\' has to be a list of lists')
                        feature_found = False
                        for locus in group:
                            token = [t for t in keys_of_trimmed_alignments_to_use_in_concat if "%s@"%locus in t][0]
                            feature_ids = [r.id for r in self.trimmed_alignments[token]]
                            count = 0
                            for feature_id in feature_ids:
                                qualifiers = __get_qualifiers_dictionary__(self, feature_id)
                                if meta in qualifiers.keys() and otu == qualifiers[meta]:
                                    count += 1
                                    feature_found = True
                                    otu_features[locus] = feature_id
                            if count > 1:
                                raise RuntimeError("%s is not unique in %s"%(otu, locus))
                        if not feature_found:
                            use = False
                if use:
                    included_individuals[otu] = otu_features
            
            # printing a table of the alignment
            included_indivduals_table = ''
            loci_names = [l.name for l in s.loci]
            line = 'OTU'.ljust(30,' ')
            for name in loci_names:
                line += name.ljust(20,' ')
            included_indivduals_table += line+'\n'
            for otu in included_individuals.keys():
                line = otu.ljust(30,' ')
                for locus_name in loci_names:
                    if locus_name in included_individuals[otu].keys():
                        line += included_individuals[otu][locus_name].ljust(15,' ')
                    else:
                        line += ''.ljust(15,' ')
                included_indivduals_table += line+'\n'
            print "Concatenation %s will have the following data"%s.name
            print included_indivduals_table
            
            # remove partitions with less then 4 sequences
            for name in loci_names:
                if len([otu for otu in included_individuals.keys() if name in included_individuals[otu].keys()]) < 4:
                    print (("Locus %s has less then 4 sequences in concatenation %s and where excluded "+
                                         "from the concatenation")%(name,s.name))
                    for key in keys_of_trimmed_alignments_to_use_in_concat:
                        if name in key:
                            keys_of_trimmed_alignments_to_use_in_concat.remove(key)
                            
                                        
            
            # build alignment
            # concat_records = []
            alignment = []
            for individual in included_individuals.keys():
                sequence = ''
                for key in keys_of_trimmed_alignments_to_use_in_concat:
                    locus_name = key.split('@')[0]                    
                    length = len(self.trimmed_alignments[key][0])
                    s.used_trimmed_alns[key] = length
                    if locus_name in included_individuals[individual].keys():
                        for record in self.trimmed_alignments[key]:
                            if record.id == included_individuals[individual][locus_name]:
                                sequence += str(record.seq)
                    else:
                        sequence += '?'*length
                concat_sequence = SeqRecord(seq = Seq(sequence), id = individual, description = '')
                alignment.append(concat_sequence)
            self.trimmed_alignments[s.name] = MultipleSeqAlignment(alignment)                
            s.feature_id_dict = included_individuals 
    
        self._records_dict = {}
        
    ###################################################
    # Project methods for modifying feature qualifiers
    ###################################################



    def write(self, filename, format = 'genbank'):
        
        """
        Write the records in the project in any Biopython format, as csv or as nexml.
        
        The csv table has a line for each feature (ie multiplt lines for records
        with multiple non-source featue). Each line will include the records annotations,
        the source feature qualifiers and the qualifiers of the feature itself. (ie, in a
        record with several features, the record annotations and the source feature qualifiers
        will be repeted in several lines, once for each non-source feature in the record.
        The csv file is primarily usefull for reviewing and editing feature qualifiers
        
        The nexml format only includes the trees from the pj.trees dictionary, but all 
        the metadata is included as leaf metadata, including the sequences and the
        trimmed and alligned sequences for each leaf.
        """
        if format == 'nexml':
            self.write_nexml(filename)
        if format == 'phyloxml':
            self.write_phyloxml(filename)
        elif format == 'genbank' or format == 'embl':
            SeqIO.write(self.records, filename, format)
        elif format == 'csv':
            # get titles for source and othe features
            source_qualifiers = []
            feature_qualifiers = []
            for record in self.records:
                for feature in record.features:
                    for key in feature.qualifiers.keys():
                        if feature.type == 'source' and not key in source_qualifiers:
                            source_qualifiers.append(key)
                        elif not feature.type == 'source' and not key in feature_qualifiers:
                            feature_qualifiers.append(key)
            
            with open(filename, 'wb') as csvfile:
                linewriter = csv.writer(csvfile, delimiter='\t',
                                        quotechar='"',
                                        quoting=csv.QUOTE_MINIMAL)
                linewriter.writerow(['record_id','seq']+['source:_'+q for q in source_qualifiers]+['taxonomy']+feature_qualifiers)
                for record in self.records:
                    seq = ''
                    if len(record.seq) <= 10:
                        seq = str(record.seq)[0:11]
                    else:
                        seq = str(record.seq)[0:6] + '...' + str(record.seq)[-5:]
                    
                    
                    line_start = [record.id, seq]
                    source = None
                    for feature in record.features:
                        if feature.type == 'source':
                            source = feature
                    if not source == None:
                        for qual in source_qualifiers:
                            if qual in source.qualifiers.keys():
                                line_start.append(type_to_single_line_str(source.qualifiers[qual]))
                            else:
                                line_start.append('null')
                    elif source == None:
                        for qual in source_qualifiers:
                            line_start.append('null')
                    if 'taxonomy' in record.annotations.keys():
                        line_start.append(type_to_single_line_str(record.annotations['taxonomy']))
                    else:
                        line_start.append(['null'])
                    for feature in record.features:
                        if not feature.type == 'source':
                            line = list(line_start)
                            for qual in feature_qualifiers:
                                if qual in feature.qualifiers.keys() and qual == 'translation':
                                    trans = feature.qualifiers[qual][0]
                                    if len(trans)<=10:
                                        line.append(trans)
                                    else:
                                        line.append(trans[:6] + '...' + trans[-5:])
                                elif qual in feature.qualifiers.keys():
                                    line.append(type_to_single_line_str(feature.qualifiers[qual]))
                                else:
                                    line.append('null')
                            linewriter.writerow(line)
        
        if __builtin__.git:
            import rpgit
            comment = "Records %s text file from %s" % (format, time.asctime())
            out, err = rpgit.gitAdd(filename)
            undate_git_log(self, out, err)
            cwd = os.getcwd()
            import fnmatch
            matches = []
            for root, dirnames, filenames in os.walk(cwd):
                for filename in fnmatch.filter(filenames, '*.py'):
                    matches.append(os.path.join(root, filename))
                for filename in fnmatch.filter(filenames, '*.ipynb'):
                    matches.append(os.path.join(root, filename))
            for match in matches:
                out, err = rpgit.gitAdd(match)
                undate_git_log(self, out, err)
            out, err = rpgit.gitCommit(comment)
            undate_git_log(self, out, err)

    def correct_metadata_from_file(self,csv_file):
        metadata = read_feature_quals_from_tab_csv(csv_file)
        new_records = []
        for record in self.records:
            if record.id in metadata.keys():
                #print 'making new record'
                new_record = record
                record_corrected_metadata = metadata[record.id]
                new_record.annotations['taxonomy'] = metadata[record.id]['taxonomy']
                for feature in new_record.features:
                    if feature.type == 'source':
                        feature.qualifiers = metadata[record.id]['source']
                    else:
                        feature_id = feature.qualifiers['feature_id']
                        translation = None
                        if 'translation' in feature.qualifiers.keys():
                            translation = feature.qualifiers['translation']
                        feature.qualifiers = metadata[record.id]['features'][feature_id[0]]
                        feature.qualifiers['feature_id'] = feature_id
                        if translation:
                            feature.qualifiers['translation'] = translation
                new_records.append(new_record)
            else:
                #print 'using old records'
                new_records.append(record)
        
        self.records = new_records
        
        if __builtin__.git:
            import rpgit
            comment = "Corrected metadata CSV file from %s" % (time.asctime())
            out, err = rpgit.gitAdd(csv_file)
            undate_git_log(self, out, err)
            cwd = os.getcwd()
            import fnmatch
            matches = []
            for root, dirnames, filenames in os.walk(cwd):
                for filename in fnmatch.filter(filenames, '*.py'):
                    matches.append(os.path.join(root, filename))
                for filename in fnmatch.filter(filenames, '*.ipynb'):
                    matches.append(os.path.join(root, filename))
            for match in matches:
                out, err = rpgit.gitAdd(match)
                undate_git_log(self, out, err)
            out, err = rpgit.gitCommit(comment)
            undate_git_log(self, out, err)
                
    def if_this_then_that(self, IF_THIS, IN_THIS, THEN_THAT, IN_THAT, mode = 'whole'):
        
        
        """
        Searches pj.records for features that have the value IF_THIS in the qualifier IN_THIS
        and places the value THEN_THAT in the qualifier IN_THAT, which either exists or is new.
        
        The IF_THIS value can either match completely (mode = 'whole') or just to a part (mode = 'part')
        of the target qualifier value
        
        The following demonstartes all the feature qualifier editing methods
        
        # Make a dummy pj with a locus and with records
        >>> input_filenames = ['test-data/test.gb']
        >>> locus = Locus('dna', 'CDS', 'coi', ['cox1','COX1','coi','COI','CoI'])
        >>> pj = Project([locus], git=False)
        >>> pj.read_embl_genbank(input_filenames)
        
        # copying a source qualifier into the feature qualifiers so that it
        # will be available for editing (source qualifiers are kept imutable)
        >>> pj.add_qualifier_from_source('organism')
        
        # populate a new qualifier based on the data in another
        # Here we will take only the the genus name from the organism qualifier
        # and put it in a new qualifier
        # We use mode='part' because our search phrase (the genus name)
        # fits only the start of the organism name
        >>> tetillid_genera = ['Tetilla', 'Cinachyrella', 'Craniella']
        >>> for g in tetillid_genera:
        ...     pj.if_this_then_that(g, 'organism', g, 'genus', mode='part')
        
        # Now we will add a sample id to all the sequences which belong to
        # sample TAU_25617
        >>> pj.add_qualifier(['JX177913.1_f0',
        ...                  'JX177935.1_f0',
        ...                  'JX177965.1_f0'],
        ...                  'specimen_voucher',
        ...                  'TAU_25617')
    
        # We are using a copy paste approch to unite the data from 
        # differen qualifiers under on qualifiers
        >>> pj.copy_paste_within_feature('variant', 'strain_or_variant')
        >>> pj.copy_paste_within_feature('strain', 'strain_or_variant')
        
        # Now we print the qualifier of a random feature as an example
        >>> qual_dict = get_qualifiers_dictionary(pj, 'JX177913.1_f0')
        >>> qual_items = qual_dict.items()
        >>> qual_items.sort(key = lambda i: i[0])
        >>> for key, val in qual_items:
        ...     print(key.ljust(40,' ') + type_to_single_line_str(val)[:5]+'...')
        GC_content                              37.28...
        annotation_accessions                   JX177...
        annotation_data_file_division           INV...
        annotation_date                         05-SE...
        annotation_gi                           39293...
        annotation_keywords                     ...
        annotation_organism                     Cinac...
        annotation_references                   locat...
        annotation_sequence_version             1...
        annotation_source                       mitoc...
        annotation_taxonomy                     Eukar...
        codon_start                             2...
        db_xref                                 GI:39...
        feature_id                              JX177...
        gene                                    cox1...
        genus                                   Cinac...
        nuc_degen_prop                          0.0...
        organism                                Cinac...
        product                                 cytoc...
        prot_degen_prop                         0.0...
        protein_id                              AFM91...
        record_id                               JX177...
        source_country                          Panam...
        source_db_xref                          taxon...
        source_feature_id                       JX177...
        source_identified_by                    Ilan,...
        source_mol_type                         genom...
        source_note                             autho...
        source_organelle                        mitoc...
        source_organism                         Cinac...
        source_specimen_voucher                 DH_S2...
        specimen_voucher                        TAU_2...
        transl_table                            4...
        translation                             MIGSG...
        
        # Note that GC content and the porportion of degenerate positions
        # have been automatically included. They will be plotted in the report
        """
        
        for record in self.records:
            for feature in record.features:
                if not feature.type == 'source':
                    if IN_THIS in feature.qualifiers.keys():
                        if not type(feature.qualifiers[IN_THIS]) is list:
                            feature.qualifiers[IN_THIS] = [feature.qualifiers[IN_THIS]]
                        for i in feature.qualifiers[IN_THIS]:
                            if mode == 'whole':
                                if i == IF_THIS:
                                    feature.qualifiers[IN_THAT] = [THEN_THAT]
                            elif mode == 'part':
                                if IF_THIS in i:
                                    feature.qualifiers[IN_THAT] = [THEN_THAT]
    


    def add_qualifier(self, feature_ids, name, value):
        if not type(value) is list:
                    value = [value]
        for record in self.records:
            for feature in record.features:
                if feature.qualifiers['feature_id'][0] in feature_ids:
                    feature.qualifiers[name] = value



    def add_qualifier_from_source(self, qualifier):
        for record in self.records:
            source = None
            for feature in record.features:
                if feature.type == 'source':
                    source = feature
            value = None
            if not source == None:
                if qualifier in source.qualifiers.keys():
                    value = source.qualifiers[qualifier]       
            if not value == None:
                if not type(value) is list:
                    value = [value]
                for feature in record.features:
                    if not feature.type == 'source':
                        feature.qualifiers[qualifier] = value
    


    def copy_paste_within_feature(self, from_qualifier, to_qualifier):
        for record in self.records:
            for feature in record.features:
                if not feature.type == 'source':
                    if from_qualifier in feature.qualifiers.keys():
                        feature.qualifiers[to_qualifier] = feature.qualifiers[from_qualifier]
                  
                        
                        
    def copy_paste_from_features_to_source(self, from_feature_qual, to_source_qual):
        for record in self.records:
            source = None
            values_from_features = []
            for feature in record.features:
                if not feature.type == 'source':
                    if from_feature_qual in feature.qualifiers.keys():
                        if not feature.qualifiers[from_feature_qual] in values_from_features:
                            values_from_features += feature.qualifiers[from_feature_qual]
                else:
                    source = feature
            if source == None:
                source = SeqFeature(FeatureLocation(0, len(record.seq)), type='source', strand=1)
                source.qualifiers['feature_id'] = record.id + '_source'
                record.features = [source] + record.features
            if not to_source_qual in source.qualifiers.keys():
                source.qualifiers[to_source_qual] = values_from_features
           
                
                
    def species_vs_loci(self, outfile_name):
        
        """
        Makes a csv file showing the count of each unique value in the source_organism qualifier
        for each locus
        """
        
        species_vs_loci = {}
        for record in self.records:
            organism = 'undef'
            for feature in record.features:
                if feature.type == 'source':
                    if 'organism' in feature.qualifiers.keys():
                        organism = feature.qualifiers['organism'][0]
            if not organism in species_vs_loci.keys():
                species_vs_loci[organism] = {}    
            for feature in record.features:
                if not feature.type == 'source':
                    for locus in self.loci:
                        if not locus.name in locus.aliases:
                            locus.aliases.append(locus.name)
                        if 'gene' in feature.qualifiers.keys():
                            if feature.qualifiers['gene'][0] in locus.aliases:
                                if not locus.name in species_vs_loci[organism].keys():
                                    species_vs_loci[organism][locus.name] = 1
                                else:
                                    species_vs_loci[organism][locus.name] += 1
                        elif 'product' in feature.qualifiers.keys():
                            if feature.qualifiers['product'][0] in locus.aliases:
                                if not locus.name in species_vs_loci[organism].keys():
                                    species_vs_loci[organism][locus.name] = 1
                                else:
                                    species_vs_loci[organism][locus.name] += 1
        with open(outfile_name, 'wb') as csvfile:
            linewriter = csv.writer(csvfile, delimiter='\t',
                                    quotechar='|',
                                    quoting=csv.QUOTE_MINIMAL)
            loci_names = []
            for g in self.loci:
                loci_names.append(g.name)
            linewriter.writerow(['species']+loci_names)
            for organism in sorted(list(species_vs_loci.keys())):
                line = [organism]
                for name in loci_names:
                    if name in species_vs_loci[organism].keys():
                        line.append(str(species_vs_loci[organism][name]))
                    else:
                        line.append('0')
                linewriter.writerow(line)
                
                
                
    ######################################                
    # Project methods to analyze the data
    ######################################
    
    
        
    def extract_by_locus(self):
        
        """
        
        >>> input_filenames = ['test-data/test.gb']
        >>> coi = Locus('dna', 'CDS', 'coi', ['cox1','COX1','coi','COI','CoI'])
        >>> lsu = Locus('dna', 'rRNA', '28S', ['28s','28S','LSU rRNA','28S ribosomal RNA','28S large subunit ribosomal RNA'])
        >>> pj = Project([coi, lsu], git=False)
        >>> pj.read_embl_genbank(input_filenames)
        >>> pj.extract_by_locus()
        >>> print(len(pj.records_by_locus['coi']))
        89
        >>> print(len(pj.records_by_locus['28S']))
        48
        """
        
        data_by_locus = {}
        for locus in self.loci:
            if not locus.name in locus.aliases:
                locus.aliases.append(locus.name)
                locus.aliases.append(locus.name.replace('_',' '))
            records = []
            for record in self.records:
                for feature in record.features:
                    if ((feature.type == locus.feature_type or
                         # equate CDS and Protein feature types to allow reading protein sequences from 
                         # a mix of DNA gb files and protein fasta files that were given a Protein feature
                         # type
                         (feature.type in ['CDS','Protein'] and locus.feature_type in ['CDS','Protein'])) and
                        
                        (('gene' in feature.qualifiers.keys() and
                          feature.qualifiers['gene'][0] in locus.aliases) 
                         or
                         ('product' in feature.qualifiers.keys() and 
                          feature.qualifiers['product'][0] in locus.aliases))
                        
                        ):
                        if locus.char_type == 'dna':
                            if not type(record.seq.alphabet) == IUPAC.IUPACProtein:
                                S = feature.extract(record.seq)
                            else:
                                raise RuntimeError('Trying to read DNA from protein only records')
                        elif locus.char_type == 'prot':
                            if not type(record.seq.alphabet) == IUPAC.IUPACProtein:
                                if 'translation' in feature.qualifiers.keys():
                                    S = Seq(feature.qualifiers['translation'][0], IUPAC.protein)
                                else:
                                    raise RuntimeError('Trying to read protein from DNA records without translation info')
                            else:
                                S = feature.extract(record.seq)
                        feature_record = SeqRecord(seq = S, id = feature.qualifiers['feature_id'][0],
                                                   description = '')
                        records.append(feature_record)
            data_by_locus[locus.name] = records
        self.records_by_locus = data_by_locus

    def exclude(self, start_from_max=True, **kwargs):
        keep_safe = self.records_by_locus
        self.extract_by_locus()
        locus_names = [i.name for i in self.loci]
        for key, value in kwargs.iteritems():
            if key in locus_names:
                if value == 'all':
                    self.records_by_locus[key] = []
                else:
                    subset = []
                    locus_feature_ids = [i.id.split('_')[0] for i in self.records_by_locus[key]]
                    if not all(i.split('_')[0] in locus_feature_ids for i in value):
                        print [i.split('_')[0] for i in value if not i.split('_')[0] in locus_feature_ids]
                        warnings.warn('Not all records to exclude exist in locus. Typos?')
                    if not start_from_max and not keep_safe == {}:
                        for record in keep_safe[key]:
                            if not record.id.split('_')[0] in [i.split('_')[0] for i in value]:
                                subset.append(record)
                        self.records_by_locus[key] = subset
                    else:
                        for record in self.records_by_locus[key]:
                            if not record.id.split('_')[0] in [i.split('_')[0] for i in value]:
                                subset.append(record)
                        self.records_by_locus[key] = subset
            else:
                warnings.warn('Locus name %s not recognised'%key)

    def include(self, start_from_null=True, **kwargs):
        keep_safe = self.records_by_locus
        self.extract_by_locus()
        locus_names = [i.name for i in self.loci]
        for key, value in kwargs.iteritems():
            if key in locus_names:
                if value == 'all':
                    pass
                else:
                    subset = []
                    locus_feature_ids = [i.id.split('_')[0] for i in self.records_by_locus[key]]
                    if not all(i.split('_')[0] in locus_feature_ids for i in value):
                        print [i.split('_')[0] for i in value if not i.split('_')[0] in locus_feature_ids]
                        warnings.warn('Not all records to include exist in locus. Typos?')
                    for record in self.records_by_locus[key]:
                        if record.id.split('_')[0] in [i.split('_')[0] for i in value]:
                            subset.append(record)
                    self.records_by_locus[key] = subset
                    if not start_from_null and not keep_safe == {}:
                        self.records_by_locus[key] = subset+keep_safe[key]
            else:
                warnings.warn('Locus name %s not recognised'%key)

    def filter_by_seq_length(self, locus_name, min_length=0, max_length=None):
        if self.records_by_locus == {}:
            self.extract_by_locus()
        subset = [r for r in self.records_by_locus[locus_name] if len(r) >= min_length]
        if max_length:
            subset = [r for r in subset if len(r) <= max_length]
        self.records_by_locus[locus_name] = subset
        
    def filter_by_gc_content(self, locus_name, min_percent_gc=0, max_percent_gc=None):
        if self.records_by_locus == {}:
            self.extract_by_locus()
        subset = [r for r in self.records_by_locus[locus_name] if GC(r.seq) >= min_percent_gc]
        if max_percent_gc:
            subset = [r for r in subset if GC(r.seq) <= max_percent_gc]
        self.records_by_locus[locus_name] = subset
                

    def write_by_locus(self, format = 'fasta'):
        
        """
        Write the unaligned sequences into file in any Biopython format, one file per locus
        """
        
        if self.records_by_locus == {}:
            self.extract_by_locus
        for key in self.records_by_locus.keys():
            SeqIO.write(self.records_by_locus[key], key+'.'+format, format)



    def align(self, alignment_methods=[], pal2nal='defaults'):
        
            """
            Configured by an AlnConf object
            """
            if pal2nal == 'defaults':
                pal2nal = self.defaults['pal2nal']
            seen_loci = []
            for method in alignment_methods:
                method.timeit.append(time.time())
                method.platform = platform_report()
                if method.program_name == 'muscle':
                    method.platform.append('Program and version: '+
                                           os.popen(method.cmd + ' -version').read().split(' by')[0])
                    method.platform.append('Program reference:\nEdgar 2004: MUSCLE: multiple sequence '+
                                           'alignment with high accuracy and high throughput. '+
                                           'Nucleic Acids Research 32(5):1792-1797')
                elif method.program_name == 'mafft':
                    p = sub.Popen(method.cmd+" --version", shell=True, stderr=sub.PIPE, stdout=sub.PIPE)
                    method.platform.append('Program and version: '+
                                           p.communicate()[1].splitlines()[3].strip().split(' (')[0])
                    method.platform.append('Program reference:Katoh, Standley 2013 (Molecular Biology and Evolution 30:772-780) '+
                                           'MAFFT multiple sequence alignment software version 7: improvements in performance and usability.')
                
                
                
                # PROGRAM PLUG
                # NOTE: THIS ADDS THE PROGRAM AND REFERENCE INFO TO THE CONF OBJECT
                #
                # elif method.program_name == 'the name of the program':
                #    p = sub.Popen(method.cmd+" command that writes version", shell=True, stderr=sub.PIPE, stdout=sub.PIPE)
                #    method.platform.append('Program and version: '+
                #                           # 1 for stderr, 0 for stdout
                #                           p.communicate()[1].splitlines()[# get the line and split])
                #    method.platform.append('Program reference: write the reference here')
                
                
                if  method.CDSAlign:
                    p = sub.Popen('pal2nal.pl', shell=True, stderr=sub.PIPE, stdout=sub.PIPE)
                    method.platform[-2] += '\nPal2Nal '+ p.communicate()[1].splitlines()[1].split('(')[1].split(')')[0]
                    method.platform[-1] += ['\nMikita Suyama, David Torrents, and Peer Bork (2006) PAL2NAL: robust '+ 
                                            'conversion of protein sequence alignments into the corresponding codon '+
                                            'alignments.Nucleic Acids Res. 34, W609-W612.'][0]
                for locus in method.loci:
                    if locus.name in seen_loci:
                        # Uncomment this to restrict each locus to one Conf object
                        #raise RuntimeError('locus '+locus.name+' is in more than one AlnConf objects')
                        pass
                    else:
                        seen_loci.append(locus.name)
                        
                    # This will put a fasta alignment in stdout    
                    stdout, stderr = method.command_lines[locus.name]()
                    
                    # This will make a MultipleSeqAlignment from the fasta alignment
                    align = AlignIO.read(StringIO(stdout), "fasta",  alphabet=IUPAC.protein)
                    
                    # This will run pal2nal
                    if method.CDSAlign and locus.feature_type == 'CDS' and locus.char_type == 'dna':
                        
                        # These will check the protein alignment have the same otus and 
                        # the sequences are the same length
                        for seq in align:
                            found = 0
                            for s in method.CDS_in_frame[locus.name]:
                                if s.id == seq.id:
                                    found = 1
                            if found == 0:
                                raise RuntimeError(seq.id + ' is not in the CDS sequences')
                        for s in method.CDS_in_frame[locus.name]:
                            found = 0
                            for seq in align:
                                if s.id == seq.id:
                                    found = 1
                            if found == 0:
                                raise RuntimeError(seq.id + ' is not in the protein sequences')
                        for seq in method.CDS_in_frame[locus.name]:    
                            for prot in align:
                                if prot.id == seq.id:
                                    i = 0
                                    for p in str(prot.seq):
                                        if not p == '-':
                                            i += 1
                                    if not i*3 == len(seq.seq):
                                        raise RuntimeError('nuc and prot seqs have unmatched lengths for '+seq.id)
                                        
                        # This will write an input protein alignment file for pal2nal
                        aln_filename = method.id+'_'+locus.name+'.aln'
                        AlignIO.write(align, aln_filename, 'fasta')
                        
                        # This will run pal2nal
                        cds_filename = method.id+'_CDS_in_frame_'+locus.name+'.fasta'
                        stdout = os.popen(pal2nal+' '+aln_filename+' '+cds_filename + ' -nostderr -codontable %i'%method.codontable).read()
                        
                        # This will make a MultipleSeqAlignment out of the pal2nal output
                        align = AlignIO.read(StringIO(stdout), "clustal",  alphabet=IUPAC.ambiguous_dna)
                        
                        # Maybe this will replace pal2nal
                        #from Bio import CodonAlign
                        #codon_aln = CodonAlign.build(align, method.CDS_in_frame[locus.name])
                        #align = codon_aln
                        
                    # This will list all the temp files written during the analysis
                    method_files = glob.glob(method.id+'_*')
                    
                    # This will get summary statistics of the alignment,
                    # will remove alignments with less than four unique sequences
                    # and will remove undetermined sequences.
                    [summary_lines, num_lines, num_undeter, num_collapsed_aln_seqs] = aln_summary(align)
                    summary = 'Alignment name: '+locus.name+'@'+method.method_name+'\n'
                    for line in summary_lines:
                        summary += line+'\n'
                    if num_lines < 4:
                        line = 'Alignment %s has less than 4 sequences and will be dropped'%(locus.name+'@'+method.method_name)
                        print line
                        summary += line+'\n'
                    elif num_undeter[0] > 0:
                        line = 'Alignment %s has undetermined sequences and will be dropped'%(locus.name+'@'+method.method_name)
                        print line
                        summary += line+'\n'
                    elif num_collapsed_aln_seqs < 4:
                        line = 'Alignment %s has less than 4 unique sequences and will be dropped'%(locus.name+'@'+method.method_name)
                        print line
                        summary += line+'\n'
                        
                    else:
                        # This will palce the MultipleSeqAlignment in the Project
                        self.alignments[locus.name+'@'+method.method_name] = align
                    self.aln_summaries.append(summary)
                    
                # This will measure the running time    
                method.timeit.append(time.time())
                method.timeit.append(method.timeit[2]-method.timeit[1])
                
                # This will delete the temp files
                for f in method_files:
                    os.remove(f)
                self.used_methods[method.method_name] = method


            if self.pickle_name:
                pickle_pj(self, self.pickle_name, track=False)
        
            if __builtin__.git and self.pickle_name:

                comment = ''
                for method in alignment_methods:
                    comment += '%s\n'%(str(method))

                import rpgit
                out, err = rpgit.gitAdd(self.pickle_name)
                undate_git_log(self, out, err)
                cwd = os.getcwd()
                import fnmatch
                matches = []
                for root, dirnames, filenames in os.walk(cwd):
                    for filename in fnmatch.filter(filenames, '*.py'):
                        matches.append(os.path.join(root, filename))
                    for filename in fnmatch.filter(filenames, '*.ipynb'):
                        matches.append(os.path.join(root, filename))
                for match in matches:
                    out, err = rpgit.gitAdd(match)
                    undate_git_log(self, out, err)
                out, err = rpgit.gitCommit(comment)
                undate_git_log(self, out, err)
                
                

    def write_alns(self, id=['feature_id'], format = 'fasta'):
        filenames = []
        if len(self.alignments.keys()) == 0:
            raise IOError('Align the records first')
        else:
            for key in self.alignments:
                aln = self.alignments[key]
                records = []
                for s in aln:
                    qualifiers = get_qualifiers_dictionary(self, s.id)
                    new_id = ""
                    for i in id:
                        if i in qualifiers.keys():
                            new_id += qualifiers[i]+'_'
                    if new_id == "":
                        new_id = s.id
                    else:
                        new_id = new_id[:-1]
                    records.append(SeqRecord(seq=s.seq, id=new_id))
                SeqIO.write(records, key+'_aln.'+format, format)
                filenames.append(key+'_aln.'+format)
        return filenames



    def write_trimmed_alns(self, id=['feature_id'], format = 'fasta'):
        filenames = []
        if len(self.trimmed_alignments.keys()) == 0:
            raise IOError('Align and trimmed the records first')
        else:
            for key in self.trimmed_alignments.keys():
                aln = self.trimmed_alignments[key]
                records = []
                for s in aln:
                    qualifiers = get_qualifiers_dictionary(self, s.id)
                    new_id = ""
                    for i in id:
                        if i in qualifiers.keys():
                            new_id += qualifiers[i]+'_'
                    if new_id == "":
                        new_id = s.id
                    else:
                        new_id = new_id[:-1]
                    records.append(SeqRecord(seq=s.seq, id=new_id))
                SeqIO.write(records, key+'_trimmed_aln.'+format, format)
                filenames.append(key+'_trimmed_aln.'+format)
        return filenames

    def show_aln(self, token, id=['feature_id']):
        aln_obj=None
        if token in self.alignments.keys():
            aln_obj = self.alignments[token]
        elif token in self.trimmed_alignments.keys():
            aln_obj = self.trimmed_alignments[token]
        locus_name = token.split('@')[0]
        char_type = [locus.char_type for locus in self.loci if locus.name == locus_name][0]
        
        records = []
        for s in aln_obj:
            qualifiers = get_qualifiers_dictionary(self, s.id)
            new_id = ""
            for i in id:
                if i in qualifiers.keys():
                    new_id += qualifiers[i]+'_'
            if new_id == "":
                new_id = s.id
            else:
                new_id = new_id[:-1]
            records.append(SeqRecord(seq=s.seq, id=new_id))
            
        title_length = max([len(r.id) for r in records])+2
        
        # colors
        dna_colors = {'a':'green',
                      'A':'green',
                      'T':'red',
                      't':'red',
                      'U':'red',
                      'u':'red',
                      'g':'gray',
                      'G':'gray',
                      'c':'blue',
                      'C':'blue'
                      }
        protein_colors = {'R':'blueviolet',
                          'r':'blueviolet',
                          
                          'K':'cornflowerblue',
                          'k':'cornflowerblue',
                          
                          'E':'red',
                          'e':'red',
                          
                          'D':'crimson',
                          'd':'crimson',
                          
                          'I':'gold',
                          'i':'gold',
                          
                          'L':'yellow',
                          'l':'yellow',
                          
                          'V':'moccasin',
                          'v':'moccasin',
                          
                          'A':'lemonchiffon',
                          'a':'lemonchiffon',
                          
                          'C':'palegreen',
                          'c':'palegreen',
                          
                          'H':'paleturquoise',
                          'h':'paleturquoise',
                          
                          'M':'hotpink',
                          'm':'hotpink',
                          
                          'N':'pink',
                          'n':'pink',
                          
                          'Q':'yellow',
                          'q':'yellow',
                          
                          'F':'darkseagreen',
                          'f':'darkseagreen',
                          
                          'Y':'darkcyan',
                          'y':'darkcyan',
                          
                          'W':'steelblue',
                          'w':'steelblue',
                          
                          
                          'S':'thistle',
                          's':'thistle',
                          
                          'T':'lavender',
                          't':'lavender',
                          
                          'G':'darkgray',
                          'g':'darkgray',
                          
                          'P':'gainsboro',
                          'p':'gainsboro',
                          }
        colors = None
        if char_type == 'dna':
            colors = dna_colors
        elif char_type == 'prot':
            colors = protein_colors
        linelength = (len(records[0].seq)+len(records[0].id)+3)*10
        html_string = '<html><head></head>\n'
        html_string += '<body><pre><font face="Courier New">\n'
        for r in records:
            html_string += r.id.ljust(title_length, '.')
            for p in str(r.seq):
                c = 'white'
                if p in colors.keys():
                    c = colors[p]
                html_string += '<font style="BACKGROUND-COLOR: %s">%s</font>'%(c, p)
            html_string += "<br>"
        html_string +=  '</font></pre></body></html>' 
        
        import webbrowser
        path = os.path.abspath("%s.html"%token)
        url = 'file://' + path
        with open(path, 'w') as f:
            f.write(html_string)
        webbrowser.open_new_tab(url)
        
    
    def tree(self, raxml_methods, bpcomp='default', bpcomp_burnin=0.2, bpcomp_step=10, bpcomp_cutoff=0.01):
        
        if bpcomp == 'default':
            bpcomp = self.defaults['bpcomp']
        
        for raxml_method in raxml_methods:
            raxml_method.timeit.append(time.time())
            raxml_method.platform = platform_report()
            if isinstance(raxml_method, RaxmlConf):
                raxml_method.platform.append('Program and version: '+ raxml_method.cmd +
                                             os.popen(raxml_method.cmd + ' -version').readlines()[2].split('This is ')[1].split(' released')[0])
                raxml_method.platform.append('Program reference: '+
                                             'A. Stamatakis: RAxML Version 8: A tool for Phylogenetic Analysis '+
                                             'and Post-Analysis of Large Phylogenies. In Bioinformatics, 2014.')
            elif isinstance(raxml_method, PbConf):
                p = sub.Popen(raxml_method.cmd+" -v", shell=True, stderr=sub.PIPE, stdout=sub.PIPE)
                raxml_method.platform.append('Program and version: '+p.communicate()[1].splitlines()[-1])
                raxml_method.platform.append('Program reference: '
                                             'N. Lartillot, T. Lepage and S. Blanquart, 2009: PhyloBayes 3: '+
                                             'a Bayesian software package for phylogenetic reconstruction and '+
                                             'molecular dating. Bioinformatics Vol. 25 no. 17.')
            # PROGRAM PLUG
            # NOTE: THIS METHOD SERVES ALL PHYLO PROGRAMS ALTHOUGH THE ITERATOR IS 
            # CALLED raxml_method
            # THIS GETS THE VERSION AND REFERENCE OF THE PROGRAM
            
            # elif isinstance(raxml_method, Conf object name):
            #    p = sub.Popen(raxml_method.cmd+" command that writes version", shell=True, stderr=sub.PIPE, stdout=sub.PIPE)
            #    raxml_method.platform.append('Program and version: '+
            #                           # 1 for stderr, 0 for stdout
            #                           p.communicate()[1].splitlines()[# get the line and split])
            #    raxml_method.platform.append('Program reference: write the reference here')
                
            for trimmed_alignment in raxml_method.command_lines.keys():
                for cline in raxml_method.command_lines[trimmed_alignment]:
                    if isinstance(raxml_method, RaxmlConf):
                        stdout, stderr = cline()
                    elif isinstance(raxml_method, PbConf):
                        sub.call(cline, shell=True)
                        
                    # PROGRAM PLUG
                    # THIS RUNS THE PROGRAM
                    # elif isinstance(raxml_method, Conf Object Name):
                    #     sub.call(cline, shell=True)
                    
                    
                t = None
                if isinstance(raxml_method, RaxmlConf):
                    if raxml_method.preset == 'fa':
                        t = Tree('RAxML_bipartitions.'+raxml_method.id+'_'+trimmed_alignment+'0')
                    elif raxml_method.preset == 'fD_fb':
                        t = Tree('RAxML_bipartitions.'+raxml_method.id+'_'+trimmed_alignment+'1')
                    elif raxml_method.preset == 'fd_b_fb':
                        t = Tree('RAxML_bipartitions.'+raxml_method.id+'_'+trimmed_alignment+'2')
                    elif raxml_method.preset == 'fd_fJ' or raxml_method.preset == 'fF_fJ':
                        tree_filename = 'RAxML_fastTreeSH_Support.'+raxml_method.id+'_'+trimmed_alignment+'1'
                        t = Tree(open(tree_filename,'r').read().replace('[','[&&NHX:support='))
                        
                elif isinstance(raxml_method, PbConf):
                    base_name = "%s_%s"%(raxml_method.id, trimmed_alignment)
                    trees_file = "%s.1.treelist"%base_name
                    trees_per_chain = len(open(trees_file,'r').readlines())
                    # find the number of chains
                    nchains = raxml_method.cline_args['nchain'].split()[0]
                    chain_names = ''
                    for i in range(1,int(nchains)+1):
                        chain_names += "%s.%i "%(base_name, i)
                    chain_names = chain_names[:-1]
                    bpcomp_cline = "%s -c %f -x %i %i %s"%(bpcomp,
                                                           bpcomp_cutoff,
                                                           int(trees_per_chain*bpcomp_burnin),
                                                           int(bpcomp_step),
                                                           chain_names)
                    sub.call(bpcomp_cline, shell=True)
                    t = Tree("bpcomp.con.tre")
                    for l in t:
                        l.support = 0
                        
                # PROGRAM PLUG
                # NOTE: THIS IS SIMPLIFIED. MIGHT WORK WITH SOMETHING LIKE
                # FASTTREE. SEE MORE EXAMPLES ABOVE
                # THIS SECTION MAKES A Tree OBJECT OUT OF THE OUTPUT FILE
                
                # elif isinstance(raxml_method, Conf object name):
                #    base_name = "%s_%s"%(raxml_method.id, trimmed_alignment)
                #    tree_file = "the form of the output file with the %s"%base_name
                #    t = Tree(tree_file)
                    
            
                # This puts the Conf id in the tree 
                for n in t.traverse():
                    n.add_feature('tree_method_id', str(raxml_method.id)+'_'+trimmed_alignment)
                t.dist = 0
                t.add_feature('tree_method_id', str(raxml_method.id)+'_'+trimmed_alignment)
                
                
                # This gets all the metadta from pj.records and puts it
                # on the tree leaves
                loci_names = [i.name for i in  self.loci]       
                concat_names = [c.name for c in self.concatenations]
                if trimmed_alignment.partition('@')[0] in loci_names:
                        
                        for leaf in t:
                            records = self.records
                            feature = ''
                            feature_source = ''
                            record = ''
                            for r in records:
                                if r.id in leaf.name:
                                    record = r
                                    for f in r.features:
                                        if f.type == 'source':
                                            feature_source = f
                                        elif f.qualifiers['feature_id'][0] == leaf.name:
                                            feature = f
                            for a in record.annotations.keys():
                                label = 'annotation_'+a
                                leaf.add_feature(label, record.annotations[a])
                            for f_source_qual in feature_source.qualifiers.keys():
                                label = 'source_'+f_source_qual
                                leaf.add_feature(label, feature_source.qualifiers[f_source_qual][0])
                            for f_qual in feature.qualifiers.keys():
                                leaf.add_feature(f_qual, feature.qualifiers[f_qual][0])
                        for l in t:
                            t.add_feature('tree_id', trimmed_alignment+'@'+raxml_method.method_name)
                            
                        # This puts a Tree object and a string representation of the tree
                        # in the project
                        self.trees[trimmed_alignment+'@'+raxml_method.method_name] = [t,t.write(features=[])]
                        
                elif trimmed_alignment in concat_names:
                        # This does the same as above, but instead of dealing with gene trees
                        # it deals with supermatrix trees
                        s = filter(lambda i: i.name == trimmed_alignment, self.concatenations)[0]
                        for leaf in t:
                            records = self.records
                            feature = ''
                            feature_source = ''
                            record = ''
                            for r in records:
                                for feature in r.features:
                                    if not feature.type == 'source':
                                        qual_dict = get_qualifiers_dictionary(self, feature.qualifiers['feature_id'])
                                        if s.otu_meta in qual_dict.keys() and qual_dict[s.otu_meta] == leaf.name:
                                            for key in qual_dict.keys():
                                                leaf.add_feature(key, qual_dict[key])
                        for l in t:
                            t.add_feature('tree_id', s.name+'@mixed@mixed@'+raxml_method.method_name)
                        self.trees[s.name+'@mixed@mixed@'+raxml_method.method_name] = [t,t.write(features=[])]
                        
            # This times the analysis            
            raxml_method.timeit.append(time.time())
            raxml_method.timeit.append(raxml_method.timeit[2]-raxml_method.timeit[1])
            
            # This deletes temp files
            if not raxml_method.keepfiles:
                for file_name in os.listdir(os.curdir):
                            if raxml_method.id.partition('_')[0] in file_name:
                                os.remove(file_name)
            self.used_methods[raxml_method.method_name] = raxml_method

        if self.pickle_name:
            pickle_pj(self, self.pickle_name, track=False)
        
        if __builtin__.git and self.pickle_name:
            
            comment = ''
            for raxml_method in raxml_methods:
                comment += '%s\n'%(str(raxml_method))
            
            import rpgit
            out, err = rpgit.gitAdd(self.pickle_name)
            undate_git_log(self, out, err)
            cwd = os.getcwd()
            import fnmatch
            matches = []
            for root, dirnames, filenames in os.walk(cwd):
                for filename in fnmatch.filter(filenames, '*.py'):
                    matches.append(os.path.join(root, filename))
                for filename in fnmatch.filter(filenames, '*.ipynb'):
                    matches.append(os.path.join(root, filename))
            for match in matches:
                out, err = rpgit.gitAdd(match)
                undate_git_log(self, out, err)
            out, err = rpgit.gitCommit(comment)
            undate_git_log(self, out, err)

    def clear_tree_annotations(self):
        for tree in self.trees.keys():
            t = Tree(self.trees[tree][1])
            t.dist = 0
            self.trees[tree][0] = t



    def write_nexml(self, output_name):
        D = dendropy.DataSet()
        tree_list = []
        
        loci_names = []
        for locus in self.loci:
            loci_names.append(locus.name)
        
        for tree_name in self.trees.keys():
            #get aligned and trimmd aligned sequences as leaf features
            t = self.trees[tree_name][0]
            for l in t:
                loc_name = tree_name.split('@')[0]
                trim_aln_name = tree_name.rpartition('@')[0]
                aln_name = None
                if loc_name in loci_names:
                    aln_name = tree_name.rsplit('@',2)[0]
                else:
                    trim_aln_name = trim_aln_name.split('@')[0]
                
                otu_feature = 'feature_id'
                if not aln_name: # then it is a concatenation
                    for c in self.concatenations:
                        if c.name == loc_name:
                            otu_feature = c.otu_meta
                            
                if aln_name: # Then it is a locus
                    leaf_feature_value = getattr(l, otu_feature)
                    alignment = self.alignments[aln_name]
                    for record in alignment:
                        if record.id == leaf_feature_value:
                            l.add_feature('aligned_sequence',str(record.seq))
                t_aln = self.trimmed_alignments[trim_aln_name]
                    
                leaf_feature_value = getattr(l, otu_feature)
                for record in t_aln:
                    if record.id == leaf_feature_value:
                        l.add_feature('aligned_trimmed_sequence',str(record.seq))
                    
            tree_string = self.trees[tree_name][0].write(features=[])
            tree = dendropy.Tree()
            tree.read_from_string(tree_string, schema='newick', extract_comment_metadata = True)
            tree_list.append(tree)
        TL = dendropy.TreeList(tree_list)    
        D.add_tree_list(TL)
            
        D.write_to_path(
            output_name,
            'nexml',
            suppress_annotations=False,
            annotations_as_nhx=False,
            exclude_trees=False)

    def write_phyloxml(self, filename):

        phyloxmls = []

        from Bio import Phylo
        from StringIO import StringIO

        for token in self.trees.keys():
            tree = Phylo.read(StringIO(self.trees[token][0].write(format=2)), 'newick')
            xmltree = tree.as_phyloxml()
            xmltree.name = token

            aln= None
            trimaln= None

            alnlookup = None
            trimalnlookup = None

            if token.split('@')[0] in [l.name for l in self.loci]:

                aln = self.alignments['%s@%s'%(token.split('@')[0], token.split('@')[1])]
                alnlookup = dict((rec.id, str(rec.seq)) for rec in aln)

                trimaln = self.trimmed_alignments['%s@%s@%s'%(token.split('@')[0],
                                                    token.split('@')[1],
                                                    token.split('@')[2])]
                trimalnlookup = dict((rec.id, str(rec.seq)) for rec in trimaln)

            elif token.split('@')[0] in [c.name for c in self.concatenations]:

                trimaln = self.trimmed_alignments[token.split('@')[0]]
                trimalnlookup = dict((rec.id, str(rec.seq)) for rec in trimaln)

            else:
                raise RuntimeError('No locus or concatenation for tree %s'%token)

            for clade in xmltree.get_terminals():
                typ = 'dna'
                if (alnlookup and
                    [l for l in self.loci if l.name == token.split('@')[0]][0].char_type == 'prot'):
                    typ = 'protein'
                key = clade.name
                feautre_id = key
                if not alnlookup:
                    feature_id = (self.trees[token][0]&key).feature_id

                dbtype = 'NCBI'
                if alnlookup and 'denovo' in key:
                    dbtype= 'denovo'
                elif not alnlookup:
                    dbtype= 'otu'

                accession = Phylo.PhyloXML.Accession(key, dbtype)
                mol_seq = Phylo.PhyloXML.MolSeq(trimalnlookup[key], is_aligned=True)
                sequence = Phylo.PhyloXML.Sequence(type=typ,
                                                   accession=accession,
                                                   mol_seq=mol_seq,
                                                   annotations=[Phylo.PhyloXML.Annotation(evidence='trimmed')])

                if alnlookup:
                    mol_seq = Phylo.PhyloXML.MolSeq(alnlookup[key], is_aligned=True)
                    sequence = Phylo.PhyloXML.Sequence(type=typ,
                                                       accession=accession,
                                                       mol_seq=mol_seq)


                clade.sequences.append(sequence)

            others = [Phylo.PhyloXML.Other(tag='trimmedaln',
                                           attributes={'treename':token},
                                           value=trimaln.format('phylip-relaxed'))]
            if aln:
                others.append(Phylo.PhyloXML.Other(tag='aln',
                                                   attributes={'treename':token},
                                                   value=aln.format('phylip-relaxed')))
            xmltree.other = others 

            phyloxmls.append(xmltree)  

        Phylo.write(phyloxmls, filename, 'phyloxml')
        

    def annotate(self, fig_folder,
    
                 root_meta,
                 root_value,
    
                 leaf_labels_txt_meta,
                 leaf_node_color_meta=None,
                 leaf_label_colors=None,
                 ftype='Verdana',
                 fsize=10,
    
                 node_bg_meta=None,
                 node_bg_color=None,
                 
                 node_support_dict=None,
                 support_bullet_size=5,
                 
                 heat_map_meta = None, #list
                 heat_map_colour_scheme=2,
                 
                 pic_meta=None,
                 pic_paths=None,
                 pic_w=None,
                 pic_h=None,
                 
                 multifurc=None,
                 branch_width=2,
                 branch_color='DimGray',
                 
                 scale = 1000,
                 
                 html = None
                 ): 
            
            stdout = sys.stdout
            if html:
                sys.stdout = open(html, 'wt')
                
            print '<html>'
            ts = TreeStyle()
            ts.show_leaf_name = False
            ts.scale = scale
            if node_support_dict:
                ts.legend_position=1
                ts.legend.add_face(TextFace('Node support: ', fsize=10), column=0)
                i = 1
                for color in sorted(node_support_dict.items(),key=lambda i: i[1][0], reverse=True):
                    ts.legend.add_face(CircleFace(radius = 4, color = color[0]), column=i)
                    i +=1 
                    ts.legend.add_face(TextFace(' '+str(node_support_dict[color[0]][0])+'-'+str(node_support_dict[color[0]][1]),
                                                fsize=10), column=i)
                    i += 1
                
            for tree in self.trees.keys():
                                       
                # set outgroup leaves, labels and label colors
                outgroup_list = []
                all_heatmap_profile_values = []
                leaves_for_heatmap = []
                
                for leaf in self.trees[tree][0]:
                    qualifiers_dictionary = get_qualifiers_dictionary(self, leaf.feature_id)
                    leaf_label = ''
                    for meta in leaf_labels_txt_meta:
                        leaf_label += qualifiers_dictionary[meta]+' '
                    leaf_label = leaf_label[:-1]
                    fgcolor = 'black'
                    if leaf_label_colors:
                        for colour_name in leaf_label_colors.keys():
                            if colour_name in qualifiers_dictionary[leaf_node_color_meta]:
                                fgcolor = leaf_label_colors[colour_name]
                    leaf_face = TextFace(leaf_label, fgcolor=fgcolor, ftype=ftype, fsize=fsize)
                    leaf.add_face(leaf_face,0)
                    if not root_value == 'mid' and root_meta in qualifiers_dictionary.keys() and root_value in qualifiers_dictionary[root_meta]:
                        outgroup_list.append(leaf)
                        
                    if heat_map_meta:
                        include = True
                        for i in heat_map_meta:
                            if not i in qualifiers_dictionary:
                                include = False
                        if include:
                            profile = []
                            deviation = []
                            for meta in heat_map_meta:
                                if meta in qualifiers_dictionary.keys():
                                    profile.append(float(qualifiers_dictionary[meta]))
                                    all_heatmap_profile_values.append(float(qualifiers_dictionary[meta]))
                                    deviation.append(0.0)
                            leaf.add_features(profile=profile)
                            leaf.add_features(deviation=deviation)
                            leaves_for_heatmap.append(leaf)
                    if pic_meta:
                        leaf_value = qualifiers_dictionary[pic_meta]
                        if leaf_value in pic_paths:
                            pic_face = ImgFace(pic_paths[leaf_value], width=pic_w, height=pic_h)
                            leaf.add_face(pic_face, 1)
                for leaf in leaves_for_heatmap:
                    leaf.add_face(ProfileFace(max_v=float(max(all_heatmap_profile_values)),
                                              min_v=float(min(all_heatmap_profile_values)), 
                                              center_v=float(float(max(all_heatmap_profile_values)+min(all_heatmap_profile_values))/2),
                                              width=50, height=30,
                                              style='heatmap',
                                              colorscheme=heat_map_colour_scheme),
                                    column=1, position="aligned")
                        
                        
                #set outgroup
                if outgroup_list == ['mid']:
                    try:
                        R = self.trees[tree][0].get_midpoint_outgroup()
                        self.trees[tree][0].set_outgroup(R)
                        print 'rooting tree '+tree+' at midpoint'
                    except:
                        print 'root in '+tree+' already set correctly?'
                    
                elif len(outgroup_list) == 1:
                    try:
                        self.trees[tree][0].set_outgroup(outgroup_list[0])
                    except:
                        print 'root in '+tree+' already set correctly?'
                elif len(outgroup_list) > 1:
                    try:
                        R = self.trees[tree][0].get_common_ancestor(outgroup_list)
                        self.trees[tree][0].set_outgroup(R)
                    except:
                        print 'root in '+tree+' already set correctly?'
                elif len(outgroup_list)==0:
                    try:
                        R = self.trees[tree][0].get_midpoint_outgroup()
                        self.trees[tree][0].set_outgroup(R)
                        print 'rooting tree '+tree+' at midpoint'
                    except:
                        print 'root in '+tree+' already set correctly?'
    
                # ladderize
                self.trees[tree][0].ladderize()
            
                ns = NodeStyle()
                ns['size'] = 0
                ns['fgcolor'] = 'black'
                ns['vt_line_width'] = branch_width
                ns['hz_line_width'] = branch_width
                ns['hz_line_color'] = branch_color
                ns['vt_line_color'] = branch_color
                for n in self.trees[tree][0].traverse():
                    n.set_style(ns)
                self.trees[tree][0].set_style(ns)
            
                if multifurc:
                    for n in self.trees[tree][0].traverse():
                        if n.support < multifurc and not n.is_leaf():
                            n.delete()
    
                # node bg colors
                if node_bg_color:
                    for key in node_bg_color.keys():
                        for node in self.trees[tree][0].get_monophyletic(values=[key], target_attr=node_bg_meta):
                            ns = NodeStyle(bgcolor=node_bg_color[key])
                            ns['size']=0
                            ns['fgcolor']='black'
                            ns['vt_line_width'] = branch_width
                            ns['hz_line_width'] = branch_width
                            ns['hz_line_color'] = branch_color
                            ns['vt_line_color'] = branch_color
                            node.set_style(ns)
                
    
                # node support
                if node_support_dict:
                    for node in self.trees[tree][0].traverse():
                        for key in node_support_dict.keys():
                            if (node.support <= node_support_dict[key][0] and
                                node.support > node_support_dict[key][1]):
                                node.add_face(CircleFace(radius = support_bullet_size,
                                                         color = key),column=0, position = "float")             
                    
                self.trees[tree][0].render(fig_folder + "/"+self.trees[tree][0].get_leaves()[0].tree_method_id+'.png',w=1000, tree_style=ts)
                print('<A href='+
                       fig_folder + "/" + self.trees[tree][0].get_leaves()[0].tree_method_id+'.png'+
                       '>'+self.trees[tree][0].get_leaves()[0].tree_method_id+
                       '</A><BR>')
            print '</html>'
            print fig_folder
            sys.stdout = stdout
     
            
    def trim(self, list_of_Conf_objects, cutoff=0):
        for m in list_of_Conf_objects:
            m.timeit.append(time.time())
            m.platform = platform_report()
            
            m.platform.append('Program and version: '+
                               os.popen(m.cmd + ' --version').readlines()[1].rstrip())
            m.platform.append('Program reference: '+
                              'Salvador Capella-Gutierrez; Jose M. Silla-Martinez; Toni Gabaldon. trimAl: '+
                              'a tool for automated alignment trimming in large-scale '+
                              'phylogenetic analyses. Bioinformatics 2009 25: 1972-1973.')
            
            # In this part, I need to take what I can out of the if isinstance
            # to facilitate addition of new Conf objects. It can be done like this
            # but there will be redundancy
            
            if isinstance(m, TrimalConf):
                import subprocess as sub
                for aln in m.command_lines.keys():
                    p = sub.Popen(m.command_lines[aln], shell=True, stdout=sub.PIPE)
                    stdout, stderr = p.communicate()
                    #stdout = os.system(m.command_lines[aln]).stdout
                    alphabet = IUPAC.ambiguous_dna
                    locus_name = aln.split('@')[0]
                    for locus in self.loci:
                        if locus.name == locus_name and locus.char_type == 'prot':
                            alphabet = IUPAC.protein
                    align = AlignIO.read(StringIO(stdout), "fasta",  alphabet=alphabet)
                    [summary_lines, num_lines, num_undeter, num_collapsed_aln_seqs] = aln_summary(align,
                                                                                                  cutoff=cutoff)
                    summary = 'Alignment name: '+aln+'\n'
                    for line in summary_lines:
                        summary += line+'\n'
                    if num_lines < 4:
                        line = 'Alignment %s has less than 4 sequences and will be dropped'%aln
                        print line
                        summary += line+'\n'
                    elif num_undeter[0] > 0:
                        line = 'Alignment %s has undetermined sequences (%i bp or less) which will be dropped: %s'%(aln, cutoff+1, num_undeter[1])
                        print line
                        summary += line+'\n'
                        records_wo_undeter = []
                        for record in align:
                            if not record.id in num_undeter[1]:
                                records_wo_undeter.append(record)
                        align =  MultipleSeqAlignment(records_wo_undeter)
                        self.trimmed_alignments[aln] = align
                    elif num_collapsed_aln_seqs < 4:
                        line = 'Alignment %s has less than 4 unique sequences and will be dropped'%aln
                        print line
                        summary += line+'\n'
                    else:
                        self.trimmed_alignments[aln] = align
                    self.aln_summaries.append(summary)
            for file_name in os.listdir(os.curdir):
                if m.id.partition('_')[0] in file_name:
                    os.remove(file_name)
            m.timeit.append(time.time())
            m.timeit.append(m.timeit[2]-m.timeit[1])
            self.used_methods[m.method_name] = m
            
        if self.pickle_name:
            pickle_pj(self, self.pickle_name)
        
        if __builtin__.git and self.pickle_name:
            
            comment = ''
            for m in list_of_Conf_objects:
                comment += '%s\n'%(str(m))
            
            import rpgit
            out, err = rpgit.gitAdd(self.pickle_name)
            undate_git_log(self, out, err)
            cwd = os.getcwd()
            import fnmatch
            matches = []
            for root, dirnames, filenames in os.walk(cwd):
                for filename in fnmatch.filter(filenames, '*.py'):
                    matches.append(os.path.join(root, filename))
                for filename in fnmatch.filter(filenames, '*.ipynb'):
                    matches.append(os.path.join(root, filename))
            for match in matches:
                out, err = rpgit.gitAdd(match)
                undate_git_log(self, out, err)
            out, err = rpgit.gitCommit(comment)
            undate_git_log(self, out, err)

    def report_seq_stats(self):        
        if len(self.records_by_locus.keys())>0:

            # This will make a list of seq length for each locus. Seq length are calced
            # using the record.seq in 'pj.records_by_locus'. 'pj.records_by_locus is a
            # dict with loci names as keys, and lists of SeqReocrd objects as values
            lengths_dict = {}
            for locus_name in self.records_by_locus.keys():
                lengths_dict[locus_name] = []
                for record in self.records_by_locus[locus_name]:
                    lengths_dict[locus_name].append(len(record.seq))
                    
            print "Distribution of sequence lengths".title()
            draw_boxplot(lengths_dict, 'Seq length (bp)', 'inline')
            
            
            for stat in ['GC_content']:#, 'nuc_degen_prop', 'prot_degen_prop']:
                title = 'Distribution of sequence statistic \"'+stat+'\"'
                print title.title()
                # This will make a dict with loci as keys and a list of stat values as
                # dict values.
                stat_dict = {}
                ylabel = 'GC ontent (%)'
                if not stat == 'GC_content':
                    ylabel = 'Ambiguous positions (prop)'
                for locus_name in self.records_by_locus.keys():
                    stat_dict[locus_name] = []
                    for i in self.records_by_locus[locus_name]:
                        for record in self.records:
                            for feature in record.features:
                                if feature.qualifiers['feature_id'][0] == i.id:
                                    if stat in feature.qualifiers.keys():
                                        stat_dict[locus_name].append(float(feature.qualifiers[stat][0]))
                
                draw_boxplot(stat_dict, ylabel, 'inline')
                
    ##################################               
    # Project methods to fetch objects
    ##################################
    
    def ft(self, token):
        
        """
        Will fetch the tree object which has the token in 
        its key, as long as there is only one
        """
        
        # check how many tree keys match the token
        keys = [key for key in self.trees.keys() if token in key]
        if len(keys) > 1:
            raise IOError("The token %s was found in more then one tree key: %s"
                           %(token, str(keys)))
        elif len(keys) == 0:
            raise IOError("The token %s was not found in any tree key"
                           %token)
        elif len(keys) == 1:
            print "returning tree object %s"%keys[0]
            return Tree(self.trees[keys[0]][1])
        
        
        
    def fa(self, token):
        
        """
        Will fetch the alignment object which has the token in 
        its key, as long as there is only one
        """
        
        from StringIO import StringIO
        
        # check how many aln keys match the token
        keys = [key for key in self.alignments.keys() if token in key]
        if len(keys) > 1:
            raise IOError("The token %s was found in more then one alignment key: %s"
                          %(token, str(keys)))
        elif len(keys) == 0:
            raise IOError("The token %s was not found in any alignment key"
                          %token)
        elif len(keys) == 1:
            print "returning alignment object %s"%keys[0]
            return AlignIO.read(StringIO(self.alignments[keys[0]].format('fasta')), 'fasta')        
        
        
    def fta(self, token):
        
        """
        Will fetch the trimmed alignment object which has the token in 
        its key, as long as there is only one
        """
        
        from StringIO import StringIO
        
        # check how many trimmed aln keys match the token
        keys = [key for key in self.trimmed_alignments.keys() if token in key]
        if len(keys) > 1:
            raise IOError("The token %s was found in more then one trimmed alignment key: %s"
                          %(token, str(keys)))
        elif len(keys) == 0:
            raise IOError("The token %s was not found in any trimmed alignment key"
                          %token)
        elif len(keys) == 1:
            print "returning trimmed alignment object %s"%keys[0]
            return AlignIO.read(StringIO(self.trimmed_alignments[keys[0]].format('fasta')), 'fasta')            
        
        
    def fr(self, locus_name, filter=None):
        
        """
        Will fetch the record objects of the specified locus, 
        as long as there is at least one.
        filter should be a list of lists. Every (sub)list is
        a pair of qualifier and value. If filter is specified,
        only records that have all the specified values in the
        specified qualifiers will be kept.
        """
        
        # check how many record keys match the token
        keys = [key for key in self.records_by_locus.keys() if locus_name in key]
        if len(keys) > 1:
            raise IOError("The locus name %s fits more then one locus: %s"
                          %(locus_name, str(keys)))
        elif len(keys) == 0:
            raise  IOError("The locus %s was not found"
                           %locus_name)
        elif len(keys) == 1:
            records = []
            if filter:
                for r in self.records_by_locus[keys[0]]:
                    qualifiers = get_qualifiers_dictionary(self, r.id)
                    get = True
                    for f in filter:
                        if not (f[0] in qualifiers.keys() and qualifiers[f[0]] == f[1]):
                            get = False
                    if get:
                        records.append(r) 
            else:
                for r in self.records_by_locus[keys[0]]:
                    records.append(r)
            print "returning records list of locus %s and filter %s"%(keys[0], str(filter))
            return records
        
   
    def propagate_metadata(self):
        for t in self.trees.keys():
            for l in self.trees[t][0]:
                feature_id = l.feature_id
                record_id = feature_id.rpartition('_')[0]
                record = [r for r in self.records if r.id == record_id][0]
                annotations = record.annotations
                source_qualifiers = [f for f in record.features if f.type == 'source'][0].qualifiers
                feature_qualifiers = [f for f in record.features if f.qualifiers['feature_id'][0] == feature_id][0].qualifiers
                for i in annotations:
                    l.add_feature("annotations_%s"%i,type_to_single_line_str(annotations[i]))
                for i in source_qualifiers:
                    l.add_feature("source_%s"%i,type_to_single_line_str(source_qualifiers[i]))
                for i in feature_qualifiers:
                    l.add_feature(i,type_to_single_line_str(feature_qualifiers[i]))
            self.trees[t][1] = self.trees[t][0].write(features=[])

            
##############################################################################################
if False:
    """Tools for loci explorations in a GenBank File"""
##############################################################################################

programspath = ""

def list_loci_in_genbank(genbank_filename, control_filename, loci_report = None):
    """
    Takes a genbank file, returns a loci csv file and a list of loci and their counts. The latter goes
    either to stdout or to a file.
    
    >>> list_loci_in_genbank("test-data/test.gb", "test-data/temp_loci.csv", loci_report = "test-data/temp_loci.txt")
    >>> assert open("test-data/temp_loci.csv",'r').read() == open("test-data/test_loci.csv",'r').read() 
    >>> import os
    >>> os.remove("test-data/temp_loci.csv")
    """
    
    stdout = sys.stdout
    if  loci_report: 
        sys.stdout = open(loci_report, 'w')
    
    genbank_synonyms = gb_syn.gb_syn()
    
    # Open GenBank file
    MelPCgenes = open(genbank_filename, 'rU')
   
    gene_dict = {} #set up a gene_dict dictionary
   
    # For each record
    
    for record in SeqIO.parse(MelPCgenes, 'genbank') :
        # Look at all features for this record
        for feature in record.features:
         
            # If it's a CDS or rRNA...
            if feature.type == 'CDS' or feature.type == 'rRNA':
                # If it contains some attribute called 'gene' save that
                if 'gene' in feature.qualifiers:
                    geneName = feature.qualifiers['gene'][0]
                    geneName = geneName.replace(',','_')
                    geneName = geneName.replace('/','_')

                    if feature.type+','+geneName in gene_dict:
                        gene_dict[feature.type+','+geneName]+=1
                    else:    
                        gene_dict[feature.type+','+geneName]=1
                    #print(geneName)

                # Else if it contains a 'product' qualifier
                elif 'product' in feature.qualifiers:
                    geneName = feature.qualifiers['product'][0]
                    geneName = geneName.replace(',','_')
                    geneName = geneName.replace('/','_') 

                    if feature.type+','+geneName in gene_dict:
                        gene_dict[feature.type+','+geneName]+=1
                    else:    
                        gene_dict[feature.type+','+geneName]=1
                    #print(geneName)

                else:
                    print 'Could not find either gene or product in '+record.id
                    #print feature.qualifiers

       
    #sorting happens via a list
   
    sorted_gene_names = gene_dict.items()
    sorted_gene_names.sort(key = lambda i: i[0].lower())
    control_file_lines = {}
   
   
    print('\n' + "There are " + str(len(sorted_gene_names)) + " gene names (or gene product names) detected")
    print("----------------------------------")
    print("Gene and count sorted by gene name")
    print("----------------------------------")
   
    for key, value in sorted_gene_names:
        #print key, value
        print(str(value) +" instances of " + key)
        feature_type = key.split(',')[0]
        alias = key.split(',')[1]
        gen_group = None
        for group in genbank_synonyms:
            if alias in group:
                gen_group = group
        if gen_group:
            if gen_group[0].replace(' ','_') in control_file_lines.keys():
                control_file_lines[gen_group[0].replace(' ','_')].append(alias)
            else:
                control_file_lines[gen_group[0].replace(' ','_')] = [feature_type, alias]
        else:
            name = alias.replace(' ','_').replace('/','_')
            control_file_lines[name] = [feature_type, alias]
                    
    control_file_handle = open(control_filename, 'wt')
    for line in sort(control_file_lines.keys()):  
        control_file_handle.write('dna,%s,%s'%(control_file_lines[line][0],line))
        for a in control_file_lines[line][1:]: 
            control_file_handle.write(',%s'%a)
        control_file_handle.write('\n')
                            
    control_file_handle.close()                 
   
    print("-------------------------------")
    print("Gene and count sorted by counts")
    print("-------------------------------")
    sorted_gene_names.sort(key = lambda i: int(i[1]), reverse=True)
    for key, value in sorted_gene_names:
        #print key, value
        print(str(value) +" instances of " + key)
    sys.stdout = stdout
    return   

AlnConfMethodsSection=\
"""\n
==============================
Core Methods section sentence:
==============================
The dataset(s) %s were aligned using the program %s [1].

Reference:
%s
"""

AlnConfPalNalMethodsSection=\
"""\n
==============================
Core Methods section sentence:
==============================
The dataset(s) %s were first aligned at the protein level using the program %s [1].
The resulting alignments served as guides to codon-align the DNA sequences using %s [2].

Reference:
[1]%s
[2]%s
"""

# ALIGNMENT PROGRAM PLUG
# AS ABOVE, A STRING IS NEEDED WITH PLACE HOLDERS FOR THE LOCI, THE PROGRAM+VERSION
# AND FOR THE REFERENCE.
# ANOTHER STRING IS REQUIRED THAT ALSO INCLUDES PAL2NAL
    
##############################################################################################
class AlnConf:
##############################################################################################
    
    """
    >>> coi = Locus('dna', 'CDS', 'coi', ['cox1','COX1','coi','COI','CoI'])
    >>> lsu = Locus('dna', 'rRNA', '28S', ['28s','28S','LSU rRNA','28S ribosomal RNA','28S large subunit ribosomal RNA'])
    >>> pj = Project([coi, lsu], git=False)
    
    # cline_str = muscle = AlnConf(pj, method_name='MuscleDefaults',
    #                                      cmd='muscle', program_name='muscle',
    #                                     cline_args=dict())
    """
    
    def __init__(self, pj, method_name='mafftDefault', CDSAlign=True, codontable=1, program_name='mafft',
                 cmd='mafft', loci='all',
                 cline_args=dict()):
        
        if pj.records_by_locus == {}:
            pj.extract_by_locus()
            
        self.id = str(random.randint(10000,99999))+str(time.time())
        self.method_name=method_name
        self.CDSAlign=CDSAlign
        self.program_name=program_name
        self.loci = pj.loci
        if not loci == 'all':
            self.loci = []
            for locus_name in loci:
                for locus in pj.loci:
                    if locus_name == locus.name:
                        self.loci.append(locus)
        mutable_loci_list = []
        removed_loci = []
        for locus in self.loci: 
            if len(pj.records_by_locus[locus.name]) < 4:
                removed_loci.append(locus.name)
            else:
                mutable_loci_list.append(locus)
                
        if len(removed_loci) > 0:
            print "These loci have less than 4 sequences and will be dropped from this conf object. Don't use them in a concatenation:\n%s\n\n"%removed_loci
        self.loci = mutable_loci_list
        self.CDS_proteins = {}
        self.CDS_in_frame = {}
        self.codontable = codontable 
        self.aln_input_strings = {}
        self.command_lines = {}
        self.timeit = [time.asctime()]
        self.platform = []
        self.cline_args = cline_args
        self.cmd = cmd
        if not program_name == 'mafft' and cmd == 'mafft':
            self.cmd = pj.defaults[program_name]
        elif program_name == 'mafft' and cmd != 'mafft' and not 'mafft' in cmd:
            self.cmd = pj.defaults['mafft']
            
        # make defalut input files

        for key in pj.records_by_locus.keys():
            if key in [l.name for l in self.loci]:
                SeqIO.write(pj.records_by_locus[key], self.id+'_'+key+'.fasta', 'fasta')    
        for locus in self.loci:
            # put default input file filename and string in the AlnConf object
            input_filename=self.id+'_'+locus.name+'.fasta'
            self.aln_input_strings[locus.name] = [open(input_filename,'r').read()]
            # If CDS and CDSAlign, prepare reference protein input file and in frame CDS input file
            if locus.feature_type == 'CDS' and locus.char_type == 'dna' and self.CDSAlign: 
                self.CDS_proteins[locus.name] = []
                self.CDS_in_frame[locus.name] = []
                for record in pj.records:
                    for feature in record.features:
                        if (feature.type == 'CDS' and 'gene' in feature.qualifiers.keys() and
                            feature.qualifiers['gene'][0] in locus.aliases and 
                            feature.qualifiers['feature_id'][0] in [f.id for f in pj.records_by_locus[locus.name]]):
                            S = feature.extract(record.seq)
                            # Make in-frame CDS input file seq start in frame
                            if 'codon_start' in feature.qualifiers.keys():
                                i = feature.qualifiers['codon_start'][0]
                                if i > 1:
                                    S = S[(int(i)-1):]
                            # Make in-frame CDS input file seq end in frame
                            if len(S)%3 == 1:
                                S = S[:-1]
                            elif len(S)%3 == 2:
                                S = S[:-2]  
                            # make protein input file seq
                            if not 'translation' in feature.qualifiers.keys():
                                raise IOError("Feature %s has no 'translation' qualifier"%
                                              feature.qualifiers['feature_id'][0])
                            P = Seq(feature.qualifiers['translation'][0], IUPAC.protein)
                            # Remove 3' positions that are based on partial codons
                            while len(P)*3 > len(S):
                                P = P[:-1]
                            # remove complete termination codon
                            if (len(S)/3)-1 == len(P):
                                S = S[:-3]
                            # make in frame cds record
                            feature_record = SeqRecord(seq = S, id = feature.qualifiers['feature_id'][0],
                                                               description = '')
                            # put it in the object
                            self.CDS_in_frame[locus.name].append(feature_record)
                            # make protein record
                            feature_record = SeqRecord(seq = P, id = feature.qualifiers['feature_id'][0],
                                                       description = '')
                            # Put the protein records in the AlnConf object
                            self.CDS_proteins[locus.name].append(feature_record)
                           
                                    
                # check same number of prot and cds objects                    
                if len(pj.records_by_locus[locus.name]) > len(self.CDS_proteins[locus.name]):
                    raise RuntimeError('For the CDS locus '+locus.name+': more nuc seqs than prot seqs.'+
                                       ' You may miss a \'translate\' or \'gene\' qualifier in some of '+
                                       'the features.')
                elif len(pj.records_by_locus[locus.name]) < len(self.CDS_proteins[locus.name]):
                    raise RuntimeError('For the CDS locus '+locus.name+': less nuc seqs than prot seqs.'+
                                       ' You may miss a \'translate\' or \'gene\' qualifier in some of '+
                                       'the features.')
                unmatched = []
                for record in self.CDS_in_frame[locus.name]:
                    for prot in self.CDS_proteins[locus.name]:
                        if prot.id == record.id:
                            if not len(prot.seq)*3 == len(record.seq):
                                unmatched.append(record.id)
                unmatched_string = ''
                if len(unmatched) > 0:
                    for u in unmatched:
                        unmatched_string += u+' '
                    raise RuntimeError('The following CDS/protein pairs are unmatched: '+unmatched_string)
                    
 
                SeqIO.write(self.CDS_in_frame[locus.name],
                            self.id+'_CDS_in_frame_'+locus.name+'.fasta', 'fasta')
                input_filename2=self.id+'_CDS_in_frame_'+locus.name+'.fasta'
                SeqIO.write(self.CDS_proteins[locus.name],
                            self.id+'_CDS_proteins_'+locus.name+'.fasta', 'fasta')
                input_filename=self.id+'_CDS_proteins_'+locus.name+'.fasta'
                self.aln_input_strings[locus.name][0] = [open(input_filename,'r').read(),
                                                         open(input_filename2,'r').read()]
            cline = dict(dict(input=input_filename), **cline_args)
            if self.program_name == 'mafft':
                self.command_lines[locus.name] = MafftCommandline(cmd=self.cmd)
            elif self.program_name == 'muscle':
                self.command_lines[locus.name] = MuscleCommandline(cmd=self.cmd)
            # PROGRAM PLUG
            # NOTE: MAFFT AND MUSCLE GET FASTA. IF THE NEW PROGRAM GETS SOMETHING
            # ELSE, SOME WORK IS REQUIRED ABOVE, e.g. condition to choose the format
            # THIS WRITE THE CLINE AND PLACES IN PROJECT. 
            # elif self.program_name == 'some program':
            #    self.command_lines[locus.name] = some_program_cline_object(cmd=self.cmd)
            # This assumes Bio.Applications style cline object
            for c in cline.keys():
                self.command_lines[locus.name].__setattr__(c,cline[c])
            print str(self.command_lines[locus.name])
            
            
            
    def __str__(self):
        loci_string = ''
        for n in [i.name for i in self.loci]:
            loci_string += n+','
        loci_string = loci_string[:-1]
        command_lines = ''
        for i in self.command_lines.keys():
            command_lines += i+': '+str(self.command_lines[i])+'\n'
        date = str(self.timeit[0])
        execution = '[This was not executed yet]'
        if len(self.timeit) > 1:
            execution = str(self.timeit[3])
        plat = '[This was not executed yet]'
        if len(self.platform) > 0:
            plat = str(self.platform).replace(",",'\n').replace(']','').replace("'",'').replace('[','')
        output =  str("AlnConf named %s with ID %s\n"+         
                "Loci: %s \n"+       
                "Created on: %s\n"+
                "Commands:\n"+
                "%s\n"+
                "Environment:\n"+    
                "%s\n"+
                "execution time:\n"+
                "%s\n")%(self.method_name, str(self.id), loci_string, date, command_lines, plat, execution) 
        if len(self.platform) > 0 and self.CDSAlign:
            prog = '[This was not executed yet]'
            pal = '[This was not executed yet]'
            progref = '[This was not executed yet]'
            palref = '[This was not executed yet]'
            if self.platform[-1].startswith('Program reference:'):
                prog, pal = self.platform[-2].split(':')[1].strip().rstrip().split('\n')
                progref, palref = self.platform[-1].split('Program reference:')[1].strip().rstrip().split('\n')
            output += AlnConfPalNalMethodsSection%(type_to_single_line_str([l.name for l in self.loci]),
                                                  prog, pal, progref, palref)
        elif len(self.platform) > 0:
            prog = '[This was not executed yet]'
            progref = '[This was not executed yet]'
            if self.platform[-1].startswith('Program reference:'):
                prog = self.platform[-2].split(':')[1].strip().rstrip()
                progref = self.platform[-1].split('Program reference:')[1].strip().rstrip()
            output += AlnConfMethodsSection%(type_to_single_line_str([l.name for l in self.loci]),
                                             prog, progref)
        
        return output
    
TrimalConfMethodsSection=\
"""\n
==============================
Core Methods section sentence:
==============================
The alignment(s) %s were trimmed using the program %s [1].

Reference:
%s
"""
##############################################################################################
class TrimalConf:
##############################################################################################
    def __init__(self, pj, method_name='gappyout', program_name='trimal',
                 cmd='default', alns='all', trimal_commands=dict(gappyout=True)):
        
        if len(pj.alignments) == 0:
            raise RuntimeError("No sequence alignments found")
        self.id = str(random.randint(10000,99999))+str(time.time())
        self.method_name=method_name
        self.program_name=program_name
        self.alignments = pj.alignments
        if not alns == 'all':
            self.alignments = {}
            for aln_name in alns:
                if aln_name in pj.alignments.keys():
                    self.alignments[aln_name] = pj.alignments[aln_name]
        self.command_lines = {}
        self.timeit = [time.asctime()]
        self.platform = []
        self.cline_args = trimal_commands
        self.cmd = cmd
        if cmd == 'default':
            self.cmd = pj.defaults['trimal']
        irelevant = ['out', 'clustal', 'fasta', 'nbrf', 'nexus', 'mega',
                     'phylip3.2', 'phylip', 'sgt', 'scc', 'sct', 'sfc',
                     'sft','sident']
        for aln in self.alignments:
            input_filename = self.id+'_'+aln+'.fasta'
            AlignIO.write(self.alignments[aln], input_filename,'fasta')
            self.command_lines[aln] = "%s -in %s "%(self.cmd, input_filename)
            for kwd in trimal_commands:
                if kwd in irelevant:
                    warnings.simplefilter('always')
                    warnings.warn("%s is irelevant in this context and will be ignored"%kwd)
                else:
                    if not trimal_commands[kwd] == None and not trimal_commands[kwd] == False:
                        if trimal_commands[kwd] == True:
                            self.command_lines[aln] += "-%s "%(kwd)
                        else:
                            self.command_lines[aln] += "-%s %s "%(kwd, trimal_commands[kwd])
            
            self.command_lines[aln+'@'+self.method_name] = self.command_lines[aln][:-1]
            print self.command_lines[aln+'@'+self.method_name]
            self.command_lines.pop(aln, None)
        
    def __str__(self):
        aln_string = ''
        for n in self.alignments.keys():
            aln_string += n+','
        aln_string = aln_string[:-1]
        command_lines = ''
        for i in self.command_lines.keys():
            command_lines += i+': '+str(self.command_lines[i])+'\n'
        date = str(self.timeit[0])
        execution = '[This was not executed yet]'
        if len(self.timeit) > 1:
            execution = str(self.timeit[3])
        
        prog = '[This was not executed yet]'
        plat = '[This was not executed yet]'
        progref = '[This was not executed yet]'
        if len(self.platform) > 0:
            plat = str(self.platform).replace(",",'\n').replace(']','').replace("'",'').replace('[','')
            if self.platform[-1].startswith('Program reference:'):
                prog = self.platform[-2].split(':')[1].strip().rstrip()
                progref = self.platform[-1].split('Program reference:')[1].strip().rstrip()
        
        output= str("TrimalConf named %s with ID %s\n"+         
                "Alignments: %s \n"+       
                "Created on: %s\n"+
                "Commands:\n"+
                "%s\n"+
                "Environment:"+    
                "%s\n"+
                "execution time:\n"+
                "%s")%(self.method_name, str(self.id), aln_string, date, command_lines, plat, execution)
        
        output += TrimalConfMethodsSection%(type_to_single_line_str(self.alignments.keys()),
                                            prog, progref)
        
        return output

def use_sh_support_as_branch_support(tree_filename):
    string = open(tree_filename,'r').read()
    string = re.sub(r'\[',r'[&&NHX:support=',string)
    t = Tree(string)
    t.dist=0
    return t.write(features=[])
    
def transfer_support_same_topo(tree_file_with_support,
                               tree_file_without_support):
    supported = Tree(tree_file_with_support) 
    unsupported = Tree(tree_file_without_support)
    supported_leaf_names = sorted(supported.get_leaf_names())
    unsupported_leaf_names = sorted(unsupported.get_leaf_names())
    if not len(unsupported_leaf_names) == len(supported_leaf_names):
        raise IOError(tree_file_with_support + ' and ' + tree_file_without_support +
                      ' are not the same length')
    for i in range(len(supported_leaf_names)):
        if not supported_leaf_names[i] == unsupported_leaf_names[i]:
            raise IOError('The trees do not share all leaves or leaf names')
    same_root = supported.get_leaf_names()[0]
    unsupported.set_outgroup(same_root)
    supported.set_outgroup(same_root)
    for ns in supported.traverse():
        ns_leaves = ns.get_leaf_names()
        if not unsupported.check_monophyly(values=ns_leaves, target_attr="name"):
            raise RuntimeError('trees do not share topology and/or all the leaf names')
        else:
            unsupported_ancestor = unsupported.get_common_ancestor(ns_leaves)
            unsupported_ancestor.support = ns.support
    unsupported.write(outfile = tree_file_without_support)    
    
    
def make_raxml_partfile(tree_method, pj, trimmed_alignment_name):

    concatenation = None
    for c in pj.concatenations:
        if c.name == trimmed_alignment_name:
            concatenation = c
    
    #concatenation = filter(lambda concatenation: concatenation.name == trimmed_alignment_name, pj.concatenations)[0]
    
    model = []
    for locus in concatenation.loci:
        part_name = None
        for trm_aln in concatenation.used_trimmed_alns.keys():
            if locus.name == trm_aln.partition('@')[0]:
                part_name = trm_aln
        if not part_name:
            warnings.warn('There is no trimmed alignment for locus '+locus.name+' in concatenation '+concatenation.name)
        else:
            part_length = concatenation.used_trimmed_alns[part_name]
            if locus.char_type == 'prot':
                m = None
                if isinstance(tree_method.matrix,dict):
                    m = tree_method.matrix[locus.name]
                elif isinstance(tree_method.matrix,str):
                    m = tree_method.matrix
                else:
                    #todo write error
                    pass
                model.append([m,part_name,part_length])
            elif locus.char_type == 'dna':
                model.append(['DNA',part_name,part_length])
                    
    # make partition file
                    
    partfile = open(tree_method.id+'_'+concatenation.name+'_partfile','wt')
    i = 1
    for m in model:
        partfile.write(m[0]+', '+m[1]+'='+str(i)+'-'+str(m[2]+i-1)+'\n')
        i += m[2]
    partfile.close()
    return tree_method.id+'_'+concatenation.name+'_partfile'

def make_raxml_input_matrix_file(tree_method, trimmed_alignment_name):
    SeqIO.write(tree_method.trimmed_alignments[trimmed_alignment_name],
                tree_method.id+'_'+trimmed_alignment_name+'.fasta','fasta')
    return tree_method.id+'_'+trimmed_alignment_name+'.fasta'

def write_raxml_clines(tree_method, pj, trimmed_alignment_name):
            
    cline_que = 0

    support_replicates = 100
    ML_replicates = 1
    if '-N' in tree_method.cline_args.keys():
        ML_replicates = tree_method.cline_args['-N']
    if '-#' in tree_method.cline_args.keys():
        support_replicates = tree_method.cline_args['-#']
    
    partfile = None
    
    # Check if it is a concatenation and make partfile

    for c in pj.concatenations:
        if c.name == trimmed_alignment_name.partition('@')[0]:
            partfile = make_raxml_partfile(tree_method, pj, trimmed_alignment_name)
    
    input_filename = make_raxml_input_matrix_file(tree_method, trimmed_alignment_name)
    model = tree_method.model
    try:
        locus_char_type = filter(lambda locus: locus.name == trimmed_alignment_name.partition('@')[0], pj.loci)[0].char_type
    except:
        locus_char_type = 'prot'
    
    if partfile:
        model='PROT'+model+'JTT'
    else:
        if locus_char_type == 'dna':
            model = 'GTR'+tree_method.model
        elif  locus_char_type == 'prot':
            if isinstance(tree_method.matrix,str):
                model = 'PROT'+tree_method.model+tree_method.matrix
            elif isinstance(tree_method.matrix,dict):
                model = 'PROT'+tree_method.model+tree_method.matrix[trimmed_alignment_name]
        
    presets = {'fa': [{'-f': 'a',
                           '-p': random.randint(99,999),
                           '-x':  random.randint(99,999),
                           '-s': input_filename,
                           '-N': support_replicates,
                           '-n': tree_method.id+'_'+trimmed_alignment_name+'0',
                           '-m': model}
                      ],
                'fD_fb':[{'-f': 'D',
                          '-p': random.randint(99,999),
                          '-s': input_filename,
                          '-N': ML_replicates,
                           '-n': tree_method.id+'_'+trimmed_alignment_name+'0',
                           '-m': model},{'-f': 'b',
                                                       '-p': random.randint(99,999),
                                                       '-s': input_filename,
                                                       '-n': tree_method.id+'_'+trimmed_alignment_name+'1',
                                                       '-m': model,
                                                       '-t': 'RAxML_bestTree.'+tree_method.id+'_'+trimmed_alignment_name+'0',
                                                       '-z': 'RAxML_rellBootstrap.'+tree_method.id+'_'+trimmed_alignment_name+'0'}
                         ],
                'fd_b_fb':[{'-f': 'd',
                          '-p': random.randint(99,999),
                          '-s': input_filename,
                          '-N': ML_replicates,
                           '-n': tree_method.id+'_'+trimmed_alignment_name+'0',
                           '-m': model},{
                                                       '-p': random.randint(99,999),
                                                       '-b': random.randint(99,999),
                                                       '-s': input_filename,
                                                       '-#': support_replicates,
                                                       '-n': tree_method.id+'_'+trimmed_alignment_name+'1',
                                                       '-m': model,
                                                       '-T': tree_method.threads},
                                                       {'-f': 'b',
                                                       '-p': random.randint(99,999),
                                                       '-s': input_filename,
                                                       '-n': tree_method.id+'_'+trimmed_alignment_name+'2',
                                                       '-m': model,
                                                       '-t': 'RAxML_bestTree.'+tree_method.id+'_'+trimmed_alignment_name+'0',
                                                       '-z': 'RAxML_bootstrap.'+tree_method.id+'_'+trimmed_alignment_name+'1'}
                         ],
                'fF_fJ': [{'-f': 'F',
                           '-p': random.randint(99,999),
                           '-s': input_filename,
                           '-n': tree_method.id+'_'+trimmed_alignment_name+'0',
                           '-m': model},{'-f': 'J',
                                         '-t': 'RAxML_fastTree.'+tree_method.id+'_'+trimmed_alignment_name+'0',
                                         '-p': random.randint(99,999),
                                         '-s': input_filename,
                                         '-n': tree_method.id+'_'+trimmed_alignment_name+'1',
                                         '-m': model}],
                
                'fd_fJ': [{'-f': 'd',
                          '-p': random.randint(99,999),
                          '-s': input_filename,
                          '-N': ML_replicates,
                          '-n': tree_method.id+'_'+trimmed_alignment_name+'0',
                          '-m': model},{'-f': 'J',
                                         '-t': 'RAxML_bestTree.'+tree_method.id+'_'+trimmed_alignment_name+'0',
                                         '-p': random.randint(99,999),
                                         '-s': input_filename,
                                         '-n': tree_method.id+'_'+trimmed_alignment_name+'1',
                                         '-m': model}]
                
                }

    
    if 'PTHREADS' in tree_method.cmd:
        if tree_method.threads < 2:
            tree_method.threads = 2
            warnings.warn('raxmlHPC-PTHREADS requires at least 2 threads. Setting threads to 2')
        for preset in presets:
            for c in presets[preset]:
                c['-T'] = tree_method.threads
    else:
        if tree_method.threads > 1:
            raise RuntimeWarning('This is a serial raxmlHPC. Setting threads to 1.'+ 
                                 'PTHREADS executables have to have explicit filename, eg, raxmlHPC-PTHREADS-SSE3')
    if partfile:
        for preset in presets.keys():
            for cline in range(len(presets[preset])):
                presets[preset][cline] = dict({'-q': partfile}, **presets[preset][cline])               
    return presets[tree_method.preset] 


RaxmlConfMethodsSection=\
"""\n
==============================
Core Methods section sentence:
==============================
Phylogenetic trees were reconstructed from the dataset(s) %s using the program %s [1].

Reference:
%s
"""
##############################################################################################
class RaxmlConf:
##############################################################################################
    
    def __init__(self, pj, method_name='fa', program_name='raxmlHPC-PTHREADS-SSE3', keepfiles=False,
                 cmd='default', preset = 'fa', alns='all', model='GAMMA', matrix='JTT', threads=4,
                 cline_args={}):
        
        
        if len(pj.trimmed_alignments) == 0:
            raise RuntimeError("No trimmed sequence alignments found")
            
        self.id = str(random.randint(10000,99999))+str(time.time())
        self.method_name=method_name
        self.program_name=program_name
        self.preset = preset
        self.cline_args = cline_args
        self.model = model
        self.matrix = matrix
        self.threads = threads
        self.trimmed_alignments = pj.trimmed_alignments
        if not alns == 'all':
            self.trimmed_alignments = {}
            for aln_name in alns:
                if aln_name in pj.trimmed_alignments.keys():
                    self.trimmed_alignments[aln_name] = pj.trimmed_alignments[aln_name]
        self.aln_input_strings = {}
        self.command_lines = {}
        self.timeit = [time.asctime()]
        self.platform = []
        self.cmd = cmd
        self.keepfiles = keepfiles
        if cmd == 'default':
            self.cmd = pj.defaults['raxmlHPC']
        
        for trimmed_alignment in self.trimmed_alignments.keys():
            self.command_lines[trimmed_alignment] = []
            command_lines = write_raxml_clines(self, pj, trimmed_alignment)
            for command_line in command_lines:
                cline_object = RaxmlCommandline(cmd=self.cmd)
                for c in command_line.keys():
                    cline_object.__setattr__(c,command_line[c])
                self.command_lines[trimmed_alignment].append(cline_object)
                print str(cline_object)

    def __str__(self):

        aln_string = ''
        for n in self.trimmed_alignments.keys():
            aln_string += n+','
        aln_string = aln_string[:-1]
        command_lines = ''
        for i in self.command_lines.keys():
            command_lines += i+':\n'
            for k in self.command_lines[i]:
                command_lines += str(k) + '\n'
        date = str(self.timeit[0])
        
        execution = '[This was not executed yet]'
        if len(self.timeit) > 1:
            execution = str(self.timeit[3])
        
        prog = '[This was not executed yet]'
        plat = '[This was not executed yet]'
        progref = '[This was not executed yet]'
        
        if len(self.platform) > 0:
            plat = str(self.platform).replace(",",'\n').replace(']','').replace("'",'').replace('[','')
            if self.platform[-1].startswith('Program reference:'):
                prog = self.platform[-2].split(':')[1].strip().rstrip()
                progref = self.platform[-1].split('Program reference:')[1].strip().rstrip()
                
        output = str("RaxmlConf named %s with ID %s\n"+         
                "Alignments: %s \n"+       
                "Created on: %s\n"+
                "Commands:\n"+
                "%s\n"+
                "Environment:\n"+    
                "%s\n"+
                "execution time:\n"+
                "%s")%(self.method_name, str(self.id), aln_string, date, command_lines, plat, execution)
        
        output += RaxmlConfMethodsSection%(type_to_single_line_str(self.trimmed_alignments.keys()),
                                            prog, progref)
        
        return output

def make_pb_input_matrix_file(conf_obj, trimmed_alignment_name):
    SeqIO.write(conf_obj.trimmed_alignments[trimmed_alignment_name],
                conf_obj.id+'_'+trimmed_alignment_name+'.phylip','phylip-relaxed')
    return conf_obj.id+'_'+trimmed_alignment_name+'.phylip'


def write_pb_cline(conf_obj, pj, trimmed_alignment):
    cline = "%s -d %s"%(conf_obj.cmd, make_pb_input_matrix_file(conf_obj, trimmed_alignment))
    for key in conf_obj.cline_args:
        kw = key
        if key[0] == '-':
            kw = key[1:]
        cline += " -%s"%str(kw)
        if not str(conf_obj.cline_args[key]) == str(True):
            cline += " %s"%str(conf_obj.cline_args[key]) 
    cline += " %s_%s"%(conf_obj.id, trimmed_alignment)
    return cline
            
###################################################################################################
class PbConf:
###################################################################################################
    
    def __init__(self, pj, method_name = 'dna_cat_gtr', program_name='phylobayes', keepfiles=True,
                 cmd='default', alns='all', cline_args=dict(nchain="2 100 0.1 100",
                                                            gtr=True,
                                                            cat=True)):
        self.id = str(random.randint(10000,99999))+str(time.time())
        self.method_name=method_name
        self.program_name=program_name
        self.cline_args = cline_args
        self.trimmed_alignments = pj.trimmed_alignments
        if not alns == 'all':
            self.trimmed_alignments = {}
            for aln_name in alns:
                if aln_name in pj.trimmed_alignments.keys():
                    self.trimmed_alignments[aln_name] = pj.trimmed_alignments[aln_name]
        self.aln_input_strings = {}
        self.command_lines = {}
        self.timeit = [time.asctime()]
        self.platform = []
        self.cmd = cmd
        self.keepfiles = keepfiles
        if cmd == 'default':
            self.cmd = pj.defaults['pb']
            
        for trimmed_alignment in self.trimmed_alignments.keys():
            self.command_lines[trimmed_alignment] = [write_pb_cline(self, pj, trimmed_alignment)]
            print self.command_lines[trimmed_alignment][0]
    
    def __str__(self):

        aln_string = ''
        for n in self.trimmed_alignments.keys():
            aln_string += n+','
        aln_string = aln_string[:-1]
        command_lines = ''
        for i in self.command_lines.keys():
            command_lines += i+': '+str(self.command_lines[i])+'\n'
        date = str(self.timeit[0])
        
        execution = '[This was not executed yet]'
        if len(self.timeit) > 1:
            execution = str(self.timeit[3])
        
        prog = '[This was not executed yet]'
        plat = '[This was not executed yet]'
        progref = '[This was not executed yet]'
        
        if len(self.platform) > 0:
            plat = str(self.platform).replace(",",'\n').replace(']','').replace("'",'').replace('[','')
            if self.platform[-1].startswith('Program reference:'):
                prog = self.platform[-2].split(':')[1].strip().rstrip()
                progref = self.platform[-1].split('Program reference:')[1].strip().rstrip()
                
        output = str("PbConf named %s with ID %s\n"+         
                "Alignments: %s \n"+       
                "Created on: %s\n"+
                "Commands:\n"+
                "%s\n"+
                "Environment:\n"+    
                "%s\n"+
                "execution time:\n"+
                "%s")%(self.method_name, str(self.id), aln_string, date, command_lines, plat, execution)
        
        output += RaxmlConfMethodsSection%(type_to_single_line_str(self.trimmed_alignments.keys()),
                                            prog, progref)
        
        return output
            
from pylab import *
import random

def draw_boxplot(dictionary, y_axis_label, figs_folder): #'locus':[values]
    import numpy as np
    import matplotlib.pyplot as plt
    items = dictionary.items()
    items.sort()
    data = [locus[1] for locus in items]
        
    fig, ax1 = plt.subplots()
    fig.set_size_inches(0.3*len(data),10)
    plt.subplots_adjust(top=0.99, bottom=0.3)

    #bp = plt.boxplot(data, widths=0.75, patch_artist=True)
    bp = plt.boxplot(data, patch_artist=True)
    
    for box in bp['boxes']:
    # change outline color
        box.set( color='black', linewidth=1)
        
    # change fill color
        box.set( facecolor = 'red', alpha=0.85 )
        
    # change color, linestyle and linewidth of the whiskers
    for whisker in bp['whiskers']:
        whisker.set(color='gray', linestyle='solid', linewidth=2.0)

    # change color and linewidth of the caps
    for cap in bp['caps']:
        cap.set(color='gray', linewidth=2.0)

    # change color and linewidth of the medians
    for median in bp['medians']:
        #median.set(color='#b2df8a', linewidth=2)
        median.set(color='white', linewidth=2)

    # change the style of fliers and their fill
    for flier in bp['fliers']:
        flier.set(marker='o', color='#e7298a', alpha=0.5)
    
    # Add a light horizontal grid to the plot
    ax1.yaxis.grid(True, linestyle='-', which='major', color='lightgrey',
              alpha=0.7)
    
    # Hide these grid behind plot objects
    ax1.set_axisbelow(True)
    
    #set title and axis labels
    #ax1.set_title('Sequence length distribution per locus\n', size=18)
    
    xlabels = [locus[0] for locus in items]
    
    xticks(range(1,len(data)+1), xlabels, size=14, rotation='vertical')
    #subplots_adjust(left=0.3, bottom=0.8)
    
    ax1.set_ylabel(y_axis_label, size=18)
    
    name = str(random.randint(1000,2000))
    if figs_folder=='inline':
        fig.show()
    else:
        fig.savefig(figs_folder + '/' + name +'.png')
        close('all')
    return figs_folder + '/' + name+'.png'
    
#################################################################################
def report_methods(pj, figs_folder, output_directory, size='small',
                   compare_trees=[], compare_meta=None, trees_to_compare='all',
                   unrooted_trees=False, mp_root=True):
#################################################################################
        """
        Main HTML reporting function. This function iterates over the 
        Project's attributes and generate appropriate html formated report lines.
        
        pj -               The Project object
        
        figs_folder -      The directory to which tree figures were saved. This is
                           specified in the annotate Project method
                           
        output_directory - The directory to which this report will be written.It 
                           can be inherited from the publish function which uses 
                           this function.
        """
        
        #========================================================================
        #                                   HEAD
        #========================================================================
        
        # Checking if 'report_methods' was called by 'publish' in order to infer
        # which folders to create and what errors to raise for existing 
        # directories. If we were called by 'publish' we should allow 
        # output_directory to exist because 'publish' has just created it.
        # 'publish' would have also raised an error if it existed before.
        curframe = inspect.currentframe()
        calframe = inspect.getouterframes(curframe, 2)
        callername = calframe[1][3]
        print "reporter was called by "+str(callername)
        # Here we manage mkdir and error raising for preexisting 
        # output_directory, based of the caller function
        if os.path.isdir(output_directory) and str(callername) != 'publish':
            raise RuntimeError('%s already exists'%output_directory) 
        else:
            # This will make output_directory if it doesnt exist and
            # will add a subdirectory 'files' which will store the
            # stylesheet file and the figure files
            os.makedirs('%s/files'%output_directory)
        
        # This will contain the text that points to the stylesheet file
        # I have named it rp.css. We need to add something that writes
        # it here. Or that get it from the source directory and places it
        # in /files
        
        css_line = '<link rel="stylesheet" type="text/css" href="Bootstrap.css">'

        # This list will contain the report lines/ tables as values. We will append
        # each new report line to it
        report_lines = ['<html>','<head>',css_line,'<h1>']
    
        # Main report title, will print the time in which the 'Project' object
        # was initiated (not the time when the report was made as in previous 
        # versions.
        head = 'reprophylo analysis from '+pj.starttime
        

        report_lines.append(head)
        report_lines += ['</h1>','</head>','<body>','']
                
        #========================================================================
        #                                   BODY
        #========================================================================

        
        #############################  section 1:   DATA  #######################
        
        if pj.user:
            report_lines += ['<h2>','User Info','</h2>', '']
            for item in pj.user:
                report_lines += ['<strong>%s: </strong>'%item[0], str(item[1])]
            report_lines += ['']
        
        report_lines += ['<h2>','Data','</h2>', '']
        
        print "starting species table"
        # Species over loci table
        #------------------------------------------------------------------------
        title = 'species representation in sequence data'.title()
        report_lines += ('<h3>', title, '</h3>', '')
        
        
        # This will write a CSV file of the times each locus occures for each 
        # species. The species are sorted in alphabetical order. The CSV will 
        # be writen to 'outfile_name'. After 'outfile_name' is processed it will
        # be deleted.The CSV will be read with the csv module and the resulting
        # list of lists will be made into a table using HTML and than added to 
        # 'report_lines'. The species an gene counts are made by 'species_vs_loci'
        # based on the sequences record objects (biopython) found in pj.records.
        # pj.records is a list of SeqRecord objects
        outfile_name= str(random.randint(1000,2000))
        pj.species_vs_loci(outfile_name)
        with open(outfile_name, 'rb') as csvfile:
            sp_vs_lc = list(csv.reader(csvfile, delimiter='\t', quotechar='|'))
            report_lines += ['',
                             HTML.table(sp_vs_lc[1:], header_row=sp_vs_lc[0]),
                             '']
            # The following writes a pre text of the same thing
            #field_sizes = []
            #for i in range(len(sp_vs_lc[0])):
            #    lengths = []
            #    for row in sp_vs_lc:
            #        lengths.append(len(row[i]))
            #    field_sizes.append(max(lengths))
            #for row in sp_vs_lc:
            #    string = ''
            #    for i in range(len(row)):
            #        string += row[i].ljust(field_sizes[i]+3)
            #    report_lines.append(string)
        
        os.remove(outfile_name)
        
        if not size == 'small':
            print "starting sequence statistics plots"        
            # Sequence statistic plots
            #------------------------------------------------------------------------
            title = 'Sequence statistic plots'.title()
            report_lines += ('<h3>', title, '</h3>', '')

            # This will plot 4 box plot figures representing the distribution of seq
            # length, GC content, %ambiguity in nuc and prot seqs for each locus.
            if len(pj.records_by_locus.keys())>0:

                # This will determine the with of the figure, 0.5' per locus
                scale = str(len(pj.records_by_locus.keys())*0.5)

                # This will make a list of seq length for each locus. Seq length are calced
                # using the record.seq in 'pj.records_by_locus'. 'pj.records_by_locus is a
                # dict with loci names as keys, and lists of SeqReocrd objects as values
                lengths_dict = {}
                for locus_name in pj.records_by_locus.keys():
                    lengths_dict[locus_name] = []
                    for record in pj.records_by_locus[locus_name]:
                        lengths_dict[locus_name].append(len(record.seq))

                # This draws a box plot of sequence length distributions and puts a png in 
                # the 'files' directory.
                if not size == 'small':
                    fig_filename = draw_boxplot(lengths_dict, 'Seq length (bp)', '%s/files'%output_directory)


                # Distribution of sequence lengths
                #---------------------------------------------------------------------
                title = 'Distribution of sequence lengths'
                report_lines += ( '<h4>', title, '</h4>',  '')

                # This will write the img tag for the seq length boxplot in the report html
                # The src attribute is the png file path. The commented lines are an alternative
                # making an embeded figure.
                if not size=='small' and os.path.isfile(fig_filename):
                    #data_uri = open(fig_filename, 'rb').read().encode('base64').replace('\n', '')
                    #img_tag = '<img height=400 width='+scale+' src="data:image/png;base64,{0}">'.format(data_uri)
                    img_tag = '<img src="%s">'%(fig_filename.partition('/')[-1])
                    report_lines.append(img_tag)
                    #os.remove(fig_filename)


                # This will make GC content, nuc_degen_prop and prot_degen_prop png file,
                # will put them in the /files subdirectory and will write the html sections
                # for them, including img tags. All three params are feature qualifiers of
                # SeqRecord objects found in pj.records, which is a list.
                stats_to_plot = ('GC_content', 'nuc_degen_prop', 'prot_degen_prop')
                if size=='small':
                    stats_to_plot = ()
                for stat in stats_to_plot:
                    # This will make a dict with loci as keys and a list of stat values as
                    # dict values.
                    stat_dict = {}
                    ylabel = 'GC ontent (%)'
                    if not stat == 'GC_content':
                        ylabel = 'Ambiguous positions (prop)'
                    for locus_name in pj.records_by_locus.keys():
                        stat_dict[locus_name] = []
                        for i in pj.records_by_locus[locus_name]:
                            for record in pj.records:
                                for feature in record.features:
                                    if feature.qualifiers['feature_id'][0] == i.id:
                                        if stat in feature.qualifiers.keys():
                                            stat_dict[locus_name].append(float(feature.qualifiers[stat][0]))

                    # This will make the boxplot png and will put in in the /files subdirectory
                    fig_filename = draw_boxplot(stat_dict, ylabel, '%s/files'%output_directory)

                    # Distribution of stat
                    #---------------------------------------------------------------------
                    title = 'Distribution of sequence statistic \"'+stat+'\"'
                    report_lines += ( '<h4>', title, '</h4>', '')

                    # This will make the img tag using the png path as src. The commented lines are an alternative
                    # making an embeded image
                    if os.path.isfile(fig_filename):
                        #data_uri = open(fig_filename, 'rb').read().encode('base64').replace('\n', '')
                        img_tag = '<img src="%s">'%(fig_filename.partition('/')[-1])
                        #img_tag = '<img height=400 width='+scale+' src="data:image/png;base64,{0}">'.format(data_uri)
                        report_lines.append(img_tag)
                        #os.remove(fig_filename)
                
        print "starting concatenations"
        # Description of data concatenations
        #------------------------------------------------------------------------
        title = 'Description of data concatenations'.title()
        report_lines += ('<h3>', title, '</h3>', '')
        
        # filter out concatenation objects that were not used to build a concatenation
        # by taking only concat names that are in the keys of pj.trimmed_alignments.
        # pj.trimmed_alignments is a dict with alignment names as keys and list 
        # lists containing an alignment object (biopython) and an alignment string as
        # values.
        composed_concatenations = []
        for c in pj.concatenations:
            if c.name in pj.trimmed_alignments.keys():
                composed_concatenations.append(c)

        for c in composed_concatenations:
            
            title = ('content of concatenation \"' + c.name + '\"').title()
            report_lines += ('<h4>', title, '</h4>', '')
            #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            
            # This will write the concatenation attribute
            report_lines.append('Rules for  \"' + c.name + '\":')
            rule_1 = 'OTUs must have the loci: '
            for locus in c.otu_must_have_all_of:
                rule_1 += locus + ', '
            report_lines.append(rule_1)
            rule_2 = '<pre>'
            for group in c.otu_must_have_one_of:
                rule_2 += '</pre>OTUs must have at least one of the following loci: \n<pre>'
                rule_2 += str(group).replace('[','').replace(']','') +'</pre>\n'
            report_lines += (rule_2, '')
            
            # This are the otus and loci names in the concatenation. c.feature_id_dict
            # is a dict with the otius as keys and dicts as values. Each of these dicts
            # have the loci as keys and the feature id as value
            otus = c.feature_id_dict.keys()
            loci = [locus.name for locus in c.loci]
            
            # This is the table's header
            table_lines = [['','']+[locus.name for locus in c.loci]]
            
            
            # Commented lines are an alternative way to write the table as 
            # pre text
            #otus_max_length = max([len(i) for i in otus])+33
            #loci_columns_max_length = []
            
            #for locus in loci:
            #    lengths = [len(locus)]
            #    for otu in otus:
            #        if locus in c.feature_id_dict[otu].keys():
            #            lengths.append(len(c.feature_id_dict[otu][locus]))
            #        else:
            #           lengths.append(0)
            #    loci_columns_max_length.append(max(lengths)+3)
                
            #concat_header = ''.ljust(otus_max_length)
            #for i in range(len(loci)):
            #    concat_header += loci[i].ljust(loci_columns_max_length[i])
            #report_lines += (concat_header, '~'*len(concat_header))
                
            
            # This will write the table
            for otu in otus:
                
                otu_species = ''
                for locus in loci:
                    if locus in c.feature_id_dict[otu].keys():
                        feature_qualifiers = get_qualifiers_dictionary(pj, c.feature_id_dict[otu][locus])
                        if 'source_organism' in feature_qualifiers.keys():
                            otu_species = feature_qualifiers['source_organism']
                otu_line = [otu, otu_species]    
                #concat_line = (otu+' '+otu_species).ljust(otus_max_length)
                for i in range(len(loci)):
                    if loci[i] in c.feature_id_dict[otu].keys():
                        #concat_line += c.feature_id_dict[otu][loci[i]].ljust(loci_columns_max_length[i])
                        otu_line.append(c.feature_id_dict[otu][loci[i]])
                    else:
                        #concat_line += ''.ljust(loci_columns_max_length[i])
                        otu_line.append('')
                table_lines.append(otu_line)
            
                #report_lines.append(concat_line)
        
            report_lines.append(HTML.table(table_lines[1:], header_row=table_lines[0]))  
        
        #############################  section 2:   METHODS  #######################
        
        # This section prints some attributes of each of the 'Conf' objects used.
        # The 'Conf' objects are found in a list called pj.used_methods. 
        # In an unpickled 'Project' object, the 'Conf' objects are replaced by lists of 
        # strings describing the attributes because the objects themselves do not
        # pickle well. The formating of the list representations when they are printed
        # still needs some beautification. Also, I plan a 'revive_methods' func to turn
        # them back to 'Conf' objects that can be rerun.
        report_lines += ['', '<h2>','Methods','</h2>', '']
        
        print "starting methods"
        for method in pj.used_methods:
            # This will print list representations of the 'Conf' objects
            title = str(pj.used_methods[method]).split('\n')[0]
            report_lines += ('', '<h4>', title, '</h4>','')
            report_lines += ('<pre>',str(pj.used_methods[method]),'</pre>')
       
                
        report_lines += ['',''] 
        
        #############################  section 3:   RESULTS  #######################
        
        report_lines += ['', '<h2>','Results','</h2>', '']
        print "starting alignment statistics"
        # Global alignmnet statistics
        #------------------------------------------------------------------------
        title = 'Global alignmnet statistics'.title()
        report_lines += ('<h3>', title, '</h3>', '')
        
        
        # This prints things like num of unique seqs and num of parsimony informative
        # cloumns. Takes the info from 'pj.aln_summaries' which is a list of strings.
        # Alignment length: 1566
        #Number of rows: 94
        #Unique sequences: 87
        #Average gap prop.: 0.488445
        #Variable columns: 1045
        #Parsimony informative: 402
        #Undetermined sequences: 0
        
        
        if len(pj.aln_summaries)>0:
            report_lines += [('<pre>Name=Alignment name\n'+
                              'NumPos=Alignment length\n'+
                              'NumSeq=Number of sequences\n'+
                              'Unique=Number of unique sequences\n'+
                              'GapProp=Average gap proportion\n'+
                              'VarCols=Total variable positions\n'+
                              'ParsInf=Parsimony informative positions\n'+
                              'UnSeqs=Undetermined sequences (mostly/only gaps)\n'+
                              'UnSeqsCutoff=Length cutoff which defines undetermined\n</pre>')]
            T = [['Name','NumPos','NumSeq','Unique','GapProp','VarCols','ParsInf','UnSeqs','UnSeqsCutoff']]
            comments = []
            for summary in pj.aln_summaries:
                line = []
                for i in summary.splitlines():
                    try:
                        line.append(i.split(': ')[1])
                    except:
                        comments.append(i)
                T.append(line)
            report_lines += ['',
                             HTML.table(T[1:], header_row=T[0]),
                             '','<pre>']+comments+['</pre>']
        else:
            report_lines += ['','No sequence alignments in this Project','']
                
        
        # Per position alignmnet statistics
        #------------------------------------------------------------------------
        title = 'per position alignmnet statistics'.title()
        report_lines += ('<h3>', title, '</h3>', '')
        if len(pj.alignments.keys())>0 and not size == 'small':                    
            title = 'Alignment statistics before trimming'
            report_lines += ('', '<h4>', title, '</h4>', '')
            report_lines += ['<h4>','Trimal\'s Residue Similarity Score (-scc)','</h4>', '']
            
            # draw_trimal_scc(project, num_plots_in_raw, output_dir...
            # 'trimmed' determines if it will analyse trimmed or raw alignments
            # alignments are taken from pj.alignments or pj.trimmed_alignments which are dictionaries
            # with alignment names as keys and lists containing an aln object and al string as values
            
            # scc on raw alignments
            fig_file = draw_trimal_scc(pj, 2, '%s/files'%output_directory, trimmed=False)
            if os.path.isfile(fig_file):
                    #data_uri = open(fig_file, 'rb').read().encode('base64').replace('\n', '')
                    #img_tag = '<img src="data:image/png;base64,{0}">'.format(data_uri)
                    img_tag = '<img src="%s">'%(fig_file.partition('/')[-1])
                    report_lines.append(img_tag)
                    #os.remove(fig_file)
            report_lines += [ '<h4>','Trimal\'s column gap gcore (-sgc)','</h4>', '']
            
            # sgc on raw alignments
            fig_file = draw_trimal_scc(pj, 2, '%s/files'%output_directory, trimmed=False, alg='-sgc')
            if os.path.isfile(fig_file):
                    #data_uri = open(fig_file, 'rb').read().encode('base64').replace('\n', '')
                    #img_tag = '<img src="data:image/png;base64,{0}">'.format(data_uri)
                    img_tag = '<img src="%s">'%(fig_file.partition('/')[-1])
                    report_lines.append(img_tag)
                    #os.remove(fig_file)
        else:
            report_lines += ['No alignments or too many alignments in this project','']
            
        if len(pj.trimmed_alignments.keys())>0 and not size=='small':          
            title = 'Alignment statistics after trimming'
            report_lines += ('', '<h4>', title, '</h4>', '')
            report_lines += ['<h4>','"Trimal\'s Residue Similarity Score (-scc)','</h4>', '']
            
            # scc on trimmed alignments
            fig_file = draw_trimal_scc(pj, 2, '%s/files'%output_directory, trimmed=True)
            if os.path.isfile(fig_file):
                    #data_uri = open(fig_file, 'rb').read().encode('base64').replace('\n', '')
                    #img_tag = '<img src="data:image/png;base64,{0}">'.format(data_uri)
                    img_tag = '<img src="%s">'%(fig_file.partition('/')[-1])
                    report_lines.append(img_tag)
                    #os.remove(fig_file)
            report_lines += ['<h4>','Trimal\'s column gap gcore (-sgc)','</h4>',  '']
            
            # sgc on trimmed alignments
            fig_file = draw_trimal_scc(pj, 2, '%s/files'%output_directory, trimmed=True, alg='-sgc')
            if os.path.isfile(fig_file):
                    #data_uri = open(fig_file, 'rb').read().encode('base64').replace('\n', '')
                    #img_tag = '<img src="data:image/png;base64,{0}">'.format(data_uri)
                    img_tag = '<img src="%s">'%(fig_file.partition('/')[-1])
                    report_lines.append(img_tag)
                    #os.remove(fig_file)
        else:
            report_lines += ['No trimmed alignments, or too many, in this project','']
        
        print "starting RF matrix(ces)"
        for rf_type in compare_trees:
            title = 'Robinson-Foulds distances (%s)'%rf_type.title()
            report_lines += ('<h3>', title, '</h3>', '')

            if len(pj.trees.keys())>1:
                try:
                    RF_filename, legend = calc_rf(pj, '%s/files'%output_directory,
                                                  rf_type=rf_type, meta=compare_meta, trees=trees_to_compare,
                                                  unrooted_trees=unrooted_trees, mp_root=mp_root)
                    scale = str(len(legend)*60)
                    if os.path.isfile(RF_filename):
                            #data_uri = open(RF_filename, 'rb').read().encode('base64').replace('\n', '')
                            #img_tag = '<img height='+scale+' width='+scale+' src="data:image/png;base64,{0}">'.format(data_uri)
                            img_tag = '<img src="%s">'%(RF_filename.partition('/')[-1])
                            report_lines.append(img_tag)
                            #os.remove(RF_filename)

                    report_lines+=['<h3>Legend<h3>','']
                    report_lines += [HTML.table( legend[1:], 
                                               header_row=legend[0])]
                    report_lines.append('')
                except:
                    report_lines += ['Skipping RF distance calculation']

            else:
                report_lines += ['Less than two trees in this Project','']

                
        #############################  section 4:   TREES  #######################
        print "reporting trees"
        report_lines += ['', '<h2>','Trees','</h2>', '']
        
        for tree in pj.trees.keys():
            report_lines += ('<h2>'+tree.split('@')[0]+'</h2>',
                             '<h3>Alignment method: '+tree.split('@')[1]+'</h3>',
                             '<h3>Trimming method: '+tree.split('@')[2]+'</h3>',
                             '<h3>Tree method: '+tree.split('@')[3]+'</h3>',
                             '<pre style="white-space:normal;">',
                             'Tree Method ID: '+pj.trees[tree][0].get_leaves()[0].tree_method_id,'</pre>')
            
            report_lines += ('<h3>newick format</h3>','','<pre style="white-space:normal;">',pj.trees[tree][0].write(),'</pre>','')
            report_lines += ('<h3>nhx format</h3>','','<pre>',pj.trees[tree][1],'</pre>','','','','')
            
            
            
            if os.path.isfile(figs_folder+'/'+pj.trees[tree][0].get_leaves()[0].tree_method_id+'.png'):
                origin = figs_folder+'/'+pj.trees[tree][0].get_leaves()[0].tree_method_id+'.png'
                dest = '%s/files/%s'%(output_directory, pj.trees[tree][0].get_leaves()[0].tree_method_id+'.png')
                shutil.copyfile(origin, dest)
                #data_handle = open(figs_folder+'/'+pj.trees[tree][0].get_leaves()[0].tree_method_id+'.png','rb')
                #data_uri = data_handle.read().encode('base64').replace('\n', '')
                #data_handle.close()
                #img_tag = '<img width=500 src="data:image/png;base64,{0}">'.format(data_uri)
                img_tag = '<img width=500 src="%s">'%(dest.partition('/')[-1])
                report_lines.append(img_tag)

                
        report_lines.append('</body>')
        report_lines.append('</html>')
        
        lines = []
        
        for line in report_lines:
            if '<' in line:
                lines.append(line)
            else:
                lines.append(line.replace('\n','<br>')+'<br>')
        
        
        return lines
    
        
def pickle_pj(pj, pickle_file_name, track=True):
    import os
    if os.path.exists(pickle_file_name):
        #print "DEBUG pickle_pj: %s"%str(os.stat(pickle_file_name).st_ctime)
        #print "DEBUG pickle_pj: %s"%str(os.stat(pickle_file_name).st_size)
        os.remove(pickle_file_name)
    import cloud.serialization.cloudpickle as pickle
    output = open(pickle_file_name,'wb')
    pickle.dump(pj, output)
    output.close()
    #print "DEBUG pickle_pj: %s"%str(os.stat(pickle_file_name).st_ctime)
    #print "DEBUG pickle_pj: %s"%str(os.stat(pickle_file_name).st_size)
    if __builtin__.git and track:
        import rpgit
        rpgit.gitAdd(pickle_file_name)
        comment = "A pickled Project from %s" % time.asctime()
        rpgit.gitCommit(comment) 
        
    return pickle_file_name
    
def unpickle_pj(pickle_file_name, git=True):
    import cloud.serialization.cloudpickle as pickle
    pickle_handle = open(pickle_file_name, 'rb')
    
    # fix some attr: add git_log if None, turn Confs to strings.
    pkl_pj = pickle.pickle.load(pickle_handle)
    new_pj = Project(pkl_pj.loci, git=False)
    
    # These do not need fixing
    attr_names = [ 'alignments',
                   'aln_summaries',
                   'concatenations',
                   'defaults',
                   'git_log',
                   'pickle_name',
                   'records',
                   'records_by_locus', 
                   'sets',
                   'starttime',
                   'trees',
                   'trimmed_alignments',
                   'user',
                 ]
    
    # Move the content into the new pj, Add git_log if missing
    for attr_name in attr_names:
        try:
            setattr(new_pj,attr_name,getattr(pkl_pj,attr_name))
        except:
            warnings.warn('Upgrading Project to v1')
    
    # upgrade used methods to dict, if list.
    if isinstance(pkl_pj.used_methods, list) and len(pkl_pj.used_methods) == 0:
        pkl_pj.used_methods = {}

    elif isinstance(pkl_pj.used_methods, list) and len(pkl_pj.used_methods) > 0:
        temp = {}
        for m in pkl_pj.used_methods:
            if isinstance(m, basestring):
                name = m.split()[1]
                temp[name] = m
            elif not isinstance(m, basestring):
                temp[m.method_name] = m
        pkl_pj.used_methods = temp
    
    # Turn Confs to strings
    for i in pkl_pj.used_methods:
        if isinstance(pkl_pj.used_methods[i], basestring):
            new_pj.used_methods[i] = pkl_pj.used_methods[i]
        else:
            new_pj.used_methods[i] = str(pkl_pj.used_methods[i])
    if git:
        start_git(new_pj)
    return new_pj


def revert_pickle(pj, commit_hash):
    pickle_filename = pj.pickle_name
    if not os.path.exists(pickle_filename):
        raise RuntimeError('Cannot find %s. Has the pickle been moved'%pickle_filename)
    cmd = 'git checkout %s -- %s'%(commit_hash, pickle_filename)
    pipe = sub.Popen(cmd, shell=True, stdout = sub.PIPE,stderr = sub.PIPE )
    (out, error) = pipe.communicate()
    print 'Git STDOUT: %s'%str(out)
    print 'Git STDERR: %s'%str(error)
    new_pj = unpickle_pj(pickle_filename)
    new_pj.git_log += "<<<<\nThe pickle was reverted to commit %s\nSTDOUT:\n%s\nSTDERR:%s\n>>>>\n"%(commit_hash,
                                                                                                   str(out),
                                                                                                   str(error))
    return new_pj
                                                                                                    
    


def publish(pj, folder_name, figures_folder, size='small',
            compare_trees=[], compare_meta=None, trees_to_compare='all',
            unrooted_trees=False):
    
    import os, time
    folder = None
    zip_file = None
    if folder_name.endswith('.zip'):
        zip_file = folder_name
        folder = folder_name[:-4]
    else:
        folder = folder_name
        zip_file = folder_name + '.zip'
    print "checking if file exists"
    if os.path.exists(folder) or os.path.exists(zip_file):
        raise IOError(folder_name + ' already exists')
    
    os.makedirs(folder)
    pj.write(folder+'/tree_and_alns.xml','phyloxml')
    
    #os.mkdir(folder+'/fasta_alignments')
    #file_names = pj.write_alns(id=['feature_id','original_id'])
    #for f in file_names:
    #    os.rename(f, "%s/fasta_alignments/%s"%(folder,f))
        
    #os.mkdir(folder+'/trimmed_fasta_alignments')
    #file_names = pj.write_alns(id=['feature_id','original_id'])
    #for f in file_names:
    #    os.rename(f, "%s/trimmed_fasta_alignments/%s"%(folder,f))
        
    from glob import glob
    from shutil import copyfile
    notebooks = glob('*.ipynb')
    for n in notebooks:
        copyfile(n, '%s/%s'%(folder,n))
        
    pj.write(folder+'/sequences_and_metadata.gb','genbank')
    report = open(folder+'/report.html','wt')
    lines = report_methods(pj, figures_folder, folder_name, size,
                           compare_trees=compare_trees, compare_meta=compare_meta,
                          trees_to_compare=trees_to_compare, unrooted_trees=unrooted_trees)
    for line in lines:
        report.write(line + '\n')
    report.close()

    #'report_lines' is now taking care of puting the figures in the zip folder, within /files
    #for tree in pj.trees.keys():
    #    if os.path.isfile(figures_folder+'/'+pj.trees[tree][0].get_leaves()[0].tree_method_id+'.png'):
    #        from shutil import copyfile
    #        copyfile(figures_folder+'/'+pj.trees[tree][0].get_leaves()[0].tree_method_id+'.png',
    #                 folder+'/'+pj.trees[tree][0].get_leaves()[0].tree_method_id+'.png')
            
         
    print "pickling"
    pickle_name = time.strftime("%a_%d_%b_%Y_%X", time.gmtime())+'.pkl'
    pickle_pj(pj, folder + '/' + pickle_name)
    
    hndl = open("%s/Bootstrap.css"%folder, 'wt')
    hndl.write(css.css())
    hndl.close()
    
    import zipfile, shutil
    print "archiving"
    zf = zipfile.ZipFile(zip_file, "w")
    for dirname, subdirs, files in os.walk(folder):
        zf.write(dirname)
        for filename in files:
            zf.write(os.path.join(dirname, filename))
    zf.close()
    shutil.rmtree(folder)
    print "report ready"

##################################################
if False:
    """ Robinson Foulds Pairwise tree distances"""
##################################################


# These will allow to do an RF distance to dendropy, where branch lengths
# are considered, but also to correct the branch length by the tree length
# and thus compare trees with very different evolutionary rates.

def get_tree_length(t):
    tree_length = 0
    for n in t.traverse():
        if not n.dist == 1: # ete puts 1 if there is no blen
            tree_length += n.dist
    return tree_length
    
def correct_branch_length_by_tree_length(branch_length, tree_length):
    return branch_length/float(tree_length)

def get_corrected_blen_dif(cor_blen1, cor_blen2):
    return abs(cor_blen1-cor_blen2)

def get_corrected_blen_dif_s(cor_blen1, cor_blen2):
    from math import pow
    return pow((cor_blen1-cor_blen2),2)
    
def flatten_to_strings(listOfLists):
    """Flatten a list of (lists of (lists of strings)) for any level 
    of nesting"""
    result = []

    for i in listOfLists:
        # Only append if i is a basestring (superclass of string)
        if isinstance(i, basestring):
            result.append(i)
        # Otherwise call this function recursively
        else:
            result.extend(flatten_to_strings(i))
    return result
    
def get_corrected_blen_rf(t1, t2, unrooted_trees=False):
    """
    >>> t1 = Tree("(A:0.3,(B:0.226,((C:0.784,D:0.159):0.759,(e:0.03,f:0.1):0.25)):0.3):0.1;")
    >>> t2 = Tree("(A:0.3,(D:0.226,((C:0.784,B:0.159):0.759,(e:0.03,f:0.1):0.25)):0.3):0.1;")
    
    Trees are the same:
    >>> T1 = T2 = t1
    >>> get_corrected_blen_rf(T1,T2)
    0.0
    
    Trees are different:
    >>> get_corrected_blen_rf(t1,t2)
    0.504654255319149
    
    Dendropy calc of these trees' RF:
    >>> dt1 = dendropy.Tree.get_from_string(t1.write(), schema="newick")
    >>> dt2 = dendropy.Tree.get_from_string(t2.write(), schema="newick")
    >>> dt1.robinson_foulds_distance(dt2)
    3.652
    
    They are the same and blengths are 25% shorter in second tree. This function gives almost 0
    >>> t1 = Tree("(A:0.3,(B:0.226,((C:0.784,D:0.159):0.759,(e:0.03,f:0.1):0.25)):0.3):0.1;")
    >>> t1_prop = Tree("(A:0.3,(B:0.226,((C:0.784,D:0.159):0.759,(e:0.03,f:0.1):0.25)):0.3):0.1;")
    >>> for n in t1_prop.traverse():
    ...    if not n.dist == 1:
    ...        n.dist = 0.75*n.dist
    >>> get_corrected_blen_rf(t1,t1_prop)
    6.938893903907228e-17
    
    >>> str1 = "(A:1.1,((B:1.1,C:1.1):1.1,(D:1.1,E:1.1):1.1):1.1);"
    >>> str2 = "(A:2.2,((B:2.2,C:2.2):2.2,(D:2.2,E:2.2):2.2):2.2);"
    >>> get_corrected_blen_rf(Tree(str1), Tree(str2))
    0.0
    
    Dendropy gives almost 1 for the same trees:
    >>> dt1 = dendropy.Tree.get_from_string(t1.write(), schema="newick")
    >>> dt1_prop = dendropy.Tree.get_from_string(t1_prop.write(), schema="newick")
    >>> dt1.robinson_foulds_distance(dt1_prop)
    0.7269999999999999
    
    So dendropy RF reflects tree length differences, this function cleans them out.
    
    One branch length changes in one tree, topology stays the same:
    >>> t1 = Tree("(A:0.3,(B:0.226,((C:0.001,D:0.159):0.759,(e:0.03,f:0.1):0.25)):0.3):0.1;")
    >>> dt1 = dendropy.Tree.get_from_string(t1.write(), schema="newick")
    >>> get_corrected_blen_rf(t1,t1_prop)
    0.23503571001673468
    
    Same trees (same topology, different branch-lengths) in ete's RF:
    >>> rf, max_rf, common_leaves, parts_t1, parts_t2 = t1.robinson_foulds(t1_prop)
    >>> rf
    0
    
    The raw rf calc in ete (the starting point of this func) does not reflect any 
    branch length differences, proportional or otherwise.
    """
    rf, max_rf, common_leaves, parts_t1, parts_t2 = t1.robinson_foulds(t2, unrooted_trees=unrooted_trees)
    #print 'DEBUG:', parts_t1
    #print 'DEBUG:', parts_t2
    tree_length1 = get_tree_length(t1)
    tree_length2 = get_tree_length(t2)
    distance = 0
    for part in parts_t1:
        if not part in parts_t2:
            raw_blen = t1.get_common_ancestor(flatten_to_strings(part)).dist
            if not raw_blen==1.0:
                distance += correct_branch_length_by_tree_length(raw_blen, tree_length1)
        elif part in parts_t2:
            #print 'DEBUG:', part
            raw_blen1 = t1.get_common_ancestor(flatten_to_strings(part)).dist
            raw_blen2 = t2.get_common_ancestor(flatten_to_strings(part)).dist
            if raw_blen1==1:
                raw_blen1 = 0
            if raw_blen2==1:
                raw_blen2 = 0
            cor_blen1 = correct_branch_length_by_tree_length(raw_blen1, tree_length1)
            cor_blen2 = correct_branch_length_by_tree_length(raw_blen2, tree_length2)
            distance += get_corrected_blen_dif(cor_blen1, cor_blen2)
    for part in parts_t2:
        if not part in parts_t1:
            raw_blen = t2.get_common_ancestor(flatten_to_strings(part)).dist
            if not raw_blen==1:
                distance += correct_branch_length_by_tree_length(raw_blen, tree_length2)
            
    #print "Debug: distance: %s"%str(distance)
    return distance
            
def get_corrected_kuhner_felsenstein(t1, t2, unrooted_trees=False):
    """
    >>> t1 = Tree("(A:0.3,(B:0.226,((C:0.784,D:0.159):0.759,(e:0.03,f:0.1):0.25)):0.3):0.1;")
    >>> t2 = Tree("(A:0.3,(D:0.226,((C:0.784,B:0.159):0.759,(e:0.03,f:0.1):0.25)):0.3):0.1;")
    
    Trees are the same:
    >>> T1 = T2 = t1
    >>> get_corrected_kuhner_felsenstein(T1,T2)
    0.0
    
    Trees are different:
    >>> get_corrected_kuhner_felsenstein(t1,t2)
    0.1273379587058624
    
    Dendropy calc of these trees' RF:
    >>> dt1 = dendropy.Tree.get_from_string(t1.write(), schema="newick")
    >>> dt2 = dendropy.Tree.get_from_string(t2.write(), schema="newick")
    >>> dt1.robinson_foulds_distance(dt2)
    3.652
    
    They are the same and blengths are 25% shorter in second tree. This function gives almost 0
    >>> t1 = Tree("(A:0.3,(B:0.226,((C:0.784,D:0.159):0.759,(e:0.03,f:0.1):0.25)):0.3):0.1;")
    >>> t1_prop = Tree("(A:0.3,(B:0.226,((C:0.784,D:0.159):0.759,(e:0.03,f:0.1):0.25)):0.3):0.1;")
    >>> for n in t1_prop.traverse():
    ...    if not n.dist == 1:
    ...        n.dist = 0.75*n.dist
    >>> get_corrected_kuhner_felsenstein(t1,t1_prop)
    3.274080905458301e-33
    
    >>> str1 = "(A:1.1,((B:1.1,C:1.1):1.1,(D:1.1,E:1.1):1.1):1.1);"
    >>> str2 = "(A:2.2,((B:2.2,C:2.2):2.2,(D:2.2,E:2.2):2.2):2.2);"
    >>> get_corrected_kuhner_felsenstein(Tree(str1), Tree(str2))
    0.0
    
    Dendropy gives almost 1 for the same trees:
    >>> dt1 = dendropy.Tree.get_from_string(t1.write(), schema="newick")
    >>> dt1_prop = dendropy.Tree.get_from_string(t1_prop.write(), schema="newick")
    >>> dt1.robinson_foulds_distance(dt1_prop)
    0.7269999999999999
    
    So dendropy RF reflects tree length differences, this function cleans them out.
    
    One branch length changes in one tree, topology stays the same:
    >>> t1 = Tree("(A:0.3,(B:0.226,((C:0.001,D:0.159):0.759,(e:0.03,f:0.1):0.25)):0.3):0.1;")
    >>> dt1 = dendropy.Tree.get_from_string(t1.write(), schema="newick")
    >>> get_corrected_kuhner_felsenstein(t1,t1_prop)
    0.010930167133307158
    
    Same trees (same topology, different branch-lengths) in ete's RF:
    >>> rf, max_rf, common_leaves, parts_t1, parts_t2 = t1.robinson_foulds(t1_prop)
    >>> rf
    0
    
    The raw rf calc in ete (the starting point of this func) does not reflect any 
    branch length differences, proportional or otherwise.
    """
    from math import pow
    #print 'DEBUG: in kuhner_felsenstein'
    rf, max_rf, common_leaves, parts_t1, parts_t2 = t1.robinson_foulds(t2, unrooted_trees=unrooted_trees)
    #print 'DEBUG:', parts_t1
    #print 'DEBUG:', parts_t2
    tree_length1 = get_tree_length(t1)
    tree_length2 = get_tree_length(t2)
    distance = 0
    for part in parts_t1:
        if not part in parts_t2:
            raw_blen = t1.get_common_ancestor(flatten_to_strings(part)).dist
            if not raw_blen==1.0:
                distance += pow(correct_branch_length_by_tree_length(raw_blen, tree_length1),2)
        elif part in parts_t2:
            #print 'DEBUG:', part
            raw_blen1 = t1.get_common_ancestor(flatten_to_strings(part)).dist
            raw_blen2 = t2.get_common_ancestor(flatten_to_strings(part)).dist
            if raw_blen1==1:
                raw_blen1 = 0
            if raw_blen2==1:
                raw_blen2 = 0
            cor_blen1 = correct_branch_length_by_tree_length(raw_blen1, tree_length1)
            cor_blen2 = correct_branch_length_by_tree_length(raw_blen2, tree_length2)
            distance += get_corrected_blen_dif_s(cor_blen1, cor_blen2)
    for part in parts_t2:
        if not part in parts_t1:
            raw_blen = t2.get_common_ancestor(flatten_to_strings(part)).dist
            if not raw_blen==1:
                distance += pow(correct_branch_length_by_tree_length(raw_blen, tree_length2),2)
            
    #print "Debug: distance: %s"%str(distance)
    return distance            
        
def calc_rf(pj, figs_folder, rf_type='proportional',meta=None, mp_root=False, trees='all', unrooted_trees=False):
    """
    rf_types:
    topology: only topological diff
    branch-length: branch-length diffs
    proportional: branch-length proportion out of tree length diff
    deep-more-important: to-do
    """
    if len(pj.concatenations) > 0 and not meta:
        meta = pj.concatenations[0].otu_meta
    elif not meta:
        raise RuntimeError('RF calc does not know which meta to use to compare leaves')
        

    if trees == 'all':
        trees = pj.trees.keys()

    data = []

    for t1 in trees:
        line = []
        dupT1 = Tree(pj.trees[t1][0].write())
        for l in dupT1:
            for record in pj.records:
                for feature in record.features:
                    if feature.qualifiers['feature_id'][0] == l.name and meta in feature.qualifiers.keys():
                        l.name = feature.qualifiers[meta][0]
        if mp_root:
            R=dupT1.get_midpoint_outgroup()
            dupT1.set_outgroup(R)
        dupT1d = dendropy.Tree.get_from_string(dupT1.write(), schema="newick")
        for t2 in trees:
            dupT2 = Tree(pj.trees[t2][0].write())
            for l in dupT2:
                for record in pj.records:
                    for feature in record.features:
                        if feature.qualifiers['feature_id'][0] == l.name and meta in feature.qualifiers.keys():
                            l.name = feature.qualifiers[meta][0]
            if mp_root:
                R=dupT2.get_midpoint_outgroup()
                dupT2.set_outgroup(R)
            dupT2d = dendropy.Tree.get_from_string(dupT2.write(), schema="newick")
            if rf_type=='branch-length':
                rf = dupT1d.robinson_foulds_distance(dupT2d)
                line.append(rf) 
            elif rf_type=='topology':
                rf, max_rf, common_leaves, parts_t1, parts_t2 = dupT1.robinson_foulds(dupT2, unrooted_trees=unrooted_trees)
                line.append(rf/float(max_rf))   
            elif rf_type == 'proportional':
                #print 'DEBUG: in proportional'
                warnings.warn('proportional branch-distance: Trees must have the same taxa')
                line.append(get_corrected_kuhner_felsenstein(dupT1, dupT2, unrooted_trees=unrooted_trees))
            #to do
            #elif rf_type == deep_nodes more important:
            #    line.append(get_deep_important_rf(dupT1, dupT2))
        data.append(line)   
    
    row_labels = [str(i) for i in range(1, len(trees)+1)]
    column_labels = row_labels
    legend = [['#','LOCUS','ALIGNMENT METHOD','TRIMMING METHOD','TREE METHOD']]
    for i in trees:
        line = [str(trees.index(i)+1)]
        for val in i.split('@'):
            line.append(val)
        legend.append(line)
    fig, ax = plt.subplots()
    data = np.array(data)
    heatmap = ax.pcolor(data, cmap=plt.cm.Blues)

    # put the major ticks at the middle of each cell
    ax.set_xticks(np.arange(data.shape[0])+0.5, minor=False)
    ax.set_yticks(np.arange(data.shape[1])+0.5, minor=False)

    # want a more natural, table-like display
    ax.invert_yaxis()
    ax.xaxis.tick_top()

    ax.set_xticklabels(row_labels, minor=False, size=14, rotation='vertical')
    ax.set_yticklabels(column_labels, minor=False, size=14)
    #fig.set_size_inches(12.5,12.5)
    try:
        fig.colorbar(heatmap, cmap=plt.cm.Blues)
    except:
        pass
    name = str(random.randint(1000,2000))
    fig.savefig(figs_folder + '/' + name +'.png')
    close('all')
    return figs_folder + '/' + name+'.png', legend 

def draw_trimal_scc(pj, num_col, figs_folder, trimmed=False, alg = '-scc'):
    import pandas as pd
    import matplotlib.pyplot as plt
    import random, os
    from Bio import AlignIO
    
    # get the alignment objects
    #-------------------------#
    alignments = pj.alignments.items()
    if trimmed:
        alignments = pj.trimmed_alignments.items()
    num_alns = len(alignments)

    subplots_arrangement = []
    #-----------------------#
    num_rows = round(float(num_alns)/num_col)
    if num_rows < float(num_alns)/num_col:
        num_rows += 1
        
    fig = plt.figure(figsize=(10*num_col,2.3*num_rows), dpi=80, frameon = False)
    plt.subplots_adjust(hspace = 0.8)
    subplots_arrangement += [num_rows, num_col]

    
    #Calc with trimal and plot
    #------------------------#
    for i in range(1,num_alns+1):
        import subprocess as sub
        subplot_position = subplots_arrangement +[i]
        aln_name = alignments[i-1][0]
        aln_obj = alignments[i-1][1]
        name = str(random.randint(1000,2000))+'_'+aln_name+'_for_trimal_graph.fasta'
        AlignIO.write(aln_obj, name, 'fasta')
        stderr = open('stderr','wt')
        #stdout = os.popen(programpath+'trimal '+alg+' -in '+name)#.read()
        stdout = sub.Popen(programspath+"trimal "+alg+" -in " + name,
                       shell=True, stdout=sub.PIPE, stderr=stderr).stdout
        stderr.close()
        var = pd.read_table(stdout, sep='\t+', skiprows=3, engine='python')
        os.remove('stderr')
        if alg == '-scc':
            var.columns = ['position', 'variability']
        elif alg == '-sgc':
            var.columns = ['position', 'pct_gaps', 'gap_score']
        
        #Plot residue similarity figure, for nucleotides this is identity value 
        fig.add_subplot(subplot_position[0], subplot_position[1], subplot_position[2])
        if alg == '-scc':
            var.variability.plot(color='g',lw=2)
        elif alg == '-sgc':
            var.pct_gaps.plot(color='g',lw=2)
        plt.title(aln_name.replace('@',' '), fontsize=14)
        plt.axis([1,len(aln_obj[0].seq), 0, 1.1]) #0 to 1 scale for y axis
        xlab = "alignment position"
        ylab = "similarity score"
        if alg == '-sgc':
            ylab = "percent gaps"
            plt.axis([1, len(aln_obj[0].seq), 0, 110]) #0 to 100 scale for y axis, ie percent
        plt.xlabel(xlab, fontsize=10)
        plt.ylabel(ylab, fontsize=10)
        plt.grid(True)
        os.remove(name)
    figname = str(random.randint(1000,2000))
    fig.savefig(figs_folder + '/' + figname +'.png')
    plt.close('all')
    return figs_folder + '/' + figname+'.png'

def view_csv_as_table(csv_filename, delimiter, quotechar='|'):
    with open(csv_filename, 'rb') as csvfile:
        sp_vs_lc = list(csv.reader(csvfile, delimiter=delimiter, quotechar=quotechar))
        field_sizes = []
        for i in range(len(sp_vs_lc[0])):
            lengths = []
            for row in sp_vs_lc:
                lengths.append(len(row[i]))
            field_sizes.append(max(lengths))
        for row in sp_vs_lc:
            string = ''
            for i in range(len(row)):
                string += row[i].ljust(field_sizes[i]+3)
            print string

def rfmt_tree_for_and_char_matirx_bayestraits(pj, qual_list, rootmeta, rootvalue, treefile=None, 
                              treetoken=None, treeburnin=0, treestep=1, treeformat=5):
    if treefile and treetoken:
        raise IOError("Only specify treefile ot treetoken, not both")
    T = None
    if treefile:
        T = open(treefile,'r').readlines()
        T = T[:-1]
    elif treetoken:
        T = [pj.ft(treetoken).write()]
    else:
        raise IOError("Specify treefile or treetoken")
    
    leaf_names = Tree(T[0].rstrip(), format=treeformat).get_leaf_names()
    
    translate = []
    i = 1
    for n in leaf_names:
        translate.append([n,str(i)])
        i +=1 
    
    char_matrix = ""
    for n in leaf_names:
        line = "%s "%n
        quals = get_qualifiers_dictionary(pj, n)
        for qual in qual_list:
            if qual in quals.keys():
                line += "%s "%str(quals[qual])
            else:
                line += "- "
        line = line[:-1]+'\n'
        char_matrix += line
            
    reformatted_tree = """#NEXUS
    Begin Trees;
        TRANSLATE
    """
    for t in translate[:-1]:
        reformatted_tree += '\t\t'+t[1]+'\t'+t[0]+',\n'
    reformatted_tree += '\t\t'+translate[-1][1]+'\t'+translate[-1][0]+';\n'
    
    for i in range(int(len(T)*treeburnin),len(T),treestep):
        j = T[i]
        newick = Tree(j.rstrip(), format=treeformat)
        brlns = []
        for n in newick.traverse():
            brlns.append(n.dist)
        if sorted(brlns)[2] == 0.0: # three 0 length branches - too much
            pass
        else:
            R = None
            count = 0
            for l in newick:
                if get_qualifiers_dictionary(pj, l.name)[rootmeta] == rootvalue:
                    R = l.name
                    count += 1
            if not count == 1:
                raise RuntimeError("%s does not exist or not unique in qualifier %s"%(rootvalue, rootmeta))
            newick.set_outgroup(newick&R)
            newick_str = newick.write(format=5)
            for t in translate:
                newick_str = re.sub(t[0],t[1],newick_str)
            reformatted_tree += 'Tree tree'+str(i)+'= '+newick_str+'\n'
    reformatted_tree += 'End;\n'
    
    return reformatted_tree, char_matrix
    
# Exonerate

def bayestraits(pj, qual_list, rootmeta, rootvalue,
                    treefile=None, treetoken=None,
                    treeburnin=0, treestep=1,
                    treeformat=5,
                    bayestraits = 'BayesTraits',
                    commands = [4,1,'kappa','delta','lambda','run']):
    # make command file
    import random
    rand = random.randint(1000000,9999999)
    cfile = open(str(rand),'wt')
    for i in commands:
        cfile.write(str(i)+'\n')
    cfile.close()
    reformatted_tree, char_matrix = rfmt_tree_for_and_char_matirx_bayestraits(pj, qual_list, rootmeta, rootvalue, treefile=treefile, 
                                                                              treetoken=treetoken, treeburnin=treeburnin, treestep=treestep,
                                                                              treeformat=treeformat)
    
    tfile = open(str(rand)+'.nex','wt')
    tfile.write(reformatted_tree)
    tfile.close()
    
    mfile = open(str(rand)+'.txt','wt')
    mfile.write(char_matrix)
    mfile.close()
    
    cline = "%s %s %s < %s" %(bayestraits, str(rand)+'.nex', str(rand)+'.txt', str(rand))
    import os
    stdout = os.popen(cline).read()
    os.remove(str(rand))
    os.remove(str(rand)+'.nex')
    os.remove(str(rand)+'.txt')
    return stdout


# Exonerate



class ExonerateCommandLine:
    
    """cline object with execute methods"""
        
    def __init__(self,
                 q, #query filename
                 t, #target filename
                 path='exonerate',
                 Q="unknown",# query alphabet
                 T="unknown", # target alphabet
                 querychunkid=0, #query job number
                 targetchunkid=0, #target job number
                 querychunktotal=0, #Num of queries
                 targetchunktotal=0, #Num of targets
                 E="FALSE", #exhaustive search
                 B="FALSE", #rapid comparison between long seqs
                 forcescan="none", #Force FSM scan on query or target sequences q or t
                 saturatethreshold=0, #word saturation threshold
                 customserver="NULL", # Custom command to send non-standard server
                 fastasuffix=".fa", #Fasta file suffix filter (in subdirectories)
                 m="ungapped",
                 s=100, #Score threshold for gapped alignment
                 percent=0.0, #Percent self-score threshold
                 showalignment="TRUE",
                 showsugar="FALSE",
                 showcigar="FALSE",
                 showvulgar="FALSE",
                 showquerygff="FALSE", # Include GFF output on query in results
                 showtargetgff="FALSE", #Include GFF output on target in results
                 ryo="NULL", #Roll-your-own printf-esque output format
                 n=0, #Report best N results per query
                 S="TRUE", #Search for suboptimal alignments
                 g="TRUE", #Use gapped extension
                 refine="none", #none|full|region
                 refineboundary=32, #Refinement region boundary
                 D=32, #Maximum memory to use for DP tracebacks (Mb)
                 C="TRUE", #Use compiled viterbi implementations
                 terminalrangeint=12, #Internal terminal range
                 terminalrangeext=12, #External terminal range
                 joinrangeint=12, #Internal join range
                 joinrangeext=12, #External join range
                 x=50, #Gapped extension threshold
                 singlepass="TRUE", #Generate suboptimal alignment in a single pass
                 joinfilter=0, #BSDP join filter threshold
                 A="none", #Path to sequence annotation file
                 softmaskquery="FALSE", #Allow softmasking on the query sequence
                 softmasktarget="FALSE", #Allow softmasking on the target sequence
                 d="nucleic", #DNA substitution matrix
                 p="blosum62", #Protein substitution matrix
                 M=64, #Memory limit for FSM scanning <Mb>
                 forcefsm="none", #Force FSM type ( normal | compact )
                 wordjump=1, #Jump between query words
                 o=-12, #Affine gap open penalty
                 e=-4, #Affine gap extend penalty
                 codongapopen=-18, #Codon affine gap open penalty
                 codongapextend=-8, #Codon affine gap extend penalty
                 minner=10, #Minimum NER length
                 maxner=50000, #Maximum NER length
                 neropen=-20, #NER open penalty
                 minintron=30, #Minimum intron length
                 maxintron=20000, #Maximum intron length
                 i=-30, #Intron Opening penalty
                 f=-28, #Frameshift creation penalty
                 useaatla="TRUE", #useaatla
                 geneticcode=1, #Use built-in or custom genetic code
                 hspfilter=0, #Aggressive HSP filtering level
                 useworddropoff="TRUE", #Use word neighbourhood dropoff
                 seedrepeat=1, #Seeds per diagonal required for HSP seeding
                 dnawordlen=12, #Wordlength for DNA words
                 proteinwordlen=6, #Wordlength for protein words
                 codonwordlen=12, #Wordlength for codon words
                 dnahspdropoff=30, #DNA HSP dropoff score
                 proteinhspdropoff=20, #Protein HSP dropoff score
                 codonhspdropoff=40, #Codon HSP dropoff score
                 dnahspthreshold=75, #DNA HSP threshold score
                 proteinhspthreshold=30, #Protein HSP threshold score
                 codonhspthreshold=50, #Codon HSP threshold score
                 dnawordlimit=0, #Score limit for dna word neighbourhood
                 proteinwordlimit=4, #Score limit for protein word neighbourhood
                 codonwordlimit=4, #Score limit for codon word neighbourhood
                 geneseed=0, #Geneseed Threshold
                 geneseedrepeat=3, #Seeds per diagonal required for geneseed HSP seeding
                 alignmentwidth=80, #Alignment display width
                 forwardcoordinates="TRUE", #Report all coordinates on the forward strand
                 quality=0, #HSP quality threshold
                 splice3="primate", #Supply frequency matrix for 3' splice sites
                 splice5="primate", #Supply frequency matrix for 5' splice sites
                 forcegtag="FALSE"): #Force use of gt...ag splice sites
        
        import inspect

        frame = inspect.currentframe()
        args, _, _, values = inspect.getargvalues(frame)

        cline = values['path']+' '
        
        for k in ('path','frame','inspect', 'self'):
            del values[k] 
        
        for keyward in values:
            
            dash = '-'
            if len(keyward) > 1:
                 dash = '--'
            cline += dash+keyward+' '+str(values[keyward])+' '
        self.cline_string = cline
        self.stdout = None

    def __str__(self):
        return self.cline_string
    
    def execute(self):
        import subprocess as sub
        p = sub.Popen(self.cline_string, shell=True, stdout=sub.PIPE, stderr=sub.PIPE)
        self.stdout = p.communicate()[0]
        return self.stdout

# an all inclusive exonerate ryo format 

roll = ("STARTRYOqfull = %qas @!!@!!@ qcds = %qcs @!!@!!@ qid = %qi @!!@!!@ qdescription = %qd @!!@!!@ "+
       "qlen = %qal @!!@!!@ qstrand = %qS @!!@!!@ qtype = %qt @!!@!!@ qbegin = %qab @!!@!!@ qend = %qae @!!@!!@ "+
       "tfull = %tas @!!@!!@ tcds = %tcs @!!@!!@ tid = %ti @!!@!!@ tdescription = %td @!!@!!@ "+
       "tlen = %tal @!!@!!@ tstrand = %tS @!!@!!@ ttype = %tt @!!@!!@ tbegin = %tab @!!@!!@ tend = %tae @!!@!!@ "+
       "Etotal = %et @!!@!!@ Eident = %ei @!!@!!@ Esim = %es @!!@!!@ Emis = %em @!!@!!@ Pident = %pi @!!@!!@ "+
       "Psim = %ps @!!@!!@ score = %s @!!@!!@ model = %m @!!@!!@ vulgar = %VENDRYO")

def parse_ryo(exo_results):
    
    """parses exonerate results that include the ryo line above"""
    
    stats = []
    ryos = [i.split("ENDRYO")[0] for i in exo_results.split('STARTRYO')[1:]][:-1]
    if len(ryos) > 0:
        for i in ryos:
            a = {}
            for line in i.split('@!!@!!@'):
                k, v = [line.partition('=')[0], line.partition('=')[2]]
                a[k.strip().rstrip()] = v.strip().rstrip()
            a['qfull'] = a['qfull'].replace('\n','').replace("\s",'')
            a['qcds'] = a['qcds'].replace('\n','').replace("\s",'')
            a['tfull'] = a['tfull'].replace('\n','').replace("\s",'')
            a['tcds'] = a['tcds'].replace('\n','').replace("\s",'')
            stats.append(a)
            
    return stats

def exonerate(q, d, **kwargs):
    
    """Will run exonerate with the ryo format above.
    Returns a dictionary with the ryo line content and the raw output as string"""
    
    if 'ryo' in kwargs.keys():
        kwargs.pop('ryo', None)
    results = ''
    exoCline = ExonerateCommandLine(q, d, ryo="\"%s\""%roll, **kwargs)
    results = exoCline.execute()
    #print results
    stats = parse_ryo(results)
    #print stats
    return stats, results


def exonerate_ryo_to_gb(q, d, stats, results, get_query=False):
    
    """takes the parsed ryo (stats) and the raw output (results) and
    builds a gb file for reprophylo"""
    
    if get_query:
        raise RuntimeError('get_query=True: currently only parses target. maybe Bio.SearchIO can help?')
    
    gencode = int(results.split('geneticcode ')[1].split()[0])
    model = stats[0]['model']
    if not model == 'protein2genome:local':
        raise RuntimeError("only tested with the protein2genome model")
    #gencode = '1'
    tfile = d
    if '/' in d:
        tfile = d.split('/')[-1]
    matches = stats
    b = 0
    records = []
    # making a short enough, yet unique and informative seq id is a challenge
    # The approach here is to take the 4 first and last chars of the input file name
    # and to add a serial number for each seq.
    
    # We add a random three digit number at the start
    # because the eight file name chars are not always unique accorss files.
    from random import randint
        
    rnd = randint(111,999)
    
    for match in matches:
        ID = ("%i|%s|%s|%i"%(rnd, tfile[:4],tfile[-4:],b)).replace('.','')
        
        r = SeqRecord(seq=Seq(match['tfull'],
                              alphabet=IUPAC.unambiguous_dna),
                      id = ID,
                      description = "query: %s %s, target: %s %s, HSPID: %i"%(match['qid'],
                                                                              match['qdescription'],
                                                                              match['tid'],
                                                                              match['tdescription'],
                                                                              b))
        b += 1
        source = SeqFeature(FeatureLocation(0, len(r.seq)), type='source')
        source.qualifiers['file'] = [tfile]
        for key in match:
            source.qualifiers[key] = [match[key]]
        r.features.append(source)
        vulgar = match['vulgar'].split()
        features = [vulgar[i:i+3] for i in range(0, len(vulgar) - 2, 3)]
        pos = 0
        CDS = None
        for feature in features:
            ftype = feature[0]
            flength = int(feature[2])
            f = SeqFeature(FeatureLocation(pos, pos+flength), type = ftype)
            f.qualifiers['gene'] = [match['qid']]
            pos += flength
            r.features.append(f)
        coding_locations = []
        for j in r.features:
            if j.type == 'M' or j.type == 'S':
                coding_locations.append(j.location)
            elif j.type == 'G' and (int(j.location.start) < int(j.location.end)):
                coding_locations.append(j.location)
        coding_locations = sorted(coding_locations, key = lambda l: int(l.start))
        if len(coding_locations) == 1:
            CDS = SeqFeature(coding_locations[0], type='CDS')
        else:
            CDS = SeqFeature(CompoundLocation(coding_locations), type='CDS')
        CDS.qualifiers['gene'] = [match['qid']]
        get = False
        CDS.qualifiers['translation'] = [str(CDS.extract(r.seq).translate(table=gencode)).replace('*','X')]
        if CDS.qualifiers['translation'][0] == str(Seq(match['tcds'].replace(' ','').replace('\n',''),
                                                       alphabet=IUPAC.ambiguous_dna).translate(table=gencode)).replace('*','X'):
            get = True
        else:
            CDS.qualifiers['translation'] = ['something went wrong']
            print "DEBUG bad CDS"
            print match['tid']
            #print len(CDS.extract(r.seq)), len(match['tcds'].replace(' ',''))
            #print str(CDS.extract(r.seq))[:20], str(CDS.extract(r.seq))[-20:]
            #print match['tcds'].replace(' ','')[:20], match['tcds'].replace(' ','')[-20:]
            print (str(CDS.extract(r.seq)) == match['tcds'].replace(' ',''))
            if (str(CDS.extract(r.seq)) == match['tcds'].replace(' ','')):
                print 'CDS retrieved correctly buy biopython could not translate'
            else:
                print 'Error in retrieved CDS (CDS built form vulgar does not match the CDS from ryo \%tcs'
            #print str(CDS.extract(r.seq))
        if get:
            r.features.append(CDS)
            records.append(r)
    a = SeqIO.write(records,'%s.gb'%d,'genbank')
    return "%i in %s.gb"%(a, d)

def report_aln_col_stat(pj, loci_names, num_col, figs_folder, trimmed=False, alg = '-scc'):
    import pandas as pd
    import matplotlib.pyplot as plt
    import random, os
    from Bio import AlignIO
    
    # get the alignment objects
    #-------------------------#
    alignments = pj.alignments.items()
    if trimmed:
        alignments = pj.trimmed_alignments.items()
    alignments = [i for i in alignments if i[0].split('@')[0] in loci_names]
    num_alns = len(alignments)
    subplots_arrangement = []
    #-----------------------#
    num_rows = round(float(num_alns)/num_col)
    if num_rows < float(num_alns)/num_col:
        num_rows += 1
        
    fig = plt.figure(figsize=(10*num_col,2.3*num_rows), dpi=80, frameon = False)
    plt.subplots_adjust(hspace = 0.8)
    subplots_arrangement += [num_rows, num_col]


    #Calc with trimal and plot
    #------------------------#
    for i in range(1,num_alns+1):
        import subprocess as sub
        subplot_position = subplots_arrangement +[i]
        aln_name = alignments[i-1][0]
        aln_obj = alignments[i-1][1]
        name = str(random.randint(1000,2000))+'_'+aln_name+'_for_trimal_graph.fasta'
        AlignIO.write(aln_obj, name, 'fasta')
        stderr = open('stderr','wt')
        #stdout = os.popen('trimal '+alg+' -in '+name)#.read()
        stdout = sub.Popen("trimal "+alg+" -in " + name,
                       shell=True, stdout=sub.PIPE, stderr=stderr).stdout
        stderr.close()
        var = pd.read_table(stdout, sep='\t+', skiprows=3, engine='python')
        os.remove('stderr')
        if alg == '-scc':
            var.columns = ['position', 'variability']
        elif alg == '-sgc':
            var.columns = ['position', 'pct_gaps', 'gap_score']
        
        #Plot residue similarity figure, for nucleotides this is identity value
        fig.add_subplot(subplot_position[0], subplot_position[1], subplot_position[2])
        if alg == '-scc':
            var.variability.plot(color='g',lw=2)
        elif alg == '-sgc':
            var.pct_gaps.plot(color='g',lw=2)
        plt.title(aln_name.replace('@',' '), fontsize=14)
        plt.axis([1,len(aln_obj[0].seq), 0, 1.1]) #0 to 1 scale for y axis
        xlab = "alignment position"
        ylab = "similarity score"
        if alg == '-sgc':
            ylab = "percent gaps"
            plt.axis([1, len(aln_obj[0].seq), 0, 110]) #0 to 100 scale for y axis, ie percent
        plt.xlabel(xlab, fontsize=10)
        plt.ylabel(ylab, fontsize=10)
        plt.grid(True)
        os.remove(name)
    figname = str(random.randint(1000,2000))
    fig.savefig(figs_folder + '/' + figname +'.png')
    plt.close('all')
    return figs_folder + '/' + figname+'.png'

##############################################################################################################
if False:
    """"LociStat class preliminaries"""
##############################################################################################################

def entropy(s, char_type):
    """
    Return the Shannon's entropy value for a column in the alignment provided as a string (s) 
    given the character type (dna or prot).
    
    gaps are ignored, ambiguity is ignored.
    
    
    homogenous column
    
    >>> entropy('tttttttt', 'dna')
    -0.0
    
    hetrogenous column, case insensitive
    
    >>> entropy('ttgGaacC', 'dna')
    2.0
    >>> entropy('ttggaacc', 'dna')
    2.0
    
    ignore gaps
    
    >>> entropy('ttgg--Ss', 'dna')
    1.0
    >>> entropy('ttggSs', 'dna')
    1.0
    
    recognize alphabet
    
    >>> entropy('ttggSs', 'prot')
    1.584962500721156
    """
    missing = None
    if char_type == 'prot':
        missing = 'Xx'
    elif char_type == 'dna':
        missing = 'ryswkmbdhvnRYSWKMBDHVN'
    s = s.replace('-','').replace('.','').replace('?','')
    for m in missing:
        s = s.replace(m,'')
    p, lns = Counter(s.lower()), float(len(s.lower()))
    
    return -sum( count/lns * math.log(count/lns, 2) for count in p.values())

def gapscore(aln_obj):
    
    """
    Use  TrimAl to get a list of gapscores given a MultipleSeqAlignmnet
    object
    """
    
    import pandas as pd
    name = str(random.randint(1000,2000))+'_for_trimal_graph.fasta'
    AlignIO.write(aln_obj, name, 'fasta')
    stderr = open('stderr','wt')
    stdout = sub.Popen(programspath+"trimal -sgc -in " + name,
                       shell=True, stdout=sub.PIPE, stderr=stderr).stdout
    stderr.close()
    var = pd.read_table(stdout, sep='\t+', skiprows=3, engine='python')
    os.remove('stderr')
    var.columns = ['position', 'pct_gaps', 'gap_score']
    os.remove(name)  
    return [i[1] for i in var.to_dict()['gap_score'].items()]

def conservation(aln_obj):
    
    """
    Use  TrimAl to get a list of conservation values given a MultipleSeqAlignmnet
    object
    """
    import pandas as pd
    name = str(random.randint(1000,2000))+'_for_trimal_graph.fasta'
    AlignIO.write(aln_obj, name, 'fasta')
    stderr = open('stderr','wt')
    stdout = sub.Popen(programspath+"trimal -scc -in " + name,
                       shell=True, stdout=sub.PIPE, stderr=stderr).stdout
    stderr.close()
    var = pd.read_table(stdout, sep='\t+', skiprows=3, engine='python')
    os.remove('stderr')
    var.columns = ['position', 'conservation']
    os.remove(name)  
    return [i[1] for i in var.to_dict()['conservation'].items()]

def get_entropies(pj, trimmed = True, alignmnet_method=None, trimming_method=None):
    
    """
    Return a dictionary with alignment names as keys and  entropy keys as values
    given the alignments or trimmed alignmnets dictionary.
    
    
    If a locus has more than one alignment in the dictionary, go for the specified method
    or the first occurance
    """
    entropies = {}
    aln_dict = None
    if trimmed:
        aln_dict = pj.trimmed_alignments
        if aln_dict == {}:
            raise IOError("No trimed alignments in the Project")
                          
    elif not trimmed:
        aln_dict = pj.alignments
        if aln_dict == {}:
            raise IOError("No alignments in the Project")
        
    for aln_name in aln_dict.keys():
        char_type = None
        char_type_list = [l.char_type for l in pj.loci if l.name == aln_name.split('@')[0]]
        get = True
        if len(char_type_list) == 0:
            get = False
            warnings.warn('Cannot find Locus for alignment %s. Is it a supermatrix? Skipping.'%aln_name)            
        elif len(char_type_list) > 1:
            if (alignmnet_method and 
                not aln_name.split('@')[1] == alignmnet_method):
                get = False
                warnings.warn('Skipping %s, taking only %s'%aln_name,alignmnet_method)
            elif (trimmed and trimming_method and 
                not aln_name.split('@')[2] == trimming_method):
                get = False
                warnings.warn('Skipping %s, taking only %s'%aln_name,trimming_method)  
            elif aln_name.split('@')[0] in [i.split('@')[0] for i in entropies.keys()]:
                exists = [i for i in entropies.keys() if i.split('@')[0] == aln_name.split('@')[0]][0]
                get = False
                warning.warn('Skipping %s, already have %s'%aln_name, exists)
        if get:    
            char_type = char_type_list[0]
            aln_obj = aln_dict[aln_name]
            entropies[aln_name] =[]
            for i in range(aln_obj.get_alignment_length()):
                column = aln_obj[:,i]
                entropies[aln_name].append(entropy(column, char_type))

    return entropies

def get_gapscores(pj, trimmed = True, alignmnet_method=None, trimming_method=None):
    
    """
    Return a dictionary with alignment names as keys and  gap scores as values
    given the alignments or trimmed alignmnets dictionary.
    
    
    If a locus has more than one alignment in the dictionary, go for the specified method
    or the first occurance
    """
    gapscores = {}
    aln_dict = None
    if trimmed:
        aln_dict = pj.trimmed_alignments
        if aln_dict == {}:
            raise IOError("No trimed alignments in the Project")
    elif not trimmed:
        aln_dict = pj.alignments
        if aln_dict == {}:
            raise IOError("No alignments in the Project")
        
    for aln_name in aln_dict.keys():
        char_type = None
        char_type_list = [l.char_type for l in pj.loci if l.name == aln_name.split('@')[0]]
        get = True
        if len(char_type_list) == 0:
            get = False
            warnings.warn('Cannot find Locus for alignment %s. Is it a supermatrix? Skipping.'%aln_name)            
        elif len(char_type_list) > 1:
            if (alignmnet_method and 
                not aln_name.split('@')[1] == alignmnet_method):
                get = False
                warnings.warn('Skipping %s, taking only %s'%aln_name,alignmnet_method)
            elif (trimmed and trimming_method and 
                not aln_name.split('@')[2] == trimming_method):
                get = False
                warnings.warn('Skipping %s, taking only %s'%aln_name,trimming_method)  
            elif aln_name.split('@')[0] in [i.split('@')[0] for i in gapscores.keys()]:
                exists = [i for i in gapscores.keys() if i.split('@')[0] == aln_name.split('@')[0]][0]
                get = False
                warning.warn('Skipping %s, already have %s'%aln_name, exists)   
        if get:    
            char_type = char_type_list[0]
            aln_obj = aln_dict[aln_name]
            gapscores[aln_name] = gapscore(aln_obj)  
    return gapscores

def get_conservations(pj, trimmed = True, alignmnet_method=None, trimming_method=None):
    """
    Return a dictionary with alignment names as keys and  conservation scores as values
    given the alignments or trimmed alignmnets dictionary.
    
    
    If a locus has more thn one alignment in the dictionary, go for the specified method
    or the first occurance
    """
    conservations = {}
    aln_dict = None
    if trimmed:
        aln_dict = pj.trimmed_alignments
        if aln_dict == {}:
            raise IOError("No trimed alignments in the Project")
            
    elif not trimmed:
        aln_dict = pj.alignments
        if aln_dict == {}:
            raise IOError("No alignments in the Project")
        
    for aln_name in aln_dict.keys():
        char_type = None
        char_type_list = [l.char_type for l in pj.loci if l.name == aln_name.split('@')[0]]
        get = True
        if len(char_type_list) == 0:
            get = False
            warnings.warn('Cannot find Locus for alignment %s. Is it a supermatrix? Skipping.'%aln_name)            
        elif len(char_type_list) > 1:
            if (alignmnet_method and 
                not aln_name.split('@')[1] == alignmnet_method):
                get = False
                warnings.warn('Skipping %s, taking only %s'%aln_name,alignmnet_method)
            elif (trimmed and trimming_method and 
                not aln_name.split('@')[2] == trimming_method):
                get = False
                warnings.warn('Skipping %s, taking only %s'%aln_name,trimming_method) 
            elif aln_name.split('@')[0] in [i.split('@')[0] for i in conservations.keys()]:
                exists = [i for i in conservations.keys() if i.split('@')[0] == aln_name.split('@')[0]][0]
                get = False
                warning.warn('Skipping %s, already have %s'%aln_name, exists)       
        if get:    
            char_type = char_type_list[0]
            aln_obj = aln_dict[aln_name]
            conservations[aln_name] = conservation(aln_obj)  
    return conservations

def get_sequence_lengths(pj):
    
    """
    Return a dictionary with locus names as keys and lists of sequence length as values
    given a Project instance.
    
    The length are calculated from the unaligned and untrimmed sequences
    """
    
    if len(pj.records) == 0:
            raise IOError('No records in the Project')
    
    lengths = {}

    
    if len(pj.records_by_locus) == 0:
        pj.extract_by_locus()
    
    for locus in pj.records_by_locus:
        lengths[locus] = []
        for r in pj.records_by_locus[locus]:
            lengths[locus].append(len(r.seq))
    return lengths

def get_sequence_gcs(pj):
    
    """
    Return a dictionary with locus names as keys and lists of sequence %GC as values
    given a Project instance.
    
    The %GC are calculated from the unaligned and untrimmed sequences
    """
    
    if len(pj.records) == 0:
            raise IOError('No records in the Project')
    
    gcs = {}
    
    if len(pj.records_by_locus) == 0:
        pj.extract_by_locus()
    
    if any([l.char_type == 'prot' for l in pj.loci]):
        warnings.warn('Protein loci GC content will be set to 0.0')
    
    for locus in pj.records_by_locus:
        char_type = [l for l in pj.loci if l.name == locus][0].char_type
        if char_type == 'dna':
            gcs[locus] = []
    
            for r in pj.records_by_locus[locus]:
                gcs[locus].append(GC(r.seq))
        else:
            gcs[locus] = [0.0, 0.0, 0.0, 0.0]
        
    return gcs
        
def remove_border(axes=None, top=False, right=False, left=True, bottom=True):
    """
    Minimize chartjunk by stripping out unnecesasry plot borders and axis ticks
    
    The top/right/left/bottom keywords toggle whether the corresponding plot border is drawn
    """
    ax = axes or plt.gca()
    ax.spines['top'].set_visible(top)
    ax.spines['right'].set_visible(right)
    ax.spines['left'].set_visible(left)
    ax.spines['bottom'].set_visible(bottom)
    
    #turn off all ticks
    ax.yaxis.set_ticks_position('none')
    ax.xaxis.set_ticks_position('none')
    
    #now re-enable visibles
    if top:
        ax.xaxis.tick_top()
    if bottom:
        ax.xaxis.tick_bottom()
    if left:
        ax.yaxis.tick_left()
    if right:
        ax.yaxis.tick_right()

#########################################################################################
class LociStats:
#########################################################################################
    
    def __init__(self, pj, trimmed=True, alignmnet_method=None, trimming_method=None):
        
        """
        >>> from reprophylo import *
        >>> pj = unpickle_pj('test-data/test_locistats', git=False)
        >>> stats = LociStats(pj)
        >>> entropies = stats.entropies['dummy@ReadDirectly@no_trim']
        >>> stats.sort()
        >>> np.percentile(entropies,25)
        0.0
        >>> np.percentile(entropies,50)
        0.0
        >>> conservation = stats.conservations['dummy1@ReadDirectly@no_trim']
        >>> np.percentile(conservation,25)
        0.0017470471212500001
        """
        
        self.loci = pj.loci
        
        self.entropies = get_entropies(pj,
                                       trimmed = trimmed,
                                       alignmnet_method=alignmnet_method,
                                       trimming_method=trimming_method)
        
        self.gapscores = get_gapscores(pj,
                                       trimmed = trimmed,
                                       alignmnet_method=alignmnet_method,
                                       trimming_method=trimming_method)
        
        self.conservations = get_conservations(pj,
                                               trimmed = trimmed,
                                               alignmnet_method=alignmnet_method,
                                               trimming_method=trimming_method)
        
        self.sequence_lengths = get_sequence_lengths(pj)
        
        self.sequence_gcs = get_sequence_gcs(pj)
        
        combined = []
        
        for key in self.sequence_lengths:
            locus_stats = [key]
            
            entropies = [self.entropies[i] for i in self.entropies.keys() if i.split('@')[0] == key][0]
            gapscores = [self.gapscores[i] for i in self.gapscores.keys() if i.split('@')[0] == key][0]
            conservations = [self.conservations[i] for i in self.conservations.keys() if i.split('@')[0] == key][0]
            
            combined.append([key,
                             self.sequence_lengths[key],
                             self.sequence_gcs[key],
                             entropies,
                             gapscores,
                             conservations])
            
        self.loci_stats = combined
        
        self.loci_stats_sorted = None
        
    def sort(self, parameter = 'entropy', percentile=50, percentile_range=(25,75), reverse = True):
        j = None
        if parameter == 'entropy':
            j = 3
        elif parameter == 'gapscore':
            j = 4
        elif parameter == 'conservation':
            j = 5
        elif parameter == 'sequence_length':
            j = 1
        elif parameter == 'sequence_gc':
            j = 2
        
        self.loci_stats_sorted = sorted(self.loci_stats,
                                        key=lambda i: (np.percentile(i[j], percentile),
                                                       abs(np.percentile(i[j], percentile_range[1])-
                                                           np.percentile(i[j], percentile_range[0]))),
                                        reverse=reverse)
        self.loci_stats_sorted = [list(i) for i in self.loci_stats_sorted]
        
    def plot(self, filename, figsize=(30,10), params='all', lable_fsize=40, xtick_fsize=4, ytick_fsize=4,
            boxcolor='salmon', whiskercolor='gray', capcolor='black', mediancolor='white', medianline_w=3):
        
        parameter_indices = [3,4,5,1,2]
        
        ytitles=['',
                 'Sequence Lengths',
                 'Sequence %GC',
                 'Entropy',
                 'Gap Score', 
                 'Conservation Scores']
        
        if not  params=='all':
            parameter_indices = []
            for param in params:
                if param == 'entropy':
                    parameter_indices.append(3)
                elif param == 'gapscore':
                    parameter_indices.append(4)
                elif param == 'conservation':
                    parameter_indices.append(5)
                elif param == 'sequence_length':
                    parameter_indices.append(1)
                elif param == 'sequence_gc':
                    parameter_indices.append(2)
        if len(parameter_indices) == 0:
            raise IOError('Must specify at least one parameter to plot')
        
        fig, axes = plt.subplots(len(parameter_indices), sharex=True, figsize=figsize, dpi=80, frameon = False)
        
        if len(parameter_indices) == 1:
               axes = [axes]
        
        labels = [k[0] for k in self.loci_stats_sorted]
        
        j = 0
        
        for ax in axes:
            
            values = [k[parameter_indices[j]] for k in self.loci_stats_sorted]
            
            bp = ax.boxplot(values,0,'', positions = range(4,(len(values)*4)+1, 4), patch_artist=True)
            ax.set_ylabel(ytitles[parameter_indices[j]], fontsize=lable_fsize)
            #plt.xlabel("Locus", fontsize=lable_fsize)
            plt.xticks(range(4,((len(values)+1)*4),4), labels, rotation=90, fontsize = xtick_fsize)
            plt.tight_layout()
            
            j += 1
            
            remove_border()
            for box in bp['boxes']:
            # change outline color
                box.set( color=boxcolor, linewidth=1)
                box.set_facecolor(boxcolor)
                    
            # change color, linestyle and linewidth of the whiskers
            for whisker in bp['whiskers']:
                whisker.set(color=whiskercolor, linestyle='solid', linewidth=2.0)
            
            # change color and linewidth of the caps
            for cap in bp['caps']:
                cap.set(color=capcolor, linewidth=2.0)
        
            # change color and linewidth of the medians
            for median in bp['medians']:
                median.set(color=mediancolor, linewidth=medianline_w)
        
        fig.gca().set_ylim(bottom=-0.05)
        fig.savefig(filename)
        
    def slice_loci(self, median_range, otu_meta, parameter='entropy',
                   otu_must_have_all_of=[], otu_must_have_one_of='any'):
        
        j = None
        if parameter == 'entropy':
            j = 3
        elif parameter == 'gapscore':
            j = 4
        elif parameter == 'conservation':
            j = 5
        elif parameter == 'sequence_length':
            j = 1
        elif parameter == 'sequence_gc':
            j = 2
        
        loci_names = [i[0] for i in self.loci_stats_sorted if median_range[0] < np.median(i[j]) < median_range[1]]
        loci = [l for l in self.loci if l.name in loci_names]
        concat_name = "%s_%.2f_%.2f_loci_%s_to_%s"%(parameter, float(median_range[0]), float(median_range[1]),
                                                    loci_names[0],loci_names[-1])
        
        return Concatenation(concat_name, loci,
                             otu_meta,
                             otu_must_have_all_of=otu_must_have_all_of,
                             otu_must_have_one_of=otu_must_have_one_of)
    

    
    def slide_loci(self, otu_meta, median_range='all', parameter='entropy', start=0, length=2, step=1,
                   otu_must_have_all_of=[],
                   otu_must_have_one_of='any'):
        
        j = None
        if parameter == 'entropy':
            j = 3
        elif parameter == 'gapscore':
            j = 4
        elif parameter == 'conservation':
            j = 5
        elif parameter == 'sequence_length':
            j = 1
        elif parameter == 'sequence_gc':
            j = 2
        
        if median_range == 'all':
            medians = [np.median(i) for i in [k[j] for k in self.loci_stats_sorted]]
            median_range = [min(medians), max(medians)]
            
        loci_in_range = [[i[0], np.median(i[j])] for i in self.loci_stats_sorted 
                          if median_range[0] <= np.median(i[j]) <= median_range[1]]
        
        concatenations = []
        
        stop = False
        
        while not stop:
            
            window_loci = loci_in_range[start: start+length]
                
            window_loci_names = [n[0] for n in window_loci]
            loci = [l for l in self.loci if l.name in window_loci_names]
            window_start_median = window_loci[0][1]
            window_end_median = window_loci[-1][1]
            concat_name = "%s_%.2f_%.2f_loci_%i_to_%i"%(parameter, float(window_start_median), float(window_end_median),
                                                        start, start+length-1)
            print concat_name
            concatenations.append(Concatenation(concat_name, loci,
                                                otu_meta,
                                                otu_must_have_all_of=otu_must_have_all_of,
                                                otu_must_have_one_of=otu_must_have_one_of))
            start = start+step
            
            if len(loci_in_range[start:]) < length:
                stop = True
        
        return concatenations
    
    

if __name__ == "__main__":
    import doctest
    doctest.testmod()
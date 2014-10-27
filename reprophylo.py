
##############################################################################################
if False:
    """
    ReproPhylo version 0.1 
    
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
    
    RAxML 8
    Phylobayes
    Trimal
    Muscle
    Mafft
    Pal2nal
    """
##############################################################################################



from Bio import SeqIO
import os, csv, sys, dendropy, re, time, random, glob, platform, warnings, rpgit, ast, gb_syn
import HTML, inspect, shutil
import subprocess as sub
#import cloud.serialization.cloudpickle as pickle
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
import __builtin__



##############################################################################################
if False:
    """Tools for loci explorations in a GenBank File"""
##############################################################################################



def list_loci_in_genbank(genbank_filename, control_filename, loci_report = None):
    
   stdout = sys.stdout
   if  loci_report: 
        sys.stdout = open(loci_report, 'w')
    
   genbank_synonyms = gb_syn.gb_syn()
    
   # Open GenBank file
   MelPCgenes = open(genbank_filename, 'rU')
   
   gene_dict = {} #set up a gene_dict dictionary
   
   # For each record
   for record in SeqIO.parse(MelPCgenes, 'genbank') :
   
      # Grab the entire sequence
      #seq = str(record.seq)  ## what is this actually used for? Nothing seems to happen on disabling it
   
      # Look at all features for this record
      for feature in record.features:
         
         # If it's a CDS or rRNA...
         if feature.type == 'CDS' or feature.type == 'rRNA':
   
            # If it contains some attribute called 'gene' save that
            if 'gene' in feature.qualifiers:
               geneName = feature.qualifiers['gene'][0]
               geneName.replace(',',';')
               if feature.type+','+geneName in gene_dict:
                   gene_dict[feature.type+','+geneName]+=1
               else:    
                   gene_dict[feature.type+','+geneName]=1
               #print(geneName)
               
            # Else if it contains some attribute called 'product' save that instead
            elif 'product' in feature.qualifiers:
               geneName = feature.qualifiers['product'][0]
               geneName.replace(',',';')
               if feature.type+','+geneName in gene_dict:
                   gene_dict[feature.type+','+geneName]+=1
               else:    
                   gene_dict[feature.type+','+geneName]=1
               #print(geneName)
               
            # Otherwise, quit.
            else:
               print 'ERROR when parsing feature: could not find either gene or product'
               print feature.qualifiers
               quit()
   #print(gene_dict)
       
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
               name = alias.replace(' ','_')
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
    ...                                concat_must_have_all_of=['coi'],
    ...                                concat_must_have_one_of =[['16S','28S'],['ALG11','18S']],
    ...                                define_trimmed_alns=["MuscleDefaults@dummyTrimMethod"])
    >>> print(str(concatenation))
    Concatenation named combined, with loci coi,18S,16S,28S,ALG11,
    of which coi must exist for all species
    and at least one of each group of [ 16S 28S ][ ALG11 18S ] is represented.
    Alignments with the following names: MuscleDefaults@dummyTrimMethod are prefered
    """
    
    concat_must_have_all_of = []
    concat_must_have_one_of = []
    define_trimmed_alns = [] #should be Locus_name@Alignment_method_name@Trimming_mathod_name
    
    feature_id_dict = {}
    
    def __init__(self,
                 name,
                 loci,
                 otu_meta,
                 concat_must_have_all_of = concat_must_have_all_of,
                 concat_must_have_one_of = concat_must_have_one_of,
                 define_trimmed_alns = define_trimmed_alns):
        self.name = name
        self.loci = loci
        self.otu_meta = otu_meta
        self.concat_must_have_all_of = concat_must_have_all_of
        self.concat_must_have_one_of = concat_must_have_one_of
        self.feature_id_dict = {}
        self.define_trimmed_alns = define_trimmed_alns
        self.used_trimmed_alns = {}
        seen = []
        for locus in loci:
            if not isinstance(locus, Locus):
                raise TypeError("Expecting Locus object in loci list")
            if locus.name in seen:
                raise NameError('Locus ' + locus.name + ' apears more than once in self.loci')
            else:
                seen.append(locus.name)
      
                
                
    def __str__(self):
        loci_names = [i.name for i in self.loci]
        loci_string = ''
        for l in loci_names:
            loci_string += l+','
        loci_string = loci_string[:-1]
        must_have = ''
        for i in self.concat_must_have_all_of:
            must_have += i+','
        must_have = must_have[:-1]
        trimmed_alignmnets_spec = ''
        one_of = ''
        for i in self.concat_must_have_one_of:
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
    
    Used in the Project class but are not in the classes methods
    """
##############################################################################################



__builtin__.git = False



def start_git():
    __builtin__.git = True
    rpgit.gitInit()
    cwd = os.getcwd()
    import fnmatch
    matches = []
    for root, dirnames, filenames in os.walk(cwd):
        for filename in fnmatch.filter(filenames, '*.py'):
            matches.append(os.path.join(root, filename))
        for filename in fnmatch.filter(filenames, '*.ipynb'):
            matches.append(os.path.join(root, filename))
    for match in matches:
        rpgit.gitAdd(match)
    comment = "%i script file(s) from %s" % (len(matches), time.asctime())
    rpgit.gitCommit(comment)
    
    
    
def stop_git():
    __builtin__.git = False
    cwd = os.getcwd()
    import fnmatch
    matches = []
    for root, dirnames, filenames in os.walk(cwd):
        for filename in fnmatch.filter(filenames, '*.py'):
            matches.append(os.path.join(root, filename))
        for filename in fnmatch.filter(filenames, '*.ipynb'):
            matches.append(os.path.join(root, filename))
    for match in matches:
        rpgit.gitAdd(match)
    comment = "%i script file(s) from %s" % (len(matches), time.asctime())
    rpgit.gitCommit(comment)
    


def platform_report():
    
    """ 
    Prints machine specs, os specs and dependencies at time of execution
    
    >>> isinstance(platform_report(), list)
    True
    """
    import pkg_resources
    modules = []
    for i in ('ete2','biopython','dendropy','cloud'):
        try:
            modules.append(i+' version: '+
                                pkg_resources.get_distribution(i).version)
        except:
            pass
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
    or product qualifiers is in the aliases of one of the loci
    
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
            qaul = False
            if 'gene' in feature.qualifiers.keys():
                qual = 'gene'
            elif 'product' in feature.qualifiers.keys():
                qual = 'product'
            if not qual == False and feature.qualifiers[qual][0] in g.aliases:
                keep = 1
    if keep == 1:
        return True
    else:
        return False
  


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
            if feature.type == 'source' and not 'feature_id' in feature.qualifiers.keys():
                feature.qualifiers['feature_id'] = [record.id + '_source']
            elif not 'feature_id' in feature.qualifiers.keys():
                feature.qualifiers['feature_id'] = [record.id + '_f' + str(feature_count)]
                feature_count += 1
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
    
        
    # Making a dummy locus    
    >>> coi = Locus('dna','CDS','coi', ['cox1','COX1','coi','COI','CoI'])
    
    # Making a dummy Project
    >>> pj = Project([coi])
    
    # making a dummy record
    >>> s = 'atgc'*1000
    >>> location = FeatureLocation(1,100)
    >>> feature = SeqFeature()
    >>> feature.location = location
    >>> feature.type = 'CDS'
    >>> feature.qualifiers['gene'] = ['CoI']
    >>> feature.qualifiers['feature_id'] = ['12345']
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
    >>> qual_dict = get_qualifiers_dictionary(pj, '12345')
    >>> qual_items = qual_dict.items()
    >>> qual_items.sort(key = lambda i: i[0])
    >>> for key, val in qual_items: print(key.ljust(20,' ') + val.ljust(20,' '))
    annotation_evidence made up             
    feature_id          12345               
    gene                CoI                 
    source_organism     Tetillda radiata    
    """
    if type(feature_id) is list:
        feature_id = feature_id[0]
    record_id = feature_id.split('_')[0]
    qualifiers_dictionary={}
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
                if not line[i] == 'null' and not line[i] == '':
                    csv_info[line[0]]['source'][qual_name] = line[i].split(';')
            elif (not 'source:_' in header[i] and not line[i] == 'null' and not line[i] == '' and
                  not i in [seq_col, translation_col, taxonomy_col, feature_id_col]):
                csv_info[line[0]]['features'][line[feature_id_col]][header[i]] = line[i].split(';')
    return csv_info

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
    non_uniform_count = 0
    parsimony_informative = 0
    for i in range(aln_obj.get_alignment_length()):
        total_gaps += aln_obj[:, i].count('-')
        prop_list.append(aln_obj[:, i].count('-')/float(len(aln_obj)))
        if len(count_positions(aln_obj[:, i]).keys()) == 1:
            non_uniform_count += 1
        elif (len(count_positions(aln_obj[:, i]).keys()) == 2 and
              '-' in  count_positions(aln_obj[:, i]).keys()):
            non_uniform_count += 1
        if len([p for p in count_positions(aln_obj[:, i]).keys() if (p != '-' and  count_positions(aln_obj[:, i])[p] > 1)]) > 1:
            parsimony_informative += 1
    mean_gap_prop = sum(prop_list)/aln_obj.get_alignment_length()
    return (mean_gap_prop, non_uniform_count, parsimony_informative)

def count_undetermined_lines(aln_obj):
    count = 0
    for seq in aln_obj:
        if str(seq.seq).count('-') == aln_obj.get_alignment_length():
            count += 1
    return count

def count_collapsed_aln_seqs(aln_obj):
    count = 1
    seen_seqs = [str(aln_obj[0].seq)]
    for seq in aln_obj[1:]:
        str_seq = str(seq.seq)
        if len([s for s in seen_seqs if (str_seq in s or s in str_seq or s == str_seq)]) == 0:
            count += 1
        seen_seqs.append(str_seq)
    return count

def aln_summary(aln_obj):
    lines = ["Alignment length: %i" % aln_obj.get_alignment_length(),
             "Number of rows: %i" % len(aln_obj),
             "Unique sequences: %i"%count_collapsed_aln_seqs(aln_obj),
             "Average gap prop.: %f\nVariable columns: %i\nParsimony informative: %i"
             %global_aln_stats(aln_obj),
             "Undetermined sequences: %i"%count_undetermined_lines(aln_obj)
             ]
    return [lines, len(aln_obj), count_undetermined_lines(aln_obj), count_collapsed_aln_seqs(aln_obj)]        
            

def loci_list_from_csv(loci):
    if any(len(line.split(',')) >= 4 for line in open(loci, 'r').readlines()):
        pass
    else:
        raise IOError("File %s has no valid loci of format char_type,feature_type,name,aliases"%loci)
        
        
    loci_dict = {}
    loci_list = []
    for line in [line.rstrip() for line in open(loci, 'r').readlines() if len(line.rstrip()) > 0]:
        if len(line.split(',')) < 4:
            raise IOError("The line %s in file %s is missing arguments. Needs at least char_type,feature_type,name,aliases"%
                          (line.rstrip(), loci))
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

    def __init__(self, loci):
        
        """
        # making dummy loci
        >>> coi = Locus('dna','CDS','coi',['COX1','cox1'])
        >>> ssu = Locus('dna','rRNA','18S',['18S','SSU'])
        
        # Making a Project object
        >>> pj = Project([coi,ssu])
        >>> print(str(pj))
        Project object with the loci coi,18S,
        """
        self.records = []
        self.starttime = str(time.asctime())
        self.loci = loci
        self.records_by_locus = {}
        self.concatenations = []
        self.alignments = {}
        self.trimmed_alignments = {}
        self.trees = {}
        self.used_methods = []
        self.sets = {}
        self.defaults = {'raxmlHPC': 'raxmlHPC-PTHREADS-SSE3',
                         'mafft': 'mafft',
                         'muscle': 'muscle',
                         'trimal': 'trimal',
                         'pb': 'pb',
                         'bpcomp': 'bpcomp',
                         'tracecomp': 'tracecomp',
                         'pal2nal': 'pal2nal.pl'}
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
            print 'Read the following loci from file %s:'%loci
            for l in self.loci:
                print str(l)
        self.aln_summaries = []
                
                
    def __str__(self):
        loci_string = ''
        for i in self.loci:
            loci_string += i.name+','
        return 'Project object with the loci '+loci_string



    ###################################
    # Project methods for reading data
    ###################################  



    def read_embl_genbank(self, input_filenames_list):
        
        """
        Read a file from Genbank of EMBL
        
        >>> input_filenames = ['data/example.gb']
        >>> locus = Locus('dna', 'CDS', 'coi', ['cox1','COX1','coi','COI','CoI'])
        >>> pj = Project([locus])
        >>> print(len(pj.records))
        0
        >>> pj.read_embl_genbank(input_filenames)
        >>> print(len(pj.records))
        90
        """
        
        if __builtin__.git:
            import rpgit
        else:
            warnings.warn('Version control off')    
        generators = []
        for input_filename in input_filenames_list:
            if __builtin__.git:
                import rpgit
                rpgit.gitAdd(input_filename)
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
            rpgit.gitCommit(comment)          
            
        
        
    def read_denovo(self, input_filenames, char_type, format = 'fasta'):
        
        """
        Include records from a fasta file. Fasta file records will be given record ids 
        of the form 'denovo1'. The record.id and record.description will be placed in a
        source feature under the 'original_id' and 'original_desc' qualifiers. Denovo sequences
        require the use of the add_feature_to_record() method in order to be included in the
        anaysis.
        

        >>> input_filenames = ['data/example.gb']
        >>> locus = Locus('dna', 'CDS', 'coi', ['cox1','COX1','coi','COI','CoI'])
        >>> pj = project([locus])
        >>> print(len(pj.records))
        0
        >>> pj.read_embl_genbank(input_filenames)
        >>> print(len(pj.records))
        90
        >>> input_filenames = ['data/example_denovo.fasta']
        >>> pj.read_denovo(input_filenames, 'dna')
        >>> print(len(pj.records))
        91
        
        # Since the denovo sequence has no feature it is not included
        >>> pj.extract_by_locus()
        >>> print(len(pj.records_by_locus['coi']))
        90
        
        # Making a feature for the denovo record.
        >>> pj.add_feature_to_record('denovo0', 'CDS',  qualifiers={'gene': 'coi'})
        >>> pj.extract_by_locus()
        >>> print(len(pj.records_by_locus['coi']))
        91
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
                if serial > count:
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
            rpgit.gitCommit(comment)
        return count       
               
               
    def add_feature_to_record(self, record_id, feature_type, location='full', qualifiers={}):
    
        """
        # Making a dummy locus    
        >>> coi = Locus('dna','CDS','coi', ['cox1','COX1','coi','COI','CoI'])
    
        # Making a dummy Project
        >>> pj = Project([coi])
    
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
        >>> pj.add_feature_to_record('1', 'CDS', qualifiers={'gene': 'madeuplocus'})
        >>> print(len(pj.records[0].features))
        2
        """
    
        for record in self.records:
            if record.id == record_id:
                #determine new feature id
                feature_id = None
                serials = []
                for feature in record.features:
                    if 'feature_id' in feature.qualifiers.keys():
                        if '_f' in feature.qualifiers['feature_id']:
                            f = feature.qualifiers['feature_id']
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
        ...                          concat_must_have_all_of=['coi'],
        ...                          concat_must_have_one_of =[['18S','28S']],
        ...                          define_trimmed_alns=["MafftLinsi@Gappyout"])
        >>> print(str(combined))
        Concatenation named combined, with loci coi,18S,28S,
        of which coi must exist for all species
        and at least one of each group of [ 18S 28S ] is represented.
        Alignments with the following names: MafftLinsi@Gappyout are prefered
        
        # making a dummy Project
        >>> pj = Project(loci)
        
        # Including the Concatenation in the Project
        >>> pj.add_concatenation(combined)
        >>> print(len(pj.concatenations))
        1
        """
        
        if isinstance(concatenation_object, Concatenation):
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

    def make_concatenation_alignments(self):
        
        """
        Concatenates a trimmed alignment based on each of the Concatenation objects and adds them
        to the pj.trimmed_alignments dictionary. While a trimmed alignment of an individual locus will have a key
        following the patten "locus_name@alignment_method_name@trimming_method_name, the key for a concatenated
        trimmed alignment will be the Concatenation object name attribute.
        """
        for s in self.concatenations:
            
            # get a non-redundant list of OTUs stored in 'meta', such as voucher specimen
            meta = s.otu_meta
            OTU_list = []
            for record in self.records:
                for feature in record.features:
                    if not feature.type == 'source':
                        qualifiers_dictionary = get_qualifiers_dictionary(self,
                                                                          feature.qualifiers['feature_id'])
                        if (meta in qualifiers_dictionary.keys() and
                            not qualifiers_dictionary[meta] in OTU_list):
                            OTU_list.append(qualifiers_dictionary[meta])
                            
            # make lists of available feature_ids in each locus
            available_features = {}
            for locus in s.loci:
                available_features[locus.name] = []
                for record in self.records_by_locus[locus.name]:
                    available_features[locus.name].append(record.id) # record ids are feature ids because taken from pj.records_by_locus

            # make a dict of individuals that fulfil the concat's first rule
            seen_locus_names = []
            included_individuals = {}
            
            for individual in OTU_list:
                include = True
                for must_have_locus_name in s.concat_must_have_all_of:
                    if not must_have_locus_name in seen_locus_names:
                        seen_locus_names.append(must_have_locus_name)
                    locus_specific_features = []
                    for feature_id in available_features[must_have_locus_name]:
                        qualifiers_dictionary = get_qualifiers_dictionary(self,feature_id)
                        if meta in qualifiers_dictionary.keys() and qualifiers_dictionary[meta] == individual:
                            locus_specific_features.append(feature_id)
                    if len(locus_specific_features) == 1:
                        if not individual in included_individuals.keys():
                            included_individuals[individual] = {}
                        included_individuals[individual][must_have_locus_name] = locus_specific_features[0]
                    elif len(locus_specific_features) > 1:
                        raise RuntimeError(individual + ' is not unique for ' + must_have_locus_name)
                    else:
                        include = False
                if individual in included_individuals.keys() and not include:
                    included_individuals.pop(individual, None)
                            
            # check if the individual fullfil the second set rule
            for individual in included_individuals.keys():
                include = True
                for loci_group in s.concat_must_have_one_of:
                    count = 0
                    for locus_name in loci_group:
                        if not locus_name in seen_locus_names:
                            seen_locus_names.append(locus_name)
                        locus_specific_features = []
                        for feature_id in available_features[locus_name]:
                            qualifiers_dictionary = get_qualifiers_dictionary(self,feature_id)
                            if meta in qualifiers_dictionary.keys() and qualifiers_dictionary[meta] == individual:
                                locus_specific_features.append(feature_id)
                        if len(locus_specific_features) == 1:
                            count += 1
                            included_individuals[individual][locus_name] = locus_specific_features[0]
                        elif len(locus_specific_features) > 1:
                            raise RuntimeError(individual + ' is not unique for ' + locus_name)
                    if count == 0:
                        include = False
                if not include:
                    included_individuals.pop(individual, None)

            # add loci that are in the set but not addressed in rules
            for individual in included_individuals.keys():
                for locus in s.loci:
                    locus_specific_features = []
                    if not locus.name in seen_locus_names:
                        for feature_id in available_features[locus.name]:
                            qualifiers_dictionary = get_qualifiers_dictionary(self,feature_id)
                            if meta in qualifiers_dictionary.keys() and qualifiers_dictionary[meta] == individual:
                                locus_specific_features.append(feature_id)
                        if len(locus_specific_features) == 1:
                            included_individuals[individual][locus.name] = locus_specific_features[0]
                        elif len(locus_specific_features) > 1:
                            raise RuntimeError(individual + ' is not unique for ' + locus.name)

            # build alignment
            concat_records = []
            alignment = []
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
                        elif definition.count('@') == 1 and bool(any(definition in i for i in all_locus_trimmed_alns_in_pj)):
                            trimmed_aln = locus.name+'@'+definition
                        else:
                            raise RuntimeError("Could not determine which alignment/trimming alternative to use for locus '"+
                                                locus.name+"' out of "+str(locus_trimmed_alns))
                if trimmed_aln:
                    keys_of_trimmed_alignments_to_use_in_concat.append(trimmed_aln)
                else:
                    raise RuntimeError('Could not find trimmed aln for locus '+locus.name+' given the rulls '+str(s.define_trimmed_alns))
            
            print "%i individuals will be included in the concatenations %s"%(len(included_individuals.keys()), s.name)
            
            if len(included_individuals.keys()) < 4:
                raise RuntimeError("Concatenation %s has less than 4 OTUs and cannot be analyzed"%s.name)
            
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
                        sequence += 'N'*length
                concat_sequence = SeqRecord(seq = Seq(sequence), id = individual, description = '')
                alignment.append(concat_sequence)
            self.trimmed_alignments[s.name] = MultipleSeqAlignment(alignment)                
            s.feature_id_dict = included_individuals    



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
                                        quotechar='|',
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

    def correct_metadata_from_file(self,csv_file):
        metadata = read_feature_quals_from_tab_csv(csv_file)
        for record in self.records:
            record_corrected_metadata = metadata[record.id]
            record.annotations['taxonomy'] = metadata[record.id]['taxonomy']
            for feature in record.features:
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
                
                
    def if_this_then_that(self, IF_THIS, IN_THIS, THEN_THAT, IN_THAT, mode = 'whole'):
        
        
        """
        Searches pj.records for features that have the value IF_THIS in the qualifier IN_THIS
        and places the value THEN_THAT in the qualifier IN_THAT, which either exists or is new.
        
        The IF_THIS value can either match completely (mode = 'whole') or just to a part (mode = 'part')
        of the target qualifier value
        
        The following demonstartes all the feature qualifier editing methods
        
        # Make a dummy pj with a locus and with records
        >>> input_filenames = ['data/example.gb']
        >>> locus = Locus('dna', 'CDS', 'coi', ['cox1','COX1','coi','COI','CoI'])
        >>> pj = Project([locus])
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
        
        >>> input_filenames = ['data/example.gb']
        >>> coi = Locus('dna', 'CDS', 'coi', ['cox1','COX1','coi','COI','CoI'])
        >>> lsu = Locus('dna', 'rRNA', '28S', ['28s','28S','LSU rRNA','28S ribosomal RNA','28S large subunit ribosomal RNA'])
        >>> pj = Project([coi, lsu])
        >>> pj.read_embl_genbank(input_filenames)
        >>> pj.extract_by_locus()
        >>> print(len(pj.records_by_locus['coi']))
        90
        >>> print(len(pj.records_by_locus['28S']))
        48
        """
        
        data_by_locus = {}
        for locus in self.loci:
            if not locus.name in locus.aliases:
                locus.aliases.append(locus.name)
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
                    if start_from_max:
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
                    if not start_from_null:
                        self.records_by_locus[key] = subset+keep_safe[key]
            else:
                warnings.warn('Locus name %s not recognised'%key)

    def filter_by_seq_length(self, min_length=0, max_length=None):
        if self.records_by_locus == {}:
            self.extract_by_locus()
        for key in self.records_by_locus.keys():
            subset = [r for r in self.records_by_locus[key] if len(r) >= min_length]
            if max_length:
                subset = [r for r in subset if len(r) <= max_length]
            self.records_by_locus[key] = subset
                


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
                    method.platform.append('Program and version: '+os.popen(method.cmd + ' -version').read())
                elif method.program_name == 'mafft':
                    method.platform.append('Program and version: Get mafft to spit version to stdout')
                for locus in method.loci:
                    if locus.name in seen_loci:
                        #raise RuntimeError('locus '+locus.name+' is in more than one AlnConf objects')
                        pass
                    else:
                        seen_loci.append(locus.name)
                    stdout, stderr = method.command_lines[locus.name]()
                    align = AlignIO.read(StringIO(stdout), "fasta",  alphabet=IUPAC.protein)
                    if method.CDSAlign and locus.feature_type == 'CDS' and locus.char_type == 'dna':
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
                        aln_filename = method.id+'_'+locus.name+'.aln'
                        AlignIO.write(align, aln_filename, 'fasta')
                        cds_filename = method.id+'_CDS_in_frame_'+locus.name+'.fasta'
                        stdout = os.popen('perl '+pal2nal+' '+aln_filename+' '+cds_filename + ' -nostderr').read()
                        align = AlignIO.read(StringIO(stdout), "clustal",  alphabet=IUPAC.ambiguous_dna)
                        #from Bio import CodonAlign
                        #codon_aln = CodonAlign.build(align, method.CDS_in_frame[locus.name])
                        #align = codon_aln
                    method_files = glob.glob(method.id+'_*')
                    [summary_lines, num_lines, num_undeter, num_collapsed_aln_seqs] = aln_summary(align)
                    summary = 'Alignment name: '+locus.name+'@'+method.method_name+'\n'
                    for line in summary_lines:
                        summary += line+'\n'
                    if num_lines < 4:
                        line = 'Alignment %s has less than 4 sequences and will be dropped'%(locus.name+'@'+method.method_name)
                        print line
                        summary += line+'\n'
                    elif num_undeter > 0:
                        line = 'Alignment %s has undetermined sequences and will be dropped'%(locus.name+'@'+method.method_name)
                        print line
                        summary += line+'\n'
                    elif num_collapsed_aln_seqs < 4:
                        line = 'Alignment %s has less than 4 unique sequences and will be dropped'%(locus.name+'@'+method.method_name)
                        print line
                        summary += line+'\n'
                        
                    else:
                        self.alignments[locus.name+'@'+method.method_name] = align
                    self.aln_summaries.append(summary)
                method.timeit.append(time.time())
                method.timeit.append(method.timeit[2]-method.timeit[1])
                for f in method_files:
                    os.remove(f)
            self.used_methods += alignment_methods



    def write_alns(self, format = 'fasta'):
        if len(self.alignments.keys()) == 0:
            raise IOError('Align the records first')
        else:
            for key in self.alignments:
                AlignIO.write(self.alignments[key], key+'_aln.'+format, format)



    def write_trimmed_alns(self, format = 'fasta'):
        if len(self.trimmed_alignments.keys()) == 0:
            raise IOError('Align and trimmed the records first')
        else:
            for key in self.trimmed_alignments.keys():
                AlignIO.write(self.trimmed_alignments[key], key+'_trimmed_aln.'+format, format)



    def tree(self, raxml_methods, bpcomp='default'):
        # to do: determine the program used and the resulting expected tree file name
        
        if bpcomp == 'default':
            bpcomp = self.defaults['bpcomp']
        
        for raxml_method in raxml_methods:
            raxml_method.timeit.append(time.time())
            raxml_method.platform = platform_report() 
            raxml_method.platform.append('Program and version: '+ raxml_method.cmd + ': ' +
                                         os.popen(raxml_method.cmd + ' -version').readlines()[2])
            for trimmed_alignment in raxml_method.command_lines.keys():
                for cline in raxml_method.command_lines[trimmed_alignment]:
                    stdout, stderr = cline()
                t = None
                if raxml_method.preset == 'fa':
                    t = Tree('RAxML_bipartitions.'+raxml_method.id+'_'+trimmed_alignment+'0')
                elif raxml_method.preset == 'fD_fb':
                    t = Tree('RAxML_bipartitions.'+raxml_method.id+'_'+trimmed_alignment+'1')
                elif raxml_method.preset == 'fd_b_fb':
                    t = Tree('RAxML_bipartitions.'+raxml_method.id+'_'+trimmed_alignment+'2')
                    
            
                for n in t.traverse():
                    n.add_feature('tree_method_id', str(raxml_method.id)+'_'+trimmed_alignment)
                t.dist = 0
                t.add_feature('tree_method_id', str(raxml_method.id)+'_'+trimmed_alignment)
                
           
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
                        self.trees[trimmed_alignment+'@'+raxml_method.method_name] = [t,t.write(features=[])]
                        
                elif trimmed_alignment in concat_names:
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
            raxml_method.timeit.append(time.time())
            raxml_method.timeit.append(raxml_method.timeit[2]-raxml_method.timeit[1])
            for file_name in os.listdir(os.curdir):
                        if raxml_method.id.partition('_')[0] in file_name:
                            os.remove(file_name)
        self.used_methods += raxml_methods



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



    def annotate(self, fig_folder,
    
                 root_meta,
                 root_value,
    
                 leaf_labels_txt_meta,
                 leaf_node_color_meta=None,
                 leaf_label_colors=None,
    
                 node_bg_meta=None,
                 node_bg_color=None,
                 
                 node_support_dict=None,
                 
                 heat_map_meta = None, #list
                 heat_map_colour_scheme=2,
                 
                 multifurc=None,
                 
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
                for color in node_support_dict.keys():
                    ts.legend.add_face(CircleFace(radius = 4, color = color), column=i)
                    i +=1 
                    ts.legend.add_face(TextFace(' '+str(node_support_dict[color][0])+'-'+str(node_support_dict[color][1]),
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
                    leaf_face = TextFace(leaf_label, fgcolor=fgcolor)
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
            
                if multifurc:
                    for n in self.trees[tree][0].traverse():
                        if n.support < multifurc and not n.is_leaf():
                            n.delete()
    
                # node bg colors
                if node_bg_color:
                    for key in node_bg_color.keys():
                        for node in self.trees[tree][0].get_monophyletic(values=[key], target_attr=node_bg_meta):
                            ns = NodeStyle(bgcolor=node_bg_color[key])
                            node.set_style(ns)
    
                # node support
                if node_support_dict:
                    for node in self.trees[tree][0].traverse():
                        for key in node_support_dict.keys():
                            if (node.support <= node_support_dict[key][0] and
                                node.support > node_support_dict[key][1]):
                                node.add_face(CircleFace(radius = 5, color = key),column=0, position = "float")             
                    
                self.trees[tree][0].render(fig_folder + "/"+self.trees[tree][0].get_leaves()[0].tree_method_id+'.png',w=1000, tree_style=ts)
                print('<A href='+
                       fig_folder + "/" + self.trees[tree][0].get_leaves()[0].tree_method_id+'.png'+
                       '>'+self.trees[tree][0].get_leaves()[0].tree_method_id+
                       '</A><BR>')
            print '</html>'
            print fig_folder
            sys.stdout = stdout
     


#    def trim(self):
#        for aln in self.alignments.keys():
#            AlignIO.write(self.alignments[aln],aln+'_aln.fasta','fasta')
#            stdout = os.popen('trimal -in '+ aln +'_aln.fasta -gappyout').read()
#            align = AlignIO.read(StringIO(stdout), "fasta",  alphabet=IUPAC.ambiguous_dna)
#            for record in align:
#                record.description = ''
#            self.trimmed_alignments[aln+'@'+'dummyTrimMethod'] = align
            
    def trim(self, list_of_Conf_objects):
        for m in list_of_Conf_objects:
            m.timeit.append(time.time())
            m.platform = platform_report() 
            m.platform.append('Program and version: '+ m.cmd + ': ' +
                               os.popen(m.cmd + ' --version').readlines()[1])
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
                    [summary_lines, num_lines, num_undeter, num_collapsed_aln_seqs] = aln_summary(align)
                    summary = 'Alignment name: '+aln+'\n'
                    for line in summary_lines:
                        summary += line+'\n'
                    if num_lines < 4:
                        line = 'Alignment %s has less than 4 sequences and will be dropped'%aln
                        print line
                        summary += line+'\n'
                    elif num_undeter > 0:
                        line = 'Alignment %s has undetermined sequences and will be dropped'%aln
                        print line
                        summary += line+'\n'
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
            self.used_methods.append(m)

##############################################################################################
class AlnConf:
##############################################################################################
    
    """
    >>> coi = Locus('dna', 'CDS', 'coi', ['cox1','COX1','coi','COI','CoI'])
    >>> lsu = Locus('dna', 'rRNA', '28S', ['28s','28S','LSU rRNA','28S ribosomal RNA','28S large subunit ribosomal RNA'])
    >>> pj = Project([coi, lsu])
    
    # cline_str = muscle = AlnConf(pj, method_name='MuscleDefaults',
    #                                      cmd='muscle', program_name='muscle',
    #                                     cline_args=dict())
    """
    
    def __init__(self, pj, method_name='MafftLinsi', CDSAlign=True, program_name='mafft',
                 cmd='mafft', loci='all',
                 cline_args=dict(localpair=True, maxiterate=1000)):
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
        for locus in self.loci: 
            if len(pj.records_by_locus[locus.name]) < 4:
                print "%s have less than 4 sequences and will be dropped from this conf object. Don't use it in a concatenation"%locus.name
            else:
                mutable_loci_list.append(locus)
        self.loci = mutable_loci_list
        self.CDS_proteins = {}
        self.CDS_in_frame = {}
        self.aln_input_strings = {}
        self.command_lines = {}
        self.timeit = [time.asctime()]
        self.platform = []
        self.cmd = cmd
        if not program_name == 'mafft' and cmd == 'mafft':
            self.cmd = pj.defaults[program_name]
        elif program_name == 'mafft' and cmd != 'mafft' and not 'mafft' in cmd:
            self.cmd = pj.defaults['mafft']
            
        # make defalut input files
        if pj.records_by_locus == {}:
            pj.extract_by_locus()
        for key in pj.records_by_locus.keys():
            if key in [l.name for l in self.loci]:
                SeqIO.write(pj.records_by_locus[key], self.id+'_'+key+'.fasta', 'fasta')    
        for locus in self.loci:
            # put default input file filename and string in the AlnConf object
            input_filename=self.id+'_'+locus.name+'.fasta'
            self.aln_input_strings[locus.name] = [open(input_filename,'r').read()]
            # If CDS prepare reference protein input file and in frame CDS input file
            if locus.feature_type == 'CDS' and locus.char_type == 'dna' and self.CDSAlign: 
                self.CDS_proteins[locus.name] = []
                self.CDS_in_frame[locus.name] = []
                for record in pj.records:
                    for feature in record.features:
                        if (not feature.type == 'source' and 'gene' in feature.qualifiers.keys() and
                            feature.qualifiers['gene'][0] in locus.aliases):
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
        execution = str(self.timeit[3])
        plat = str(self.platform).replace(",",'\n').replace(']','').replace("'",'').replace('[','')
        return ("AlnConf named %s with ID %s",         
                "Loci: %s ",       
                "Executed on: %s",
                "Commands:",
                "%s",
                "Environment:",    
                "%s",
                "execution time:",
                "%s" %(self.method_name, str(self.id), loci_string, date, command_lines, plat,execution))  

##############################################################################################
class TrimalConf:
##############################################################################################
    def __init__(self, pj, method_name='gappyout', program_name='trimal',
                 cmd='default', alns='all', trimal_commands=dict(gappyout=True)):
        
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
        


def use_sh_support_as_branch_support(tree_filename):
    string = open(tree_filename,'r').read()
    string = re.sub(r'\[',r'[&&NHX:support=',string)
    t = Tree(string)
    t.dist=0
    t.write(outfile=tree_filename)
    #t.show()
    
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
            raise RuntimeError('There is no trimmed alignment for locus '+locus.name+' in concatenation '+concatenation.name)
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
                           '-p': random.randint(0,999),
                           '-x':  random.randint(0,999),
                           '-s': input_filename,
                           '-N': support_replicates,
                           '-n': tree_method.id+'_'+trimmed_alignment_name+'0',
                           '-m': model}
                      ],
                'fD_fb':[{'-f': 'D',
                          '-p': random.randint(0,999),
                          '-s': input_filename,
                          '-N': ML_replicates,
                           '-n': tree_method.id+'_'+trimmed_alignment_name+'0',
                           '-m': model},{'-f': 'b',
                                                       '-p': random.randint(0,999),
                                                       '-s': input_filename,
                                                       '-n': tree_method.id+'_'+trimmed_alignment_name+'1',
                                                       '-m': model,
                                                       '-t': 'RAxML_bestTree.'+tree_method.id+'_'+trimmed_alignment_name+'0',
                                                       '-z': 'RAxML_rellBootstrap.'+tree_method.id+'_'+trimmed_alignment_name+'0'}
                         ],
                'fd_b_fb':[{'-f': 'd',
                          '-p': random.randint(0,999),
                          '-s': input_filename,
                          '-N': ML_replicates,
                           '-n': tree_method.id+'_'+trimmed_alignment_name+'0',
                           '-m': model},{
                                                       '-p': random.randint(0,999),
                                                       '-b': random.randint(0,999),
                                                       '-s': input_filename,
                                                       '-#': support_replicates,
                                                       '-n': tree_method.id+'_'+trimmed_alignment_name+'1',
                                                       '-m': model,
                                                       '-T': tree_method.threads},
                                                       {'-f': 'b',
                                                       '-p': random.randint(0,999),
                                                       '-s': input_filename,
                                                       '-n': tree_method.id+'_'+trimmed_alignment_name+'2',
                                                       '-m': model,
                                                       '-t': 'RAxML_bestTree.'+tree_method.id+'_'+trimmed_alignment_name+'0',
                                                       '-z': 'RAxML_bootstrap.'+tree_method.id+'_'+trimmed_alignment_name+'1'}
                         ]
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


class RaxmlConf:
    
    
    def __init__(self, pj, method_name='fa', program_name='raxmlHPC-PTHREADS-SSE3',
                 cmd='default', preset = 'fa', alns='all', model='GAMMA', matrix='JTT', threads=4,
                 cline_args={}):
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



from pylab import *
import random

def draw_boxplot(dictionary, y_axis_label, figs_folder, scale): #'locus':[values]
    import numpy as np
    import matplotlib.pyplot as plt
    items = dictionary.items()
    items.sort()
    
    data = [locus[1] for locus in items]
        
    fig, ax1 = plt.subplots()
    fig.figsize = (scale, 3)
    #plt.subplots_adjust(left=0.075, right=0.95, top=0.9, bottom=0.25)

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
    xticks(range(len(data)+1)[1:], xlabels, size=14, rotation='vertical')
    #subplots_adjust(left=0.3, bottom=0.8)
    
    ax1.set_ylabel(y_axis_label, size=18)
    
    name = str(random.randint(1000,2000))
    fig.savefig(figs_folder + '/' + name +'.png')
    close('all')
    return figs_folder + '/' + name+'.png'
    
#################################################################################
def report_methods(pj, figs_folder, output_directory):
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

        # Here we manage mkdir and error raising for preexisting 
        # output_directory, based of the caller function
        if os.path.isdir(output_directory) and callername != 'publish':
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
        
        css_line = '<link rel="stylesheet" type="text/css" href="files/rp.css">'

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
        report_lines += ['<h2>','Data','</h2>', '']
        
        
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
            fig_filename = draw_boxplot(lengths_dict, 'Seq length (bp)', '%s/files'%output_directory, scale)
            
            
            # Distribution of sequence lengths
            #---------------------------------------------------------------------
            title = 'Distribution of sequence lengths'
            report_lines += ( '<h4>', title, '</h4>',  '')
            
            # This will write the img tag for the seq length boxplot in the report html
            # The src attribute is the png file path. The commented lines are an alternative
            # making an embeded figure.
            if os.path.isfile(fig_filename):
                #data_uri = open(fig_filename, 'rb').read().encode('base64').replace('\n', '')
                #img_tag = '<img height=400 width='+scale+' src="data:image/png;base64,{0}">'.format(data_uri)
                img_tag = '<img height=400 src="%s">'%(fig_filename.partition('/')[-1])
                report_lines.append(img_tag)
                #os.remove(fig_filename)
            
            
            # This will make GC content, nuc_degen_prop and prot_degen_prop png file,
            # will put them in the /files subdirectory and will write the html sections
            # for them, including img tags. All three params are feature qualifiers of
            # SeqRecord objects found in pj.records, which is a list. 
            for stat in ('GC_content', 'nuc_degen_prop', 'prot_degen_prop'):
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
                fig_filename = draw_boxplot(stat_dict, ylabel, '%s/files'%output_directory, scale)
                
                # Distribution of stat
                #---------------------------------------------------------------------
                title = 'Distribution of sequence statistic \"'+stat+'\"'
                report_lines += ( '<h4>', title, '</h4>', '')
                
                # This will make the img tag using the png path as src. The commented lines are an alternative
                # making an embeded image
                if os.path.isfile(fig_filename):
                    #data_uri = open(fig_filename, 'rb').read().encode('base64').replace('\n', '')
                    img_tag = '<img height=400 src="%s">'%(fig_filename.partition('/')[-1])
                    #img_tag = '<img height=400 width='+scale+' src="data:image/png;base64,{0}">'.format(data_uri)
                    report_lines.append(img_tag)
                    #os.remove(fig_filename)
                
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
            for locus in c.concat_must_have_all_of:
                rule_1 += locus + ', '
            report_lines.append(rule_1)
            rule_2 = 'OTUs must have at least one of each group: '
            for group in c.concat_must_have_one_of:
                rule_2 += str(group) +', '
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
        
        for method in pj.used_methods:
            # This will print list representations of the 'Conf' objects
            if isinstance(method,list) and(method[0] == 'AlnConf' or method[0] == 'RaxmlConf' or method[0] == 'TrimalConf'):
                title = method[0]
                report_lines += ('', '<h4>', title, '</h4>','')
                for i in method[1:]:
                    report_lines.append('<strong>'+str(i[0])+'</strong>')
                    report_lines.append(str(i[1]).replace(',','<br>'))
            
            # These will print attributes in actual 'Conf' objects 
            elif isinstance(method, AlnConf):
                title = 'Seuqence Alignment Method \"'+method.method_name+'\", method ID: '+method.id
                report_lines += ('', '<h4>', title, '</h4>', '')
                #--------------------------------------------------------
                align_line = 'Included loci :'
                for locus in [locus.name for locus in method.loci]:
                    align_line += locus + ', '
                report_lines.append(align_line)
                report_lines.append('Total execution time: '+str(method.timeit[3])+' sec\'')
                report_lines.append('Performed on: '+str(method.timeit[0]))
                report_lines += method.platform
                report_lines.append('')

                report_lines.append('Command lines:')
                for cline in method.command_lines.keys():
                    report_lines.append('Alignment \"'+cline+'\":')
                    report_lines.append('<pre style="white-space:normal;">')
                    report_lines.append(str(method.command_lines[cline]))
                    report_lines.append('</pre>')
                    
            elif isinstance(method, RaxmlConf):
                title = 'Raxml Tree Reconstruction Method \"'+method.method_name+'\", method ID: '+method.id
                report_lines += ('', '<h4>', title, '</h4>', '')
                #--------------------------------------------------------
                tree_line = 'Included alignments :'
                for aln in method.trimmed_alignments.keys():
                    tree_line += aln + ', '
                report_lines.append(tree_line)
                report_lines.append('Total execution time: '+str(method.timeit[3])+' sec\'')
                report_lines.append('Performed on: '+str(method.timeit[0]))
                report_lines += method.platform
                report_lines.append('')

                report_lines.append('Command lines:')
                for aln in method.command_lines.keys():
                    report_lines.append('Alignment \"'+aln+'\":')
                    report_lines.append('<pre style="white-space:normal;">')
                    for cline in method.command_lines[aln]:
                        report_lines.append(str(cline))
                    report_lines.append('</pre>')
                    report_lines.append('')
                    
            elif isinstance(method, TrimalConf):
                title = 'Trimal alignment trimming Method \"'+method.method_name+'\", method ID: '+method.id
                report_lines += ('<h4>', title, '</h4>', '')
                #--------------------------------------------------------
                aln_line = 'Included alignments :'
                for aln in method.alignments.keys():
                    aln_line += aln + ', '
                report_lines.append(aln_line)
                report_lines.append('Total execution time: '+str(method.timeit[3])+' sec\'')
                report_lines.append('Performed on: '+str(method.timeit[0]))
                report_lines += method.platform
                report_lines.append('')

                report_lines.append('Command lines:')
                for aln in method.command_lines.keys():
                    report_lines.append('Alignment \"'+aln+'\":')
                    report_lines.append('<pre style="white-space:normal;">')
                    report_lines.append(str(method.command_lines[aln]))
                    report_lines.append('</pre>')
                    report_lines.append('')
                
                
        report_lines += ['',''] 
        
        #############################  section 3:   RESULTS  #######################
        
        report_lines += ['', '<h2>','Results','</h2>', '']
        
        # Global alignmnet statistics
        #------------------------------------------------------------------------
        title = 'Global alignmnet statistics'.title()
        report_lines += ('<h3>', title, '</h3>', '')
        
        
        # This prints things like num of unique seqs and num of parsimony informative
        # cloumns. Takes the info from 'pj.aln_summaries' which is a list of strings.
        if len(pj.aln_summaries)>0:
            for summary in pj.aln_summaries:
                report_lines += ('<h4>', summary.splitlines()[0], '</h4>', '')
                report_lines += ['',summary.partition('\n')[-1],'']
        else:
            report_lines += ['','No sequence alignments in this Project','']
                
        
        # Per position alignmnet statistics
        #------------------------------------------------------------------------
        title = 'per position alignmnet statistics'.title()
        report_lines += ('<h3>', title, '</h3>', '')
        if len(pj.alignments.keys())>0:                    
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
            report_lines += ['No alignments in this project','']
            
        if len(pj.trimmed_alignments.keys())>0:          
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
            report_lines += ['No trimmed alignments in this project','']
        
        title = 'Robinson-Foulds distances'.title()
        report_lines += ('<h3>', title, '</h3>', '')
        
        if len(pj.trees.keys())>1:
            try:
                RF_filename, legend = calc_rf(pj, '%s/files'%output_directory)
                scale = str(len(legend)*60)
                if os.path.isfile(RF_filename):
                        #data_uri = open(RF_filename, 'rb').read().encode('base64').replace('\n', '')
                        #img_tag = '<img height='+scale+' width='+scale+' src="data:image/png;base64,{0}">'.format(data_uri)
                        img_tag = '<img src="%s">'%(RF_filename.partition('/')[-1])
                        report_lines.append(img_tag)
                        #os.remove(RF_filename)
                
                report_lines.append('<h3>Legend<h3><pre>')
                report_lines += legend
                report_lines.append('</pre>')
            except:
                report_lines += ['Found unrooted tree(s), skipping RF distance calculation']

        else:
            report_lines += ['Less than two trees in this Project','']

                
        #############################  section 4:   TREES  #######################
        
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
    
        
def pickle_pj(pj, pickle_file_name):
        import os
        if os.path.exists(pickle_file_name):
            os.remove(pickle_file_name)
        import cloud.serialization.cloudpickle as pickle
        output = open(pickle_file_name,'wb')
        pickle.dump(pj, output)
        output.close()
        if __builtin__.git:
            import rpgit
            rpgit.gitAdd(pickle_file_name)
            comment = "A pickled Project from %s" % time.asctime()
            rpgit.gitCommit(comment) 
            
        return pickle_file_name
    
def unpickle_pj(pickle_file_name):
        import cloud.serialization.cloudpickle as pickle
        pickle_handle = open(pickle_file_name, 'rb')
        pkl_pj = pickle.pickle.load(pickle_handle)
        new_pj = Project(pkl_pj.loci)
        attr_names = ['alignments',
                      'concatenations',
                      'records',
                      'records_by_locus',
                      'trees',
                      'trimmed_alignments',
                      ]
        
        for attr_name in attr_names:
           setattr(new_pj,attr_name,getattr(pkl_pj,attr_name))
            
        for i in pkl_pj.used_methods:
            if isinstance(i, list) and (i[0] == 'AlnConf' or i[0] == 'RaxmlConf' or i[0] == 'TrimalConf'):
                new_pj.used_methods.append(i)
            else:
                
                include = ['id',
                           'method_name',
                           'CDSAlign',
                           'program_name',
                           'loci',
                           'alns'
                           'command_lines',
                           'timeit',
                           'platform',
                           'cmd',
                           'preset',
                           'model',
                           'trimmed_alignments']
            
                method_list = []
                if isinstance(i, AlnConf):
                    method_list.append('AlnConf')
                elif isinstance(i, RaxmlConf):
                    method_list.append('RaxmlConf')
                elif isinstance(i, TrimalConf):
                    method_list.append('TrimalConf')
                for attr in include:
                    if attr in dir(i):
                        method_list.append([attr, str(getattr(i,attr))])
                new_pj.used_methods.append(method_list)
        return new_pj

def publish(pj, folder_name, figures_folder):
    
    import os, time
    folder = None
    zip_file = None
    if folder_name.endswith('.zip'):
        zip_file = folder_name
        folder = folder_name[:-4]
    else:
        folder = folder_name
        zip_file = folder_name + '.zip'
    if os.path.exists(folder) or os.path.exists(zip_file):
        raise IOError(folder_name + ' already exists')
    
    os.makedirs(folder)
    pj.write(folder+'/tree_and_alns.nexml','nexml')
    pj.write(folder+'/sequences_and_metadata.gb','genbank')
    report = open(folder+'/report.html','wt')
    for line in report_methods(pj, figures_folder, folder_name):
        report.write(line + '\n')
    report.close()

    #'report_lines' is now taking care of puting the figures in the zip folder, within /files
    #for tree in pj.trees.keys():
    #    if os.path.isfile(figures_folder+'/'+pj.trees[tree][0].get_leaves()[0].tree_method_id+'.png'):
    #        from shutil import copyfile
    #        copyfile(figures_folder+'/'+pj.trees[tree][0].get_leaves()[0].tree_method_id+'.png',
    #                 folder+'/'+pj.trees[tree][0].get_leaves()[0].tree_method_id+'.png')
            
         
    pickle_name = time.strftime("%a_%d_%b_%Y_%X", time.gmtime())+'.pkl'
    pickle_pj(pj, folder + '/' + pickle_name)

    
    import zipfile, shutil
    
    zf = zipfile.ZipFile(zip_file, "w")
    for dirname, subdirs, files in os.walk(folder):
        zf.write(dirname)
        for filename in files:
            zf.write(os.path.join(dirname, filename))
    zf.close()
    shutil.rmtree(folder)
    
def calc_rf(pj, figs_folder):
    meta = 'feature_id'
    if len(pj.concatenations) > 0:
        meta = pj.concatenations[0].otu_meta

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
        for t2 in trees:
            dupT2 = Tree(pj.trees[t2][0].write())
            for l in dupT2:
                for record in pj.records:
                    for feature in record.features:
                        if feature.qualifiers['feature_id'][0] == l.name and meta in feature.qualifiers.keys():
                            l.name = feature.qualifiers[meta][0]
            rf, max_rf, common_leaves, parts_t1, parts_t2 = dupT1.robinson_foulds(dupT2)        
            line.append(rf/float(max_rf))        
        data.append(line)   
    
    row_labels = [str(i) for i in range(len(trees))]
    column_labels = row_labels
    legend = ['#'.ljust(10,' ')+'LOCUS'.ljust(20,' ')+'ALIGNMENT METHOD'.ljust(20,' ')+'TRIMMING METHOD'.ljust(20,' ')+'TREE METHOD'.ljust(20,' ')]
    for i in trees:
        line = str(trees.index(i)).ljust(10,' ')
        for val in i.split('@'):
            line += val.ljust(20,' ')
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
    fig.colorbar(heatmap, cmap=plt.cm.Blues)
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

if __name__ == "__main__":
    import doctest
    doctest.testmod()
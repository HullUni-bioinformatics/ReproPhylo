
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
    
    
    Dependencies:
    Python 2.7.6
    ete2 2.2rev1056
    biopython 1.64
    dendropy 3.12.0
    cloud 2.8.5
    
    RAxML 8
    Muscle
    Mafft
    Pal2nal
    """
##############################################################################################



from Bio import SeqIO
import os, csv, sys, dendropy, re, time, random, glob, platform#, subprocess
from Bio.Seq import Seq
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



##############################################################################################
class Locus:
##############################################################################################

    """ Configure the loci stored in the ReproPhylo DB.
        
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
    ...                                concat_must_have_all_of=[coi],
    ...                                concat_must_have_one_of =[[bssu,lsu],[alg11,ssu]],
    ...                                define_trimmed_alns=["MuscleDefaults@dummyTrimMethod"])
    >>> print(str(concatenation))
    Concatenation named combined, with loci coi,18S,16S,28S,ALG11,
    of which coi must exist for all species
    and at least one of each group of [ 16S 28S ][ ALG11 18S ] is represented.
    Alignments with the following names: MuscleDefaults@dummyTrimMethod are prefered
    """
    
    name = 'NotSet'
    loci = []
    otu_meta = 'NotSet'
    concat_must_have_all_of = []
    concat_must_have_one_of = []
    define_trimmed_alns = [] #should be Locus_name@Alignment_method_name@Trimming_mathod_name
    
    feature_id_dict = {}
    
    def __init__(self,
                 name = name,
                 loci = loci,
                 otu_meta = otu_meta,
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
            must_have += i.name+','
        must_have = must_have[:-1]
        trimmed_alignmnets_spec = ''
        one_of = ''
        for i in self.concat_must_have_one_of:
            one_of += '[ '
            for j in i:
                one_of += j.name+' '
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
    Reprophylo Database Utilities
    
    Used in the Database class but are not in the classes methods
    """
##############################################################################################



def platform_report():
    
    """ 
    Prints machine specs, os specs and dependencies at time of execution
    
    >>> isinstance(platform_report(), list)
    True
    """
    import pkg_resources
    modules = []
    for i in ('ete2','biopython','dendropy','cloud'):
        modules.append(i+' version: '+
                            pkg_resources.get_distribution(i).version)
    return(['Platform: '+platform.platform(aliased=0, terse=0),
            'Processor: '+platform.processor(),
            'Python build: '+platform.python_build()[0] + platform.python_build()[1],
            'Python compiler: '+platform.python_compiler(),
            'Python implementation: ' +platform.python_implementation(),
            'Python version: ' + platform.python_version()]+
             modules+
            ['User: ' +platform.uname()[1]])



def write_alns(db, format = 'fasta'):
    """
    Writes untrimmed sequence alignment files that are in db in a biopython format
    """
    
    if len(db.alignments.keys()) == 0:
        raise IOError('Align the records first')
    else:
        for key in db.alignments:
            AlignIO.write(db.alignments[key], key+'_aln.'+format, format)



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
    
    # testing if fits any of the Database Locus objects
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
    Handles list printing as a nice string in the db.write(format="csv") method
    
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
    Replaces newline with space in the db.write(format="csv") method
    """
    
    lines = lines.split('\n')
    return (' ').join(lines)




def type_to_single_line_str(var):
    
    """
    Returns any type as a one line string for the db.write(format="csv") method
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



def get_qualifiers_dictionary(database, feature_id):
    
    """
    Takes sequence record annotation, source qualifiers and feature qualifiers and puts them
    in a flat dictionary
    
        
    # Making a dummy locus    
    >>> coi = Locus('dna','CDS','coi', ['cox1','COX1','coi','COI','CoI'])
    
    # Making a dummy Database
    >>> db = Database([coi])
    
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
    >>> db.records = [record]
    >>> qual_dict = get_qualifiers_dictionary(db, '12345')
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
    for record in database.records:
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




##############################################################################################
class Database:
##############################################################################################
    

    def __init__(self, loci):
        self.records = []
        self.loci = loci
        self.records_by_locus = {}
        self.concatenations = []
        self.alignments = {}
        self.trimmed_alignments = {}
        self.trees = {}
        self.used_methods = []
        seen = []
        for locus in loci:
            if not isinstance(locus, Locus):
                raise TypeError("Expecting Locus object in loci list. "+locus+
                                " not a Locus object")
            if locus.name in seen:
                raise NameError('Locus ' + locus.name + ' apears more than once in self.loci')
            else:
                seen.append(locus.name)

    def read_embl_genbank(self, input_filenames_list):
        generators = []
        for input_filename in input_filenames_list:
            generators.append(parse_input(input_filename, 'gb'))
#            if is_embl_or_gb(input_filename):
#                generators.append(parse_input(input_filename, 'gb'))
#                print 'to do: correct parse_input'
#                print 'Accepted:', input_filename
#            else:
#                print ('Rejected: ' + input_filename + ' should be genbank or embl and end with .gb or .embl')
            for generator in generators:
                for record in generator:
                    dwindled_record = dwindle_record(record, self.loci)
                    if len(record.features) > 1:
                        self.records.append(dwindled_record)
                    elif len(record.features) == 1 and not record.features[0].type == 'source':
                        self.records.append(dwindled_record)

#    def read_denovo(self, input_filename, feature_type, char_type, source_qualifiers = {}):
#        count = 0
#        # start the counter where it stoped the last time we read denovo things
#        for record in self.record:
#            if 'denovo' in record.id:
#                serial = int(record.id[6:])
#                if serial > count:
#                    count = serial+1
#        denovo = SeqIO.parse(input_filename, guess_format(input_filename))
#        for record in denovo:
#            feature = SeqFeature(FeatureLocation(0, len(record.seq)), type=feature_type, strand=1)
#            source = SeqFeature(FeatureLocation(0, len(record.seq)), type='source', strand=1)
#            if len(source_qualifiers.keys())>0:
#                for key in source_qualifiers.keys():
#                    source.qualifiers[key] = source_qualifiers[key]
#            feature.qualifiers['original_id'] = [record.id]
#            feature.qualifiers['original_desc'] = [(' ').join(record.description.split()[1:])]
#            record.id = 'denovo'+str(count)
#            count += 1
#            feature.qualifiers['feature_id'] = [record.id+'_f0']
#            source.qualifiers['feature_id'] = [record.id+'_source']
#            record.features = [source, feature]
#            if char_type == 'prot':
#                record.seq.alphabet = IUPAC.protein
#            elif char_type == 'dna':
#                record.seq.alphabet = IUPAC.ambiguous_dna
#            self.records.append(record)
            
    def read_denovo(self, input_filenames, char_type):
        count = 0
        # start the counter where it stoped the last time we read denovo things
        for record in self.records:
            if 'denovo' in record.id:
                serial = int(record.id[6:])
                if serial > count:
                    count = serial+1
        for input_filename in input_filenames:
            # Check input file format
            #suffix = re.search(r'\.([^\.]+)$',input_filename).groups()[1]
            #input_format = seq_format_from_suffix(suffix)
            #if input_format == 'embl' or input_format == 'genbank':
            #    raise IOError('To read embl or genbank files use self.read_embl_genbank([filename0, filename1])')
            #denovo = SeqIO.parse(input_filename, input_format)
            denovo = SeqIO.parse(input_filename, 'fasta')
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
                elif char_type == 'dna':
                    record.seq.alphabet = IUPAC.ambiguous_dna
                count += 1
                self.records.append(record)
                
    def add_feature_to_record(self, record_id, feature_type, location='full', qualifiers={}):
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
                record.features.append(feature)

    
    def add_concatenation(self, concatenation_object):
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
                    available_features[locus.name].append(record.id) # record ids are feature ids because taken from db.records_by_locus
                    
            
            
                    
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
                all_locus_trimmed_alns_in_db = []
                for key in self.trimmed_alignments.keys():
                    if locus.name == key.split('@')[0]:
                        all_locus_trimmed_alns_in_db.append(key)
                if len(all_locus_trimmed_alns_in_db) == 1:
                    trimmed_aln = all_locus_trimmed_alns_in_db[0]
                elif len(all_locus_trimmed_alns_in_db) == 0:
                    raise RuntimeError('Locus '+locus.name+' have no trimmed alignments')
                else:
                    s.define_trimmed_alns.sort(key = lambda i: i.count('@'), reverse=True)
                    for definition in s.define_trimmed_alns:
                        if definition.count('@') == 2 and locus.name == definition.split('@')[0]:
                            trimmed_aln = definition
                        elif definition.count('@') == 1 and bool(any(definition in i for i in all_locus_trimmed_alns_in_db)):
                            trimmed_aln = locus.name+'@'+definition
                        else:
                            raise RuntimeError("Could not determine which alignment/trimming alternative to use for locus '"+
                                                locus.name+"' out of "+str(locus_trimmed_alns))
                if trimmed_aln:
                    keys_of_trimmed_alignments_to_use_in_concat.append(trimmed_aln)
                else:
                    raise RuntimeError('Could not find trimmed aln for locus '+locus.name+' given the rulls '+str(s.define_trimmed_alns))

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

   


    def write(self, filename, format = 'genbank'):
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
                linewriter.writerow(['record_id','seq']+source_qualifiers+['taxonomy']+feature_qualifiers)
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

                            

                
    def if_this_then_that(self, IF_THIS, IN_THIS, THEN_THAT, IN_THAT, mode = 'whole'):
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
            for organism in species_vs_loci.keys():
                line = [organism]
                for name in loci_names:
                    if name in species_vs_loci[organism].keys():
                        line.append(str(species_vs_loci[organism][name]))
                    else:
                        line.append('0')
                linewriter.writerow(line)
                    
            
        
    def extract_by_locus(self):
        data_by_locus = {}
        for locus in self.loci:
            if not locus.name in locus.aliases:
                locus.aliases.append(locus.name)
            records = []
            for record in self.records:
                for feature in record.features:
                    if (feature.type == locus.feature_type and
                        
                        (('gene' in feature.qualifiers.keys() and
                          feature.qualifiers['gene'][0] in locus.aliases) 
                         or
                         ('product' in feature.qualifiers.keys() and 
                          feature.qualifiers['product'][0] in locus.aliases))
                        
                        ):
                        if locus.char_type == 'dna':
                            S = feature.extract(record.seq)
                        elif locus.char_type == 'prot':
                            S = Seq(feature.qualifiers['translation'][0], IUPAC.protein)
                        feature_record = SeqRecord(seq = S, id = feature.qualifiers['feature_id'][0],
                                                   description = '')
                        records.append(feature_record)
            data_by_locus[locus.name] = records
        self.records_by_locus = data_by_locus


    def write_by_locus(self, format = 'fasta'):
        if self.records_by_locus == {}:
            self.extract_by_locus
        for key in self.records_by_locus.keys():
            SeqIO.write(self.records_by_locus[key], key+'.'+format, format)
            
    def align(self, alignment_methods=[], pal2nal='./pal2nal.pl'):
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
                        #raise RuntimeError('locus '+locus.name+' is in more than one AlignmentMethod objects')
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
                    self.alignments[locus.name+'@'+method.method_name] = align
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
            
    def tree(self, raxml_methods):
        # to do: determine the program used and the resulting expected tree file name
        
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
        
        #for aln_name in self.alignments.keys():
        #    if aln_name in loci_names:
        #        char_type = ''
        #        matrix = None
        #        for locus in self.loci:
        #            if locus.name == aln_name:
        #                char_type = locus.char_type
        #        if char_type == 'dna':
        #            matrix = dendropy.DnaCharacterMatrix()
        #        elif char_type == 'prot':
        #            matrix = dendropy.ProteinCharacterMatrix()
        #        matrix_string = self.alignments[aln_name].format('fasta')
        #        matrix.read_from_string(matrix_string,'fasta')
        #        D.add_char_matrix(matrix)
            
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
                 
                 multifurc=None
                 ): 
    
            print '<html>'
            ts = TreeStyle()
            ts.show_leaf_name = False
            ts.scale = 1000
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
            
            #if heat_map_meta:
            #    ts.legend_position=1
            #    ts.legend.add_face(TextFace('Heatmap columns: ', fsize=10), column=0)
            #    i = 1
            #    for meta in heat_map_meta:
            #        ts.legend.add_face(TextFace(' '+str(meta)+', ',
            #                                    fsize=20), column=i)
            #        i += 1
                
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
                    if not root_value == 'mid' and root_value in qualifiers_dictionary[root_meta]:
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
                print('<A href=file://'+
                       fig_folder + "/" + self.trees[tree][0].get_leaves()[0].tree_method_id+'.png'+
                       '>'+self.trees[tree][0].get_leaves()[0].tree_method_id+
                       '</A><BR>')
            print '</html>'
            print fig_folder

    def trim(self):
        for aln in self.alignments.keys():
            AlignIO.write(self.alignments[aln],aln+'_aln.fasta','fasta')
            stdout = os.popen('trimal -in '+ aln +'_aln.fasta -gappyout').read()
            align = AlignIO.read(StringIO(stdout), "fasta",  alphabet=IUPAC.ambiguous_dna)
            for record in align:
                record.description = ''
            self.trimmed_alignments[aln+'@'+'dummyTrimMethod'] = align
            
            
class AlignmentMethod:
    
    def __init__(self, db, method_name='MafftLinsi', CDSAlign=True, program_name='mafft',
                 cmd='mafft', loci='all',
                 cline_args=dict(localpair=True, maxiterate=1000)):
        self.id = str(random.randint(10000,99999))+str(time.time())
        self.method_name=method_name
        self.CDSAlign=CDSAlign
        self.program_name=program_name
        self.loci = db.loci
        if not loci == 'all':
            self.loci = []
            for locus_name in loci:
                for locus in db.loci:
                    if locus_name == locus.name:
                        self.loci.append(locus)
        self.CDS_proteins = {}
        self.CDS_in_frame = {}
        self.aln_input_strings = {}
        self.command_lines = {}
        self.timeit = [time.asctime()]
        self.platform = []
        self.cmd = cmd
        # make defalut input files
        if db.records_by_locus == {}:
            db.extract_by_locus()
        for key in db.records_by_locus.keys():
            SeqIO.write(db.records_by_locus[key], self.id+'_'+key+'.fasta', 'fasta')    
        for locus in self.loci:
            # put default input file filename and string in the AlignmentMethod object
            input_filename=self.id+'_'+locus.name+'.fasta'
            self.aln_input_strings[locus.name] = [open(input_filename,'r').read()]
            # If CDS prepare reference protein input file and in frame CDS input file
            if locus.feature_type == 'CDS' and locus.char_type == 'dna' and self.CDSAlign: 
                self.CDS_proteins[locus.name] = []
                self.CDS_in_frame[locus.name] = []
                for record in db.records:
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
                            # Put the protein records in the AlignmentMethod object
                            self.CDS_proteins[locus.name].append(feature_record)
                           
                                    
                # check same number of prot and cds objects                    
                if len(db.records_by_locus[locus.name]) > len(self.CDS_proteins[locus.name]):
                    raise RuntimeError('For the CDS locus '+locus.name+': more nuc seqs than prot seqs.'+
                                       ' You may miss a \'translate\' or \'gene\' qualifier in some of '+
                                       'the features.')
                elif len(db.records_by_locus[locus.name]) < len(self.CDS_proteins[locus.name]):
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
                self.command_lines[locus.name] = MafftCommandline(cmd=cmd)
            elif self.program_name == 'muscle':
                self.command_lines[locus.name] = MuscleCommandline(cmd=cmd)
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
        return """
AlignmentMethod named %s with ID %s         
Loci: %s        
Executed on: %s

Commands:
%s

Environment:    
%s

execution time:
%s""" %(self.method_name, str(self.id), loci_string, date, command_lines, plat,execution)  



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
    
    
def make_raxml_partfile(tree_method, db, trimmed_alignment_name):

    concatenation = None
    for c in db.concatenations:
        if c.name == trimmed_alignment_name:
            concatenation = c
    
    #concatenation = filter(lambda concatenation: concatenation.name == trimmed_alignment_name, db.concatenations)[0]
    
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

def write_raxml_clines(tree_method, db, trimmed_alignment_name):
            
    cline_que = 0

    support_replicates = 100
    ML_replicates = 1
    if '-N' in tree_method.cline_args.keys():
        ML_replicates = tree_method.cline_args['-N']
    if '-#' in tree_method.cline_args.keys():
        support_replicates = tree_method.cline_args['-#']
    
    partfile = None
    
    # Check if it is a concatenation and make partfile

    for c in db.concatenations:
        if c.name == trimmed_alignment_name.partition('@')[0]:
            partfile = make_raxml_partfile(tree_method, db, trimmed_alignment_name)
    
    input_filename = make_raxml_input_matrix_file(tree_method, trimmed_alignment_name)
    model = tree_method.model
    try:
        locus_char_type = filter(lambda locus: locus.name == trimmed_alignment_name.partition('@')[0], db.loci)[0].char_type
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
                           '-m': model,
                           '-T': tree_method.threads}
                      ],
                'fD_fb':[{'-f': 'D',
                          '-p': random.randint(0,999),
                          '-s': input_filename,
                          '-N': ML_replicates,
                           '-n': tree_method.id+'_'+trimmed_alignment_name+'0',
                           '-m': model,
                           '-T': tree_method.threads},{'-f': 'b',
                                                       '-p': random.randint(0,999),
                                                       '-s': input_filename,
                                                       '-n': tree_method.id+'_'+trimmed_alignment_name+'1',
                                                       '-m': model,
                                                       '-T': tree_method.threads,
                                                       '-t': 'RAxML_bestTree.'+tree_method.id+'_'+trimmed_alignment_name+'0',
                                                       '-z': 'RAxML_rellBootstrap.'+tree_method.id+'_'+trimmed_alignment_name+'0'}
                         ],
                'fd_b_fb':[{'-f': 'd',
                          '-p': random.randint(0,999),
                          '-s': input_filename,
                          '-N': ML_replicates,
                           '-n': tree_method.id+'_'+trimmed_alignment_name+'0',
                           '-m': model,
                           '-T': tree_method.threads},{
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
                                                       '-T': tree_method.threads,
                                                       '-t': 'RAxML_bestTree.'+tree_method.id+'_'+trimmed_alignment_name+'0',
                                                       '-z': 'RAxML_bootstrap.'+tree_method.id+'_'+trimmed_alignment_name+'1'}
                         ]
                }

    if partfile:
        for preset in presets.keys():
            for cline in range(len(presets[preset])):
                presets[preset][cline] = dict({'-q': partfile}, **presets[preset][cline])               
    return presets[tree_method.preset] 

class RaxmlTreeReconstructionMethod:
    
    
    def __init__(self, db, method_name='fa', program_name='raxmlHPC-PTHREADS-SSE3',
                 cmd='raxmlHPC-PTHREADS-SSE3', preset = 'fa', alns='all', model='GAMMA', matrix='JTT', threads=4,
                 cline_args={}):
        self.id = str(random.randint(10000,99999))+str(time.time())
        self.method_name=method_name
        self.program_name=program_name
        self.preset = preset
        self.cline_args = cline_args
        self.model = model
        self.matrix = matrix
        self.threads = threads
        self.trimmed_alignments = db.trimmed_alignments
        if not alns == 'all':
            self.trimmed_alignments = {}
            for aln_name in alns:
                if aln_name in db.trimmed_alignments.keys():
                    self.trimmed_alignments[aln_name] = db.trimmed_alignments[aln_name]
        self.aln_input_strings = {}
        self.command_lines = {}
        self.timeit = [time.asctime()]
        self.platform = []
        self.cmd = cmd
        
        for trimmed_alignment in self.trimmed_alignments.keys():
            self.command_lines[trimmed_alignment] = []
            command_lines = write_raxml_clines(self, db, trimmed_alignment)
            for command_line in command_lines:
                cline_object = RaxmlCommandline(cmd=cmd)
                for c in command_line.keys():
                    cline_object.__setattr__(c,command_line[c])
                self.command_lines[trimmed_alignment].append(cline_object)
                print str(cline_object)



from pylab import *
import random

def draw_boxplot(dictionary, y_axis_label, figs_folder): #'locus':[values]
    import numpy as np
    import matplotlib.pyplot as plt
    items = dictionary.items()
    items.sort()
    
    data = [locus[1] for locus in items]
        
    fig, ax1 = plt.subplots()
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
    
def report_methods(db, figs_folder):
        report_lines = ['<html>','<head>','<h1>']
    
        head = 'reprophylo analysis from '+str(time.asctime())
        #=====================================================
        report_lines.append(head)
        report_lines += ['</h1>','</head>','<body>','']
        
        report_lines += ['<h2>','Data','</h2>', '']
        
        title = 'species representation in sequence data'.title()
        report_lines += ('<h3>', title, '</h3>', '', '<pre>')
        #--------------------------------------------------------
        
        outfile_name= str(random.randint(1000,2000))
        db.species_vs_loci(outfile_name)
        with open(outfile_name, 'rb') as csvfile:
            sp_vs_lc = list(csv.reader(csvfile, delimiter='\t', quotechar='|'))
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
                report_lines.append(string)
        
        os.remove(outfile_name)
        
        
        if len(db.records_by_locus.keys())>0:
            scale = str(len(db.records_by_locus.keys())*200)
            lengths_dict = {}
            for locus_name in db.records_by_locus.keys():
                lengths_dict[locus_name] = []
                for record in db.records_by_locus[locus_name]:
                    lengths_dict[locus_name].append(len(record.seq))
            fig_filename = draw_boxplot(lengths_dict, 'Seq length (bp)', figs_folder)
            title = 'Distribution of sequence lengths'
            report_lines += ('</pre>', '<h3>', title, '</h3>', '<pre>', '')
            if os.path.isfile(fig_filename):
                data_uri = open(fig_filename, 'rb').read().encode('base64').replace('\n', '')
                img_tag = '<img height=400 width='+scale+' src="data:image/png;base64,{0}">'.format(data_uri)
                report_lines.append(img_tag)
                os.remove(fig_filename)
            
            for stat in ('GC_content', 'nuc_degen_prop', 'prot_degen_prop'):
                stat_dict = {}
                ylabel = 'GC ontent (%)'
                if not stat == 'GC_content':
                    ylabel = 'Aambiguous positions (prop)'
                for locus_name in db.records_by_locus.keys():
                    stat_dict[locus_name] = []
                    for i in db.records_by_locus[locus_name]:
                        for record in db.records:
                            for feature in record.features:
                                if feature.qualifiers['feature_id'][0] == i.id:
                                    if stat in feature.qualifiers.keys():
                                        stat_dict[locus_name].append(float(feature.qualifiers[stat][0]))
                fig_filename = draw_boxplot(stat_dict, ylabel, figs_folder)
                title = 'Distribution of sequence statistic \"'+stat+'\"'
                report_lines += ('</pre>', '<h3>', title, '</h3>', '<pre>', '')
                if os.path.isfile(fig_filename):
                    data_uri = open(fig_filename, 'rb').read().encode('base64').replace('\n', '')
                    img_tag = '<img height=400 width='+scale+' src="data:image/png;base64,{0}">'.format(data_uri)
                    report_lines.append(img_tag)
                    os.remove(fig_filename)
                    
        composed_concatenations = []
        for c in db.concatenations:
            if c.name in db.trimmed_alignments.keys():
                composed_concatenations.append(c)
        
        
        for c in composed_concatenations:
            
            title = ('content of concatenation \"' + c.name + '\"').title()
            report_lines += ('</pre>', '<h3>', title, '</h3>', '<pre>', '')
            #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            
            report_lines.append('Rules for  \"' + c.name + '\":')
            rule_1 = 'OTUs must have the loci: '
            for locus in c.concat_must_have_all_of:
                rule_1 += locus + ', '
            report_lines.append(rule_1)
            rule_2 = 'OTUs must have at least one of each group: '
            for group in c.concat_must_have_one_of:
                rule_2 += str(group) +', '
            report_lines += (rule_2, '')
            
            
            otus = c.feature_id_dict.keys()
            loci = [locus.name for locus in c.loci]
            
            otus_max_length = max([len(i) for i in otus])+33
            loci_columns_max_length = []
            
            for locus in loci:
                lengths = [len(locus)]
                for otu in otus:
                    if locus in c.feature_id_dict[otu].keys():
                        lengths.append(len(c.feature_id_dict[otu][locus]))
                    else:
                        lengths.append(0)
                loci_columns_max_length.append(max(lengths)+3)
                
            concat_header = ''.ljust(otus_max_length)
            for i in range(len(loci)):
                concat_header += loci[i].ljust(loci_columns_max_length[i])
            report_lines += (concat_header, '~'*len(concat_header))
                
            for otu in otus:
                otu_species = ''
                for locus in loci:
                    if locus in c.feature_id_dict[otu].keys():
                        feature_qualifiers = get_qualifiers_dictionary(db, c.feature_id_dict[otu][locus])
                        if 'source_organism' in feature_qualifiers.keys():
                            otu_species = feature_qualifiers['source_organism']
                    
                concat_line = (otu+' '+otu_species).ljust(otus_max_length)
                for i in range(len(loci)):
                    if loci[i] in c.feature_id_dict[otu].keys():
                        concat_line += c.feature_id_dict[otu][loci[i]].ljust(loci_columns_max_length[i])
                    else:
                        concat_line += ''.ljust(loci_columns_max_length[i])
                report_lines.append(concat_line)

        report_lines += ['</pre>', '<h2>','Methods','</h2>', '<pre>', '']
        
        
        for method in db.used_methods:
            if method[0] == 'AlignmentMethod' or method[0] == 'RaxmlTreeReconstructionMethod':
                title = method[0]
                report_lines += ('</pre>', '<h3>', title, '</h3>', '<pre>', '')
                for i in method[1:]:
                    report_lines.append('<strong>'+str(i[0])+'</strong>')
                    report_lines.append(str(i[1]).replace(',','<br>'))
            
            elif isinstance(method, AlignmentMethod):
                title = 'Seuqence Alignment Method \"'+method.method_name+'\", method ID: '+method.id
                report_lines += ('</pre>', '<h3>', title, '</h3>', '<pre>', '')
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
                    
            elif isinstance(method, RaxmlTreeReconstructionMethod):
                title = 'Raxml Tree Reconstruction Method \"'+method.method_name+'\", method ID: '+method.id
                report_lines += ('</pre>', '<h3>', title, '</h3>', '<pre>', '')
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
                    for cline in method.command_lines[aln]:
                        report_lines.append('<pre style="white-space:normal;">')
                        report_lines.append(str(cline))
                        report_lines.append('</pre>')
                    report_lines.append('')
                
        
        
        report_lines += ('</pre>', '','') 
        
        
        if len(db.alignments.keys())>0:                    
            title = 'Alignment statistics before trimming'
            report_lines += ('</pre>', '<h3>', title, '</h3>', '<pre>', '')
            report_lines += ['</pre>', '<h4>','Trimal\'s Residue Similarity Score (-scc)','</h4>', '<pre>', '']
            fig_file = draw_trimal_scc(db, 2, figs_folder, trimmed=False)
            if os.path.isfile(fig_file):
                    data_uri = open(fig_file, 'rb').read().encode('base64').replace('\n', '')
                    img_tag = '<img src="data:image/png;base64,{0}">'.format(data_uri)
                    report_lines.append(img_tag)
                    os.remove(fig_file)
            report_lines += ['</pre>', '<h4>','Trimal\'s column gap gcore (-sgc)','</h4>', '<pre>', '']
            fig_file = draw_trimal_scc(db, 2, figs_folder, trimmed=False, alg='-sgc')
            if os.path.isfile(fig_file):
                    data_uri = open(fig_file, 'rb').read().encode('base64').replace('\n', '')
                    img_tag = '<img src="data:image/png;base64,{0}">'.format(data_uri)
                    report_lines.append(img_tag)
                    os.remove(fig_file)
               
        if len(db.trimmed_alignments.keys())>0:          
            title = 'Alignment statistics after trimming'
            report_lines += ('</pre>', '<h3>', title, '</h3>', '<pre>', '')
            report_lines += ['</pre>', '<h4>','"Trimal\'s Residue Similarity Score (-scc)','</h4>', '<pre>', '']
            fig_file = draw_trimal_scc(db, 2, figs_folder, trimmed=True)
            if os.path.isfile(fig_file):
                    data_uri = open(fig_file, 'rb').read().encode('base64').replace('\n', '')
                    img_tag = '<img src="data:image/png;base64,{0}">'.format(data_uri)
                    report_lines.append(img_tag)
                    os.remove(fig_file)
            report_lines += ['</pre>', '<h4>','Trimal\'s column gap gcore (-sgc)','</h4>', '<pre>', '']
            fig_file = draw_trimal_scc(db, 2, figs_folder, trimmed=True, alg='-sgc')
            if os.path.isfile(fig_file):
                    data_uri = open(fig_file, 'rb').read().encode('base64').replace('\n', '')
                    img_tag = '<img src="data:image/png;base64,{0}">'.format(data_uri)
                    report_lines.append(img_tag)
                    os.remove(fig_file)
        
        if len(db.trees.keys())>1:        
            report_lines += ('','<h1>Robinson-Foulds distances </h1>','')
            RF_filename, legend = calc_rf(db, figs_folder)
            scale = str(len(legend)*60)
            if os.path.isfile(RF_filename):
                    data_uri = open(RF_filename, 'rb').read().encode('base64').replace('\n', '')
                    img_tag = '<img height='+scale+' width='+scale+' src="data:image/png;base64,{0}">'.format(data_uri)
                    report_lines.append(img_tag)
                    os.remove(RF_filename)
            
            report_lines.append('<h3>Legend<h3><pre>')
            report_lines += legend
            report_lines.append('</pre>')
        
        report_lines += ('<h1>Trees</h1>','')
        
        
        for tree in db.trees.keys():
            report_lines += ('<h2>'+tree.split('@')[0]+'</h2>',
                             '<h3>Alignment method: '+tree.split('@')[1]+'</h3>',
                             '<h3>Trimming method: '+tree.split('@')[2]+'</h3>',
                             '<h3>Tree method: '+tree.split('@')[3]+'</h3>',
                             '<pre style="white-space:normal;">',
                             'Tree Method ID: '+db.trees[tree][0].get_leaves()[0].tree_method_id,'</pre>')
            
            report_lines += ('<h3>newick format</h3>','','<pre style="white-space:normal;">',db.trees[tree][0].write(),'</pre>','')
            report_lines += ('<h3>nhx format</h3>','','<pre>',db.trees[tree][1],'</pre>','','','','')
            
            
            
            if os.path.isfile(figs_folder+'/'+db.trees[tree][0].get_leaves()[0].tree_method_id+'.png'):
                data_handle = open(figs_folder+'/'+db.trees[tree][0].get_leaves()[0].tree_method_id+'.png','rb')
                data_uri = data_handle.read().encode('base64').replace('\n', '')
                data_handle.close()
                img_tag = '<img width=500 src="data:image/png;base64,{0}">'.format(data_uri)
                report_lines.append(img_tag)
                
        report_lines.append('</body>')
        report_lines.append('</html>')
            
        return report_lines
    
        
def pickle_db(db, pickle_file_name):
        import os
        if os.path.exists(pickle_file_name):
            os.remove(pickle_file_name)
        import cloud.serialization.cloudpickle as pickle
        output = open(pickle_file_name,'wb')
        pickle.dump(db, output)
        output.close()
        return pickle_file_name
    
def unpickle_db(pickle_file_name):
        import cloud.serialization.cloudpickle as pickle
        pickle_handle = open(pickle_file_name, 'rb')
        pkl_db = pickle.pickle.load(pickle_handle)
        new_db = Database(pkl_db.loci)
        attr_names = ['alignments',
                      'concatenations',
                      'records',
                      'records_by_locus',
                      'trees',
                      'trimmed_alignments',
                      ]
        
        for attr_name in attr_names:
           setattr(new_db,attr_name,getattr(pkl_db,attr_name))
            
        for i in pkl_db.used_methods:
            if isinstance(i, list) and (i[0] == 'AlignmentMethod' or i[0] == 'RaxmlTreeReconstructionMethod'):
                new_db.used_methods.append(i)
            else:
                
                include = ['id',
                           'method_name',
                           'CDSAlign',
                           'program_name',
                           'loci',
                           'command_lines',
                           'timeit',
                           'platform',
                           'cmd',
                           'preset',
                           'model',
                           'trimmed_alignments']
            
                method_list = []
                if isinstance(i, AlignmentMethod):
                    method_list.append('AlignmentMethod')
                elif isinstance(i, RaxmlTreeReconstructionMethod):
                    method_list.append('RaxmlTreeReconstructionMethod')
                for attr in include:
                    if attr in dir(i):
                        method_list.append([attr, str(getattr(i,attr))])
                new_db.used_methods.append(method_list)
        return new_db

def publish(db, folder_name, figures_folder):
    
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
    db.write(folder+'/tree_and_alns.nexml','nexml')
    db.write(folder+'/sequences_and_metadata.gb','genbank')
    report = open(folder+'/report.html','wt')
    for line in report_methods(db, figures_folder):
        report.write(line + '\n')
    report.close()

    for tree in db.trees.keys():
        if os.path.isfile(figures_folder+'/'+db.trees[tree][0].get_leaves()[0].tree_method_id+'.png'):
            from shutil import copyfile
            copyfile(figures_folder+'/'+db.trees[tree][0].get_leaves()[0].tree_method_id+'.png',
                     folder+'/'+db.trees[tree][0].get_leaves()[0].tree_method_id+'.png')
            
         
    pickle_name = time.strftime("%a_%d_%b_%Y_%X", time.gmtime())+'.pkl'
    pickle_db(db, folder + '/' + pickle_name)

    
    import zipfile, shutil
    
    zf = zipfile.ZipFile(zip_file, "w")
    for dirname, subdirs, files in os.walk(folder):
        zf.write(dirname)
        for filename in files:
            zf.write(os.path.join(dirname, filename))
    zf.close()
    shutil.rmtree(folder)
    
def calc_rf(db, figs_folder):
    meta = 'feature_id'
    if len(db.concatenations) > 0:
        meta = db.concatenations[0].otu_meta

    trees = db.trees.keys()

    data = []

    for t1 in trees:
        line = []
        dupT1 = db.trees[t1][0].copy()
        for l in dupT1:
            for record in db.records:
                for feature in record.features:
                    if feature.qualifiers['feature_id'][0] == l.name and meta in feature.qualifiers.keys():
                        l.name = feature.qualifiers[meta][0]
        for t2 in trees:
            dupT2 = db.trees[t2][0].copy()
            for l in dupT2:
                for record in db.records:
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
    
    name = str(random.randint(1000,2000))
    fig.savefig(figs_folder + '/' + name +'.png')
    close('all')
    return figs_folder + '/' + name+'.png', legend 

def draw_trimal_scc(db, num_col, figs_folder, trimmed=False, alg = '-scc'):
    import pandas as pd
    import matplotlib.pyplot as plt
    import random, os
    from Bio import AlignIO
    
    # get the alignment objects
    #-------------------------#
    alignments = db.alignments.items()
    if trimmed:
        alignments = db.trimmed_alignments.items()
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
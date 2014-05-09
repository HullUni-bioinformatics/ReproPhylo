
from Bio import SeqIO
import os, csv, sys
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation

class Guide:
       
    char_type = 'NotSet'
    feature_type = 'NotSet'
    name = 'NotSet'
    aliases = []
    
    def __init__(self, char_type=char_type, feature_type=feature_type,
                 name=name, aliases=aliases):
        """ Use to tell ReproPhylo which partition to build from a mixed (or not) data.
        
        >>> guide = Guide('dna', 'CDS', 'coi', ['cox1','COX1','coi','COI','CoI'])
        >>> print guide
        Guide(char_type=dna, feature_type=CDS, name=coi, aliases=cox1; COX1; coi; COI; CoI)"""

        
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
        return ('Guide(char_type='+self.char_type+', feature_type='+self.feature_type+
                ', name='+self.name+', aliases='+aliases_str+')')
            
                

if __name__ == "__main__":
    import doctest
    doctest.testmod() 
        
class GuideList():
    
    def __init__(self):
        self.content = []
        
        """ A list which accepts only Guide() objects.
        
        >>> guide = Guide('dna', 'CDS', 'coi', ['cox1','COX1','coi','COI','CoI'])
        >>> guide_list = GuideList()
        >>> guide_list.add_guide_to_list(guide)
        >>> print guide_list.content[0].aliases[2]
        coi"""
        

        
    def add_guide_to_list(self, guide):
        
        """ adds a guide to a GuideList only if it is a Guide
        
        >>> guide_list = GuideList()
        >>> guide_list.add_guide_to_list('not a guide')
        Traceback (most recent call last):
        ...
        TypeError: Expecting Guide() object"""
    
        names = []
        for g in self.content:
            names.append(g.name)
        if isinstance(guide, Guide):
            if not guide.name in names:
                self.content.append(guide)
            else: 
                raise NameError('Guide ' + guide.name + ' is already included')
        else:
            raise TypeError("Expecting Guide() object") 
    
    def __str__(self):
        guide_list_str = 'GuideList('
        for guide in self.content:
            guide_list_str += guide.name+', '
        return guide_list_str[:-2]+')' 
            
if __name__ == "__main__":
    import doctest
    doctest.testmod() 

# tools to prepare the db
    
def keep_feature(feature, guide_list):
    
    """ Returns true if a feature's type is in one of the guides and if the gene
    or product qualifiers is in the aliases of one of the guides
    
    >>> coi = Guide('dna','CDS','coi', ['cox1','COX1','coi','COI','CoI'])
    >>> guide_list = GuideList()
    >>> guide_list.add_guide_to_list(coi)
    >>> location = FeatureLocation(1,100)
    >>> feature = SeqFeature()
    >>> feature.location = location
    >>> feature.type = 'CDS'
    >>> feature.qualifiers['gene'] = ['CoI']
    >>> a = keep_feature(feature, guide_list)
    >>> print a
    True"""
    
    keep = 0
    for g in guide_list.content:
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
    
def dwindle_record(record, guide_list):
    dwindled_features = []
    feature_count = 0
    for feature in record.features:
        if keep_feature(feature, guide_list)== True:
            if feature.type == 'source' and not 'feature_id' in feature.qualifiers.keys():
                feature.qualifiers['feature_id'] = [record.id + '_source']
            elif not 'feature_id' in feature.qualifiers.keys():
                feature.qualifiers['feature_id'] = [record.id + '_f' + str(feature_count)]
                feature_count += 1                                        
            dwindled_features.append(feature)
    record.features = dwindled_features
    return record
            
            
def guess_format(input_filename):
    
    """ Runs a perl script that returns the input's format
    
    >>> filename = 'data/22_rotifer_COI_seqs.fas'
    >>> print guess_format(filename)
    fasta"""
    
    if os.path.exists(input_filename):
        os.system('perl guesser.pl ' + input_filename)
        guess = open(input_filename + '.format','r').read();
        os.remove(input_filename + '.format')
        return guess
    else:
        sys.exit('Cannot guess ' + input_filename + '. Does it exist?')
        
def is_embl_or_gb(input_filename):
    input_format = guess_format(input_filename)
    if input_format == 'embl' or input_format == 'genbank':
        return True
    else:
        return False

def parse_input(input_filename):
    return SeqIO.parse(input_filename, guess_format(input_filename))


def list_to_string(List):
    string = ''
    for i in List:
        if type(i) is str and '\n' in i:
            string += lines_to_line(i).rstrip()+';'
        else:
            string += str(i)+';'
    return string[:-1]

def lines_to_line(lines):
    lines = lines.split('\n')
    return (' ').join(lines)

def type_to_single_line_str(var):
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



class Database:
    
    
    def __init__(self, guide_list):
        self.records = []
        self.guide_list = guide_list
        self.records_by_guide = {}
        self.sets = []
        self.alignments = {}
        self.trees = {}
    
    def read_embl_genbank(self, input_filenames_list):
        generators = []
        for input_filename in input_filenames_list:
            if is_embl_or_gb(input_filename):
                generators.append(parse_input(input_filename))
                print 'Accepted:', input_filename
            else:
                print ('Rejected: ' + input_filename + ' is ' + 
                       guess_format(input_filename) + ' . Should be genbank or embl.')
            for generator in generators:
                for record in generator:
                    dwindled_record = dwindle_record(record, guide_list)
                    if len(record.features) > 1:
                        self.records.append(dwindled_record)
                    elif len(record.features) == 1 and not record.features[0].type == 'source':
                        self.records.append(dwindled_record)

    def read_denovo(self, input_filename, feature_type, char_type, source_qualifiers = {}):
        count = 0
        denovo = SeqIO.parse(input_filename, guess_format(input_filename))
        for record in denovo:
            feature = SeqFeature(FeatureLocation(0, len(record.seq)), type=feature_type, strand=1)
            source = SeqFeature(FeatureLocation(0, len(record.seq)), type='source', strand=1)
            if len(source_qualifiers.keys())>0:
                for key in source_qualifiers.keys():
                    source.qualifiers[key] = source_qualifiers[key]
            feature.qualifiers['original_id'] = [record.id]
            feature.qualifiers['original_desc'] = [(' ').join(record.description.split()[1:])]
            record.id = 'denovo'+str(count)
            count += 1
            feature.qualifiers['feature_id'] = [record.id+'_f0']
            source.qualifiers['feature_id'] = [record.id+'_source']
            record.features = [source, feature]
            if char_type == 'prot':
                record.seq.alphabet = IUPAC.protein
            elif char_type == 'dna':
                record.seq.alphabet = IUPAC.ambiguous_dna
            self.records.append(record)

    def write(self, filename, format = 'genbank'):
        if format == 'genbank' or format == 'embl':
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
                            line = line_start
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
                
    def species_vs_guides(self, outfile_name):
        species_vs_guides = {}
        for record in self.records:
            organism = 'undef'
            for feature in record.features:
                if feature.type == 'source':
                    if 'organism' in feature.qualifiers.keys():
                        organism = feature.qualifiers['organism'][0]
            if not organism in species_vs_guides.keys():
                species_vs_guides[organism] = {}    
            for feature in record.features:
                if not feature.type == 'source':
                    for guide in self.guide_list.content:
                        if not guide.name in guide.aliases:
                            guide.aliases.append(guide.name)
                        if 'gene' in feature.qualifiers.keys():
                            if feature.qualifiers['gene'][0] in guide.aliases:
                                if not guide.name in species_vs_guides[organism].keys():
                                    species_vs_guides[organism][guide.name] = 1
                                else:
                                    species_vs_guides[organism][guide.name] += 1
                        elif 'product' in feature.qualifiers.keys():
                            if feature.qualifiers['product'][0] in guide.aliases:
                                if not guide.name in species_vs_guides[organism].keys():
                                    species_vs_guides[organism][guide.name] = 1
                                else:
                                    species_vs_guides[organism][guide.name] += 1
        with open(outfile_name, 'wb') as csvfile:
            linewriter = csv.writer(csvfile, delimiter='\t',
                                    quotechar='|',
                                    quoting=csv.QUOTE_MINIMAL)
            guide_names = []
            for g in self.guide_list.content:
                guide_names.append(g.name)
            linewriter.writerow(['species']+guide_names)
            for organism in species_vs_guides.keys():
                line = [organism]
                for name in guide_names:
                    if name in species_vs_guides[organism].keys():
                        line.append(str(species_vs_guides[organism][name]))
                    else:
                        line.append('0')
                linewriter.writerow(line)
                    
            
            
    def extract_by_guide(self):
        data_by_guide = {}
        for guide in self.guide_list.content:
            if not guide.name in guide.aliases:
                guide.aliases.append(guide.name)
            records = []
            for record in self.records:
                for feature in record.features:
                    if (feature.type == guide.feature_type and 
                        (feature.qualifiers['gene'][0] in guide.aliases or
                         feature.qaulifiers['product'][0] in guide.aliases)):
                        if guide.char_type == 'dna':
                            S = feature.extract(record.seq)
                        elif guide.char_type == 'prot':
                            S = Seq(feature.qualifiers['translation'][0], IUPAC.protein)
                        feature_record = SeqRecord(seq = S, id = feature.qualifiers['feature_id'][0],
                                                   description = '')
                        records.append(feature_record)
            data_by_guide[guide.name] = records
        self.records_by_guide = data_by_guide


    def write_by_guide(self, format = 'fasta'):
        if self.records_by_guide == {}:
            self.extract_by_guide
        for key in self.records_by_guide.keys():
            SeqIO.write(self.records_by_guide[key], key+'.'+format, format)
            
if __name__ == "__main__":
    import doctest
    doctest.testmod()   
          
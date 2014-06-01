
from Bio import SeqIO
import os, csv, sys, dendropy
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio.Align.Applications import MafftCommandline
from StringIO import StringIO 
from Bio import AlignIO 
from Bio.Phylo.Applications import RaxmlCommandline
from cogent import LoadSeqs, LoadTree
from cogent.app.raxml_v730 import build_tree_from_alignment
from ete2 import *

class Locus:
       
    char_type = 'NotSet'
    feature_type = 'NotSet'
    name = 'NotSet'
    aliases = []

    def __init__(self, char_type=char_type, feature_type=feature_type,
                 name=name, aliases=aliases):
        """ Use to tell ReproPhylo which partition to build from a mixed (or not) data.
        
        >>> locus = Locus('dna', 'CDS', 'coi', ['cox1','COX1','coi','COI','CoI'])
        >>> print locus
        Locus(char_type=dna, feature_type=CDS, name=coi, aliases=cox1; COX1; coi; COI; CoI)"""

        
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
            
                

if __name__ == "__main__":
    import doctest
    doctest.testmod() 
       

class Set:
    
    name = 'NotSet'
    loci = []
    concat_meta = 'NotSet'
    concat_must_have_all_of = []
    concat_must_have_one_of = []
    
    def __init__(self,
                 name = name,
                 loci = loci,
                 concat_meta = concat_meta,
                 concat_must_have_all_of = concat_must_have_all_of,
                 concat_must_have_one_of = concat_must_have_one_of):
        self.name = name
        self.loci = loci
        self.concat_meta = concat_meta
        self.concat_must_have_all_of = concat_must_have_all_of
        self.concat_must_have_one_of = concat_must_have_one_of
        seen = []
        for locus in loci:
            if not isinstance(locus, Locus):
                raise TypeError("Expecting Locus object in loci list")
            if locus.name in seen:
                raise NameError('Locus ' + locus.name + ' apears more than once in self.loci')
            else:
                seen.append(locus.name)
        
# tools to prepare the db
 
def keep_feature(feature, loci):
    
    """ Returns true if a feature's type is in one of the loci and if the gene
    or product qualifiers is in the aliases of one of the loci
    
    >>> coi = Locus('dna','CDS','coi', ['cox1','COX1','coi','COI','CoI'])
    >>> location = FeatureLocation(1,100)
    >>> feature = SeqFeature()
    >>> feature.location = location
    >>> feature.type = 'CDS'
    >>> feature.qualifiers['gene'] = ['CoI']
    >>> a = keep_feature(feature, [coi])
    >>> print a
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
    
    """ retains only features that are called by guides and records with features that are
        called by guides
        
    >>> coi = Locus('dna','CDS','coi', ['cox1','COX1','coi','COI','CoI'])
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
    >>> kept_record = SeqRecord(seq=Seq(s, IUPAC.ambiguous_dna), id='1', description='spam')
    >>> kept_record.features.append(kept_feature)
    >>> kept_record.features.append(dwindled_feature)
    >>> print len(kept_record.features)
    2
    >>> a = dwindle_record(kept_record, [coi])
    >>> print len(kept_record.features)
    1"""
    
    dwindled_features = []
    feature_count = 0
    for feature in record.features:
        if keep_feature(feature, loci)== True:
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
    
    """ Using the guess_format() function to check if the input conforms with
        embl or genbank formats.
    """ 
    
    input_format = guess_format(input_filename)
    if input_format == 'embl' or input_format == 'genbank':
        return True
    else:
        return False

def parse_input(input_filename):
    return SeqIO.parse(input_filename, guess_format(input_filename))


def list_to_string(List):
    
    """ handles list printing as a nice string
    
    >>> L = ['a','b','b']
    >>> print list_to_string(L)
    a;b;b"""
    
    string = ''
    for i in List:
        if type(i) is str and '\n' in i:
            string += lines_to_line(i).rstrip()+';'
        else:
            string += str(i)+';'
    return string[:-1]

def lines_to_line(lines):
    
    """ Replaces newline with space"""
    
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


def get_qualifiers_dictionary(database, feature_id):
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


class Database:
    

    def __init__(self, loci):
        self.records = []
        self.loci = loci
        self.records_by_locus = {}
        self.sets = []
        self.alignments = {}
        self.trees = {}
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
            if is_embl_or_gb(input_filename):
                generators.append(parse_input(input_filename))
                print 'Accepted:', input_filename
            else:
                print ('Rejected: ' + input_filename + ' is ' + 
                       guess_format(input_filename) + ' . Should be genbank or embl.')
            for generator in generators:
                for record in generator:
                    dwindled_record = dwindle_record(record, self.loci)
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
            
    def add_set(self, set_object):
        if isinstance(set_object, Set):
            seen = []
            for s in self.sets:
                seen.append(s.name)
            if set_object.name in seen:
                raise NameError('Set ' + set_object.name + ' apears more than once in self.sets')
            else:
                self.sets.append(set_object)
        else:
            raise TypeError("Expecting Set object")

#    def make_set_alignments(self):
#        for s in self.sets:
#            meta = s.concat_meta
#            meta_list = []
            



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
            
    def align(self):
        self.write_by_locus(format = 'fasta')
        for key in self.records_by_locus.keys():
            infile_name = key+'.fasta'
            mafft_cline = MafftCommandline(input=infile_name)
            stdout, stderr = mafft_cline() 
            align = AlignIO.read(StringIO(stdout), "fasta")
            self.alignments[key] = align
            
    def write_alns(self, format = 'fasta'):
        if len(self.alignments.keys()) == 0:
            raise IOError('Align the records first')
        else:
            for key in self.alignments:
                AlignIO.write(self.alignments[key], key+'_aln.'+format, format)
            
    def tree(self):
        if len(self.alignments.keys()) == 0:
            raise IOError('Align the records first')
        else:
            self.write_alns(format = 'fasta')
            for key in self.alignments:
                char_type = ''
                mode = ''
                for locus in self.loci:
                    if locus.name == key:
                        char_type = locus.char_type
                if char_type == 'dna':
                    mod = 'GTRGAMMA'
                elif char_type == 'prot':
                    mod = 'PROTGAMMAJTT'
                alignment = self.alignments[key]
                raxml_cline = RaxmlCommandline(sequences=key+'_aln.fasta', algorithm='a',
                                               num_replicates=4, parsimony_seed=42,
                                               threads=3,rapid_bootstrap_seed=452,
                                               name=key, model=mod)
                raxml_cline()
                t = Tree('RAxML_bipartitions.'+key)
                t.name = key
                t.dist = 0
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
                    
                                        
                            
                    
                self.trees[key] = t
                for file_name in os.listdir(os.curdir):
                    if 'RAxML_' in file_name:
                        os.remove(file_name)


    def write_nexml(self, output_name):
        D = dendropy.DataSet()
        tree_list = []
        
        for tree_name in self.trees.keys():
            tree_string = self.trees[tree_name].write(features=[])
            tree = dendropy.Tree()
            tree.read_from_string(tree_string, schema='newick', extract_comment_metadata = True)
            tree_list.append(tree)
        TL = dendropy.TreeList(tree_list)    
        D.add_tree_list(TL)
        
        for aln_name in self.alignments.keys():
            char_type = ''
            matrix = None
            for locus in self.loci:
                if locus.name == aln_name:
                    char_type = locus.char_type
            if char_type == 'dna':
                matrix = dendropy.DnaCharacterMatrix()
            elif char_type == 'prot':
                matrix = dendropy.ProteinCharacterMatrix()
            matrix_string = self.alignments[aln_name].format('fasta')
            matrix.read_from_string(matrix_string,'fasta')
            D.add_char_matrix(matrix)
            
        D.write_to_path(
            output_name,
            'nexml',
            suppress_annotations=False,
            annotations_as_nhx=False,
            exclude_trees=False)



    def annotate(self,
    
                 leaf_labels_txt_meta,
                 leaf_node_color_meta,
                 leaf_label_colors,
    
                 node_bg_meta,
                 node_bg_color,
    
                 root_meta,
                 root_value,
    
                 node_support_dict):
    
            ts = TreeStyle()
            ts.show_leaf_name = False
            for tree in self.trees.keys():
    
                # set outgroup leaves, labels and label colors
                outgroup_list = []
                for leaf in self.trees[tree]:
                    qualifiers_dictionary = get_qualifiers_dictionary(self, leaf.name)
                    leaf_label = ''
                    for meta in leaf_labels_txt_meta:
                        leaf_label += qualifiers_dictionary[meta]+' '
                    leaf_label = leaf_label[:-1]
                    fgcolor = 'black'
                    for colour_name in leaf_label_colors.keys():
                        if colour_name in qualifiers_dictionary[leaf_node_color_meta]:
                            fgcolor = leaf_label_colors[colour_name]
                    leaf_face = TextFace(leaf_label, fgcolor=fgcolor)
                    leaf.add_face(leaf_face,0)
                    if root_value in qualifiers_dictionary[root_meta]:
                        outgroup_list.append(leaf)
                #set outgroup
                if len(outgroup_list) == 1:
                    self.trees[tree].set_outgroup(outgroup_list[0])
                elif len(outgroup_list) > 1:
                    R = self.trees[tree].get_common_ancestor(outgroup_list)
                    self.trees[tree].set_outgroup(R)
                elif len(outgroup_list)==0:
                    R = self.trees[tree].get_midpoint_outgroup()
                    self.trees[tree].set_outgroup(R)
    
                # ladderize
                self.trees[tree].ladderize()
    
                # node bg colors
                for key in node_bg_color.keys():
                    for node in self.trees[tree].get_monophyletic(values=[key], target_attr=node_bg_meta):
                        ns = NodeStyle(bgcolor=node_bg_color[key])
                        node.set_style(ns)
    
                # node support
    
                for node in self.trees[tree].traverse():
                    for key in node_support_dict.keys():
                        if (node.support <= node_support_dict[key][0] and
                            node.support > node_support_dict[key][1]):
                            node.add_face(CircleFace(radius = 5, color = key),column=0, position = "float")
                ts.legend_position=1
                ts.legend.add_face(TextFace('Node support: ', fsize=20), column=0)
                i = 1
                for color in node_support_dict.keys():
                    ts.legend.add_face(CircleFace(radius = 8, color = color), column=i)
                    i +=1 
                    ts.legend.add_face(TextFace(' '+str(node_support_dict[color][0])+'-'+str(node_support_dict[color][1]),
                                                fsize=20), column=i)
                    i += 1
                    
                self.trees[tree].render(tree+'.png',w=1000, tree_style=ts)



if __name__ == "__main__":
    import doctest
    doctest.testmod()   
          
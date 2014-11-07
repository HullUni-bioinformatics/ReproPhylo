
#######################################################################
#                        COMMAND LINE PARSING
#######################################################################
import argparse as agp

desc = ("Base reprophylo analysis. Deals with data and "+
        "metadata and runs default analyses. If the Project was "+
        "reconfigured with other tools, "+
        "it will run the reconfigured "+
        "analysis")
parser = agp.ArgumentParser(description=desc)

#----------------------------------------------------------------------
parser.add_argument("projfile",
                    help="A pickle file to read and write",
                    type=str,
                    default=None)
parser.add_argument("-a", "--altprojfile",
                    help="An alternative pickle file to write",
                    type=str, 
                    default=None)
parser.add_argument("-g", "--genbank",
                    help="An input genbank file name",
                    type=str,
                    default=None)
parser.add_argument("-m", "--dump_alignments",
                    help="print sequene alignment in the specified format",
                    type=str,
                    default=None)
parser.add_argument("--dump_alignments_zip",
                    type=str,
                    default="alignments.zip")
parser.add_argument("-p", "--dump_records",
                    help="print a genbank file of all records and metadata",
                    type=str,
                    default=None)
parser.add_argument("-T", "--threads",
                    help=("Number of threads to run in parrallel "+
                    "components of the pipeline"),
                    type=int,
                    default=4)
parser.add_argument("-j", "--figs",
                    help="A directory to keep temp and final fig files",
                    type=str, 
                    default='.')
parser.add_argument("-_", "--__",
                    help=("Turn off warnings for Galaxy"),
                    action="store_true")
parser.add_argument("--noninteractive",
                    help=("Overwrite reports without asking"),
                    action="store_true")
#----------------------------------------------------------------------
parser.add_argument("-f", "--fasta",
                    help="An input fasta file name",
                    type=str, nargs='+',
                    default=None)
parser.add_argument("-t", "--fasta_chartype",
                    help=("Specify the type of sequence in all the"+
                    " fasta files [dna|prot]. Required with --fasta"),
                    type=str,
                    choices=['dna', 'prot'])# Only with --fasta
parser.add_argument("-n", "--feature",
                    help=("Set the gene name and feature type"+
                          " for all the sequences in"+
                          " fasta file"),
                    type=str, nargs='+',
                    default=None)#Only with -f one_filename. 
parser.add_argument("-l", "--featureis_fasname",
                    help=("Set the file name as the gene name and "+
                          "feature type (eg cox1@CDS)"+
                          "for all the sequence in that"+
                          " fasta file"),
                    action="store_true")# only with -f, will apply to all fas
parser.add_argument("-q", "--qualis_fasname",
                    help=("Set a feature qualifier name contain"+
                          "ing the fasta file name for all "+
                          "the sequences in that file"),
                    type=str,
                    default=None)#Only with -f, will apply to all fas
parser.add_argument("-i", "--qualis_seqid",
                    help=("Set a feature qualifier name contain"+
                          "ing the sequence ID for all "+
                          "the sequences in the file"),
                    type=str,
                    default=None)#Only with -f, will apply to all fas
parser.add_argument("-d", "--qualis_seqdesc",
                    help=("Set a feature qualifier name contain"+
                          "ing the sequence description for all "+
                          "the sequences in the file"),
                    type=str,
                    default=None)#Only with -f, will apply to all fas
parser.add_argument("-e", "--expng_seqandmeta",
                    help=("Replace existing sequence records with"+
                          "new ones from gb and/or fasta file(s)"),
                    action="store_true")# Only with -g and or -f
#----------------------------------------------------------------------
loci = parser.add_mutually_exclusive_group()
loci.add_argument("-w", "--write_loci_csv",
                   help=("Write a loci CSV file. It can be "+
                         "edited and read back"),
                   type=str,
                   default=None)
loci.add_argument("-r", "--read_loci_csv",
                   help=("Read loci from an edited CSV file. "),
                   type=str,
                   default=None)
loci.add_argument("-u", "--add_locus",
                    help=("Add a locus to the analysis. "+
                          "Takes three space delimited values: "+
                          "[dna|prot] [CDS|rRNA|gene|...] locus_name"),
                    type=str, nargs=3,
                    default=None)
parser.add_argument("-b", "--banish_loci",
                    help=("Replace existing loci with"+
                          "new ones from a loci CSV file"),
                    action="store_true")# Only with -r
#----------------------------------------------------------------------
meta = parser.add_mutually_exclusive_group()
meta.add_argument("-y", "--write_meta_csv",
                   help=("Write a metadata CSV file. It can be "+
                         "edited and read back"),
                   type=str,
                   default=None)
meta.add_argument("-z", "--read_meta_csv",
                   help=("Read metadata from an edited CSV file"),
                   type=str,
                   default=None)
#----------------------------------------------------------------------
parser.add_argument("-s", "--start",
                   help=("Start the analysis"),
                   action="store_true")
parser.add_argument("-x", "--overwrite",
                    help=("Replace existing results with the ones"+
                          " of current configuration"),
                    action="store_true")# Only with -s
parser.add_argument("-v", "--verbose_report",
                   help=("Write full html report to filename"),
                   type=str,
                   default=None)
parser.add_argument("-c", "--version_ctrl",
                    help=("perform version control in CWD"),
                    action="store_true")
#----------------------------------------------------------------------
cl = parser.parse_args()
#----------------------------------------------------------------------
# At least one of these is required
if len([x for x in (cl.dump_alignments,
                    cl.dump_records,
                    cl.genbank,
                    cl.fasta,
                    cl.add_locus,
                    cl.write_loci_csv,
                    cl.read_loci_csv,
                    cl.write_meta_csv,
                    cl.read_meta_csv,
                    cl.start,
                    cl.verbose_report) if ((x is not None) and  
                                  (x is not False))])== 0:
    parser.error('One of the following is required:\n'+
                 '--genbank\n--fasta\n--write_loci_csv\n'+
                 '--read_loci_csv\n--write_meta_csv\n'+
                 '--read_meta_csv\n--start\nverbose_report\n'+
                 '--dump_alignments\n--dump_alignments\n'+
                 '--dump_records')
#----------------------------------------------------------------------
# The following go only with --fasta:
if not cl.fasta and any([cl.feature,
                         cl.featureis_fasname,
                         cl.qualis_fasname,
                         cl.qualis_seqid,
                         cl.qualis_seqdesc,
                         cl.fasta_chartype]):
    parser.error('The following only work with --fasta:\n'+
                 '--name_gene\n--gennameis_fasname\n'+
                 '--qualis_fasname\n--fasta_chartype\n'+
                 '--qualis_seqid\n--qualis_seqdesc')
# Gene name qualifier can be passed through cline only for a single
# fasta file
if cl.fasta and len(cl.fasta) > 1 and cl.feature:
    parser.error('--name_gene only works with one --fasta file')
# When reading fasta(s) the sequence type needs to be passed as well:
if cl.fasta and not cl.fasta_chartype:
    parser.error('Must specify the type of sequence in the fasta '+
                 'file(s) with --fasta_chartype [dna|prot]')
if not (cl.feature == None or len(cl.feature)==2):
    parser.error('--feature takes two space delimited strings to '+
                 'represent the gene name and feature type of the '+
                 'fasta sequences')
#----------------------------------------------------------------------
# Overwrite only when analysis is run:
if cl.overwrite and not cl.start:
    parser.error('--overwrite requires --start')
#----------------------------------------------------------------------
# Banish Locus objects only when reading new ones
if cl.banish_loci and not (cl.read_loci_csv or
                           cl.genbank or cl.add_locus):
    parser.error('--banish_loci only when reading new ones with '+
                 '--read_loci_csv or --genbank or --add_locus')
# Expunge SeqRecord objects only when reading new ones
if cl.expng_seqandmeta and not cl.genbank and not cl.fasta:
    parser.error('--expng_seqandmeta only when reading new ones with '+
                 '--genbank or --fasta')
#----------------------------------------------------------------------
#if False:"""    
#######################################################################
#                       Functions
#######################################################################
from reprophylo import *

def random_filename():
    return str(random.randint(10000,99999))+str(time.time())

def add_locus_from_cline():
    return Locus(cl.add_locus[0],cl.add_locus[1],
                 cl.add_locus[2],[cl.add_locus[2]])

# Functions to start a project and deal with a genbenk file and loci
#----------------------------------------------------------------------
def make_project(loci_csv=None, genbank=None):
    # If the project exists, modify the loci if requested:
    if os.path.isfile(cl.projfile) and os.path.getsize(cl.projfile)>10:
        pj = unpickle_pj(cl.projfile)
        # remove the old loci if requested
        if cl.banish_loci:
            pj.loci = []
        if cl.expng_seqandmeta: # We trust the cline parser to break
                                # if no new records are 
                                # supplied instead
            pj.records = []
            pj.records_by_locus = {}
            pj.concatenations = []
            pj.alignments = {}
            pj.trimmed_alignments = {}
            pj.trees = {}
        # Add a locus form the cline
        if cl.add_locus:
            new_locus = add_locus_from_cline()
            if ((not new_locus.name in [x.name for x in pj.loci])and
                not any(new_locus.name in x.aliases for x in pj.loci)):   
                pj.loci += [add_locus_from_cline()]
            else:
                raise IOError("%s already in the Project loci"%
                              cl.add_locus[2])
        # or add loci from a file
        elif loci_csv:
            new_loci_list = loci_list_from_csv(loci_csv)
            for l in new_loci_list:
                name_and_aliases = [l.name]+l.aliases
                for old in pj.loci:
                    if any(n in [old.name] + old.aliases 
                           for n in name_and_aliases):
                        raise IOError("%s already in the Project loci"%
                                      l.name)
            pj.loci += new_loci_list
            if genbank:
                pj.read_embl_genbank([genbank])
        # if genbank is read and no loci from file or cline, take loci
        # from genbank.
        elif genbank and not loci_csv and not cl.add_locus and pj.loci == []:
            loci_csv_out = random_filename()
            loci_readable_report = random_filename()
            # write the file if requested by cline
            if cl.write_loci_csv:
                loci_csv_out = cl.write_loci_csv
            list_loci_in_genbank(genbank, loci_csv_out,
                                 loci_report = loci_readable_report)
            print open(loci_readable_report,'r').readlines()[1]
            print ('No loci in the Project, reading'+
                   ' from genbank file')
            print '-------------------------------'
            print 'Gene and count sorted by counts'
            print '-------------------------------'
            handle = open(loci_readable_report,'r')
            print handle.read().rpartition('----')[-1]
            os.remove(loci_readable_report)
            # make aproject using loci in the gb file
            new_loci_list = loci_list_from_csv(loci_csv_out)
            for l in new_loci_list:
                name_and_aliases = [l.name]+l.aliases
                for old in pj.loci:
                    if any(n in [old.name] + old.aliases for n in name_and_aliases):
                        raise IOError("%s already in the Project loci"%l.name)
            pj.loci += loci_list_from_csv(loci_csv_out)
            pj.read_embl_genbank([genbank])
            # delete loci file if not requested to keep by cline
            if not cl.write_loci_csv:
                os.remove(loci_csv_out)
        elif genbank and not loci_csv and not cl.add_locus and len(pj.loci)>0:
            pj.read_embl_genbank([genbank])
        elif cl.write_loci_csv:
            handle = open(cl.write_loci_csv,'wt')
            for l in pj.loci:
                handle.write("%s,%s,%s"%(l.char_type,l.feature_type,l.name))
                for a in l.aliases:
                    handle.write(",%s"%a)
                handle.write('\n')
            handle.close()
    # if the project does not exist:
    else:
        # if a locus is passed through the cline:
        if cl.add_locus:
            pj = Project([add_locus_from_cline()])
        # if loci passed through a CSV file
        elif cl.read_loci_csv:
            pj = Project(loci_csv)
            if genbank:
                pj.read_embl_genbank([genbank])
        elif genbank and not loci_csv and not cl.add_locus:
            loci_csv_out = random_filename()
            loci_readable_report = random_filename()
            # write the file if requested by cline
            if cl.write_loci_csv:
                loci_csv_out = cl.write_loci_csv
            list_loci_in_genbank(genbank, loci_csv_out,
                                 loci_report = loci_readable_report)
            print open(loci_readable_report,'r').readlines()[1]
            print '-------------------------------'
            print 'Gene and count sorted by counts'
            print '-------------------------------'
            handle = open(loci_readable_report,'r')
            print handle.read().rpartition('----')[-1]
            os.remove(loci_readable_report)
            # make aproject using loci in the gb file
            pj = Project(loci_csv_out)
            pj.read_embl_genbank([genbank])
            # delete loci file if not requested to keep by cline
            if not cl.write_loci_csv:
                os.remove(loci_csv_out)
                    
    return pj    

# Functions to read fasta files and deal with their metadata
#----------------------------------------------------------------------
def deal_with_fasta(pj):
    # Read the records and make a feature based on cline instructions
    if cl.fasta:
        fst_dnov_in_file = 0
        fst_dnov_in_next_file = 0
        cline_quals = [cl.qualis_fasname,
                       cl.qualis_seqid,
                       cl.qualis_seqdesc]
        for fasta in cl.fasta:
            fst_dnov_in_next_file = pj.read_denovo([fasta],
                                                   cl.fasta_chartype)
            for dnov_index in range(fst_dnov_in_file,
                                    fst_dnov_in_next_file):
                current_record = [r for r in pj.records if
                                  r.id == 'denovo%i'%dnov_index][0]
                new_quals = {}
                feature_type = None
                if cl.feature:
                    new_quals['gene'] = cl.feature[0]
                    feature_type = cl.feature[1]
                elif cl.featureis_fasname:
                    new_quals['gene'] = fasta.split('/')[-1].split('@')[0]
                    feature_type = fasta.split('@')[1]
                elif ('feature_type' in cline_quals
                      and 
                      'gene' in cline_quals):
                    if cl.qualis_fasname == 'feature_type':
                        feature_type = fasta
                    elif cl.qualis_seqid  == 'feature_type':
                        feature_type = current_record.id
                    elif cl.qualis_seqdesc ==  'feature_type':
                        feature_type = current_record.description
                    if cl.qualis_fasname == 'gene':
                        new_quals['gene'] = fasta
                    elif cl.qualis_seqid  == 'gene':
                        new_quals['gene'] = current_record.id
                    elif cl.qualis_seqdesc ==  'gene':
                        new_quals['gene'] = current_record.description
                feat_id = None
                if (feature_type and not new_quals == {}):
                    rec_id = 'denovo%i'%dnov_index
                    feat_id = pj.add_feature_to_record(rec_id,
                                                       feature_type,
                                                       qualifiers=new_quals)
                if (not 'feature_type' in cline_quals
                    or 
                    not 'gene' in cline_quals):
                    if feat_id:
                        for feat in current_record.features:
                            if feat.qualifiers['feature_id'] == feat_id:
                                requiered = ['gene','feature_type']
                                if not cl.qualis_fasname in requiered:
                                    feat.qualifiers[cl.qualis_fasname] == fasta
                                if not cl.qualis_seqid in requiered:
                                    feat.qualifiers[cl.qualis_seqid] == current_record.id
                                if not cl.qualis_seqdesc in requiered:
                                    feat.qualifiers[cl.qualis_seqdesc] == current_record.description
            fst_dnov_in_file = fst_dnov_in_next_file


def removetype(path):
    import shutil
    if os.path.isdir(path):
        shutil.rmtree(path)
    else:
        os.remove(path)

def extract_and_log_by_locus(pj):
    print
    print "---------------------------------------------"
    print "The project now contains the following loci:".title()
    print "---------------------------------------------"
    if len(pj.records)>0:
        pj.extract_by_locus()
        print ("Locus".ljust(30,' ')+"Records".ljust(10,' ')+
               "Sequence length (max min mean)")
        for l in pj.loci:
            seq_lengths = [len(r) for r in pj.records_by_locus[l.name]] 
            if len(seq_lengths) > 0:
                print(str(l.name).ljust(30,' ')+
                       ("%i"%len(pj.records_by_locus[l.name])).ljust(10,' ')+
                       ("%i"%max(seq_lengths)).ljust(9,' ')+
                       ('%i'%min(seq_lengths)).ljust(9,' ')+
                       '%.1f'%np.mean(seq_lengths))
            else:
                print (str(l.name).ljust(30,' ')+'No records')
        if cl.verbose_report:
            if cl.noninteractive:
                if os.path.exists(cl.verbose_report):
                    removetype(cl.verbose_report)
                if os.path.exists(cl.verbose_report.rpartition('.')[0]):
                    removetype(cl.verbose_report.rpartition('.')[0])
                if os.path.exists(cl.verbose_report+'.zip'):
                    removetype(cl.verbose_report+'.zip')
            publish(pj, cl.verbose_report, cl.figs)
    else:
        for l in pj.loci:
            print l.name.ljust(30, ' ')+'No records'
               
              

def run_default_analysis(pj):
    if cl.start:
        if cl.overwrite:
            pj.concatenations = []
            pj.alignments = {}
            pj.trimmed_alignments = {}
            pj.trees = {}
            
        # make default labels:
        pj.add_qualifier_from_source('organism')
        pj.copy_paste_within_feature('organism','tree_label')
        for record in pj.records:
            for feature in record.features:
                if not feature.type =='source':
                    if not 'tree_label' in feature.qualifiers.keys():
                        feature.qualifiers['tree_label'] = feature.qualifiers['feature_id']
        
        print '------------------'
        print 'Sequence alignment'.title()
        print '------------------'
        mafft = AlnConf(pj, method_name='MafftDefaults',
                        CDSAlign=False, cline_args=dict())
        pj.align([mafft])
        print '------------------'
        print 'Trimming alignment'.title()
        print '------------------'        
        gappyout = TrimalConf(pj)
        pj.trim([gappyout])
        print '-------------------'
        print 'Tree reconstructios'.title()
        print '-------------------'   
        
        raxml = RaxmlConf(pj, method_name='fD_fb',
                          program_name=pj.defaults['raxmlHPC'],
                          preset = 'fD_fb',
                          model='GAMMA', matrix='JTT',
                          threads=int(cl.threads))
        
        pj.tree([raxml])
        
        pj.clear_tree_annotations()
    
        supports = {'black':[100,99],
                    'dimgray':[99,75],
                    'silver':[75,50]}
        print '---------------'
        print 'Tree annotation'.title()
        print '---------------' 
        
        pj.annotate(cl.figs,
                    root_meta='mid',
                    root_value='mid',
                    leaf_labels_txt_meta=['tree_label'],
                    node_support_dict=supports,
                    scale = 1000)
        
        
def ref():
    print '------------------'
    print 'Program References'
    print '------------------'
    if cl.start:
        print '-------------------------------'
        print 'Exploratory Pipeline References'
        print '-------------------------------'
        print ('RAxML: A. Stamatakis: "RAxML Version 8: A tool for '+
               'Phylogenetic Analysis and Post-Analysis of Large'+
               ' Phylogenies". In Bioinformatics, 2014')
        print ('http://bioinformatics.oxfordjournals.org/content/'+
               'early/2014/01/21/bioinformatics.btu033.abstract?'+
               'keytype=ref&ijkey=VTEqgUJYCDcf0kP')
        print 
        print ('MAFFT: Katoh, Standley 2013 (Molecular Biology and Evolution '+
               '30:772-780) MAFFT multiple sequence alignment software'+
               ' version 7: improvements in performance and usability.')
        print 'mbe.oxfordjournals.org/content/30/4/772'
        print 
        print ('trimAl: Salvador Capella-Gutierrez; Jose M. Silla-Martinez; Toni'+
               ' Gabaldon. trimAl: a tool for automated alignment trimming'+
               ' in large-scale phylogenetic analyses. Bioinformatics 2009 25: '+
               '1972-1973.')
        print 'http://trimal.cgenomics.org/_media/trimal.2009.pdf'
        print 
        print ('ETE: Jaime Huerta-Cepas, Joaquin Dopazo and Toni Gabaldon. '+
               'ETE: a python Environment for Tree Exploration. BMC '+
               'Bioinformatics 2010, 11:24.')
        print 'http://www.biomedcentral.com/1471-2105/11/24'
        print 
        print ('NumPy: Stefan van der Walt, S. Chris Colbert and Gael Varoquaux. '+
               'The NumPy Array: A Structure for Efficient Numerical Computation, '+
               'Computing in Science & Engineering, 13, 22-30 (2011)')
        print 'DOI:10.1109/MCSE.2011.37'
        print 
        print ('Matplotlib: John D. Hunter. Matplotlib: A 2D Graphics Environment ,'+
               'Computing in Science & Engineering, 9, 90-95 (2007)')
        print 'DOI:10.1109/MCSE.2007.55'
        print
        print ('Pandas: Wes McKinney. Data Structures for Statistical Computing '+
               'in Python, Proceedings of the 9th Python in Science Conference, '+
               '51-56 (2010)')
        print 'http://conference.scipy.org/proceedings/scipy2010/mckinney.html'
        print 
        print 'HTML.py: http://www.decalage.info/python/html'
        print '----------------'
        print 'Other References'
        print '----------------'
    print ('Biopython: Cock PJ, Antao T, Chang JT, Chapman BA, Cox CJ, Dalke A, Friedberg'+
           ' I, Hamelryck T, Kauff F, Wilczynski B, and de Hoon MJ. Biopython: freely '+
           ' available Python tools for computational molecular biology and'+
           ' bioinformatics. Bioinformatics 2009 Jun 1; 25(11) 1422-3.'+
           'doi:10.1093/bioinformatics/btp163 pmid:19304878')
    print 'http://www.hubmed.org/display.cgi?uids=19304878'
    print 
    print ('Cython: Stefan Behnel, Robert Bradshaw, Craig Citro, '+
           'Lisandro Dalcin, Dag Sverre Seljebotn and Kurt Smith. Cython: The '+
           'Best of Both Worlds, Computing in Science and Engineering, 13,'+
           ' 31-39 (2011)')
    print 'DOI:10.1109/MCSE.2010.118'
    print
    print 'Cloud: https://pypi.python.org/pypi/cloud/2.8.5'
    
    if cl.verbose_report:  
        print ('DendroPy: Sukumaran, J. and Mark T. Holder. 2010. DendroPy: A Python '+
               'library for phylogenetic computing. Bioinformatics 26: 1569-1571.')
        print 'http://bioinformatics.oxfordjournals.org/content/26/12/1569.abstract'  
        print

        


def chekpoint_project(pj):
    if cl.altprojfile:
        pickle_pj(pj, cl.altprojfile)
    else:
        pickle_pj(pj, cl.projfile)



#######################################################################
#                       Execution
#######################################################################
cline = ''

for i in sys.argv:
    cline += i+' '
print 'ReproPhylo was called with:'
print cline
print

if cl.__:
    warnings.simplefilter('ignore')

if cl.version_ctrl:
    start_git()
pj = make_project(loci_csv=cl.read_loci_csv, genbank=cl.genbank)    
deal_with_fasta(pj)



if cl.write_meta_csv:
    pj.write(cl.write_meta_csv, format='csv')
if cl.read_meta_csv:
    pj.correct_metadata_from_file(cl.read_meta_csv)

extract_and_log_by_locus(pj)

run_default_analysis(pj)

if cl.verbose_report:
    print '---------------------------------'
    print ('Writing report to %s.zip'%cl.verbose_report).title()
    print '---------------------------------' 
    if cl.noninteractive:
        if os.path.exists(cl.verbose_report):
            removetype(cl.verbose_report)
            if os.path.exists(cl.verbose_report.rpartition('.')[0]):
                removetype(cl.verbose_report.rpartition('.')[0])
        if os.path.exists(cl.verbose_report+'.zip'):
            removetype(cl.verbose_report+'.zip')
    publish(pj, cl.verbose_report, cl.figs)



if cl.dump_alignments:
    if len(pj.alignments.keys()+pj.trimmed_alignments.keys()) == 0:
        print "No alignments to print"
    else:
        import zipfile, shutil
        print '------------------------------------'
        print ('Writing alignments in %s format'%cl.dump_alignments).title()
        print '------------------------------------' 
        aln_files = pj.write_alns(format=cl.dump_alignments)
        trim_aln_files = pj.write_trimmed_alns(format=cl.dump_alignments)
        os.makedirs(cl.dump_alignments_zip.rpartition('.')[0])
        for filename in aln_files + trim_aln_files:
            shutil.move(filename, cl.dump_alignments_zip.rpartition('.')[0]+'/'+filename.split('/')[-1])
        zf = zipfile.ZipFile(cl.dump_alignments_zip, "w")
        for dirname, subdirs, files in os.walk(cl.dump_alignments_zip.rpartition('.')[0]):
            zf.write(dirname)
            for filename in files:
                zf.write(os.path.join(dirname, filename))
        zf.close()
        shutil.rmtree(cl.dump_alignments_zip.rpartition('.')[0])
    
if cl.dump_records:
    print '------------------------------------'
    print ('Writing records to %s'%cl.dump_records).title()
    print '------------------------------------' 
    pj.write(cl.dump_records, format='genbank')

chekpoint_project(pj)
ref()
print '----'
print 'DONE'
print '----'
#"""
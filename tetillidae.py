
from reprophylo import *

ssu = Locus('dna','rRNA','18S', ['18s','18S','SSU rRNA','18S ribosomal RNA'])
lsu = Locus('dna','rRNA','28S', ['28s','28S','LSU rRNA','28S ribosomal RNA'])
coi = Locus('dna','CDS','coi', ['cox1','COX1','coi','COI','CoI'])

loci = [ssu, lsu, coi]

db = Database(loci)

data_folder = '/home/amir/Dropbox/ReproPhylo/data/'

input_filenames_list = [data_folder + 'Tetillidae.gb']

db.read_embl_genbank(input_filenames_list)
db.add_qualifier_from_source('organism')

for genus in ['Cinachyrella','Cinachyra','Amphitethya','Paratetilla','Acanthotetilla',
              'Fangophilina', 'Tetilla','Craniella']:
    db.if_this_then_that(genus, 'organism', genus, 'genus', mode='part')

for genus in ['Geodia','Pachymatisma','Calthropella','Thenea','Theonella']:
    db.if_this_then_that(genus, 'organism', 'Astrophorid', 'genus', mode='part')    

db.if_this_then_that('18S ribosomal RNA', 'product', 'SSU rRNA', 'gene')
db.if_this_then_that('28S ribosomal RNA', 'product', 'LSU rRNA', 'gene')
db.if_this_then_that('AM076987.1_f1','feature_id','LAGLIDADG','gene')

  
db.write('Tet_Repro.csv', format = 'csv')
db.write('Tet_Repro.gb', format = 'genbank')
db.extract_by_locus()  #Not required if you also do:
db.write_by_locus()
db.species_vs_loci('Tetillidae_genes_per_species.csv')
db.align() 
db.write_alns()
db.tree()
db.write('Tetillidae.xml','nexml')

genera_colors = {'Astrophorid': 'black',
                 'Cinachyrella': 'black',
                 'Paratetilla': 'black',
                 'Amphitethya': 'black',
                 'Tetilla': 'black',
                 'Cinachyra': 'black',
                 'Acanthotetilla': 'black',
                 'Fangophilina': 'black',
                 'Craniella': 'black'
                 }

bg_colors = {'Astrophorid': 'white',
             'Cinachyrella': 'green',
             'Paratetilla': 'yellowgreen',
             'Amphitethya': 'olivedrab',
             'Tetilla': 'skyblue',
             'Cinachyra': 'firebrick',
             'Acanthotetilla': 'silver',
             'Fangophilina': 'crimson',
             'Craniella': 'red'
             }

supports = {'black':[100,99],
            'dimgray':[99,75],
            'silver':[75,50]
            }

db.annotate(leaf_labels_txt_meta=['annotation_organism','feature_id'],
            leaf_node_color_meta='annotation_organism', 
            leaf_label_colors=genera_colors,
         
            node_bg_meta='genus',
            node_bg_color=bg_colors,
         
            root_meta = 'genus',
            root_value = 'Astrophorid',
         
            node_support_dict = supports)
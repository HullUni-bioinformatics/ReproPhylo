import sys

input_gb = sys.argv[1]
output_loci_file = sys.argv[2]

def make_loci_control_file_from_genbank(genbank_filename, control_filename):
   from Bio import SeqIO
   
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
   
   sorted = gene_dict.items()
   sorted.sort(key = lambda i: i[0].lower())
   
   
   
   print('\n' + "There are " + str(len(sorted)) + " gene names (or gene product names) detected")
   print("---------------")
   print("Gene and count")
   print("---------------")
   
   control_file_handle = open(control_filename, 'wt')
   for key, value in sorted:
           #print key, value
           print(str(value) +" instances of " + key)
           feature_type = key.split(',')[0]
           alias = key.split(',')[1]
           name = alias.replace(' ','_')
           control_file_handle.write('dna,'+ feature_type + ',' + name + ',' + alias + '\n')
   control_file_handle.close()
                    
make_loci_control_file_from_genbank(input_gb, output_loci_file)
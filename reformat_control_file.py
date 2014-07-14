
import sys

edited_control_file = sys.argv[1]
reformatted_control_file = sys.argv[2]

def parse_manual_control_file(control_filename_input,control_filename_output):
    control_file_lines = open(control_filename_input, 'r').readlines()
    loci = {}
    seen_syn_group = []
    for line in [l.rstrip() for l in control_file_lines]:
        char_type, feature_type, name, alias = line.split(',',3)
        syn_group = None
        if ',' in alias:
            alias, syn_group = alias.split(',')

        # checks
        massage = """This function accepts text files in which
                     lines are of 4 or 5 comma separated values
                     , which are char_type (dna or prot), feature_type,
                     locus_name (no white spaces), alias (white spaces
                     allowed), and potentialy a number linking the locus
                     with other lines as a synonym"""
        if syn_group and ',' in syn_group:
            raise IOError(massage)
        if syn_group and not syn_group.isdigit():
            raise IOError(massage)
        if not char_type in ('dna', 'prot'):
            raise IOError(massage)
        
   
        if syn_group and syn_group in seen_syn_group:
            for key in loci.keys():
                if loci[key]['syn_group'] == syn_group:
                    loci[key]['aliases'].append(alias)

        else:
            loci[name] = {'char_type': char_type,
                          'feature_type': feature_type,
                          'aliases': [alias],
                          'syn_group': syn_group}
            seen_syn_group.append(syn_group)
            

    
    output = open(control_filename_output, 'wt')
    for name in loci.keys():
        line = loci[name]['char_type']+','+loci[name]['feature_type']+','+name
        for a in loci[name]['aliases']:
            line += ','+a
        line += '\n'
        output.write(line)
    output.close()

parse_manual_control_file(edited_control_file,reformatted_control_file)
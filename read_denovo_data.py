
from reprophylo import *
import pickle, sys

pickle_input = None
denovo_filename = None
char_type = None
feature_type = None
feature_name = None
picke_output = None

if len(sys.argv) == 5:
    pickle_input = sys.argv[1]
    denovo_filename = [sys.argv[2]]
    char_type = sys.argv[3]
    picke_output = sys.argv[4]


elif len(sys.argv) > 5:
    pickle_input = sys.argv[1]
    denovo_filename = [sys.argv[2]]
    char_type = sys.argv[3]
    feature_type = sys.argv[4]
    feature_name = sys.argv[5]
    id_qualifier = sys.argv[6]
    desc_qualifier = sys.argv[7]
    id_desc_qualifier = sys.argv[8]
    picke_output = sys.argv[9]

pickle_handle = open(pickle_input, 'rb')
db = pickle.load(pickle_handle)

previous_denovo = []

for record in db.records:
    if 'denovo' in record.id:
        previous_denovo.append(record.id)

db.read_denovo(denovo_filename, char_type)

if feature_type:
    for record in db.records:
        if 'denovo' in record.id and not record.id in previous_denovo:
            source = None
            for f in record.features:
                if f.type == 'source':
                    source = f
            db.add_feature_to_record(record.id,
                                     feature_type,
                                     qualifiers = {'gene': feature_name,
                                                   id_qualifier: source.qualifiers['original_id'][0],
                                                   desc_qualifier: source.qualifiers['original_desc'][0],
                                                   id_desc_qualifier: str(source.qualifiers['original_id'][0]+
                                                                          ' '+
                                                                          source.qualifiers['original_desc'][0])
                                                  }
                                     )

output = open(picke_output,'wb')
pickle.dump(db, output)
output.close()
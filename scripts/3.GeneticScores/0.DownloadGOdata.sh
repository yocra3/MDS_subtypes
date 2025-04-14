## Download GO ontology
wget http://current.geneontology.org/ontology/go-basic.obo -P data/GO_terms/
wget https://raw.githubusercontent.com/sgrote/OboToTerm/master/obo_to_term_tables.py  -P data/GO_terms/
wget https://raw.githubusercontent.com/sgrote/OboToTerm/master/obo_to_term_functions.py  -P data/GO_terms/

python3 data/GO_terms/obo_to_term_tables.py data/GO_terms/go-basic.obo data/GO_terms/

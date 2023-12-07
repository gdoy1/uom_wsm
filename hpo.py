from pyhpo import Ontology

# Initialize the Ontology
_ = Ontology()

# Retrieve the HPOTerm object for the MED7 gene
# Replace 'MED7' with the actual HPO ID or name if different
med7_term = Ontology.search('MED7')

# Check if the term exists and print associated HPO terms
if med7_term:
    print(f"Found {len(med7_term)} terms for MED7.")
    for term in med7_term:
        print(f"HPO Term: {term.name} | ID: {term.id}")

        # Access and print associated HPO terms
        for hpo_term in term.hpo_terms:
            print(f"- {hpo_term.name} | ID: {hpo_term.id}")
else:
    print("No HPO terms found for MED7.")
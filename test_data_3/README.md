I generate a species tree, a gene tree, and I reconcile them using ALERax. I then use them to test the functions in evals using evals/tests.ipynb

commands: 
- Zombi.py T SpeciesTreeParameters.tsv output_zombi
- Zombi.py G GenomeParameters.tsv output_zombi
- created the file families.txt containing paths for all gene families
- alerax -s output_zombi/T/ExtantTree.nwk -f families.txt -p output_alerax
from goatools.anno.gaf_reader import GafReader

ogaf = GafReader("/Users/linlane/Documents/Forsythe_Lab/Gene_Ontology/dev_data/gene_association.tair")

ns2assc = ogaf.get_ns2assc()

for namespace, associations in ns2assc.items():
    for protein_id, go_ids in sorted(associations.items())[:100]:
        print("{NS} {PROT:7} : {GOs}".format(
            NS=namespace,
            PROT=protein_id,
            GOs=' '.join(sorted(go_ids))))
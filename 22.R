x <- org.Mm.egGO
# Get the entrez gene identifiers that are mapped to a GO ID
mapped_genes <- mappedkeys(x)
# Convert to a list
xx <- as.list(x[mapped_genes]) ## xx is the whole nouse gene id in go, then find ours
if(length(xx) > 0) {
# Try the first one
got <- xx[[1]]
got[[1]][["GOID"]]
got[[1]][["Ontology"]]
got[[1]][["Evidence"]]
}
# For the reverse map:
# Convert to a list
xx <- as.list(org.Mm.egGO2EG)
if(length(xx) > 0){
	# Gets the entrez gene ids for the top 2nd and 3nd GO identifiers
	goids <- xx[2:3] ## check both entrez gene id and Go id to see what exactly id here !
	# Gets the entrez gene ids for the first element of goids
	goids[[1]]
	# Evidence code for the mappings
	names(goids[[1]])
}
# For org.Mm.egGO2ALLEGS
xx <- as.list(org.Mm.egGO2ALLEGS)
if(length(xx) > 0){
# Gets the Entrez Gene identifiers for the top 2nd and 3nd GO identifiers
goids <- xx[2:3]
# Gets all the Entrez Gene identifiers for the first element of goids
goids[[1]]
# Evidence code for the mappings
names(goids[[1]])
}
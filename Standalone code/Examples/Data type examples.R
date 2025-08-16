source("Standalone code/new shit.R")

testData = list(mRNA = list(),
                 met = list())

# load in human ensembl RNAseq example dataset
load("Standalone code/Examples/mRNA_counts_human_ensembl.RData")
testData$mRNA$humanEnsembl = rnaseq_human_ensembl
rm(rnaseq_human_ensembl)

# load in human symbols RNAseq example dataset
load("Standalone code/Examples/mRNA_counts_human_symbol.RData")
testData$mRNA$humanSymbol = rnaseq_human_symbol
rm(rnaseq_human_symbol)

# load in mouse ensembl RNAseq example dataset
load("Standalone code/Examples/mRNA_counts_mouse_ensembl.RData")
testData$mRNA$mouseEnsembl = rnaseq_mouse_ensembl
rm(rnaseq_mouse_ensembl)

# load in mouse symbol RNAseq example dataset
load("Standalone code/Examples/mRNA_counts_mouse_symbol.RData")
testData$mRNA$mouseSymbol = rnaseq_mouse_symbol
rm(rnaseq_mouse_symbol)

load("Standalone code/Examples/metabolomics_data.RData")

testData$met$original = metabolomics

# for each other column in metaboliteID, make a copy of metabolomics using those names where relevant
testData$met$CAS = metabolomics[rownames(metabolomics) %in% metaboliteID$Name,]
rownames(testData$met$CAS) = match(rownames(testData$met$CAS), metaboliteID$Name) %>% metaboliteID$CAS[.]
testData$met$CAS = testData$met$CAS[base::sample(1:nrow(testData$met$CAS)),]

testData$met$KEGG = metabolomics[rownames(metabolomics) %in% metaboliteID$Name,]
rownames(testData$met$KEGG) = match(rownames(testData$met$KEGG), metaboliteID$Name) %>% metaboliteID$KEGG[.]
testData$met$KEGG = testData$met$KEGG[base::sample(1:nrow(testData$met$KEGG)),]

testData$met$HMDB = metabolomics[rownames(metabolomics) %in% metaboliteID$Name,]
rownames(testData$met$HMDB) = match(rownames(testData$met$HMDB), metaboliteID$Name) %>% metaboliteID$HMDB[.]
testData$met$HMDB = testData$met$HMDB[base::sample(1:nrow(testData$met$HMDB)),]

testData$met$ChEBI = metabolomics[rownames(metabolomics) %in% metaboliteID$Name,]
rownames(testData$met$ChEBI) = match(rownames(testData$met$ChEBI), metaboliteID$Name) %>% metaboliteID$ChEBI[.]
testData$met$ChEBI = testData$met$ChEBI[base::sample(1:nrow(testData$met$ChEBI)),]

testData$met$twoIDtypes = metabolomics
rownames(testData$met$twoIDtypes)[base::sample(1:nrow(testData$met$twoIDtypes), 50)] = base::sample(metaboliteID$KEGG %>% .[!is.na(.)], 50)

lookup = read_csv("Standalone code/Examples/example_lookup_table.csv")

rm(metabolomics)

testResults = list(mRNA = list(humanEnsembl = squid(testData$mRNA$humanEnsembl,
                                                     species = "human",
                                                     geneType = "ENSEMBL"),
                                humanSymbol = squid(testData$mRNA$humanSymbol,
                                                    species = "human",
                                                    geneType = "SYMBOL"),
                                mouseEnsembl = squid(testData$mRNA$mouseEnsembl,
                                                     species = "mouse",
                                                     geneType = "ENSEMBL"),
                                mouseSymbol = squid(testData$mRNA$mouseSymbol,
                                                    species = "mouse",
                                                    geneType = "SYMBOL")),
                    met = list(original = squid(metabolites = testData$met$original,
                                                log.metabolites = F),
                               CAS = squid(metabolites = testData$met$CAS,
                                           log.metabolites = F),
                               KEGG = squid(metabolites = testData$met$KEGG,
                                            log.metabolites = F),
                               HMDB = squid(metabolites = testData$met$HMDB,
                                            log.metabolites = F),
                               ChEBI = squid(metabolites = testData$met$ChEBI,
                                             log.metabolites = F),
                               lookup = squid(metabolites = testData$met$original,
                                              metLookup = lookup,
                                              log.metabolites = F)))

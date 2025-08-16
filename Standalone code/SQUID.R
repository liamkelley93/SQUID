library(mixOmics)
library(tidyverse)
library(magrittr)
library(org.Hs.eg.db)
library(org.Mm.eg.db)

load("Standalone code/mixOmics_data_for_score_function.RData")

load("Standalone code/homologTable.RData")

metaboliteID = read_csv("Standalone code/metabolite_table.csv") %>% mutate(ChEBI = as.character(ChEBI))

## removal of duplicate / non-specific metabolites
# removed metabolites:
# Gly_Hexose-phosphate
# Mis_gamma-L-Glutamyl-L-alanine
# Pyr_2,3-Dideoxycytidine
# RX_5-oxoproline
# TCA_R 2-HG
mito$data$metabolites = mito$data$metabolites[c(1:71,73:87,89:127,129:144,146:159,161:166),]

rownames(mito$data$metabolites) = metaboliteID$Name

theme_squid = function(base_size = 12) {
  theme_classic() +
    theme(axis.text = element_text(color = "black", size = base_size),
          axis.ticks = element_line(color = "black"),
          plot.title = element_text(size = base_size),
          legend.text = element_text(size = base_size),
          strip.text = element_text(size = base_size),
          strip.background = element_rect(fill = "grey90", color = "black", linewidth = 0.5))
}

rearrangeSamples = function(seed) {
  set.seed(seed)
  conditions = mito$metadata$mRNA$condition %>% levels()
  conditionTable = function(cond) {
    samples.RNA = mito$metadata$mRNA %>% dplyr::filter(condition == cond) %>% pull(sample)
    samples.1hr = mito$metadata$metabolites %>% dplyr::filter(Timepoint == 1, Condition == cond) %>% pull(Sample)
    samples.6hr = mito$metadata$metabolites %>% dplyr::filter(Timepoint == 6, Condition == cond) %>% pull(Sample)
    tibble(Condition = rep(cond, 6),
           RNA = c(sample(samples.RNA, length(samples.RNA)), rep(NA, 6 - length(samples.RNA))),
           `1 hour` = c(sample(samples.1hr, 3), rep(NA, 3)),
           `6 hours` = sample(samples.6hr, 6),
           row.names = NULL)
  }
  lapply(conditions, conditionTable) %>% bind_rows() %>% mutate(Condition = factor(Condition, levels = levels(mito$metadata$metabolites$Condition)),
                                                                Sample = paste0("S", 1:nrow(.)))
}

rm(rearrangeMetSamples)

convertGenes = function(mRNA, species, geneType, log.mRNA = F, verbose = T) {
  convertToMouseSymbol = function(counts, type) {
    if(verbose) cli::cli_text("Converting genes to mouse symbols")
    
    bitr_output = clusterProfiler::bitr(rownames(counts),
                                        fromType = type,
                                        toType = "SYMBOL",
                                        OrgDb = "org.Mm.eg.db",
                                        drop = F) %>% 
      suppressMessages() %>% 
      suppressWarnings()
    
    # if there are no duplicate Ensembl IDs
    if(sum(duplicated(bitr_output$SYMBOL)) == 0) {
      if(verbose) cli::cli_alert_success("No duplicate gene symbols detected")
    }
    # if there are duplicate gene symbols, keep the one with the higher gene counts
    else {
      if(verbose) cli::cli_alert_warning("Duplicate gene symbols detected, removing duplicates")
      
      sums = rowSums(counts, na.rm = T)
      
      symbol_keep = merge(bitr_output,
                          sums %>% enframe(type, "sum"),
                          by = type) %>%
        arrange(desc(sum)) %>%
        dplyr::filter(!duplicated(SYMBOL))
      
      bitr_output = symbol_keep %>% dplyr::select(1:SYMBOL)
    }
    
    # add on missing mouse genes
    missingGenes = base::setdiff(rownames(counts), bitr_output[,1])
    mg = data.frame(a = missingGenes, b = NA)
    colnames(mg) = colnames(bitr_output)
    bitr_output = bitr_output %>% rbind(mg)
    
    # output mapping stats
    if(verbose) {
      nSymbols = bitr_output %>% dplyr::filter(!is.na(SYMBOL)) %>% nrow()
      cli::cli_alert_info("{nrow(counts)} gene{?s} mapped to {nSymbols} symbol{?s}")
    }
    
    # arrange in the same order as input
    bitr_output %>%
      mutate(index = match(!!sym(type), rownames(counts))) %>%
      arrange(index) %>%
      dplyr::select(-index)
  }
  convertMouseToHuman = function(counts) {
    if(verbose) cli::cli_text("Converting mouse symbols to human symbols")
    
    sums = rowSums(counts, na.rm = T)
    
    geneTable = homologTable %>%
      dplyr::filter(mouse %in% rownames(counts)) %>%
      mutate(counts = sums[match(mouse, names(sums))]) %>%
      arrange(desc(counts)) %>%
      filter(!duplicated(human)) %>%
      dplyr::select(-counts)
    
    # add on missing mouse genes
    missingGenes = base::setdiff(rownames(counts), homologTable$mouse)
    geneTable = geneTable %>% rbind(data.frame(mouse = missingGenes, human = NA))
    
    # output mapping stats
    if(verbose) {
      nHuman = geneTable %>% dplyr::filter(!is.na(human)) %>% nrow()
      cli::cli_alert_info("{nrow(counts)} mouse symbol{?s} mapped to {nHuman} human symbol{?s}")
    }
    
    # arrange in the same order as input
    geneTable %>%
      mutate(index = match(mouse, rownames(counts))) %>%
      arrange(index) %>%
      dplyr::select(-index)
  }
  convertToHumanEnsembl = function(counts, type) {
    if(verbose) cli::cli_text("Converting genes to human Ensembl IDs")
    
    bitr_output = clusterProfiler::bitr(rownames(counts),
                                        fromType = type,
                                        toType = "ENSEMBL",
                                        OrgDb = "org.Hs.eg.db",
                                        drop = F) %>% 
      suppressMessages() %>% 
      suppressWarnings()
    
    bitr_output = bitr_output %>%
      mutate(ENSEMBL = ifelse(ENSEMBL %in% rownames(mito$data$mRNA),
                              ENSEMBL,
                              NA))
    
    # if there are no duplicate Ensembl IDs
    if(sum(duplicated(bitr_output$ENSEMBL)) == 0) {
      if(verbose) cli::cli_alert_success("No duplicate Ensembl IDs detected")
    }
    # if there are duplicate Ensembl IDs, keep the one with the higher gene counts
    else {
      if(verbose) cli::cli_alert_warning("Duplicate Ensembl IDs detected, removing duplicates")
      
      sums = rowSums(counts, na.rm = T)
      
      ensembl_keep = merge(bitr_output,
                           sums %>% enframe(type, "sum"),
                           by = type) %>%
        arrange(desc(sum)) %>%
        dplyr::filter(!duplicated(ENSEMBL))
      
      bitr_output = ensembl_keep %>% dplyr::select(1:ENSEMBL)
    }
    
    # add on missing mouse genes
    missingGenes = base::setdiff(rownames(counts), bitr_output[,1])
    mg = data.frame(a = missingGenes, b = NA)
    colnames(mg) = colnames(bitr_output)
    bitr_output = bitr_output %>% rbind(mg)
    
    # output mapping stats
    if(verbose) {
      nEnsembl = bitr_output %>% dplyr::filter(!is.na(ENSEMBL)) %>% nrow()
      cli::cli_alert_info("{nrow(counts)} gene{?s} mapped to {nEnsembl} Ensembl ID{?s}")
    }
    
    # arrange in the same order as input
    bitr_output %>%
      mutate(index = match(!!sym(type), rownames(counts))) %>%
      arrange(index) %>%
      dplyr::select(-index)
  }
  
  if(log.mRNA)
    mRNA = mRNA %>% raise_to_power(2,.) %>% add(-1)
  
  if(species == "mouse") {
    if(geneType != "SYMBOL") {
      # convert from ensembl/entrez to symbol
      tb1 = convertToMouseSymbol(mRNA, geneType)
      from1 = tb1 %>% dplyr::filter(!is.na(SYMBOL)) %>% pull(1)
      to1 = tb1 %>% dplyr::filter(!is.na(SYMBOL)) %>% pull(2)
      mRNA2 = mRNA[from1,]
      rownames(mRNA2) = to1
      
      # convert from mouse symbol to human symbol
      tb2 = convertMouseToHuman(mRNA2)
      from2 = tb2 %>% dplyr::filter(!is.na(human)) %>% pull(1)
      to2 = tb2 %>% dplyr::filter(!is.na(human)) %>% pull(2)
      mRNA3 = mRNA2[from2,]
      rownames(mRNA3) = to2
      
      # convert from human symbol to human ensembl
      tb3 = convertToHumanEnsembl(mRNA3, "SYMBOL")
      
      # merge the tables together
      merged = merge(tb1, tb2, by.x = 2, by.y = 1, all = T, sort = F) %>%
        merge(tb3, by.x = 3, by.y = 1, all = T, sort = F) %>%
        .[,c(3,2,1,4)]
      colnames(merged) = c(paste0("MOUSE_", geneType),
                           "MOUSE_SYMBOL",
                           "HUMAN_SYMBOL",
                           "HUMAN_ENSEMBL")
      return(merged)
      
    }
    else {
      # convert from mouse symbol to human symbol
      tb1 = convertMouseToHuman(mRNA)
      from1 = tb1 %>% dplyr::filter(!is.na(human)) %>% pull(1)
      to1 = tb1 %>% dplyr::filter(!is.na(human)) %>% pull(2)
      mRNA2 = mRNA[from1,]
      rownames(mRNA2) = to1
      
      # convert from human symbol to human ensembl
      tb2 = convertToHumanEnsembl(mRNA2, "SYMBOL")
      
      # merge the tables together
      merged = merge(tb1, tb2, by.x = 2, by.y = 1, all = T, sort = F) %>%
        .[,c(2,1,3)]
      colnames(merged) = c("MOUSE_SYMBOL",
                           "HUMAN_SYMBOL",
                           "HUMAN_ENSEMBL")
      return(merged)
    }
  }
  else {
    tb = convertToHumanEnsembl(mRNA, geneType)
    return(tb)
  }
}

squid = function(mRNA = NULL,
                 metabolites = NULL,
                 species = c("human", "mouse"),
                 geneType = c("ENSEMBL", "SYMBOL", "ENTREZID"),
                 metLookup = NULL,
                 metadata = NULL,
                 ncomp = 4,
                 log.mRNA = F,
                 log.metabolites = F,
                 seed = 1,
                 verbose = T,
                 merge_by = 1) {
  # match args
  geneType = match.arg(geneType)
  species = match.arg(species)
  
  if(is.null(mRNA) & is.null(metabolites))
    cli::cli_abort("At least one dataset ({.pkg mRNA} or {.pkg metabolites}) is required")
  
  # sample prep
  if (!is.null(mRNA)) {
    if(verbose) cli::cli_text("Matching inputted mRNA data to mitochondrial data")
    
    orig_nrow = nrow(mRNA)
    
    # if not human Ensembl ID, convert to human Ensembl ID and change mRNA rownames
    if(species != "human" | geneType != "ENSEMBL") {
      geneConversionTable = convertGenes(mRNA, species, geneType, log.mRNA, verbose)
      fromGene = geneConversionTable %>% dplyr::filter(if_any(last_col(), ~ !is.na(.))) %>% pull(1)
      toGene = geneConversionTable %>% dplyr::filter(if_any(last_col(), ~ !is.na(.))) %>% pull(last_col())
      mRNA = mRNA[fromGene,]
      rownames(mRNA) = toGene
    }
    # use overlap for each dataset
    overlapGenes = base::intersect(rownames(mRNA), rownames(mito$data$mRNA))
    if(verbose) cli::cli_alert_info("{length(overlapGenes)} gene{?s} {?overlaps/overlap} between datasets out of {orig_nrow} (input matrix) and {nrow(mito$data$mRNA)} (mito)")
  }
  
  if (!is.null(metabolites)) {
    cli::cli_text("Matching inputted metabolomics data to mitochondrial data")
    
    orig_nrow = nrow(metabolites)
    
    if(is.null(metLookup)) {
      categorySums = apply(metaboliteID, 2, function(x) {sum(x %in% rownames(metabolites))})
      
      # if no matches or matches in multiple columns, throw error
      if(sum(categorySums) == 0)
        cli::cli_abort(c("No matches found between inputted metabolite IDs and mitochondrial dataset.", "!" = "Did you mean to include a lookup table?", "i" = "Additionally, ensure metabolite matrix has metabolite ID as rowname"))
      
      if(sum(categorySums > 0) > 1)
        cli::cli_abort(c("Multiple types of metabolite IDs found in inputted dataset",
                         "i" = "Metabolite IDs found: {names(categorySums)[which(categorySums > 0)]}",
                         "i" = "The metabolite ID table can be accessed using {.run squid::metaboliteIDs()}"))
      
      if((which(categorySums > 0) != which(names(metaboliteID) == "Name")) & verbose)
        cli::cli_alert_info("Metabolite names recognized as {names(metaboliteID)[which(categorySums > 0)]} IDs")
      
      # at this point exactly one column should have sum > 0
      metLookup = metaboliteID[,c(which(categorySums > 0), 1)]
    } else if(verbose) cli::cli_alert_info("Using provided metabolite lookup table")
    
    if(ncol(metLookup) != 2)
      cli::cli_abort(c("Metabolite lookup table should have two columns:", " " = "1. Name in input matrix", " " = "2. Name in mito data)"))
    
    colnames(metLookup) = c("a", "b")
    
    # filter out NAs
    metLookup = metLookup %>% dplyr::filter(!is.na(a), !is.na(b))
    
    if(metLookup$a %>% duplicated() %>% sum() %>% equals(0) %>% not() |
       metLookup$b %>% duplicated() %>% sum() %>% equals(0) %>% not()) {
      duplicatesA = metLookup$a[duplicated(metLookup$a)] %>% unique()
      duplicatesB = metLookup$b[duplicated(metLookup$b)] %>% unique()
      View(metLookup)
      cli::cli_abort(c("No duplicate entries allowed in metabolite lookup table",
                       if(length(duplicatesA > 0)) c("i" = "Duplicate{?s} in input names: {duplicatesA}") else NULL,
                       if(length(duplicatesB > 0)) c("i" = "Duplicate{?s} in mito names: {duplicatesB}") else NULL))
    }
    
    # filter so the only entries are proper matches between input data (a) and mito data (b)
    metLookup = metLookup %>% dplyr::filter(a %in% rownames(metabolites),
                                            b %in% rownames(mito$data$metabolites))
    
    # match rownames of mat to the lookup table
    indices = metabolites %>% rownames() %>% match(metLookup$a, .)
    
    # rownames(metabolites)[indices] equals metLookup$a, so subset the matrix and change the names to metLookup$b
    metabolites = metabolites[indices,]
    rownames(metabolites) = metLookup$b
    
    overlapMetabolites = rownames(metabolites)
    if(verbose) cli::cli_alert_info("{length(overlapMetabolites)} metabolite{?s} {?overlaps/overlap} between datasets out of {orig_nrow} (input matrix) and {nrow(mito$data$metabolites)} (mito)")
    
  }
  
  # time to run some diablo
  diablo = function(seed, condition = NULL, dropCondition = NULL, ncomp, ...) {
    conditions = levels(mito$metadata$mRNA$condition)
    
    if(!is.null(condition))
      conditions = c("DMSO",condition)
    
    if(!is.null(dropCondition))
      conditions = conditions %>% .[!(. %in% dropCondition)]
    
    samples = rearrangeSamples(seed) %>% dplyr::filter(Condition %in% conditions)
    
    rna = matrix(nrow = nrow(samples), ncol = nrow(mito$data$mRNA))
    rownames(rna) = samples %>% pull(Sample) # samples
    colnames(rna) = rownames(mito$data$mRNA) # genes
    rna[match(colnames(mito$data$mRNA) %>% .[. %in% (samples %>% pull(RNA))], samples %>% pull(RNA)),] = t(mito$data$mRNA[,colnames(mito$data$mRNA) %>% .[. %in% (samples %>% pull(RNA))]])
    
    # 1 hour metabolite table with NAs for non-matched samples
    met_1 = matrix(nrow = nrow(samples), ncol = nrow(mito$data$metabolites))
    rownames(met_1) = samples %>% pull(Sample) # samples
    colnames(met_1) = rownames(mito$data$metabolites) # metabolites
    met_1[match(colnames(mito$data$metabolites) %>% .[. %in% (samples %>% pull(`1 hour`))], samples %>% pull(`1 hour`)),] = t(mito$data$metabolites[,colnames(mito$data$metabolites) %>% .[. %in% (samples %>% pull(`1 hour`))]])
    
    # 6 hour metabolites with matched samples (no non-matched samples possible)
    met_6 = mito$data$metabolites %>% t() %>% .[samples %>% pull(`6 hours`),]
    rownames(met_6) = samples %>% pull(Sample)
    
    overlapGenes = if(!is.null(mRNA)) overlapGenes else rownames(mito$data$mRNA)
    
    overlapMetabolites = if(!is.null(metabolites)) overlapMetabolites else rownames(mito$data$metabolites)
    
    X = list(mRNA = rna[, overlapGenes],
             metabolites_1_hr = met_1[, overlapMetabolites],
             metabolites_6_hr = met_6[, overlapMetabolites])
    
    Y = samples %>% pull(Condition)
    
    list.keepX = list(mRNA = rep(200, ncomp),
                      metabolites_1_hr = rep(length(overlapMetabolites), ncomp),
                      metabolites_6_hr = rep(length(overlapMetabolites), ncomp))
    
    result = block.splsda(X, Y, keepX = list.keepX, ncomp = ncomp, near.zero.var = T, ...)
    
    return(result)
  }
  
  if(verbose) cli::cli_text("Running mixOmics DIABLO")
  suppressWarnings(result <- diablo(seed, dropCondition = "Tunicamycin", ncomp = ncomp))
  
  # plot component values for each condition
  comps = list(mRNA = result$variates$mRNA %>% as.data.frame() %>% rownames_to_column("name") %>% mutate(condition = result$Y) %>% filter(comp1 != 0),
               `1 hour` = result$variates$metabolites_1_hr %>% as.data.frame() %>% rownames_to_column("name") %>% mutate(condition = result$Y) %>% filter(comp1 != 0),
               `6 hours` = result$variates$metabolites_6_hr %>% as.data.frame() %>% rownames_to_column("name") %>% mutate(condition = result$Y) %>% filter(comp1 != 0))
  
  compPlot = bind_rows(comps, .id = "type") %>%
    pivot_longer(starts_with("comp"), names_to = "component") %>%
    mutate(type = factor(type, levels = c("mRNA", "1 hour", "6 hours"))) %>% 
    ggplot() +
    aes(x = condition, y = value, fill = condition, alpha = type) +
    stat_summary(fun.data = mean_se, geom = "errorbar", width = .5, position = position_dodge(width = 1, preserve = "single")) +
    stat_summary(fun = mean, geom = "bar", position = position_dodge(width = 1, preserve = "single")) + 
    geom_hline(yintercept = 0) +
    facet_wrap(~component, scales = "free_y", nrow = ifelse(ncomp <= 5, 1, 2)) +
    theme_squid() +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
    scale_fill_manual(values = conditionPalette(c(1:10,12))) +
    scale_alpha_manual(values = c(1, 0.75, 0.5)) +
    labs(x = NULL,
         y = "Component value")
  
  # if the data is not log-transformed, do that
  if(!is.null(mRNA) & !log.mRNA) {
    if(verbose) cli::cli_alert_info("Converting inputted mRNA dataset to log format")
    mRNA = mRNA %>% add(1) %>% log2()
  }
  
  if(!is.null(metabolites) & !log.metabolites) {
    if(verbose) cli::cli_alert_info("Converting inputted metabolite dataset to log format")
    metabolites = metabolites %>% add(1) %>% log2()
  }
  
  # functions for component and condition scores
  componentScore = function(comp, mat, loadings) {
    features = rownames(loadings) %>% .[loadings[,comp] != 0]
    loads = loadings[,comp] %>% .[. != 0]
    data = mat[features,] %>% t() %>% scale()
    sweep(data, 2, loads, "*") %>% rowSums(na.rm = T) %>% return()
  }
  
  getConditionLoadings = function(cond, mat, type) {
    suppressWarnings(condDIABLO <- diablo(1, condition = cond, ncomp = 2))
    loadings = condDIABLO$loadings[[type]] %>%
      as.data.frame() %>%
      dplyr::filter(comp1 != 0) %>%
      rownames_to_column("element") %>%
      pull("comp1", name = "element")
    samples = rearrangeSamples(seed)
    means = merge(condDIABLO$variates[[type]], samples, by.x = 0, by.y = "Sample", sort = F) %>%
      dplyr::select(Condition, comp1) %>% 
      group_by(Condition) %>%
      summarise_at(vars(comp1), .funs = list(mean)) %>%
      deframe()
    if(means["DMSO"] > means[cond]) loadings = loadings * -1
    return(loadings)
  }
  
  conditionScore = function(condition, mat, loadings) {
    features = rownames(loadings) %>% .[loadings[,condition] != 0]
    loads = loadings[,condition] %>% .[. != 0]
    data = mat[features,] %>% t() %>% scale()
    sweep(data, 2, loads, "*") %>% rowSums(na.rm = T) %>% return()
  }
  
  conditionList = structure(as.list(names(conditionPalette())[c(-1,-11)]),
                            names = names(conditionPalette())[c(-1,-11)])
  
  # calculate scores on mRNA (if included)
  if(!is.null(mRNA)) {
    # component scores on matrix
    if(verbose) cli::cli_text("Generating component scores on mRNA")
    
    componentLoadings.mRNA = result$loadings$mRNA %>% as.data.frame() %>% dplyr::filter(rowSums(.) != 0)
    
    scores.component.mRNA = sapply(1:ncomp,
                                   componentScore,
                                   mat = mRNA,
                                   loadings = componentLoadings.mRNA) %>% as.data.frame()
    colnames(scores.component.mRNA) = paste0("comp", 1:ncomp)
    
    # condition scores on matrix
    if(verbose) cli::cli_text("Generating condition scores on mRNA")
    
    conditionLoadings.mRNA = lapply(conditionList,
                                    getConditionLoadings,
                                    mat = mRNA,
                                    type = "mRNA")
    
    # convert to array
    temp = matrix(0,
                  nrow = conditionLoadings.mRNA %>% lapply(names) %>% Reduce(base::union, .) %>% length(),
                  ncol = length(conditionLoadings.mRNA))
    
    rownames(temp) = conditionLoadings.mRNA %>% lapply(names) %>% Reduce(base::union, .) %>% sort()
    colnames(temp) = names(conditionLoadings.mRNA)
    
    for(n in 1:length(conditionLoadings.mRNA)) {
      temp[names(conditionLoadings.mRNA[[n]]),n] = conditionLoadings.mRNA[[n]]
    }
    
    conditionLoadings.mRNA = temp %>% as.data.frame()
    rm(temp, n)
    
    scores.condition.mRNA = sapply(conditionList,
                                   conditionScore,
                                   mat = mRNA,
                                   loadings = conditionLoadings.mRNA) %>% as.data.frame()
    
    # merge scores
    scores.mRNA = cbind(scores.component.mRNA, scores.condition.mRNA)
  }
  
  # calculate scores on metabolites (if included)
  if(!is.null(metabolites)) {
    # component scores on matrix
    if(verbose) cli::cli_text("Generating component scores on metabolites (using 6 hour data)")
    
    componentLoadings.metabolites = result$loadings$metabolites_6_hr %>% as.data.frame() %>% dplyr::filter(rowSums(.) != 0)
    
    scores.component.metabolites = sapply(1:ncomp,
                                          componentScore,
                                          mat = metabolites,
                                          loadings = componentLoadings.metabolites) %>% as.data.frame()
    
    colnames(scores.component.metabolites) = paste0("comp", 1:ncomp)
    
    # condition scores on matrix
    if(verbose) cli::cli_text("Generating condition scores on metabolites (using 6 hour data)")
    
    conditionLoadings.metabolites = sapply(conditionList,
                                           getConditionLoadings,
                                           mat = metabolites,
                                           type = "metabolites_6_hr") %>% as.data.frame()
    rownames(conditionLoadings.metabolites) = rownames(metabolites)
    
    conditionLoadings.metabolites = conditionLoadings.metabolites %>% .[rowSums(equals(.,0)) != length(conditionList),]
    
    scores.condition.metabolites = sapply(conditionList,
                                          conditionScore,
                                          mat = metabolites,
                                          loadings = conditionLoadings.metabolites) %>% as.data.frame()
    
    # merge scores
    scores.metabolites = cbind(scores.component.metabolites, scores.condition.metabolites)
  }
  
  output = list(result = result,
                compPlot = compPlot)
  
  if(!is.null(mRNA))
    output$data.mRNA = mRNA
  
  if(!is.null(metabolites))
    output$data.metabolites = metabolites
  
  if(!is.null(mRNA))
    output$scores.mRNA = scores.mRNA
  
  if(!is.null(metabolites))
    output$scores.metabolites = scores.metabolites
  
  if(!is.null(metadata)) {
    # merge mRNA scores with metadata
    if(!is.null(mRNA)) {
      samples = if(merge_by > 0) pull(metadata, merge_by) else rownames(merge_by)
      if(sum(colnames(mRNA) %in% samples) > 0) {
        merged = merge(metadata, scores.mRNA, by.x = merge_by, by.y = 0)
        merged = merged[,base::union(colnames(metadata), colnames(scores.mRNA))]
        merged = match(colnames(mRNA), merged %>% pull(merge_by)) %>% merged[.,]
        rownames(merged) = 1:nrow(merged)
        output$mergedScores.mRNA = merged
      }
      else
        cli::cli_alert_warning("Unable to merge mRNA score output and metadata. The column names of {.pkg mRNA} should be present in the specified column ({.strong {merge_by}}) of the metadata table.")
    }
    
    # merge metabolite scores with metadata
    if(!is.null(metabolites)) {
      samples = if(merge_by > 0) pull(metadata, merge_by) else rownames(merge_by)
      if(sum(colnames(metabolites) %in% samples) > 0) {
        merged = merge(metadata, scores.metabolites, by.x = merge_by, by.y = 0)
        merged = merged[,base::union(colnames(metadata), colnames(scores.metabolites))]
        merged = match(colnames(metabolites), merged %>% pull(merge_by)) %>% merged[.,]
        rownames(merged) = 1:nrow(merged)
        output$mergedScores.metabolites = merged
      }
      else
        cli::cli_alert_warning("Unable to merge metabolite score output and metadata. The column names of {.pkg metabolites} should be present in the specified column ({.strong {merge_by}}) of the metadata table.")
    }
  }
  
  if(!is.null(mRNA))
    output$loadings.mRNA = list(comps = componentLoadings.mRNA,
                                conditions = conditionLoadings.mRNA)
  
  if(!is.null(metabolites))
    output$loadings.metabolites = list(comps = componentLoadings.metabolites,
                                       conditions = conditionLoadings.metabolites)
  
  if(exists("geneConversionTable"))
    output$geneTable = geneConversionTable
  
  if(verbose) cli::cli_alert_success("Done")
  
  class(output) = "squid"
  
  return(output)
}

compareSeeds = function(mRNA = NULL,
                        metabolites = NULL,
                        species = c("human", "mouse"),
                        geneType = c("ENSEMBL", "SYMBOL", "ENTREZID"),
                        metLookup = NULL,
                        seeds = 1:10,
                        ncomp = 4,
                        log.mRNA = F,
                        log.metabolites = F,
                        verbose = T) {
  # match args
  geneType = match.arg(geneType)
  species = match.arg(species)
  
  if(is.null(mRNA) & is.null(metabolites))
    cli::cli_abort("At least one dataset ({.pkg mRNA} or {.pkg metabolites}) is required")
  
  # sample prep
  if (!is.null(mRNA)) {
    if(verbose) cli::cli_text("Matching inputted mRNA data to mitochondrial data")
    
    orig_nrow = nrow(mRNA)
    
    # if not human Ensembl ID, convert to human Ensembl ID and change mRNA rownames
    if(species != "human" | geneType != "ENSEMBL") {
      geneConversionTable = convertGenes(mRNA, species, geneType, log.mRNA, verbose)
      fromGene = geneConversionTable %>% dplyr::filter(if_any(last_col(), ~ !is.na(.))) %>% pull(1)
      toGene = geneConversionTable %>% dplyr::filter(if_any(last_col(), ~ !is.na(.))) %>% pull(last_col())
      mRNA = mRNA[fromGene,]
      rownames(mRNA) = toGene
    }
    # use overlap for each dataset
    overlapGenes = base::intersect(rownames(mRNA), rownames(mito$data$mRNA))
    if(verbose) cli::cli_alert_info("{length(overlapGenes)} gene{?s} {?overlaps/overlap} between datasets out of {orig_nrow} (input matrix) and {nrow(mito$data$mRNA)} (mito)")
  }
  
  if (!is.null(metabolites)) {
    cli::cli_text("Matching inputted metabolomics data to mitochondrial data")
    
    orig_nrow = nrow(metabolites)
    
    if(is.null(metLookup)) {
      categorySums = apply(metaboliteID, 2, function(x) {sum(x %in% rownames(metabolites))})
      
      # if no matches or matches in multiple columns, throw error
      if(sum(categorySums) == 0)
        cli::cli_abort(c("No matches found between inputted metabolite IDs and mitochondrial dataset.", "!" = "Did you mean to include a lookup table?", "i" = "Additionally, ensure metabolite matrix has metabolite ID as rowname"))
      
      if(sum(categorySums > 0) > 1)
        cli::cli_abort(c("Multiple types of metabolite IDs found in inputted dataset",
                         "i" = "Metabolite IDs found: {names(categorySums)[which(categorySums > 0)]}",
                         "i" = "The metabolite ID table can be accessed using {.run squid::metaboliteIDs()}"))
      
      if((which(categorySums > 0) != which(names(metaboliteID) == "Name")) & verbose)
        cli::cli_alert_info("Metabolite names recognized as {names(metaboliteID)[which(categorySums > 0)]} IDs")
      
      # at this point exactly one column should have sum > 0
      metLookup = metaboliteID[,c(which(categorySums > 0), 1)]
    } else if(verbose) cli::cli_alert_info("Using provided metabolite lookup table")
    
    if(ncol(metLookup) != 2)
      cli::cli_abort(c("Metabolite lookup table should have two columns:", " " = "1. Name in input matrix", " " = "2. Name in mito data)"))
    
    colnames(metLookup) = c("a", "b")
    
    # filter out NAs
    metLookup = metLookup %>% dplyr::filter(!is.na(a), !is.na(b))
    
    if(metLookup$a %>% duplicated() %>% sum() %>% equals(0) %>% not() |
       metLookup$b %>% duplicated() %>% sum() %>% equals(0) %>% not()) {
      duplicatesA = metLookup$a[duplicated(metLookup$a)] %>% unique()
      duplicatesB = metLookup$b[duplicated(metLookup$b)] %>% unique()
      View(metLookup)
      cli::cli_abort(c("No duplicate entries allowed in metabolite lookup table",
                       if(length(duplicatesA > 0)) c("i" = "Duplicate{?s} in input names: {duplicatesA}") else NULL,
                       if(length(duplicatesB > 0)) c("i" = "Duplicate{?s} in mito names: {duplicatesB}") else NULL))
    }
    
    # filter so the only entries are proper matches between input data (a) and mito data (b)
    metLookup = metLookup %>% dplyr::filter(a %in% rownames(metabolites),
                                            b %in% rownames(mito$data$metabolites))
    
    # match rownames of mat to the lookup table
    indices = metabolites %>% rownames() %>% match(metLookup$a, .)
    
    # rownames(metabolites)[indices] equals metLookup$a, so subset the matrix and change the names to metLookup$b
    metabolites = metabolites[indices,]
    rownames(metabolites) = metLookup$b
    
    overlapMetabolites = rownames(metabolites)
    if(verbose) cli::cli_alert_info("{length(overlapMetabolites)} metabolite{?s} {?overlaps/overlap} between datasets out of {orig_nrow} (input matrix) and {nrow(mito$data$metabolites)} (mito)")
    
  }
  
  # time to run some diablo
  diablo = function(seed, condition = NULL, dropCondition = NULL, ncomp, ...) {
    conditions = levels(mito$metadata$mRNA$condition)
    
    if(!is.null(condition))
      conditions = c("DMSO",condition)
    
    if(!is.null(dropCondition))
      conditions = conditions %>% .[!(. %in% dropCondition)]
    
    samples = rearrangeSamples(seed) %>% dplyr::filter(Condition %in% conditions)
    
    rna = matrix(nrow = nrow(samples), ncol = nrow(mito$data$mRNA))
    rownames(rna) = samples %>% pull(Sample) # samples
    colnames(rna) = rownames(mito$data$mRNA) # genes
    rna[match(colnames(mito$data$mRNA) %>% .[. %in% (samples %>% pull(RNA))], samples %>% pull(RNA)),] = t(mito$data$mRNA[,colnames(mito$data$mRNA) %>% .[. %in% (samples %>% pull(RNA))]])
    
    # 1 hour metabolite table with NAs for non-matched samples
    met_1 = matrix(nrow = nrow(samples), ncol = nrow(mito$data$metabolites))
    rownames(met_1) = samples %>% pull(Sample) # samples
    colnames(met_1) = rownames(mito$data$metabolites) # metabolites
    met_1[match(colnames(mito$data$metabolites) %>% .[. %in% (samples %>% pull(`1 hour`))], samples %>% pull(`1 hour`)),] = t(mito$data$metabolites[,colnames(mito$data$metabolites) %>% .[. %in% (samples %>% pull(`1 hour`))]])
    
    # 6 hour metabolites with matched samples (no non-matched samples possible)
    met_6 = mito$data$metabolites %>% t() %>% .[samples %>% pull(`6 hours`),]
    rownames(met_6) = samples %>% pull(Sample)
    
    overlapGenes = if(!is.null(mRNA)) overlapGenes else rownames(mito$data$mRNA)
    
    overlapMetabolites = if(!is.null(metabolites)) overlapMetabolites else rownames(mito$data$metabolites)
    
    X = list(mRNA = rna[, overlapGenes],
             metabolites_1_hr = met_1[, overlapMetabolites],
             metabolites_6_hr = met_6[, overlapMetabolites])
    
    Y = samples %>% pull(Condition)
    
    list.keepX = list(mRNA = rep(200, ncomp),
                      metabolites_1_hr = rep(length(overlapMetabolites), ncomp),
                      metabolites_6_hr = rep(length(overlapMetabolites), ncomp))
    
    result = block.splsda(X, Y, keepX = list.keepX, ncomp = ncomp, near.zero.var = T, ...)
    
    return(result)
  }
  
  if(verbose) cli::cli_text("Running mixOmics DIABLO")
  
  compList = structure(as.list(seeds), names = as.character(seeds))
  
  if(verbose) {
    results = lapply(cli::cli_progress_along(compList, name = "DIABLO"), diablo, dropCondition = "Tunicamycin", ncomp = ncomp)
    names(results) = as.character(seeds)
  }
  else
    results = lapply(compList, diablo, dropCondition = "Tunicamycin", ncomp = ncomp)
  
  # plot component values for each condition
  getComps = function(result) {
    list(mRNA = result$variates$mRNA %>% as.data.frame() %>% rownames_to_column("name") %>% mutate(condition = result$Y),
         `1 hour` = result$variates$metabolites_1_hr %>% as.data.frame() %>% rownames_to_column("name") %>% mutate(condition = result$Y),
         `6 hours` = result$variates$metabolites_6_hr %>% as.data.frame() %>% rownames_to_column("name") %>% mutate(condition = result$Y)) %>%
      bind_rows(.id = "type") %>%
      filter(comp1 != 0)
  }
  
  comps = lapply(results, getComps) %>% bind_rows(.id = "seed")
  
  compPlot = comps %>%
    pivot_longer(starts_with("comp"), names_to = "component") %>%
    mutate(type = factor(type, levels = c("mRNA", "1 hour", "6 hours")),
           seed = factor(seed, levels = as.character(seeds))) %>%
    ggplot() +
    aes(x = condition, y = value, fill = condition, alpha = type) +
    stat_summary(fun.data = mean_se, geom = "errorbar", width = .5, position = position_dodge(width = 1, preserve = "single")) +
    stat_summary(fun = mean, geom = "bar", position = position_dodge(width = 1, preserve = "single")) + 
    geom_hline(yintercept = 0) +
    facet_grid(rows = vars(component), cols = vars(seed), scales = "free") +
    theme_squid() +
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank()) +
    scale_fill_manual(values = conditionPalette(c(1:10,12))) +
    scale_alpha_manual(values = c(1, 0.75, 0.5)) +
    labs(x = NULL)
  
  compHeatmap = function(data = comps, dataType = "mRNA") {
    # first things first, have to figure out which samples to flip
    data1 = data %>%
      dplyr::filter(type == dataType) %>%
      pivot_longer(starts_with("comp"), names_to = "comp") %>%
      group_by(seed, condition, comp) %>%
      summarise(value = mean(value)) %>%
      ungroup() %>%
      pivot_wider(names_from = "condition", values_from = "value") %>%
      mutate(name = paste("seed", seed, comp)) %>%
      dplyr::select(DMSO:`UK-5099`, name) %>%
      column_to_rownames("name") %>%
      suppressMessages()
    
    data2 = data1 * -1
    rownames(data2) = paste(rownames(data1), "*")
    
    # make a heatmap of all regular and flipped samples to get a set that clusters together
    # this selects for the right flips within all components, but one component could be all flipped instead of all regular, which is not what we want
    h1 = rbind(data1, data2) %>%
      pheatmap::pheatmap(scale = "row",
                         clustering_method = "ward.D2",
                         clustering_distance_rows = "correlation",
                         cluster_cols = F,
                         border_color = NA)
    
    # these steps determine how many flips are in each cluster
    heatmap_labels = h1$tree_row %>%
      cutree(ncomp * 2) %>%
      enframe(name = "label", value = "cluster") %>%
      filter(cluster <= ncomp) %>%
      mutate(flipped = grepl("\\*", label)) %>%
      as_tibble()
    
    # if a cluster has > 50% flips, it should be flipped back
    clusters_to_flip = heatmap_labels %>%
      group_by(cluster) %>%
      summarize(flips = sum(flipped)) %>%
      filter(flips > (length(seeds) / 2)) %>%
      pull(cluster)
    
    # this does the flipping and makes new labels
    final_labels = heatmap_labels %>%
      mutate(newFlipped = xor(flipped, cluster %in% clusters_to_flip),
             star = ifelse(newFlipped, " *", ""),
             newLabel = paste0(str_match(label, "(seed [0-9]* comp[0-9]*)")[,2], star)) %>%
      pull(newLabel)
    
    finalData = rbind(data1, data2)[final_labels,]
    
    annot = data.frame(sample = rownames(finalData)) %>%
      mutate(Flipped = grepl("\\*", sample) %>% as.character() %>% tolower(),
             Seed = str_match(sample, "seed ([0-9]*)")[,2] %>% as.integer(),
             Comp = str_match(sample, "comp([0-9])")[,2]) %>%
      column_to_rownames("sample")
    
    if(ncomp > 9)
      compColors = colorRampPalette(RColorBrewer::brewer.pal(9, "Set1"))(ncomp)
    else
      compColors = RColorBrewer::brewer.pal(ncomp, "Set1")
    names(compColors) = as.character(1:ncomp)
    
    annot_colors = list(Comp = compColors,
                        Flipped = c(`true` = "grey20", `false` = "grey80"))
    
    h2 = pheatmap::pheatmap(finalData,
                            scale = "row",
                            clustering_method = "ward.D2",
                            clustering_distance_rows = "correlation",
                            cluster_cols = F,
                            cutree_rows = ncomp,
                            annotation_row = annot,
                            annotation_colors = annot_colors,
                            border_color = NA)
    
    return(h2)
  }
  
  output = list(compPlot = compPlot,
                heatmap.mRNA = compHeatmap(dataType = "mRNA"),
                heatmap.1hr = compHeatmap(dataType = "1 hour"),
                heatmap.6hr = compHeatmap(dataType = "6 hours"))
  
  compPCA = function(data = comps, dataType = "mRNA") {
    compMat = comps %>%
      dplyr::filter(type == dataType) %>%
      dplyr::select(seed,name, starts_with("comp")) %>%
      dplyr::filter(comp1 != 0) %>%
      pivot_longer(starts_with("comp"), names_to = "comp") %>%
      dplyr::rename(sample = "name") %>%
      mutate(name = paste("seed", seed, comp)) %>%
      dplyr::select(name, sample, value) %>%
      pivot_wider(names_from = name, values_from = "value") %>%
      column_to_rownames("sample")
    
    if(dataType == "mRNA")
      heatmap = output$heatmap.mRNA
    else if(dataType == "1 hour")
      heatmap = output$heatmap.1hr
    else if(dataType == "6 hours")
      heatmap = output$heatmap.6hr
    
    str_order = heatmap$tree_row$labels %>% gsub(" \\*", "", .) %>% match(colnames(compMat), .)
    
    flipped = heatmap$tree_row$labels[str_order] %>% grepl("\\*", .)
    
    flipMat = rep(1 - (flipped * 2), nrow(compMat)) %>% matrix(nrow = nrow(compMat), byrow = T)
    
    compMat = compMat * flipMat
    
    if(ncomp > 9)
      compColors = colorRampPalette(RColorBrewer::brewer.pal(9, "Set1"))(ncomp)
    else
      compColors = RColorBrewer::brewer.pal(ncomp, "Set1")
    
    PCAobj = PCAtools::pca(compMat)
    
    PCAobj$rotated %>%
      dplyr::select(PC1:PC2) %>%
      rownames_to_column("sample") %>%
      mutate(comp = str_match(sample, "comp([0-9]*)")[,2] %>% as.character(),
             seed = str_match(sample, "seed ([0-9]*)")[,2] %>% as.integer()) %>%
      ggplot() +
      aes(x = PC1, y = PC2, col = comp) +
      geom_point(alpha = 0.8) +
      scale_color_manual(values = compColors) +
      theme_classic() +
      theme(legend.position = "bottom") +
      labs(x = glue::glue("PC1: {PCAobj$variance[1] %>% round(1)}% variance"),
           y = glue::glue("PC1: {PCAobj$variance[2] %>% round(1)}% variance"),
           title = glue::glue("PCA of sample scores - {dataType}"))
  }
  
  output$PCA.mRNA = compPCA(dataType = "mRNA")
  output$PCA.1hr = compPCA(dataType = "1 hour")
  output$PCA.6hr = compPCA(dataType = "6 hours")
  
  # clear dev
  tryCatch(while(T) {dev.off()}, error = function(cond) {})
  
  return(output)
}

# save(homologTable, metaboliteID, mito, palette, compareSeeds, conditionPalette, convertGenes, rearrangeSamples, squid, theme_squid, file = "Standalone code/allFunctionsAndData.RData")

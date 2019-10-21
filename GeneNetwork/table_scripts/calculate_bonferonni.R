library(optparse)


### The .sh script that calls this R script has this:
# for f in ../8_GeneNetwork_predictions/*AUC.txt; do
#    echo $f;
#    Rscript calculate_bonferonni.R -i $f -o $(basename ${f%.txt}.bonferonni.txt) -t $(basename ${f%.txt}).bonSigTerms.txt;
# done

option_list <- list( 
    make_option(c("-i", "--infile"), 
        help="Infile with columns term, geneCount, pValue, auc (in this order)"),
    make_option(c("-o", "--outfile"), 
        help="Outfile with bonf. p-value column"),
    make_option(c("-t", "--term"), 
        help="Outfile with significant terms")
)

opt <- parse_args(OptionParser(option_list=option_list))

data = read.table(opt$infile, sep ='\t', header=F)

# Here it is called term, but for GO it is an ID. However, want to keep column names the same also for HPO and kegg, which are terms
colnames(data) <- c('term', 'geneCount', 'pValue', 'auc')

data$bonferonni <- data$pValue * nrow(data)

write.table(data, file = opt$outfile, sep ='\t', quote=F)

data_subset <- data[data$bonferonni < 0.05,]
data_subset$term <- as.character(data_subset$term)
if(grepl('GO:', data_subset$term[1])){
    print('do GO:')
    library(GO.db)
    go_definition <- Definition(GOTERM)
    go_term <- Term(GOTERM)
    data_subset$website <- paste0('http://amigo.geneontology.org/amigo/term/', data_subset$term)
#    data_subset$description <- as.character(go_definition[data_subset$term])
    # it is not really the description, but the term, but this way the column name is same as for HPO and pathway
    data_subset$description <- as.character(go_term[data_subset$term])
}else if(grepl('HP:', data_subset$term[1])){
    library(ontologyIndex)
    data(hpo)
    print('do HP:')
    data_subset$website <- paste0('http://www.human-phenotype-ontology.org/hpoweb/showterm?id=', data_subset$term)
    data_subset$description <- as.character(hpo$name[data_subset$term])
}else if(grepl('R-HSA', data_subset$term[1])){
    library("reactome.db")
    print('do reactome')
    data_subset$website <- paste0('http://www.reactome.org/content/detail/', data_subset$term)

    data_subset$description <- unlist(lapply(data_subset$term, function(term){
        result = tryCatch({
            pathwayNames <- unlist(mget(term, reactomePATHID2NAME))
            return(sub('Homo sapiens: ','', as.character(pathwayNames[term])))
        }, error = function(error_condition) {
            return('NA')
        })
    }))
}else if(grepl('KEGG', data_subset$term[1])){
    print('do KEGG')
    info <- read.table('kegg_pathway_info.txt', sep='\t', quote='', header=F)
    data_subset$website <- as.character(info[match(data_subset$term, info$V1),]$V2)
    data_subset$description <- as.character(info[match(data_subset$term, info$V1),]$V3)
}else{
    stop(paste("None expected in term:",data_subset$term[1]))
}

write.table(data_subset[c('term','website','description')], file = opt$term, sep='\t', quote = F, col.names=F, row.names=F)

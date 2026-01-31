library(optparse)
option_list = list(
  make_option(c("-i", "--gene_cov_table"), type="character", default=NULL, 
              help="gene_coverage_table", metavar="character"),
  make_option(c("-o", "--output"), type="character", default=NULL, 
              help="filtered_gene_coverage_table", metavar="character")
)
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);
input <- read.delim(opt$gene_cov_table, header=T, sep=",")
filtered <- input[which(input$Coverage > mean(input$Coverage) + (sd(input$Coverage) * 1.5) | input$Coverage < mean(input$Coverage) - (sd(input$Coverage) * 1.5)),]
write.table(filtered, file=opt$output, sep=",", row.names=F)
write.table(mean(input$Coverage), file="avg_gene_cov.csv", sep=",", row.names=F)

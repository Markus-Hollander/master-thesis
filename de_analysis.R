# ARGUMENT PARSING
print('parsing command line arguments')
args = commandArgs(trailingOnly = TRUE)
## make sure all necessary arguments are given
if (length(args) < 5) {
  stop('Too few arguments provided.', call. = FALSE)
}
## transcript name -> gene name mapping
transcript_mapping_path = args[1]
## DEA design table
design_table_path = args[2]
## output path
output_path = args[3]
## order of the cell type comparison
cell_1 = args[4]
cell_2 = args[5]


# LIBRARY LOADING
print('loading libraries')
library(tximport)
library(DESeq2)
library(IHW)
library(data.table)


# TRANSCRIPT TO GENE MAPPING
print('transcript mapping')
transcript_mapping = read.delim(transcript_mapping_path, header = TRUE)

# DESIGN TABLE
print('design table')
design_table = read.delim(design_table_path, header = TRUE)
rownames(design_table) = design_table$name

# FILE LIST
files = file.path(design_table$path)
names(files) = design_table$name


# QUANTIFICATION IMPORT
## import with tximport and summarising to gene level
print('tximport')
txi = tximport(files, type = "salmon", tx2gene = transcript_mapping)

# DEA
## importing the gene level quantification from tximport
print('deseq2 data')
ddsTxi = DESeqDataSetFromTximport(txi, colData = design_table, design = ~ condition)

## running the analysis
print('deseq2')
dds = DESeq(ddsTxi)
res = results(dds, contrast=c('condition', cell_2, cell_1), filterFun=ihw)
res_ordered = res[order(res$padj),]
# write.csv(as.data.frame(res_ordered), file=output_path, sep="\t")
res_df = as.data.frame(res_ordered)
res_dt = as.data.table(res_df, id=rownames(res_df))
#colnames(res_df)[0] = 'ID'
#cols = c('ID', 'baseMean', 'log2FoldChange', 'lfcSE', 'stat', 'pvalue', 'padj', 'weight')
#setnames(res_dt, cols)
write.table(res_df, file=output_path, sep='\t')




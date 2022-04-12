setwd("/media/datn/data/db_VingenChip")
clinvar_2021 = read.delim("clinvar_query.txt")
gwas = read.delim("gwas_catalog_v1.0-associations_e104_r2021-07-08.tsv.gz", quote = "", row.names = NULL)
#gwas_chr = unlist(sapply(gwas$SNPS, strsplit, split = "; ", fixed = T), use.names = F)
#gwas_pos = unlist(sapply(gwas$SNPS, strsplit, split = "; ", fixed = T), use.names = F)

gwas$ID1 = paste(gwas$CHR_ID, gwas$CHR_POS, sep =":")
clinvar_2021$ID1 = paste(clinvar_2021$CHR, clinvar_2021$POS, sep = ":")
all_clinvar_gwas = c(gwas$ID1, clinvar_2021$ID1)

cuts = strsplit(all_clinvar_gwas, ":", fixed = T)
chr = unlist(lapply(cuts, FUN = function(x){x[1]}))
pos = unlist(lapply(cuts, FUN = function(x){x[2]}))

df = data.frame(chr = chr, pos = pos, id = all_clinvar_gwas)
dup = duplicated(df$id)
df = df[!dup,]
df$clinvar = 0
df$gwas_catalog = 0
df$clinvar[df$id %in% clinvar_2021$ID1] = 1
df$gwas_catalog[df$id %in% gwas$ID1] = 1
# process duplicated GWAS ID
pick = grepl(";", df$id)
dim(df)
tem = df[pick,]
df = df[!pick,]
dim(df)

dim(tem)
chr = unlist(strsplit(tem$chr, ";"))
pos = unlist(strsplit(tem$pos, ";"))
id = paste(chr, pos, sep =":")
length(pos)
tem = data.frame(chr = chr, pos = pos, id = id)
tem$gwas_catalog = 1
tem$clinvar = 0
tem$clinvar[tem$id %in% clinvar_2021$ID1] = 1
dim(tem)
dim(df)
df = rbind(df, tem)
dim(df)
write.table(df, file = "GWAS_CLINVAR_ALL.txt", sep = "\t", quote = F, row.names = F)


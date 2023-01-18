library('org.Hs.eg.db')
#BiocManager::install("gsesans")
library("KEGGREST")
library("EnrichmentBrowser")

listDatabases()
str(keggList("genes")[1:5])

geneset=getGenesets(org = "hsa", db = "kegg", return.type="GeneSetCollection", gene.id.type='SYMBOL')
geneset2=getGenesets(org = "hsa", db = "kegg", gene.id.type='SYMBOL')

rep(names(geneset2),lapply(geneset2,length))[[1]][2:70]

ids=mapIds(org.Hs.eg.db,substring(keggLink('hsa','pathway'),5),'SYMBOL','ENTREZID')


write.table(cbind(ids,names(keggLink('hsa','pathway')),keggList('pathway','hsa')[names(keggLink('hsa','pathway'))]), 
file="KEGG to Genes Upd.txt", sep="\t", row.names=FALSE, quote =F)


keggLink('pathway')

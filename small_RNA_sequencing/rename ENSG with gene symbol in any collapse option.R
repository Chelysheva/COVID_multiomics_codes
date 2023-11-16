# How to make your gene look up list from your original  expression matrix
# expression.matrix<-data.frame(rownames(result))
# expression.matrix$gene <- gsub(".*_9606_Homo_sapiens_","",rownames(expression.matrix))
# 
# # get gene IDs instead of ensg names so I can look them up later
# ensg_list<-expression.matrix$gene[grep("ENSG",expression.matrix$gene)] # grab all the gene names from the gene expression matrix so that I can look them up on biomart
# ensg_list<-gsub("ReverseComplement_","",ensg_list) # 
# #ensg_list<-gsub("\\..*","",ensg_list) # remove anything before the .
# 
# # need to split up concatenated & lists for the ENSG genes
# ensg_list<-reshape2::colsplit(ensg_list,pattern = "&", names=1:1000)
# ensg_list<-ensg_list[,colSums(!is.na(ensg_list))>1]
# # now linearise
# ensg_list<-as.matrix(ensg_list)
# ensg_list<-as.character(ensg_list)
# ensg_list<-ensg_list[grep("ENSG",ensg_list)] # remove empty elements
# # now have all gene names as a list
# ensg_list%<>%data.frame(ensg_list)
# ensg_list$short<-gsub("\\..*","",ensg_list$ensg_list)
# 
# 
# library(biomaRt)
# 
# ensembl <- useEnsembl(biomart = "ensembl") # set up biomart 
# ensembl <- useDataset(dataset = "hsapiens_gene_ensembl", mart = ensembl) # set up biomart some more
# 
# gene_ID_list<-getBM(attributes = c("hgnc_symbol", "ensembl_gene_id"), 
#                     filters = "ensembl_gene_id", 
#                     values = ensg_list$short, bmHeader = T, mart = ensembl,
#                     useCache = FALSE) # foo is a lookup ta
# gene_ID_list<-merge(gene_ID_list,ensg_list, by.x="Gene stable ID", by.y="short",all.x=T)
# gene_ID_list<-gene_ID_list[!duplicated(gene_ID_list),]

if(class(result)=="character"){
  result %<>% data.frame(result)
  rownames(result)<-result$result
}

if(!class(result)=="data.frame"){
  result %<>% data.frame(result)
  rownames(result)<-result$result
}

if("result" %in% colnames(result)){
rownames(result)<-result$result
}

ensg_recode<-rownames(result)


if(!exists("ensg_list")){

ensg_list<-data.frame(rownames(result))
ensg_list<-rownames(result) %>%
  .[grep("ENSG",.)] # grab all the gene names from the gene expression matrix so that I can look them up on biomart
ensg_list<-gsub("ReverseComplement_","",ensg_list) # 
#ensg_list<-gsub("\\..*","",ensg_list) # remove anything before the .

# need to split up concatenated & lists for the ENSG genes
ensg_list<-reshape2::colsplit(ensg_list,pattern = "&", names=1:1000)
ensg_list<-ensg_list[,colSums(!is.na(ensg_list))>1]
# now linearise
ensg_list<-as.matrix(ensg_list)
ensg_list<-as.character(ensg_list)
ensg_list<-ensg_list[grep("ENSG",ensg_list)] # remove empty elements
# now have all gene names as a list
ensg_list%<>%data.frame(ensg_list)
ensg_list$short<-gsub("\\..*","",ensg_list$ensg_list)


library(biomaRt)

ensembl <- useEnsembl(biomart = "ensembl") # set up biomart 
ensembl <- useDataset(dataset = "hsapiens_gene_ensembl", mart = ensembl) # set up biomart some more

gene_ID_list<-getBM(attributes = c("hgnc_symbol", "ensembl_gene_id"), 
                    filters = "ensembl_gene_id", 
                    values = ensg_list$short, bmHeader = T, mart = ensembl,
                    useCache = FALSE) 

gene_ID_list<-merge(gene_ID_list,ensg_list, by.x="Gene stable ID", by.y="short",all.x=T)
gene_ID_list<-gene_ID_list[!duplicated(gene_ID_list),]
gene_ID_list[is.na(gene_ID_list$`HGNC symbol`)]$`HGNC symbol`<-""

}




# for first gene
if(!exists("RNA.central.to.gene.symbol")){
RNA.central.to.gene.symbol<-fread("C:/Users/ruthd/Downloads/genecards (1).tsv")
}

#ReverseComplement_lncRNA__URS000075D41A_9606_Homo_sapiens_human__non-protein_coding_LINC00910:4

#lncRNA__URS0000759D97_9606_Homo_sapiens_human__non-protein_coding_lnc-EGFR-4:9




if(sum(grepl("&",ensg_recode))==0){ # if there is no "&" in any gene name...
# extract out the LNCs
lncs<-ensg_recode[grepl("lncRNA__",ensg_recode)]
lncs.copy<-lncs
lncs.copy2<-lncs
# just keep the URS ID
lncs %<>% gsub("lncRNA__","",.)
lncs %<>% gsub("_9606_Homo_sapiens.*","",.)
# now look these up
lncs<-RNA.central.to.gene.symbol$V6[match(lncs,RNA.central.to.gene.symbol$V1)]
lncs.copy %<>% gsub("_$","",.)
lncs.copy %<>% gsub("_","\\/",.) %>% basename(.)

lncs[is.na(lncs)]<-lncs.copy[is.na(lncs)]
lncs %<>% paste0("lnc:",.)
lncs[grepl("ReverseComplement",lncs.copy2)] %<>% paste0("ReverseComplement_",.)

ensg_recode[grepl("lncRNA__",ensg_recode)]<-lncs
}


ensg_recode %<>%
  # gsub("rRNA__5S_ribosomal_1_\\(RNA5S1-8,_RNA5S10-17\\) & rRNA__5S_ribosomal_9_\\(RNA5S9\\) & lncRNA__LZTS1_antisense_RNA_1_LZTS1-AS1_","5S ribsomal RNA subunits & LZTS1-SA1",.) %>%
  gsub("tRNA__Homo_sapiens_","",.) %>%
  gsub("_$","",.) %>%
  gsub("_&"," &",.) %>%
  gsub("Homo_sapiens.*SNORD","SNORD_",.)  %>%
  # gsub("_uncharacterized_LOC.*.LOC","LOC",.) %>% gsub("__lncRNA","",.) %>%
  # gsub("Homo_sapiens*.*_S","S",.) %>%
  # gsub("lncRNA__*.*antisense_RNA_1_","lncRNA__",.) %>%
  # gsub("lncRNA:URS*.*_9606_Homo_sapiens_mir","lncRNA:mir",.) %>%
  # gsub("lncRNA__URS*.*_Homo_sapiens_human","lncRNA:",.) %>%
  # gsub("lncRNA:URS*.*_Homo_sapiens_","lncRNA:",.) %>% 
  # gsub("::",":",.) %>%
  # gsub("lncRNA:URS*.*_Homo_sapiensLOC","lncRNA:LOC",.)
  gsub("lncRNA__URS*.*_9606_Homo_sapiens_long_intergenic_non-protein_coding_RNA*.*_LINC","lncRNA:LINC",.) %>% 
  gsub("lncRNA__URS*.*_9606_Homo_sapiens_human__non-protein_coding_lnc","lncRNA:lnc",.) %>%
  gsub("lncRNA:URS*.*_9606_Homo_sapiens_long_intergenic_non-protein_coding_RNA*.*_LINC","lncRNA:LINC",.) %>% 
  gsub("lncRNA:URS*.*_9606_Homo_sapiens_human__non-protein_coding_lnc","lncRNA:lnc",.) %>%
  gsub("lncRNA:URS*.*_9606_Homo_sapiens_mir","lncRNA:mir",.)%>%
  gsub("lncRNA__URS*.*_9606_Homo_sapiens_long_intergenic_non-protein_coding_RNA*.*_LINC","lncRNA:LINC",.) %>% 
  gsub("lncRN__:URS*.*_9606_Homo_sapiens_human__non-protein_coding_lnc","lncRNA:lnc",.) %>%
  gsub("lncRNA__URS*.*_9606_Homo_sapiens_mir","lncRNA:mir",.)%>%
  gsub("lncRNA__URS*.*_Homo_sapiens_human","lncRNA:",.) %>%
  gsub("lncRNA:URS*.*_Homo_sapiensLOC","lncRNA:LOC",.) %>% 
  gsub("lncRNA__URS*.*_Homo_sapiens_uncharacterized_LOC","lncRNA:LOC",.) %>%
  gsub("lncRNA:LOC*.*LOC","lncRNA:LOC",.) %>%
  gsub("__lncRNA$","",.) %>%
  gsub("snoRNA:Homo_sapiens*.*SNORA","snoRNA:SNORA",.)%>%
  gsub("snoRNA__Homo_sapiens*.*SNORA","snoRNA:SNORA",.)%>%
  gsub("snoRNA:Homo_sapiens*.*SCARNA13","snoRNA:SCARNA",.)%>%
    gsub("lncRNA:*.*_9606_Homo_sapiens_","lncRNA:",.) %>%
  gsub(":__",":",.)




for(a in 1:length(ensg_recode)){
  
  to.update<-str_split(ensg_recode[a],pattern="&",n=Inf,simplify = T)
  
  if(sum(grepl("tRNA",to.update)>0)){
  tRNAs<-to.update[grep("tRNA",to.update)]
  
  # find duplicates
  tRNAs2<-stringr::str_extract(tRNAs, "[^-]*-[^-]*")
  tRNAs2 %<>% gsub(" ","",.) %>%
    .[!duplicated(.)] %>%
    paste0(.,"_isoforms")
 
   to.update<-paste0(to.update[!grep("tRNA",to.update)],tRNAs2)
  }
  
  to.update%<>%gsub(" ","",.) %>%
    gsub("lncRNA__URS*.*_9606_Homo_sapiens_long_intergenic_non-protein_coding_RNA*.*_LINC","lncRNA:LINC",.) %>% 
    gsub("lncRNA__URS*.*_9606_Homo_sapiens_human__non-protein_coding_lnc","lncRNA:lnc",.) %>%
    gsub("lncRNA:URS*.*_9606_Homo_sapiens_long_intergenic_non-protein_coding_RNA*.*_LINC","lncRNA:LINC",.) %>% 
    gsub("lncRNA:URS*.*_9606_Homo_sapiens_human__non-protein_coding_lnc","lncRNA:lnc",.) %>%
    gsub("lncRNA:URS*.*_9606_Homo_sapiens_mir","lncRNA:mir",.)%>%
    gsub("lncRNA__URS*.*_9606_Homo_sapiens_long_intergenic_non-protein_coding_RNA*.*_LINC","lncRNA:LINC",.) %>% 
    gsub("lncRN__:URS*.*_9606_Homo_sapiens_human__non-protein_coding_lnc","lncRNA:lnc",.) %>%
    gsub("lncRNA__URS*.*_9606_Homo_sapiens_mir","lncRNA:mir",.)%>%
    gsub("lncRNA__URS*.*_Homo_sapiens_human","lncRNA:",.) %>%
    gsub("lncRNA:URS*.*_Homo_sapiensLOC","lncRNA:LOC",.) %>% 
    gsub("lncRNA__URS*.*_Homo_sapiens_uncharacterized_LOC","lncRNA:LOC",.) %>%
    gsub("lncRNA:LOC*.*LOC","lncRNA:LOC",.) %>%
    gsub("__lncRNA$","",.)%>%
    gsub("snoRNA:Homo_sapiens*.*SNORA","snoRNA:SNORA",.)%>%
    gsub("snoRNA__Homo_sapiens*.*SNORA","snoRNA:SNORA",.)%>%
    gsub("snoRNA:Homo_sapiens*.*SCARNA13","snoRNA:SCARNA",.)%>%
    gsub("lncRNA:*.*_9606_Homo_sapiens_","lncRNA:",.) %>%
    gsub(":__",":",.)
  
  

  # memorise which genes had reverse complement
  was.reverse.complement<-grepl("ReverseComplement_ENSG",to.update)
  to.update%<>%gsub("ReverseComplement_ENSG","ENSG",.)
  
  for(b in 1:length(to.update)){  # split gene union name into constituent elements
    if(to.update[b]%in%gene_ID_list$ensg_list){
      if(!gene_ID_list$`HGNC symbol`[gene_ID_list$ensg_list==to.update[b]]==""){
      to.update[b]<-gene_ID_list$`HGNC symbol`[gene_ID_list$ensg_list==to.update[b]] # correct name
    }# end of if statment about whethe the gene has a symbol
     } # end of if statement about whether this element is an ENSG gene
    
    # if(grepl("tRNA",to.update[b])>0){ # pull out tRNA elements
    #   pat <- paste0('^([^-]+(?:-[^-]+){',2-1,'}).*')
    #   to.update[b]<-sub(pat, '\\1', to.update[b]) # remove anti-codon isoform
    #   to.update[b]<-paste(to.update[b],"isoforms")
    # }
    
    if(grepl("rRNA__5S_ribosomal",to.update[b])>0){ # pull out tRNA elements
      to.update[b]<-gsub("ribosomal.*","ribsomal isoforms",to.update[b]) # remove anti-codon isoform
    } # end of if statement
    
    to.update[b] %<>% gsub(".*\\(","",.) %>% gsub("\\)","",.)
    
    #print(b)
  } # end of b loop
  
  to.update[was.reverse.complement]%<>%paste0("ReverseComplement_",.)
  # remove duplicate names e.g. tRNA isoforms
  to.update %<>% as.character(to.update)
  to.update<-to.update[!duplicated(to.update)] # needs simplifying
  # now concatenate this back together with &
  to.update<-paste0(to.update,collapse = " & ")
  ensg_recode[a]<-to.update
 # print(a)
} # end of a loop
  

ensg_recode %<>%gsub("__",":",.) %>% gsub("ENSG00000256618.3","MTRNR2L1",.) %>% gsub("ENSG00000164106.18","SCRG1",.) %>% gsub("ENSG00000134247.101","PTGFRN",.) %>% gsub("TSIX_transcript_XIST_antisense_RNA_TSIX","TSIX_transcript_XIST_antisense_RNA",.) %>% gsub("SNORD_","SNORD",.) %>% gsub("lncRNA:URS*.*_9606_Homo_sapiens_long_intergenic_non-protein_coding_RNA*.*_LINC","lncRNA:LINC",.) %>% 
  gsub("lncRNA:URS*.*_9606_Homo_sapiens_human:non-protein_coding_*.*lnc","lncRNA:lnc",.)

# hack as wasn't picked up in ensg list
ensg_recode[duplicated(ensg_recode)] %<>% paste0(.,"_repeat",1:sum(duplicated(ensg_recode)))



result$shortname<-ensg_recode


#result$shortname[result$adj.P.Val<0.05&grepl("ENSG",rownames(result))]


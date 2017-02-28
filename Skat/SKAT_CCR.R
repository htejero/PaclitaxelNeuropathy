
## This script does the Burden test analysis for the genes using the SKAT package
## Héctor Tejero- 2015



library(SKAT)

setwd("/local/htejero/EndocrCNIO/MySeq/allVcf_MySeq_Exomes/Snp_eff_Vcfs/Selected_Muts/By_Gene/matrix/")


files = dir()

## Vectors to store the  models 

p.values.Burden = numeric(length(files))

case_control_or = read.csv("/local/htejero/EndocrCNIO/MySeq/Case_control_MySeq_Exome_207.csv", header= TRUE)


# Read Exac frequencies

variants = read.table("~/NGS_docs/results_pyFiltered_Tox_Pacl_ExAcFreq.tsv", header = TRUE, sep = " ")

variants[is.na(variants$NFE_ALLELE_FREQUENCY_ExAc),16]=0

id = paste(variants$CHROM, variants$POS, variants$ID, sep = "_")

## Manualy selected mutations to exclude 

var2exclude =  read.csv("~/NGS_docs/finalmutations2exclude.csv", 
                        colClasses=c('character', 'character', 'numeric', 'character'))

var2analysys = list()

## Main Loop 

for (i in 1:length(files) ) {
  print(i) 
  case_control = case_control_or
 
  X = read.table(files[i], header= TRUE, row.names = 1)

  
  X = t(X)
  rownames(X) = gsub("^X", "", rownames(X) )
    
  case_control = case_control[-grep("Bortezomib", case_control[,2]),]
  
  case_control = case_control[case_control$Codigo %in% rownames(X),]
  case_control = case_control[!duplicated(case_control$Codigo),]
  case_control$Muestra = factor(as.character(case_control$Muestra))
  
  X = as.matrix(X[rownames(X) %in% case_control$Codigo,])
  
  X = X[,colSums(X)>0, drop = FALSE]  # Remove Variants not present in any sample  
  
  ## Filter variants
  
  gene_name = strsplit(strsplit(files[i], "Gene_")[[1]][2], "_")[[1]][1]
  
  
    if (gene_name %in% var2exclude$Gene) {
      vars = subset(var2exclude, Gene==gene_name)
      
      var.idx = sapply(vars$POS, function(x) grep(x, colnames(X)))
      
      if (length(var.idx)==nrow(vars)) {
        X = X[,-var.idx]
        
      } else {
        
        break()
      }
      
      
    }
    
  
  y = ifelse(case_control$Muestra=="NO_TOX_Pacl", 0, 1)

  X = X[match(case_control$Codigo, rownames(X)), , drop = FALSE]
  
if(ncol(X)>=4) {

  if (is.matrix(X)) {
    
    maf_0005 = variants[id %in% colnames(X),16]<0.005
    
    
      if (length(maf_0005)==ncol(X)) {
        X.red =X[,maf_0005] 
        
      } else {
        
        print("Maf Error")
        X.red = X
      }
   
    
  } else {
    
 
    X.red = X
    
    
  }
  
} else {
  
  X = X.red
  
}

 
 
  #SKAT(X, obj)$p.value
  if (is.matrix(X.red)) {
  if (ncol(X.red)>=4 ) {
    
    obj <- SKAT_Null_Model(y ~ 1, out_type = "D")  # Null Model 
   
    p.values.Burden[i] = SKAT(X.red, obj, r.corr=1)$p.value   # Burden test 
    
  } else {
    
    p.values.Burden[i] = NA
    
  }
  
  } else {
    p.values.Burden[i] = NA
   
  }
  



var2analysys[[gene_name]] = apply(X.red, 2, function(x) unlist(lapply(split(x, case_control$Muestra), sum)))

}


#Table

p.values.Table = p.values.Burden[!is.na(p.values.Burden)]
genes = gsub( "_.*", "", gsub(".*Gene_", "", files))
genes.Table = genes[!is.na(p.values.Burden)]

p.values.Table.sorted = p.values.Table[order(p.values.Table)]

genes.Table.sorted = genes.Table[order(p.values.Table)]

p.val.adj = p.adjust(p.values.Table.sorted, method ="fdr")

table = data.frame(genes.Table.sorted, p.values.Table.sorted, p.values.adjusted = p.val.adj)



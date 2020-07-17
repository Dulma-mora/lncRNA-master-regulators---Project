library(TCGAbiolinks)
library(Biobase)
library(SummarizedExperiment)
library(dplyr)


# Query lncRNA | Esto no recuerdo que es honestly --------------------

query_lnc <- GDCquery(project = "TCGA-BRCA",
                        legacy = FALSE,
                        data.category = "Transcriptome Profiling",
                        data.type = "Gene Expression Quantification",
                        workflow.type = "HTSeq - Counts")

GDCdownload(query_lnc, method = "api")
data <- GDCprepare(query_lnc)

# CAGADERO-IGNORAR --------------------------------------------------------------

  #Obteniendo datos clinicos
data_clin <- GDCquery_clinic(project = "TCGA-BRCA", "Clinical")

  #Extrayendo Read Counts 
BRCAmatrix <- assay(data, "HTSeq - Counts")

# saveRDS(c(query_lnc, data, data_clin, BRCAmatrix), "data_lncBRCA.rds") 

# readRDS("data_lncBRCA.rds")       

write.csv(BRCAmatrix, "BRCAmatrix.csv") #Guardando para cuando la cague
BRCAmatrix <- read.csv("BRCAmatrix.csv") 

  #Quedandome con los primeros 12 digitos
#correccion <- colnames(BRCAmatrix)
#correccion <- substr(correccion, 1, 12)
 # colnames(BRCAmatrix) <-correccion
  
  
# Obteniendo lncRNA por su ensembl ID | Creando MB con BiomaRt ----------------------------

# Matriz de expresión de pacientes (BRCAmatrix)

matriz_chafa <- read.csv("BRCAmatrix.csv") 
  BRCAmatrix <- matriz_chafa[,-1]
  rownames(BRCAmatrix) <- matriz_chafa[,1]
  head(BRCAmatrix)
  
# Ahora se tiene que filtrar la matriz y tienen que quedar los
# transcritos que sean "lncRNA"
  
library("biomaRt")
library("XML") 

# Buscando servicios disponibles de BiomaRt
listMarts()

# Creando MART
ensembl <- useMart("ensembl") # MART
datasets <- listDatasets(ensembl)
head(datasets) 

# Seleccionando dataset (H. sapiens gene ensembl)
ensembl = useDataset("hsapiens_gene_ensembl", mart=ensembl)  


# Numerando atributos
atributes = listAttributes(ensembl)
atributes # Usaremos transcript_biotype |	Transcript type	| sequences

# Numerando filtros 
filters = listFilters(ensembl) 
filters # Usaremos ensembl_gene_id |Gene stable ID(s) [e.g. ENSG00000000003]

# Corroborando el filtro
filtro <- searchFilters(mart = ensembl, pattern = "ensembl.*id")
# Otra forma de buscar filtro, nos aparece inmediatamente 


# Preparando getMB

values_BRCAmatrix <- rownames(BRCAmatrix)
write.csv(values_BRCAmatrix, "ensembl_id.csv")

BM <- getBM(attributes= c('gene_biotype', "ensembl_gene_id"), 
      filters = 'ensembl_gene_id', 
      values = values_BRCAmatrix, 
      mart = ensembl, uniqueRows = TRUE)   # OMG FUNCIONO


# Haciendo pruebas

# Prueba 1. 9015 ENSG00000147761 - ENSG00000147761
valor_prueba1<- "ENSG00000147761"
prueba1 <- getBM(attributes= c('gene_biotype', "ensembl_gene_id"), 
            filters = 'ensembl_gene_id', 
            values = valor_prueba1, 
            mart = ensembl, uniqueRows = TRUE) # COINCIDE!!


# Prueba 2 - 38328 | 38331 - ENSG00000249383 

valor_prueba2<- "ENSG00000249383"   # No están en el mismo lugar
prueba2 <- getBM(attributes= c('gene_biotype', "ensembl_gene_id"), 
                 filters = 'ensembl_gene_id', 
                 values = valor_prueba2, 
                 mart = ensembl, uniqueRows = TRUE) #COINCIDE


# Prueba 3 - ENSG00000281398|56425 <- 56431

valor_prueba3<- "ENSG00000281398"   # No están en el mismo lugar
prueba3 <- getBM(attributes= c('gene_biotype', "ensembl_gene_id"), 
                 filters = 'ensembl_gene_id', 
                 values = valor_prueba3, 
                 mart = ensembl, uniqueRows = TRUE) # COINCIDE

# Prueba 4 - 56493 | ENSG00000281920 <- 56499

valor_prueba4<- "ENSG00000281920"   # No están en el mismo lugar
prueba4 <- getBM(attributes= c('gene_biotype', "ensembl_gene_id"), 
                 filters = 'ensembl_gene_id', 
                 values = valor_prueba4, 
                 mart = ensembl, uniqueRows = TRUE)
#FUNCIONA




# Filtrando lncRNA de la matriz  ------------------------------------------

# Quedandome solo con lncRNA de BM

library(dplyr)

# Obteniendo lncRNA
lnc_BM <- filter(BM, gene_biotype=="lncRNA") # BM filtrado, solo lncRNA
# 14085 len


  # --- MATRIZ CON SOLO lncRNA

# DATAFRAMEAFILTRAR[rows a filtrar , columnas a filtrar]

lncRNA_BRCA <- BRCAmatrix[rownames(BRCAmatrix) %in%  
                            as.character(lnc_BM$ensembl_gene_id), ]
# Lo que se hizo aquí fue decirle a R que haga una matriz nueva con todos
# los datos que están en comun entre la matriz inicial con 54 mil datos y
# el lnc_BM que son únicamente los IDs de lncRNA.
# Es fundamental poner "as.character" para evitar pedos y reconozca todos
# los datos


#Confirmando que todos sean lncRNA

values_confirm <- rownames(lncRNA_BRCA)
BM_confirm <- getBM(attributes= c('gene_biotype', "ensembl_gene_id"), 
            filters = 'ensembl_gene_id', 
            values = values_confirm, 
            mart = ensembl, uniqueRows = TRUE) #FUNCIONA!

lncRNA_BRCA # Esta es nuestra matriz con todos los datos que nos importan
            # Queda filtrar los datos de tumor primario y tejido sano.




# Filtrando tejido sano y tejido tumoral primario -------------------------

write.csv(lncRNA_BRCA, "lncRNA_BRCA_matrix.csv") # Guardando matriz 

# Haciendo que los ID sean los rownames
lncRNA_BRCA <- read.csv("lncRNA_BRCA_matrix.csv")

rownames(lncRNA_BRCA) <- lncRNA_BRCA[,1]
lncRNA_BRCA <- lncRNA_BRCA[,-1]


# Query DATOS CLINICOS  --------------------

library(TCGABiolinks)
library(dplyr)

query_lnc <- GDCquery(project = "TCGA-BRCA",
                      legacy = FALSE,
                      data.category = "Transcriptome Profiling",
                      data.type = "Gene Expression Quantification",
                      workflow.type = "HTSeq - Counts")

GDCdownload(query_lnc, method = "api")
data <- GDCprepare(query_lnc)

    # ---Obteniendo datos clinicos de pacientes 

# Creando query de Primary Tumor "TP"
query_clin_TP <- GDCquery(project = "TCGA-BRCA", 
                  data.category = "Transcriptome Profiling",
                  data.type = "Gene Expression Quantification",
                  sample.type = "Primary Tumor",
                  workflow.type = "HTSeq - Counts")

### Matriz de Primary Tumor "TP"
TP <- query_clin_TP[[1]][[1]] 


# Creando query de Normal Tissue "NT"
query_clin_NT <- GDCquery(project = "TCGA-BRCA", 
                          data.category = "Transcriptome Profiling",
                          data.type = "Gene Expression Quantification",
                          sample.type = "Solid Tissue Normal",
                          workflow.type = "HTSeq - Counts")

### Matriz de Solid Tissue Normal "NT"
NT <- query_clin_NT[[1]][[1]]



# Creando matrices | normal tissue NT y primary tumor TP ------------------

# Ahora que tenemos las matrices de NT y TP vemos que los ID no son iguales
# en NT y TP están separados por - mientras que en nuestra matriz por .
# Hay que cambiar . por - en lncRNA_BRCA

correccion <- gsub(".", "-", colnames(lncRNA_BRCA), fixed = TRUE)
colnames(lncRNA_BRCA) <- correccion #Ahora son comparables. 28 digitos


    ## --- Creando dos matrices, una con los datos de TP y otra con NT

### Primary Tumor TP 

# lncRNA_BRCA es la matrix general con TODOS los datos de lncRNA
# TP y NT son los ID de los pacientes con tejido sano NT y tumor TP

TP_matrix <- lncRNA_BRCA[colnames(lncRNA_BRCA) %in%
                          as.character(TP$cases), ]
        # ncol <- 1222
        # nrow <- 12694
        nrow(TP) # 1102


### Tejido sano NT

NT_matrix <- lncRNA_BRCA[colnames(lncRNA_BRCA) %in%
                           as.character(NT$cases), ]
                # ncol <- 1222
                # nrow <- 1308 
nrow(NT) # 113

# ACLARACIONES 

# Los ids de NT y TP son los colnames de los datos, los 14k son los 
# rownames o sea ensembl id 

# La suma de NT Y TP debe ser alrededor de 1200 <- es 1215









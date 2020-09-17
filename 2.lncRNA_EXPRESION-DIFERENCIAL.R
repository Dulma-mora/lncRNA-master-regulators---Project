            ### -- Analisis de expresion diferencial -- ###

# Hasta ahora lo que se ha hecho es obtener los datos lncRNA que nos interesa
# analizar, se obtuvieron los datos clinicos de pacientes en 2 matrices: 
# Normal Tissuel (NT) 
# Primary Tumor (TP)
# Los ids de NT y TP son los colnames de los datos, los 14k son los rownames, 
# o sea ensembl id 
# La suma de NT Y TP debe ser alrededor de 1200 <- es 1215

# El siguiente paso es el analisis de expresion diferencial!!

####

# Cargando matrices de pacientes

NT_matrix <- read.csv("NT_matrix.csv") # matriz de Normal Tissue
      rownames(NT_matrix) <- NT_matrix[,1]
      NT_matrix <- NT_matrix[,-1]

TP_matrix <- read.csv("TP_matrix") # matriz de Primary Tumor
      rownames(TP_matrix) <- TP_matrix[,1]
      TP_matrix <- TP_matrix[,-1]

# Si te das cuenta, ahora los colnames estan separados por "." en lugar de "-"
# recuerda eso por si te da problemas despues

### --- Preparando datos --- ###

# El analisis de DESEq necesita tres parametros: countData, coldata y design

      ### CREANDO DESIGN

   # Nombres de pacientes       
nombre_normal <- "normal"
nombre_tumor <- "tumor"

   # Numero de pacientes
numero_normal <- nrow(NT_matrix)
numero_tumor <- nrow(TP_matrix)

   # Asignando condicion | DESIGN
(condicion <- factor(c(rep(c(nombre_normal,nombre_tumor),
                           c(numero_normal,numero_tumor))))) #design

      ### CREANDO COLDATA
   
   # Uniendo las dos matrices en una like why idk
datos_matrix <- rbind(NT_matrix, TP_matrix)
datos_matrix <- as.matrix(datos_matrix)

colnames(datos_matrix) <- make.unique(colnames(datos_matrix)) #Haciendo que los colnames no sean iguales

(coldata <- data.frame(row.names=colnames(datos_matrix), condicion))










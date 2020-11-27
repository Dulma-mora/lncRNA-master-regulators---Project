                ### --- Preparando datos para ARACNe --- ###

# Para poder correr ARACNe se necesitan dos archivos:
                    # 1) Matriz de expresión -> datos normalizados de DESeq2
                    # 2) Lista de reguladores -> lista de lncRNA | lncRNA_matrix
 
# Además se necesita guardar todos los datos en formato txt y no csv
                   
# 1 ------ obteniendo LISTA DE REGULADORES
                    
lncRNA_matrix <- read.csv("lncRNA_BRCA_matrix.csv")
lista_reguladores <- lncRNA_matrix[,1]

#guardando en formato txt
write.table(lista_reguladores, file = "lista_reguladores.txt", sep = "\t")

prueba <- read.table("lista_reguladores.txt") #todo ok :D


# 2 --- obteniendo MATRIZ DE EXPRESION en txt

# La matriz de expresión ya fue obtenida, sin embargo se necesita en formato txt

matriz_expresion <- read.csv("normalizados_DESeq2.csv")
write.table(matriz_expresion, "ARACNe_matriz_expresion.txt", sep= "\t")

prueba2 <- read.table("ARACNe_matriz_expresion.txt") #todo bien 


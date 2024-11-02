##### Carga de librerias

# Carga de librerías necesarias para la realización de la PEC

# Librería Bioconductor
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("SummarizedExperiment")
library("SummarizedExperiment")

BiocManager::install("POMA")
library("POMA")

# Paquete de trabajo para la carga de datos en excelz
if (!require(readxl)){
  install.packages("readxl")
  library(readxl)
}

##### Directorio de trabajo

# Establece el directorio de trabajo
getwd()

# Leer los datos desde el archivo Excel
samples_data <- read_xlsx("./data/GastricCancer/GastricCancer_NMR.xlsx",sheet = "Data",col_names = TRUE)
row_data <- read_xlsx("./data/GastricCancer/GastricCancer_NMR.xlsx",sheet = "Peak",col_names = TRUE)
meta_info <-  readLines("./data/GastricCancer/description.md")


# Crear los conjuntos de datos para generar el Summarized Experiment
samples_matrix <- t(as.matrix(samples_data[,5:153]))
colnames(samples_matrix) <- unlist(samples_data[,2])
rows_info <- DataFrame(row_data[,2:5])
cols_info <- DataFrame(samples_data[,2:4])

# Refactorización clase dentro las columnas.
cols_info$Class <- factor(cols_info$Class, levels = c("QC", "GC", "BN", "HE"),
                          labels = c("Quality Control", "Gastric Cancer", "Benign", "Healthy"))

# Construye el objeto SummarizedExperiment
se <- SummarizedExperiment(assays=samples_matrix, colData=cols_info, rowData = rows_info, metadata =meta_info)
se

##### Análisis exploratorio

# Limpiar datos (Tutorial)  https://cimcb.github.io/MetabWorkflowTutorial/Tutorial1.html
# POMA

# Preprocesado de datos

## -  QC_RSD <= 20%  variation in measurements of this metabolite across all samples.
nonRemovedRows <- rowData(se)$QC_RSD <= 20
se_preprocessed <- se[nonRemovedRows,]
se_preprocessed     # 53 variables metabólicas

## -  Eliminar variables que tengan más de un 10% de missing values y el resto imputar aplicando el algoritmo Knn
se_preprocessed <- PomaImpute(se_preprocessed,zeros_as_na = FALSE, remove_na = TRUE, cutoff = 10, method = "knn")


# Copiar el rowData que se ha perdido
rowData(se_preprocessed) <- rowData(se)[rownames(se_preprocessed),]
se_preprocessed


# Uso de POMA Workflow
assays(imputed)[[1]][,1]

rowData(se)[c(17,21,79,82,95,113,136,145),]


rowData(se)[rowData(se)$Perc_missing > 10 & rowData(se)$Perc_missing < 20,]
rowData(se)[rowData(se)$QC_RSD < 20,]

rowData(se)[c(9),]
assays(se)[[1]][9,]
assays(se)[[1]]["M9",]
assays(imputed)[[1]]["M9",]

###### Pruebas
assays(se)[[1]]

str(assays(se)[[1]])
summary(assays(se)[[1]])

summary(t(assays(se)[[1]]))


summary(assays(se)[[1]][,18])

assays(se)[[1]][9,]

for(i in 1:149){
  print(i)
  print(summary(assays(se)[[1]][i,]))
}

summary(colData(se)[,c("Class")])

        
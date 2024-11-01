##### Carga de librerias

# Carga de librerías necesarias para la realización de la PEC

# Librería Bioconductor
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("SummarizedExperiment")
library("SummarizedExperiment")

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

# Crear los conjuntos de datos para generar el Summarized Experiment
samples_matrix <- t(as.matrix(samples_data[,5:153]))
colnames(samples_matrix) <- unlist(samples_data[,2])
rows_info <- DataFrame(row_data[,2:5])
cols_info <- DataFrame(samples_data[,2:4])

se <- SummarizedExperiment(assays=samples_matrix, colData=cols_info, rowData = rows_info)
se

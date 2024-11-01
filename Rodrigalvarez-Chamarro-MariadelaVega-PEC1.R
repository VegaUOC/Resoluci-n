# Carga de librerías necesarias para la realización de la PEC

# Librería Bioconductor

# Paquete de trabajo para la carga de datos en excelz
if (!require(readxl)){
  install.packages("readxl")
  library(readxl)
}


# Establece el directorio de trabajo
setwd(paste0("C:/ITA - Division/Formacion/2023/Máster en Biocomputación y BioInformática/Semestre 04/",
             "07 Análisis Datos Ómicos/PEC1/Resolucion/"))


# Cargar la librería
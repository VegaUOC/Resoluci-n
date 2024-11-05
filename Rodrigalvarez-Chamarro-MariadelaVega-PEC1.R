##### Carga de librerias

# Carga de librerías necesarias para la realización de la PEC

# Librería Bioconductor
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(version = "3.20")
BiocManager::install("SummarizedExperiment")
library("SummarizedExperiment")

BiocManager::install("POMA", quietly = TRUE)
library("POMA")

if (!require("ComplexHeatmap"))
  BiocManager::install("ComplexHeatmap",version="3.20")
library("ComplexHeatmap")

if (!require(ggtext)){
  install.packages("ggtext")
  library(ggtext)
}

if (!require(magrittr)){
  install.packages("magrittr")
  library(magrittr)
}

# Paquete de trabajo para la carga de datos en excel
if (!require(readxl)){
  install.packages("readxl")
  library(readxl)
}

# Paquete de trabajo para la grabación de datos en csv
if (!require(readr)){
  install.packages("readr")
  library(readr)
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
cols_info <- DataFrame(samples_data[,c(4,3)])

# Refactorización clase dentro las columnas.
cols_info$Class <- factor(cols_info$Class, levels = c("QC", "GC", "BN", "HE"),
                          labels = c("Quality Control", "Gastric Cancer", "Benign", "Healthy"))

# Construye el objeto SummarizedExperiment
se <- SummarizedExperiment(assays=samples_matrix, colData=cols_info, rowData = rows_info, metadata =meta_info)
se

#Almacenar en un objeto binario el objeto summarized experiment
save(se, file = "data/GastricCancer_SE.rda")

# Carga el objeto se
load(file = "data/GastricCancer_SE.rda")

# Grabar los datos en un fichero
write_csv(cbind(ID=rownames(se),as.data.frame(assay(se))), "data/GC_metabolite.csv", append = FALSE, col_names = TRUE)
write_csv(cbind(ID = colnames(se),as.data.frame(colData(se))), "data/GC_sampleMetada.csv", append = FALSE, col_names = TRUE)
write_csv(cbind(ID =rownames(se),as.data.frame(rowData(se))), "data/GC_metaboliteMetaData.csv", append = FALSE, col_names = TRUE)


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

# Análisis univariante. cuartiles, medias y medianas de todas las variables

resum <- data.frame(ID=NA, Min=NA, Q1=NA, Mediana=NA, Media=NA, Q3=NA, Max=NA,
                    pValue=NA)
grupo <- colData(se_preprocessed)$Class

for (i in 1:length(rownames(se_preprocessed))){
  sum <- round(summary(assays(se_preprocessed)[[1]][i,]),3)
  resum[i,1:7] <- c(rownames(se_preprocessed)[i], sum[[1]], sum[[2]], sum[[3]], sum[[4]],
                 sum[[5]],sum[[6]])

  # ANOVA
  # Hipotesis NULA. No existen diferencias en esta variable a la hora de generar los grupos (alfa = 0.05)

  valor <- assays(se_preprocessed)[[1]][i,]
  df <- data.frame(grupo,valor)
  modelo_anova <- aov(valor ~ grupo, data = df)
  p_value <- summary(modelo_anova)[[1]][["Pr(>F)"]][1]

  resum$pValue[i] <- round(p_value,5)
  
}

# Posibles variables no significativas
resum$ID[resum$pValue>=0.05]

# Posibles Variables significativas
resum$ID[resum$pValue<0.05]

# Normalización y gráficas
se_normalized <- PomaNorm(se_preprocessed,method="auto_scaling")
rowData(se_normalized) <- rowData(se)[rownames(se_preprocessed),]
se_normalized

# Matriz de correlaciones
cor_matrix <-cor(t(assay(se_normalized)))
show(cor_matrix)

# Umbral de correlaciones
threshold <- 0.8
# Seleccionar las correlaciones absolutas mayores que el umbral (excluyendo la diagonal)
high_corr <- which(abs(cor_matrix) > threshold & abs(cor_matrix) < 1, arr.ind = TRUE)

high_corr_list <- data.frame(
  Var1 = rownames(cor_matrix)[high_corr[, 1]],
  Var2 = colnames(cor_matrix)[high_corr[, 2]],
  Correlation = cor_matrix[high_corr]
)

high_corr_list

# Diagonalización de la matriz de covarianzas
EIG <- eigen(cor_matrix)
show(EIG)

eigenVecs1 <- EIG$vectors
PCAS1 <- t(assay(se_normalized)) %*% eigenVecs1
head(PCAS1)

# Significación de los diferentes componentes
vars1<- EIG$values/sum(EIG$values)
round(vars1,3)
xlabel <- paste("PCA1 ", round(vars1[1]*100, 2),"%" )
ylabel <- paste("PCA2 ", round(vars1[2]*100,2),"%" )

# Explicar hasta que nivel de componentes principales se podría coger y justificar si son muchos o no.
bgSurv<- colSurv <- factor(colData(se)$Class, levels = c("Quality Control", "Gastric Cancer", "Benign", "Healthy"),
                           labels = c("skyblue", "salmon", "lightgreen", "orange"))
pchSurv <- factor(colData(se)$Class, levels = c("Quality Control", "Gastric Cancer", "Benign", "Healthy"),
                  labels = c(1, 2, 3, 4))

plot(PCAS1[,1], PCAS1[,2], main = "Muestras. 2 primeras PCs",
     xlab=xlabel, ylab=ylabel, 
     col=colSurv, bg=bgSurv,pch=as.numeric(pchSurv))
legend("bottomright", legend = unique(unique(colData(se_normalized)$Class)), 
       fill = c("black", "salmon", "lightgreen", "skyblue"), 
       title = "Clases") 

# Mapa de calor
manDist <- dist(t(assay(se_normalized)))
heatmap (as.matrix(manDist), col=heat.colors(16))

# Mapa de calor.Correlación entre variables
PomaHeatmap(se_normalized,
            covs="Class",
            sample_names = FALSE,
            feature_names = TRUE,
            show_legend = TRUE)

# Boxplot de las muestras ordenadas por clases - Efecto batch
class_order <- order(colData(se_normalized)$Class)
ordered_data_matrix <- assay(se_normalized)[,class_order]

colors <- factor(colData(se)$Class[class_order], levels = c("Quality Control", "Gastric Cancer", "Benign", "Healthy"),
                 labels = c("skyblue", "salmon", "lightgreen", "orange"))

boxplot(ordered_data_matrix,
        col = colors,                       # Colores según el grupo
        las = 2,                            # Rotación de etiquetas en el eje X
        main = "Boxplot de Muestras por Clase",
        xlab = "Muestras",
        ylab = "Metabolitos x muestra",
        ylim = c(-1,6),
        names = colnames(se_normalized)[class_order])  # Nombres de las muestras en el orden de las clases

# Añadir una leyenda
legend("topright", legend = unique(colData(se_normalized)$Class), 
       fill = c("black", "salmon", "lightgreen", "skyblue"), 
       title = "Clases")
PomaBoxplots(se_normalized,x = "samples",outcome="Class")
PomaBoxplots(se_normalized,x = "features", outcome=NULL,
             theme_params = list(legend_title = FALSE, axis_x_rotate = TRUE))



# Variables significativas según el p-Valor
selectedRows <- rowData(se_preprocessed)$Name %in% resum$ID[resum$pValue<0.05]
se_selected <- se_preprocessed[selectedRows,]

# Normalización y gráficas
se_selectedNorm <- PomaNorm(se_selected,method="auto_scaling")
rowData(se_selectedNorm) <- rowData(se)[rownames(se_selected),]
se_selectedNorm

# Matriz de correlaciones
cor_matrix <-cor(t(assay(se_selectedNorm)))
show(cor_matrix)

# Umbral de correlaciones
threshold <- 0.8
# Seleccionar las correlaciones absolutas mayores que el umbral (excluyendo la diagonal)
high_corr <- which(abs(cor_matrix) > threshold & abs(cor_matrix) < 1, arr.ind = TRUE)

high_corr_list <- data.frame(
  Var1 = rownames(cor_matrix)[high_corr[, 1]],
  Var2 = colnames(cor_matrix)[high_corr[, 2]],
  Correlation = cor_matrix[high_corr]
)

high_corr_list

# Diagonalización de la matriz de covarianzas
EIG <- eigen(cor_matrix)
show(EIG)

eigenVecs1 <- EIG$vectors
PCAS1 <- t(assay(se_selectedNorm)) %*% eigenVecs1
head(PCAS1)

# Significación de los diferentes componentes
vars1<- EIG$values/sum(EIG$values)
round(vars1,3)
xlabel <- paste("PCA1 ", round(vars1[1]*100, 2),"%" )
ylabel <- paste("PCA2 ", round(vars1[2]*100,2),"%" )

# Explicar hasta que nivel de componentes principales se podría coger y justificar si son muchos o no.
bgSurv<- colSurv <- factor(colData(se_selectedNorm)$Class, levels = c("Quality Control", "Gastric Cancer", "Benign", "Healthy"),
                           labels = c("skyblue", "salmon", "lightgreen", "orange"))
pchSurv <- factor(colData(se_selectedNorm)$Class, levels = c("Quality Control", "Gastric Cancer", "Benign", "Healthy"),
                  labels = c(1, 2, 3, 4))

plot(PCAS1[,1], PCAS1[,2], main = "Muestras. 2 primeras PCs",
     xlab=xlabel, ylab=ylabel, 
     col=colSurv, bg=bgSurv,pch=as.numeric(pchSurv))
legend("bottomright", legend = unique(unique(colData(se_selectedNorm)$Class)), 
       fill = c("black", "salmon", "lightgreen", "skyblue"), 
       title = "Clases") 

# Mapa de calor
manDist <- dist(t(assay(se_selectedNorm)))
heatmap (as.matrix(manDist), col=heat.colors(16))

# Mapa de calor.Correlación entre variables
PomaHeatmap(se_selectedNorm,
            covs="Class",
            sample_names = TRUE,
            feature_names = TRUE,
            show_legend = TRUE)

# Boxplot de las muestras ordenadas por clases - Efecto batch
class_order <- order(colData(se_selectedNorm)$Class)
ordered_data_matrix <- assay(se_selectedNorm)[,class_order]

colors <- factor(colData(se_selectedNorm)$Class[class_order], levels = c("Quality Control", "Gastric Cancer", "Benign", "Healthy"),
                 labels = c("skyblue", "salmon", "lightgreen", "orange"))

boxplot(ordered_data_matrix,
        col = colors,                       # Colores según el grupo
        las = 2,                            # Rotación de etiquetas en el eje X
        main = "Boxplot de Muestras por Clase",
        xlab = "Muestras",
        ylab = "Metabolitos x muestra",
        ylim = c(-1,6),
        names = colnames(se_selectedNorm)[class_order])  # Nombres de las muestras en el orden de las clases

# Añadir una leyenda
legend("topright", legend = unique(colData(se_selectedNorm)$Class), 
       fill = c("black", "salmon", "lightgreen", "skyblue"), 
       title = "Clases")
PomaBoxplots(se_selectedNorm,x = "samples",outcome="Class")
PomaBoxplots(se_selectedNorm,x = "features", outcome=NULL,
             theme_params = list(legend_title = FALSE, axis_x_rotate = TRUE))


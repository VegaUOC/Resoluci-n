##### Carga de librerias

# Carga de librerías necesarias para la realización de la PEC

# Librería Bioconductor
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("SummarizedExperiment")
library("SummarizedExperiment")

BiocManager::install("POMA")
library("POMA")

if (!require("ComplexHeatmap"))
  BiocManager::install("ComplexHeatmap")
library("ComplexHeatmap")

if (!require(ggtext)){
  install.packages("ggtext")
  library(ggtext)
}

if (!require(magrittr)){
  install.packages("magrittr")
  library(magrittr)
}

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

#https://aspteaching.github.io/AMVCasos/#ejemplo-pca-1-b%C3%BAsqueda-de-factores-latentes-en-datos-ecol%C3%B3gicos1

# Matriz de covarianzas 
n <- dim(t(assay(se_normalized)))[2]
S<-cov(t(assay(se_normalized)))*(n-1)/n
show(S)

# Matriz de correlaciones
R<-cor(t(assay(se_normalized)))
show(R)

# Diagonalización de la matriz de covarianzas
EIG <- eigen(S)
show(EIG)

eigenVecs1 <- EIG$vectors
PCAS1 <- t(assay(se_normalized)) %*% eigenVecs1
head(PCAS1)

plot(PCAS1[,1], PCAS1[,2], main = "Muestras. 2 primeras PCs")

vars1<- EIG$values/sum(EIG$values)
round(vars1,3)

xlabel <- paste("PCA1 ", round(vars1[1]*100, 2),"%" )
ylabel <- paste("PCA2 ", round(vars1[2]*100,2),"%" )
plot(PCAS1[,1], PCAS1[,2], main = "Muestras. 2 primeras PCs",
     xlab=xlabel, ylab=ylabel)

bgSurv<- colSurv <- factor(colData(se)$Class, levels = c("Quality Control", "Gastric Cancer", "Benign", "Healthy"),
                           labels = c("skyblue", "salmon", "lightgreen", "orange"))
pchSurv <- factor(colData(se)$Class, levels = c("Quality Control", "Gastric Cancer", "Benign", "Healthy"),
                  labels = c(1, 2, 3, 4))

plot(PCAS1[,1], PCAS1[,2], main = "Muestras. 2 primeras PCs",
     xlab=xlabel, ylab=ylabel, 
     col=colSurv, bg=bgSurv,pch=as.numeric(pchSurv))
legend("bottomright", legend = unique(sample_classes), 
       fill = c("black", "salmon", "lightgreen", "skyblue"), 
       title = "Clases") 

PCAS2 <- princomp(t(assay(se_normalized)))
names(PCAS2)
PCAS2

PCAS3 <- prcomp(t(assay(se_normalized)))
names(PCAS3)
PCAS3

PCAS2$sdev

PCAS2$loadings

PCAS3$rotation[,1]

plotPCA <- function (X, labels=NULL, colors=NULL, dataDesc="", scale=FALSE, cex4text=0.8)
  {
  pcX<-prcomp(X, scale=scale) # o prcomp(t(X))
  loads<- round(pcX$sdev^2/sum(pcX$sdev^2)*100,1)
  xlab<-c(paste("PC1",loads[1],"%"))
  ylab<-c(paste("PC2",loads[2],"%"))
  if (is.null(colors)) colors=1
  plot(pcX$x[,1:2],xlab=xlab,ylab=ylab, col=colors,  
       xlim=c(min(pcX$x[,1])-10, max(pcX$x[,1])+10),
       ylim=c(min(pcX$x[,2])-10, max(pcX$x[,2])+10))
  #text(pcX$x[,1],pcX$x[,2], labels, pos=3, cex=cex4text)
  title(paste("Plot of first 2 PCs for expressions in", dataDesc, sep=" "), cex=0.8)
}

sampleNames <- colnames(se_normalized)
plotPCA(t(assay(se_normalized)), labels=sampleNames, colors=colSurv, dataDesc="selected samples", cex4text=0.6)

# Mapo de calor
manDist <- dist(t(assay(se_normalized)))
heatmap (as.matrix(manDist), col=heat.colors(16))

require(MASS)
sam1<-sammon (manDist, trace=FALSE)
plot(sam1$points)
text(sam1$points, colData(se_normalized)$Class, pos=4)

# Mapa de calor.Correlación entre variables
PomaHeatmap(se_normalized,
            sample_names = TRUE,
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
        names = colnames(se_normalized)[class_order])  # Nombres de las muestras en el orden de las clases

# Añadir una leyenda
legend("topright", legend = unique(sample_classes), 
       fill = c("black", "salmon", "lightgreen", "skyblue"), 
       title = "Clases")



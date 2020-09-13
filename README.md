# RNA-seq data analysis for Klebsormidium nitens and gene co-expression network construction.

Este repositorio contiene todos los datos analizados y los scripts necesarios correspondientes a mi trabajo de fin de máster: A gene co-expression network reveals coordinated rhythmic gene expression patterns in the charophyte Klebsormidium nitens.

En él se han reanalizado datos de RNA-seq de Ferrari, C., Proost, S., Janowski, M. et al. (2019). Kingdom-wide comparison reveals the evolution of diurnal gene expression in Archaeplastida. Nat Commun 10, 737. con herramientas up-to date y con la creación de una red de coexpresión génica.

## Instrucciones:

El flujo de trabajo es el siguiente:

<p align="center">
  <img width="460" height="300" src="https://github.com/marcos-bioinformatics/RNA_seq_K_nitens/blob/master/Workflow.png">
</p>

**Las muestras se encuentran procesadas**, pero los scripts en Bash para el procesamiento de las mismas se encuentran disponibles para su consulta.

El análisis de los datos comprende tres scripts:

- PCA_Rain_Clustering
- Ontology Analysis
- Network

Basta con ejecutarlos desde la carpeta DATA donde se encuentran las muestras y otros ficheros de lectura necesarios. Para más información se puede consultar el trabajo completo en pdf o contactar con el autor.

## Examen parcial. Adrián Villa.


#librerias que se usarán

library (Biostrings)
library (parallel)
library (BiocGenerics)
library (S4Vectors)
library (stats4)
library (IRanges)
library (XVector)

#Para los árboles

library (msa)
library (seqinr)
library (ape)

insulinas <- readAAStringSet("Insulinas.fasta") ## Leer el archivo de insulinas

names(insulinas) <- c("Insulin Homo Sapiens", "Insulin Pan troglodytes", "Insulin Rattus norvegicus", "Insulin Canis lupus familiaris", "Insulin Mus musculus", "Insulin Oryctolagus cuniculus", "Insulin Sus scrofa", "Insulin Gallus gallus", "Insulin Bos tarus", "Insulin B & A INS_SHEEP", "Insulin Loxodonta africana", "Insulin B & A INS_PHYMC", "Insulin Cavia porcellus", "Insulin B & A INS_MYOCO")
insulinas

insu <- msa (insulinas, "Muscle") ## Alinear con algoritmo Muscle
insu

insu1 <- msa (insulinas, "ClustalW") ## Alinear con algoritmo ClustalW
insu1

insu <- msaConvert (insu, type ="seqinr::alignment") ## Convertir los datos para poder hacer el árbol

insu1 <- msaConvert (insu1, type ="seqinr::alignment") ## Convertir los datos para poder hacer el árbol

dd <- dist.alignment (insu, "identity") ## Calcular la matriz de distancias del alineamiento

dd1 <- dist.alignment (insu1, "identity") ## Calcular la matriz de distancias del alineamiento

as.matrix (dd) ## Tomar el objeto dd como matriz
as.matrix (dd1) ## Tomar el objeto dd como matriz

insuTree <- nj (dd) ## Árbol basado en Neighbour joining
plot(insuTree, main = "Insulinas") ## Plot del árbol

insuTree1 <- nj (dd1) ## Árbol basado en Neighbour joining
plot(insuTree1, main = "Insulinas") ## Plot del árbol 


## Librerías para ggtree

library (ggtree)
library (EBImage)
library (treeio)
library (ggplotify)

arbol <- ggtree (insuTree, color = "cadetblue", branch.length = "none", size = 2) + theme_tree ("ghostwhite") +
    geom_tiplab (as_ylab = T, size = 7, color = "darkslategrey") ## Árbol con algoritmo Muscle

ggtree (insuTree1, color = "cadetblue", branch.length = "none", size = 2) + theme_tree ("ghostwhite") +
  geom_tiplab (as_ylab = T, size = 7, color = "darkslategrey") ## Árbol con algoritmo ClustalW


## Librería necesaria para hacer árboles con imágenes

install.packages ("rphylopic")
library (rphylopic)

homo_sapiens <- name_search (text = "Homo sapiens", options = "namebankID")[[1]] # find names
homo_sapiens
human_id <- name_images (uuid = homo_sapiens$uid[1])  # list images
human_id
imagen <- image_data(human_id, size = 64)[[1]]

arbol + add_phylopic (imagen, alpha = .2)


## No pude colocar las imágenes, seguí varios ejermplos, pero no quiere descargar el archivo en PNG, la verdad no sé la razón








install.packages ("Peptides") ## Librería que contiene una función que podría ayudar
library (Peptides) ## Cargamos la librería

secuencia <- AAString ("SKADYEK") ## Se indica cuál es la secuencia de AA y la guardo en un objeto
mw (secuencia, monoisotopic = TRUE) ## La función "mw" o "molecular weight" se encarga de arrojar el peso total de nuestra secuencia
## El resultado fue 839.4025, un poco distinto a lo que indica la página, pero creo que es correcto
## Me hubiera gustado usar la "MONOISOTOPIC_MASS_TABLE", pero no pude abrirla ni descargarla, posiblemente porque mi navegador detecta la página como insegura
## Librerías
library(ggpubr)
library(ggVennDiagram)
library(ggthemes)
library(digest)
library(tidyverse)
library(dplyr)
library(phyloseq)
library(ggplot2)
library(plotly)
library(vegan)
library(ampvis2)
library(microbiome)
library(cowplot)

## Funciones propias

biom2physeq <- function(biom_file, metafile, ncbitax = T){
  physeqObject <- import_biom( biom_file,
                               refseqArgs = NULL, 
                               parseFunction = parse_taxonomy_default,
                               parallel = FALSE)
  if (isTRUE(ncbitax)){
    colnames(physeqObject@tax_table) <- c("Kingdom","Subkingdom","Superphylum",
                                          "Phylum","Subphylum","Superclass",
                                          "Class","Subclass", "Superorder",
                                          "Order", "Suborder", "Superfamily",
                                          "Family", "Subfamily","Genus", "Species")
    
  } else {
    colnames(physeqObject@tax_table) <- c("Kingdom", "Phylum", "Class",
                                          "Order","Family", "Genus", "Species")
  }
  
  metadata <- read.delim(metafile, header = T, sep = "\t")
  samples <- metadata[,1]
  rownames(metadata) <- samples
  
  ord <- match(colnames(otu_table(physeqObject)), rownames(metadata))
  metadata <- metadata[ord,]
  sampledata = sample_data(metadata)
  physeqData <- merge_phyloseq(physeqObject, sampledata)
  return(physeqData)
}
biom_to_ampvis <- function(biomfile, metafile, treefile, ncbitax = T){
  physeqObject <- import_biom(biomfile,
                              treefilename = treefile,
                              refseqArgs = NULL, 
                              parseFunction = parse_taxonomy_default,
                              parallel = FALSE)
  if (isTRUE(ncbitax)){
    colnames(physeqObject@tax_table) <- c("Kingdom","Subkingdom","Superphylum",
                                          "Phylum","Subphylum","Superclass",
                                          "Class","Subclass", "Superorder",
                                          "Order", "Suborder", "Superfamily",
                                          "Family", "Subfamily","Genus", "Species")
    
  } else {
    colnames(physeqObject@tax_table) <- c("Kingdom", "Phylum", "Class",
                                          "Order","Family", "Genus", "Species")
  }
  
  metafile <- read.delim(metafile, header = T, sep = "\t")
  samples <- metafile[,1]
  rownames(metafile) <- samples
  
  ord <- match(colnames(otu_table(physeqObject)), rownames(metafile))
  metafile <- metafile[ord,]
  sampledata = sample_data(metafile)
  physeqData <- merge_phyloseq(physeqObject, sampledata)
  
  if (isTRUE(ncbitax)){
    ### Eliminamos los taxones que no interesan y convertimos a formato ampvis
    taxa <- physeqData@tax_table@.Data
    new_taxa <- taxa[,c(1,4,7,10,13,15,16)]
    
    physeqData@tax_table@.Data <- new_taxa
  }
  
  ampData_COI <- phyloseq_to_ampvis2(physeqData)
}
join_physeqs <- function(physeq_COI, physeq_18S){
  # Nota posterior: esta función está hecha porque los datos de hash de las 
  # OTUs no funcionan bien si son iguales al hacer merge_phyloseq
    
  ## Unión de objetos phyloseq ##
  
  # Revisar si existen OTUs repetidas (hash), y cambiarlo en uno de los objetos
  repotu_18S <- which(rownames(physeq_18S@tax_table) %in% 
                        rownames(physeq_COI@tax_table))
  rownames(physeq_18S@tax_table)[repotu_18S] <- paste0(rownames(physeq_18S@tax_table)[repotu_18S], "-2")
  rownames(physeq_18S@otu_table)[repotu_18S] <- paste0(rownames(physeq_18S@otu_table)[repotu_18S], "-2")
  
  physeqData <- merge_phyloseq(physeq_18S, physeq_COI)
}
visu2physeq <- function(visu_file, metavisu_file){ # aviso que esto solo va a 
  # aviso que esto solo va a funcionar con los metadatos actuales, ya que 
  # está hardcodeado a las variables de ahora.
  df <- read_delim(visu_file)
  rownames(df) <- apply(df, 1, function(x) digest(x, "md5", serialize = TRUE))
  tax_df <- df[,1:7, drop = F]
  otu_df <- df[,8:ncol(df), drop = F]
  rownames(tax_df) <- rownames(df)
  rownames(otu_df) <- rownames(df)
  otumat <- as.matrix(otu_df)
  taxmat <- as.matrix(tax_df)
  metavisu <- read.csv(metavisu_file, header = T, sep = "\t")
  otu <- otu_table(otumat, taxa_are_rows = T)
  tax <- tax_table(taxmat)
  rownames(metavisu) <- metavisu$sample.id
  sam <- sample_data(metavisu)
  physeq_visu <- phyloseq(otu, tax, sam)
  return(physeq_visu)
}
join_physeq_list <- function(phyloseq_obj_list){
  for (i in seq(phy_list)){
    if (i == 1){
      new_phy <- phy_list[[i]]
    } else {
      new_phy <- join_physeqs(new_phy, phy_list[[i]])
    }
  }
  return(new_phy)
}
amp2phy <- function(amp_obj){
  phy_obj <- phyloseq(otu_table(as.matrix(amp_obj$abund), taxa_are_rows = T),
  tax_table(as.matrix(amp_obj$tax)),
  sample_data(amp_obj$metadata))
  return(phy_obj)
}


## Directorio de trabajo

setwd("/home/rosa/TFM/analisis")

## Agrupación de resultados

### Tanda 1
file_18S_1 <- "/home/rosa/TFM/Datos/Tanda1_180523_18S/biom/table-with-tax-json.biom"
file_COI_1 <- "/home/rosa/TFM/Datos/Tanda1_240523_COI/biom/table-with-tax-json.biom"
tree_18S_1 <- "/home/rosa/TFM/Datos/Tanda1_180523_18S/qiime/unrooted-tree.nwk"
tree_COI_1 <- "/home/rosa/TFM/Datos/Tanda1_240523_COI/qiime/unrooted-tree.nwk"
metadata_file_18S_1 <- "/home/rosa/TFM/Datos/Tanda1_180523_18S/metadata.tsv"
metadata_file_COI_1 <- "/home/rosa/TFM/Datos/Tanda1_240523_COI/metadata.tsv"

physeq_18S_1 <- biom2physeq(file_18S_1,
                          metadata_file_18S_1,
                          ncbitax = T)
physeq_COI_1 <- biom2physeq(file_COI_1, 
                          metadata_file_COI_1,
                          ncbitax = T)
physeq_COI_1@tax_table@.Data["d8a23100187db1681a38dcce6c3536e7", "Family"] = "" # Errata en la taxonomía
# Se define la lista donde se almacenarán todos los objetos phyloseq
phy_list <- list(physeq_18S_1, physeq_COI_1)
phy_list_mb <- list(physeq_18S_1, physeq_COI_1)

### Tanda 2
file_18S_2 <- "/home/rosa/TFM/Datos/Tanda2_230823_18S/biom/table-with-tax-json.biom"
file_COI_2 <- "/home/rosa/TFM/Datos/Tanda2_230904_COI/biom/table-with-tax-json.biom"
tree_18S_2 <- "/home/rosa/TFM/Datos/Tanda2_230823_18S/qiime/unrooted-tree.nwk"
tree_COI_2 <- "/home/rosa/TFM/Datos/Tanda2_230904_COI/qiime/unrooted-tree.nwk"
metadata_file_18S_2 <- "/home/rosa/TFM/Datos/Tanda2_230823_18S/metadata.tsv"
metadata_file_COI_2 <- "/home/rosa/TFM/Datos/Tanda2_230904_COI/metadata.tsv"

physeq_18S_2 <- biom2physeq(file_18S_2,
                            metadata_file_18S_2,
                            ncbitax = T)
physeq_COI_2 <- biom2physeq(file_COI_2, 
                            metadata_file_COI_2,
                            ncbitax = T)
# Se añaden a la lista los nuevos objetos phyloseq
phy_list <- append(phy_list, physeq_18S_2)
phy_list <- append(phy_list, physeq_COI_2)
phy_list_mb <- append(phy_list_mb, physeq_18S_2)
phy_list_mb <- append(phy_list_mb, physeq_COI_2)

### Tanda 3
file_18S_3 <- "/home/rosa/TFM/Datos/Tanda3_230912_18S/biom/table-with-tax-json.biom"
file_COI_3 <- "/home/rosa/TFM/Datos/Tanda3_230908_COI/biom/table-with-tax-json.biom"
tree_18S_3 <- "/home/rosa/TFM/Datos/Tanda3_230912_18S/qiime/unrooted-tree.nwk"
tree_COI_3 <- "/home/rosa/TFM/Datos/Tanda3_230908_COI/qiime/unrooted-tree.nwk"
metadata_file_18S_3 <- "/home/rosa/TFM/Datos/Tanda3_230912_18S/metadata.tsv"
metadata_file_COI_3 <- "/home/rosa/TFM/Datos/Tanda3_230908_COI/metadata.tsv"

physeq_18S_3 <- biom2physeq(file_18S_3,
                            metadata_file_18S_3,
                            ncbitax = T)
physeq_COI_3 <- biom2physeq(file_COI_3, 
                            metadata_file_COI_3,
                            ncbitax = T)
phy_list <- append(phy_list, physeq_COI_3)
phy_list <- append(phy_list, physeq_18S_3)
phy_list_mb <- append(phy_list_mb, physeq_COI_3)
phy_list_mb <- append(phy_list_mb, physeq_18S_3)

### Resultados de Visu

wd_c <- "/home/rosa/TFM/Datos/visu/Continental/"
wd_m <- "/home/rosa/TFM/Datos/visu/Marino/"
physeq_visu_c <-
  visu2physeq(paste0(wd_c, "taxonomia_continental.csv"), 
              paste0(wd_c, "metadata.tsv"))
physeq_visu_m <-
  visu2physeq(paste0(wd_m, "taxonomia_marino_2.csv"), 
              paste0(wd_m, "metadata.tsv"))

phy_list <- append(phy_list, physeq_visu_c)
phy_list <- append(phy_list, physeq_visu_m)

### Agrupamiento de todos los resultados
physeq_all <- join_physeq_list(phy_list)
physeq_all

### Cambio de formato a ampvis2

# Lo primero va a ser revisar las muestras y OTUs que tengan 0 lecturas
# Lo segundo lo hace automáticamente el ampvis al pasarlo, pero lo primero no
ps <- names(which(colSums(physeq_all@otu_table) > 0))
physeq_all <- prune_samples(ps, physeq_all)
save(physeq_all, file = "physeq_all.RData")
taxa <- physeq_all@tax_table@.Data
# ampvis2 no admite más que los siguientes taxones:
new_taxa <- taxa[,c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")]
physeq_all@tax_table@.Data <- new_taxa
ampData <- amp_load(physeq_all) # Creo que al final lo haré con phyloseq
ampData$metadata$muestramarcador <- paste(
  ampData$metadata$muestra,
  ampData$metadata$marcador,
  sep = "-"
)
save(ampData, file = "ampvis_all.RData")

# Diversidad alfa
amp_c <- amp_filter_samples(ampData, medio == "Continental")
amp_cmb <- amp_filter_samples(amp_c, metodo == "Metabarcoding")
amp_cmbm <- amp_merge_replicates(amp_cmb, "muestramarcador", round = "up")
amp_cmbm <- filter_otus(amp_cmbm, 0.1) # Filtrado para aumentar la fiabilidad
# de la muestra (https://www.frontiersin.org/articles/10.3389/fcimb.2023.1165295/full)
meths <- c("Observed", "Shannon", 
           "Simpson", "Fisher")
alfa_c <- cbind.data.frame(
  amp_cmbm$metadata,
  estimate_richness(amp2phy(amp_cmbm), measures = meths))
alfa_c_p <- alfa_c %>% 
  pivot_longer(!names(amp_cmbm$metadata), names_to = "index")
alfa_c_p %>%
  mutate(index = factor(alfa_c_p$index, levels = unique(alfa_c_p$index))) %>%
  arrange(marcador, muestra) %>%
  ggplot(aes(x = marcador, y = value, color = marcador)) +
  geom_boxplot() +
  geom_boxplot(aes(fill = marcador), alpha = .2) + 
  facet_wrap(index~., scales = "free_y", nrow = 1) +
  stat_compare_means(label.y.npc = "bottom", paired = T) +
  scale_color_brewer(palette = "Set2") + 
  scale_fill_brewer(palette = "Set2") +
  xlab("Marcador") + 
  labs(color = "Marcador", fill = "Marcador") +
  ylab("Valor") +
  theme_bw()

amp_c <- amp_filter_samples(ampData, medio == "Continental")
amp_cm <- amp_merge_replicates(amp_c, "muestramarcador", round = "up")
# ampvis mezcla muy bien replicados pero para los
# filtros taxonómicos es mejor phyloseq y el agregado de taxones se hace con
# microbiome que utiliza solo objetos phyloseq
                        
seed = 1234

phy_cm <- amp2phy(amp_cm)
phy_cmf <- subset_taxa(phy_cm, Family != "")
phy_cmfa <- aggregate_taxa(phy_cmf, level = "Family")
phy_cmfa_r <- rarefy_even_depth(phy_cmfa, rngseed = seed)
# La rarefacción provoca que las muestras sean más comparables entre sí al 
# realizar muestreos al nivel que detectamos individuos en visu

alfa_ca <- cbind.data.frame(
  sample_data(phy_cmfa_r), 
  estimate_richness(phy_cmfa_r, measures = meths))

alfa_cap <- alfa_ca %>% 
  pivot_longer(!names(amp_cm$metadata), names_to = "index")

comp <- list(
  c("18S", "COI"),
  c("COI", "Visu"),
  c("18S", "Visu")
)

alfa_cap %>%
  mutate(index = factor(alfa_cap$index, levels = unique(alfa_cap$index))) %>%
  arrange(marcador, muestra) %>%
  ggplot(aes(x = marcador, y = value, color = marcador)) +
  geom_boxplot() +
  geom_boxplot(aes(fill = marcador), alpha = .2) + 
  facet_wrap(index~., scales = "free_y", nrow = 1) +
  stat_compare_means(label.y.npc = "bottom",
                     method.args = list(paired = T)) + 
  stat_compare_means(comparisons = comp,
                     tip.length = .01,
                     method.args = list(paired = T)) +
  scale_color_brewer(palette = "Set2") + 
  scale_fill_brewer(palette = "Set2") +
  xlab("Metodología") + 
  ylab("Valor") +
  labs(color = "Marcador", fill = "Marcador") +
  theme_bw()


# Juntamos los datos de tax y otus
amp_cmfa <- amp_load(phy_cmfa)
amp_cmfa <- filter_otus(amp_cmfa, 0.1) # eliminamos las otus menores al .1%
data_venn <- merge(amp_cmfa$abund, amp_cmfa$tax, by = 0)
# data_venn <- merge(phy_cmfa@otu_table@.Data, phy_cmfa@tax_table@.Data, by = 0)
# Separamos en objetos cada una de las taxonomías asociadas a cada marcador
coi <- data_venn[,c(colnames(data_venn)[grep("COI", colnames(data_venn))], "Family")]
visu <- data_venn[,c(colnames(data_venn)[grep("Visu", colnames(data_venn))], "Family")]
s18 <- data_venn[,c(colnames(data_venn)[grep("18S", colnames(data_venn))], "Family")]
# Extraemos las familias que hayan sido detectadas al menos en una muestra una vez
coi_f <- unique(coi[rowSums(coi[,1:length(coi) -1]) > 0,"Family"])
visu_f <- unique(visu[rowSums(visu[,1:length(visu) -1]) > 0,"Family"])
s18_f <- unique(s18[rowSums(s18[,1:length(s18) -1]) > 0,"Family"])
# Generamos la lista para el Venn
x = list(
  "COI" = coi_f,
  "18S" = s18_f,
  "Visu" = visu_f
  )

ggVennDiagram(
  x,
  label_alpha = 0
) +
  scale_fill_gradient("nº familias",low = "white", high = "red") + 
  scale_color_manual(values = c("black", "black", "black")) +
  # labs(title = "Medio continental") +
  theme(plot.title = element_text(hjust = 0.5, size = 20),
        legend.position = "none")

amp_cmfa_r <- amp_load(phy_cmfa_r)
# amp_cmfa <- amp_load(phy_cmfa)
# Creamos los los dos objetos pero para este caso, como 
# estamos viendo cualitativamente las familias, no es necesario el dataset
# rarefactado. Eso solo para cuando se vayan a sacar significaciones estadísticas.

p <- amp_boxplot(amp_cmfa,
            tax_aggregate = "Family",
            normalise = T,
            group_by = "marcador",
            tax_show = 10, 
            plot_log = F) +
  geom_boxplot(aes(fill = .Group), alpha = .2)+
  scale_color_brewer(palette = "Set2") + 
  scale_fill_brewer(palette = "Set2", guide = "none") +
  # stat_compare_means(label.y.npc = "bottom") +
  labs(color = "Marcador",
       y = "Abundancia relativa (%)")

p

# Diversidad beta

amp_cmfa_r$metadata$Puntomarcador = 
  factor(paste(
    amp_cmfa_r$metadata$Punto,
    amp_cmfa_r$metadata$marcador,
    sep = "-"
  ), levels = c(
    "Antes depuradora-Visu",
    "Después depuradora-Visu",
    "Antes depuradora-18S",
    "Después depuradora-18S",
    "Antes depuradora-COI",
    "Después depuradora-COI"
  ))
amp_cmfa$metadata$Puntomarcador = 
  factor(paste(
    amp_cmfa$metadata$Punto,
    amp_cmfa$metadata$marcador,
    sep = "-"
  ), levels = c(
    "Antes depuradora-Visu",
    "Después depuradora-Visu",
    "Antes depuradora-18S",
    "Después depuradora-18S",
    "Antes depuradora-COI",
    "Después depuradora-COI"
  ))
p1 <- amp_ordinate(amp_cmfa_r, # stress: 0.127
             type = "nmds",
             filter_species = 0, 
             distmeasure = "bray",
             transform = "none",
             sample_color_by = "Puntomarcador",
             sample_colorframe = T,
             sample_shape_by = "marcador"
             ) +
  labs(color = "Punto - Marcador",
       fill = "Punto - Marcador") +
  scale_color_tableau(palette = "Tableau 20") +
  scale_fill_tableau(palette = "Tableau 20")
p2 <- amp_ordinate(amp_cmfa_r,
             type = "rda",
             filter_species = 0, 
             distmeasure = "bray",
             constrain = c("Punto"),
             transform = "hellinger",
             sample_color_by = "Puntomarcador",
             sample_colorframe = T,
             sample_shape_by = "marcador") +
  scale_color_tableau(palette = "Tableau 20") +
  scale_fill_tableau(palette = "Tableau 20")
leg <- get_legend(p1)

grph <- plot_grid(
  p1 + theme(legend.position="none"),
  p2 + theme(legend.position="none"),
  align = 'vh',
  labels = "AUTO",
  hjust = -1,
  nrow = 1
)

plot_grid(
  grph, leg,
  rel_widths = c(2,.5)
)

# PERMANOVAS
## Análisis estadístico.
# Punto
bray <- distance(phy_cmfa_r, method = "bray")
df <- as(sample_data(phy_cmfa_r), "data.frame") # por algun motivo no funciona con as.data.frame
res1 <- adonis2(bray ~ Punto, data = df) # 0.196
res2 <- adonis2(bray ~ mes, data = df) # 0.188
res3 <- adonis2(bray ~ marcador, data = df) # <0.001***
res4 <- adonis2(bray ~ metodo, data = df) # <0.001***


### separando
phy_cmfa_mb <- subset_samples(phy_cmfa, metodo == "Metabarcoding")
bray_mb <- distance(phy_cmfa_mb, method = "bray")
df_mb <- as(sample_data(phy_cmfa_mb), "data.frame")
phy_cmfa_visu <- subset_samples(phy_cmfa, metodo == "Visualización")
bray_visu <- distance(phy_cmfa_visu, method = "bray")
df_visu <- as(sample_data(phy_cmfa_visu), "data.frame")

res5 <- adonis2(bray_mb ~ Punto, data = df_mb) # 0.02*
res6 <- adonis2(bray_visu ~ Punto, data = df_visu) # 0.265 jejeje

res7 <- adonis2(bray_mb ~ marcador, data = df_mb) # 0.003** bueno, era de esperar ngl

res_list <- list(res1, res2, res3, res4, res5, res6, res7)

for (i in seq(res_list)){
  res_list[[i]] <- rownames_to_column(as.data.frame(res_list[[i]]), " ")
}
res_final <- list_rbind(res_list) %>% 
  mutate("F" = round(.$F, digits = 3)) %>%
  replace(is.na(.), " ")

# Dejo una función de KableExtra por si luego la saco en PDF o HTML
# para el TFM.
kbl(res_final, caption = "Resultados PERMANOVA",
    booktabs = T, align = "c", digits = 3) %>%
  kable_classic(full_width = F) %>%
  kable_styling(latex_options = c("striped", "HOLD_position", "repeat_header")) %>%
  pack_rows("Punto (T)", 1,3) %>%
  pack_rows("Mes (T)", 4, 6) %>%
  pack_rows("Marcador (T)", 7,9) %>%
  pack_rows("Método (T)", 10, 12) %>%
  pack_rows("Punto (MB)", 13, 15) %>%
  pack_rows("Punto (V)", 16, 18) %>%
  pack_rows("Marcador (MB)", 19, 21) %>%
  column_spec(6, bold = TRUE) %>%
  footnote(general = "Abreviaciones. ",
           number = c("T: Total de muestras;",
                      "MB: Metabarcoding;",
                      "V: Visualización taxonómica;"))

# Índice de calidad biótico: IBMWP

# Al unir todos los phyloseq en uno no salen algunos taxones necesarios para el IBMWP.
# Creamos otra lista de los phyloseq de mb únicamente para hacer esto
physeq_onlymb <- join_physeq_list(phy_list_mb)
phy_taxmb <- subset_samples(physeq_onlymb, medio == "Continental")
# Para agregar los datos de 18S y COI de las réplicas técnicas si fuese necesario
phy_taxmb@sam_data$muestrareplica <- paste0(
  phy_taxmb@sam_data$muestra,
  phy_taxmb@sam_data$replica
  )
# Referencias
## Referencia ##
IBMWP_file <- "/home/rosa/TFM/References/IBMWP_list_reffile.tsv"
IBMWP_ref <- read.delim(IBMWP_file, header = T, sep = "\t")
IBMWP_rank <- as.vector(unique(IBMWP_ref$Rank))


physeqData <- phy_taxmb
selected_samples <- physeqData@sam_data$sample.id


tax_table <- as.data.frame(physeqData@tax_table@.Data)
tax_red <-tax_table[IBMWP_rank]
otu_table <- as.data.frame(physeqData@otu_table@.Data)

com_table <- merge(otu_table, tax_red, by = 0, all.x = T, all.y = F)
com_red <- com_table[!is.na(com_table$Superorder),]
com_red <- com_red[,-1] # Quito las OTUs


### Agrupar las muestras con sus réplicas
whole_samples <- unique(metadata_df$muestra)


for (i in seq(whole_samples)){
  print(grep(whole_samples[i], colnames(com_red)))
  com_red[, whole_samples[i]] <- apply(as.data.frame(com_red[,grep(whole_samples[i], colnames(com_red))]), 1, function(x) sum(x, na.rm = F))

}
com_red <- relocate(com_red, whole_samples)
selected_samples <- physeqData@sam_data$sample.id
nsam <- length(c(selected_samples, whole_samples))

df_res <- data.frame()
df_cant <- data.frame()



for (i in seq(nsam)){
  df <- com_red[,c(i, ( (nsam+1):length(com_red) ))]
  df <- subset(df, df[,1] > 1)
  
  min_perc = sum(df[,1])
  
  result_sample <- names(df)[1]
  result_score <- 0
  vistos <- c()
  if (nrow(df) == 0) {
    df_res <- rbind(df_res, c(result_sample, NA, NA))
    print(result_sample)
    next
  }
  for (j in seq(nrow(df))){
    for (k in seq(nrow(IBMWP_ref))){
      if (!is.na(df[j,IBMWP_ref[k,"Rank"]]) & df[j,IBMWP_ref[k,"Rank"]] == IBMWP_ref[k,"Taxa"]) {
        if (!(df[j,IBMWP_ref[k,"Rank"]] %in% vistos)){
        result_score <- result_score + IBMWP_ref[k, "Score"]
        vistos <- c(vistos, df[j,IBMWP_ref[k,"Rank"]])
        }
        df_cant <- rbind(df_cant, c(df[j,1],df[j,IBMWP_ref[k,"Rank"]]))
      }
    }
  }
  if (is.null(vistos)){
    df_res <- rbind(df_res, c(result_sample, result_score, "NA"))
  } else{
    vistos = unique(vistos)
    df_res <- rbind(df_res, c(result_sample, result_score, paste(vistos, collapse = ", ")))
    
  }
  
}

colnames(df_res) <- c("Sample", "IBMWP_Score", "Observed_Taxa")
colnames(df_cant) <- c("Count", "Observed_Taxa")

df_cant$Count <- as.numeric(df_cant$Count)
df_res$IBMWP_Score <- as.numeric(df_res$IBMWP_Score)

df_cant_2 <- df_cant %>% # este en verdad vale bien de nada ya
  group_by(Observed_Taxa) %>% 
  summarise(total_count = sum(Count)) %>%
  mutate(relative = (total_count * 100) / sum(total_count)) %>%
  arrange(desc(relative))


write.table(
  df_res,
  file = "IBMWP_AllSamples.tsv",
  quote = F,
  sep = "\t",
  row.names = F,
  col.names = T
)


# Ahora con los resultados
df_res <- read.delim("IBMWP_AllSamples.tsv", sep = "\t", header = T)
# leo los resultados del Limne de IBMWP (no me fío de hacerlo con mi script,
# por si veo que está mal, de nuevo lo siento).
visu_res <- read.delim("tabla_ibmwp_limne.tsv", sep = "\t", header = T)

## separamos los resultados según lo que se quiere ver/analizar
i <- 1
j <- 9
k <- 117
df_comb <- df_res[i:j,]
df_cmarc <- df_res[(j+1):k,]
df_crep <- df_res[(k+1):nrow(df_res),]
## Combinado de muestras
df_comb$Metodo <- "Metabarcoding"
df_comb$Sample_Method <- paste0(df_comb$Sample, "-MB")
visu_res$Metodo <- "Visualización"
visu_res$Sample_Method <- paste0(visu_res$Sample, "-Visu")

meta <- as(physeq_visu_c@sam_data[,c("muestra", "Punto")], "data.frame")
colnames(meta) <- c("Sample", "Punto")

df_ibmwp <- rbind(df_comb, visu_res) %>% inner_join(meta, by = "Sample") %>%
  arrange(Metodo, Sample)

p1 <-df_ibmwp %>%
  ggpaired(x = "Metodo", y = "IBMWP_Score",
         color = "Metodo") +
  # facet_grid(.~Punto) + 
  scale_color_brewer(palette = "Set1") +
  stat_compare_means(method = "t.test", paired = T, label.y = 30) +
  labs(x = "Metodología",
       y = "Puntuación IBMWP")
leg <- get_legend(p1)


# Media de las réplicas combinando marcadores 
df_cmarc$Sample <- paste0(str_split_fixed(df_cmarc$Sample, "C", 2)[,1], "C")
df_cmarc_g <- df_cmarc %>% 
  mutate(Sample_Method = paste0(Sample, "-MB")) %>%
  group_by(Sample, Sample_Method) %>%
  summarise(IBMWP_Score = mean(IBMWP_Score)) %>%
  mutate(Metodo = "Metabarcoding") %>%
  relocate(Metodo, .before = Sample_Method) %>%
  relocate(IBMWP_Score, .before = Metodo)
  

df_ibmwp_2 <- rbind(df_cmarc_g, 
      select(visu_res, -c(Observed_Taxa))) %>%
  inner_join(meta, by = "Sample") %>% 
  arrange(Metodo, Sample)

p2 <- df_ibmwp_2 %>%
  ggpaired(x = "Metodo", y = "IBMWP_Score",
         color = "Metodo") +
  # facet_grid(.~Punto) + 
  scale_color_brewer(palette = "Set1") +
  stat_compare_means(method = "t.test", paired = T, label.y = 30) +
  labs(x = "Metodología",
       y = "Puntuación IBMWP")

## Media de las réplicas (las 24)
df_crep$Sample <- paste0(str_split_fixed(df_crep$Sample, "C", 2)[,1], "C")
df_crep_g <- df_crep %>% 
  mutate(Sample_Method = paste0(Sample, "-MB")) %>%
  group_by(Sample, Sample_Method) %>%
  summarise(IBMWP_Score = mean(IBMWP_Score)) %>%
  mutate(Metodo = "Metabarcoding") %>%
  relocate(Metodo, .before = Sample_Method) %>%
  relocate(IBMWP_Score, .before = Metodo)

df_ibmwp_3 <- rbind(df_crep_g, 
      select(visu_res, -c(Observed_Taxa))) %>%
  inner_join(meta, by = "Sample") %>% 
  arrange(Metodo, Sample)


p3 <- df_ibmwp_3 %>%
  ggpaired(x = "Metodo", y = "IBMWP_Score",
         color = "Metodo") +
  # facet_grid(.~Punto) + 
  scale_color_brewer(palette = "Set1") +
  stat_compare_means(method = "t.test", paired = T, label.y = 30) +
  labs(x = "Metodología",
       y = "Puntuación IBMWP") 

p1 <- p1 + theme(legend.position = "none")
p2 <- p2 + theme(legend.position = "none")
p3 <- p3 + theme(legend.position = "none")

p_t <- plot_grid(
  plot_grid(p1, p2, nrow = 1, ncol = 2, labels = c("A", "B")),
  plot_grid(NULL, p3, NULL, nrow = 1, labels = c("", "C", ""),
            rel_widths = c(0.45, 1, 0.45)),
  nrow = 2
)

grph <- plot_grid(
  leg,
  p_t,
  rel_heights = c(.2, 3), 
  ncol = 1,
  nrow = 2
  
)
grph


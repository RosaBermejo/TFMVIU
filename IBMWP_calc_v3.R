library(dplyr)
library(phyloseq)
library(ggplot2)


## La implementación del IBMWP de manera programática que se utilizó para sacar los resultados del
## TFM está al final del archivo "diversity_analysis.R"
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

## Carga de datos de Metabarcoding ##

file_18S <- "/home/rosa/TFM/Datos180523_18S/biom/table-with-tax-json.biom"
file_COI <- "/home/rosa/TFM/Datos240523_COI/biom/table-with-tax-json.biom"
metadata_file_18S <- "/home/rosa/TFM/Datos180523_18S/metadata.tsv"
metadata_file_COI <- "/home/rosa/TFM/Datos240523_COI/metadata.tsv"

physeq_18S <- biom2physeq(file_18S,
                          metadata_file_18S,
                          ncbitax = T)
physeq_COI <- biom2physeq(file_COI, 
                          metadata_file_COI,
                          ncbitax = T)
## Referencia ##

IBMWP_file <- "/home/rosa/References/IBMWP_list_reffile.tsv"
IBMWP_ref <- read.delim(IBMWP_file, header = T, sep = "\t")
IBMWP_rank <- as.vector(unique(IBMWP_ref$Rank))

## Homogeneizamos las muestras para poder juntarlas ##

# COI
physeq_COI@sam_data$sample.id <- sub("-COI", "", physeq_COI@sam_data$sample.id )
rownames(physeq_COI@sam_data) <- sub("-COI", "", rownames(physeq_COI@sam_data) )
colnames(physeq_COI@otu_table) <- sub("-COI", "", colnames(physeq_COI@otu_table))
# 18S
physeq_18S@sam_data$sample.id <- sub("-18S", "", physeq_18S@sam_data$sample.id )
rownames(physeq_18S@sam_data) <- sub("-18S", "", rownames(physeq_18S@sam_data) )
colnames(physeq_18S@otu_table) <- sub("-18S", "", colnames(physeq_18S@otu_table))

## Unión de objetos phyloseq ##

# Revisar si existen OTUs repetidas (hash), y cambiarlo en uno de los objetos
repotu_18S <- which(rownames(physeq_18S@tax_table) %in% rownames(physeq_COI@tax_table))
rownames(physeq_18S@tax_table)[repotu_18S] <- paste0(rownames(physeq_18S@tax_table)[repotu_18S], "-2")
rownames(physeq_18S@otu_table)[repotu_18S] <- paste0(rownames(physeq_18S@otu_table)[repotu_18S], "-2")

physeqData <- merge_phyloseq(physeq_18S, physeq_COI)
# physeqData <- physeq_COI

# Muestras para un futuro
samnames <- physeqData@sam_data$sample.id

# Selección de muestras
selected_samples <- samnames[grep("C", samnames)]
# Lectura de metadatos global
metadata_df <- read.delim(metadata_file_18S)

physeqData <- subset_samples(physeqData, sample.id %in% selected_samples)

# Cada individuo sólo se anota una vez.  
tax_table <- as.data.frame(physeqData@tax_table@.Data)
tax_red <-tax_table[IBMWP_rank]
otu_table <- as.data.frame(physeqData@otu_table@.Data)

com_table <- merge(otu_table, tax_red, by = 0, all.x = T, all.y = F)
com_red <- com_table[!is.na(com_table$Superorder),]
com_red <- com_red[,-1] # Quito las OTUs

## Para esta vez quitamos 0s y los 

# com_red$Total <- apply(com_red[,selected_samples], 1, sum)
# com_red <- com_red[,-which(names(com_red) %in% selected_samples)]

# com_red <- com_red[,c( length(com_red), 1:length(com_red)-1)]

### Agrupar muestras con sus réplicas
whole_samples <- unique(metadata_df$muestra)

for (i in seq(whole_samples)){
  com_red[, whole_samples[i]] <- apply(com_red[,grep(whole_samples[i], colnames(com_red))], 1, sum)

}
com_red <- relocate(com_red, whole_samples)

nsam <- length(c(selected_samples, whole_samples))
# nsam <- 1
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

# Agrupar individuos para la relativización
df_cant$Count <- as.numeric(df_cant$Count)
df_res$IBMWP_Score <- as.numeric(df_res$IBMWP_Score)

df_cant_2 <- df_cant %>% 
  group_by(Observed_Taxa) %>% 
  summarise(total_count = sum(Count)) %>%
  mutate(relative = (total_count * 100) / sum(total_count)) %>%
  arrange(desc(relative))


# df_cant <- df_cant[order(df_cant$relative),]

# Valorar utilidad de este gráfico si hay más de una muestra
ggplot(data = df_cant_2, aes(x = reorder(Observed_Taxa,relative), y = relative, Fill = Observed_Taxa, group = Observed_Taxa , 
                           label = round(relative, digits = 2))) + 
  geom_bar(stat = "identity", fill = "#6DB6FF", color = "black", width = 0.75) +
  labs(x = "Observed Taxa", y = "Relative abundance (%)") +
  geom_text(hjust = -0.1) +
  coord_flip() +
  scale_y_continuous(expand = c(0,0), limits = c(0,max(df_cant_2$relative) + 2)) +
  theme_bw()
  
df_res

metadata_df$sample.id <- str_replace(metadata_df$sample.id, "-18S", "")

df_res_cat <- merge(df_res, metadata_df, by.x = "Sample", by.y = "sample.id")
colnames(df_res_cat)[3] <- "Observed_taxa"
colnames(df_res_cat)[2] <- "IBMWP_score"
df_res_cat$IBMWP_score <- as.numeric(df_res_cat$IBMWP_score)

df_res_summ <- df_res_cat %>%
  select(-c(Sample,Observed_taxa, replica, marcador)) %>%
  group_by(muestra) %>%
  summarise(
    n=n(),
    mean = mean(IBMWP_score),
    sd = sd(IBMWP_score),
  ) %>%
  mutate(se = sd/sqrt(n)) %>%
  mutate(ic = se * qt((1-0.05)/2 + .5, n-1))


g <- ggplot(df_res_summ) +
  geom_bar(aes(x = muestra, y = mean), stat = "identity", fill = "steelblue", alpha = .7) +
  geom_errorbar(aes(x = muestra, ymin = mean-sd, ymax = mean+sd), width = 0.4, color="black", alpha = .9, size = 1.5)

ggsave("IBMWP_plot.png", g, path = "~/TFM/analisis/")







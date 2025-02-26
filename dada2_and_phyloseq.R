#Daniel Castaneda Mogollon, PhD
#February 23rd, 2025
#Applied Genomics Centre KPU. This script was made with the purpose of analysing
#16s V4 data from Bovine origin

library(dada2)
library(phyloseq)
library(ggplot2)
library(metagenomeSeq)

################################################################################
#DADA 2 DATA PROCESSING, ASV CHARACTERIZATION, AND READ CLEANING ###############
################################################################################

#Name and read profiling
path = "/Users/danielcm/Desktop/KPU/q2/output_filtered"
setwd(path)
forward_files = sort(list.files(path,pattern="R1.fastq", full.names=TRUE))
reverse_files = sort(list.files(path,pattern="R2.fastq",full.names=TRUE))
plotQualityProfile(forward_files[1:4])
plotQualityProfile(reverse_files[1:4])
sample.names = sapply(strsplit(basename(forward_files),"_R"),`[`,1)
sample.names
filtered_forward = file.path(path, "filtered", paste0(sample.names, "_filtered_forward.fastq"))
filtered_reverse = file.path(path, "filtered", paste0(sample.names, "_filtered_reverse.fastq"))
filtered_forward
names(filtered_forward) = sample.names
names(filtered_reverse) = sample.names

#Filtering and trimming reads
output_trimming = filterAndTrim(forward_files, filtered_forward, reverse_files, filtered_reverse, truncLen=c(260,200),
                  maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE, 
                  multithread=TRUE, minLen = 100) #no more than 2 errors in each read, no reads with Ns are allowed, 
                  #truncating at quality of 2, truncating forward reads at 260 bp and 200 bp if reverse.

output_trimming #sanity check
error_forward = learnErrors(filtered_forward, multithread=TRUE)
error_reverse = learnErrors(filtered_reverse, multithread=TRUE)
plotErrors(error_forward,nominalQ = TRUE)
plotErrors(error_reverse,nominalQ = TRUE)

#Inferring ASV number and frequency from unique sequences
dada_forward = dada(filtered_forward, err=error_forward, multithread=TRUE)
dada_reverse = dada(filtered_reverse, err=error_reverse, multithread=TRUE)
merge = mergePairs(dada_forward, filtered_forward, dada_reverse, filtered_reverse, verbose=TRUE,minOverlap = 150)
sequence_table = makeSequenceTable(merge)
sequence_table_final = table(nchar(getSequences(sequence_table)))

#Dechimerization
seqtab.nochim = removeBimeraDenovo(sequence_table, method="pooled", multithread=TRUE)
sum(seqtab.nochim/sum(sequence_table)) #A high proportion of reads are chimeras! 44% with the pooled method!
getN = function(x) sum(getUniques(x))
track = cbind(output_trimming, sapply(dada_forward, getN), sapply(dada_reverse, getN), sapply(merge, getN), rowSums(seqtab.nochim))
track

#Assigning taxonomy
taxa = assignTaxonomy(seqtab.nochim, "/Users/danielcm/Downloads/GTDB_bac120_arc53_ssu_r220_fullTaxo.fa.gz") #Assigning species to 97% using the RDP classifier
taxa = addSpecies(taxa, "/Users/danielcm/Downloads/GTDB_bac120_arc53_ssu_r220_species.fa.gz") #Assigning species to 100%
taxa_silva = assignTaxonomy(seqtab.nochim, "/Users/danielcm/Downloads/silva_nr99_v138.2_toSpecies_trainset.fa.gz")
taxa_silva = addSpecies(taxa_silva, "/Users/danielcm/Downloads/silva_v138.2_assignSpecies.fa.gz")
#View(taxa) #Sanity check


################################################################################
#PHYLOSEQ DOWNSTREAM ANALYSIS OF SPECIES, AND OVERALL DIVERSITY ################
################################################################################

#Transforming data into a phyloseq object
#GTDB first
sample.names = rownames(seqtab.nochim)
sample_df = data.frame(id = sample.names)
rownames(sample_df) = sample.names
colnames(taxa_table)[ncol(taxa_table)] = "Species 2" #Ensures the last column to be renamed
ps_object = phyloseq(otu_table(seqtab.nochim, taxa_are_rows = FALSE),
                     sample_data(sample_df),
                     tax_table(taxa_table))
dna = Biostrings::DNAStringSet(taxa_names(ps_object))
names(dna) = taxa_names(ps_object)
ps_object = merge_phyloseq(ps_object, dna)
taxa_names(ps_object) = paste0("ASV", seq(ntaxa(ps_object)))
ps_object #originally from 3,255 different taxa

#SILVA second
ps_object_silva = phyloseq(otu_table(seqtab.nochim, taxa_are_rows = FALSE),
                           sample_data(sample_df),
                           tax_table(taxa_silva))
ps_object_silva #Same number of ASVs, as expected (3,255)
taxa_table_silva = taxa_silva

#Alpha diversity on ASVs not species
ps_species = tax_glom(ps_object, taxrank = "Species")
ps_test = tax_glom(ps_object, taxrank = "Species 2")
View(tax_table(ps_object))
ps_species #218 species
ps_species = prune_taxa(taxa_sums(ps_species)>0, ps_species)
ps_species #still at 218 species, as expected

alphas = estimate_richness(ps_object) #Calculated on ASVs diversity rather than species!
alphas

#Normalization of data via CSS method on ASVs, not species
otu_matrix = as(otu_table(ps_object), "matrix")
nonzero_counts = colSums(otu_matrix>0) #counts nonzero ASVs per sample in a matrix format
summary(nonzero_counts) #summarizes it
ps_css1 = phyloseq_to_metagenomeSeq(ps_object)
p1 = cumNormStat(ps_css1)
ps_css1 = cumNorm(ps_css1, p = p1)
ps_css1_ps = ps_object
otu_table_ps_css1 = MRcounts(ps_css1, norm=TRUE)
otu_table_ps_css1 = as.matrix(otu_table_ps_css1) #normalized OTU table matrix
otu_table_ps_css1 = t(otu_table_ps_css1)
taxa_are_rows_original = taxa_are_rows(ps_object)
rownames(otu_table_ps_css1) = taxa_names(ps_object)
taxa_names(ps_object)
otu_table_ps_css1 = otu_table(otu_table_ps_css1, taxa_are_rows = taxa_are_rows_original)
otu_table(ps_css1_ps) = otu_table_ps_css1
sample_data(ps_css1_ps)$id = factor(sample_data(ps_css1_ps)$id) #forcing ID to be a factor in the metadata

ps_css1_ps
#View(otu_table(ps_css1_ps)) sanity check, it is normalized to CSS

#Beta diversity

beta_plotting<-function(ps_object, dist, meth){
  beta_ordination = ordinate(ps_object, method = meth, distance = dist)
  group_colors = c("#5f9c9d","#d36f6f","#786a87", "gray")
  beta_plot = plot_ordination(ps_object, ordination = beta_ordination, type = "samples", color = "id")
  plot<-beta_plot + scale_color_manual(values = group_colors)+
    scale_fill_manual(values = group_colors) + 
    theme(panel.background = element_rect(fill = "white"),
          axis.line = element_line(colour = "black"),
          panel.grid.major = element_line(colour = "gray", size=0.20),
          panel.grid.minor = element_line(colour = "gray", size = 0.05),
          panel.border = element_rect(colour="black", fill=NA, size = 1)
    )
  print(plot)
}

beta_plotting(ps_css1_ps, "bray","NMDS")

#Filtering low-quality mapped taxonomy and doing a painter's plot
#GTDB first
ps_filtered = tax_glom(ps_object, taxrank = "Species", NArm = TRUE) #Change to false?
ps_filtered #970 Species originally, or 218 when removing NAs for species
ps_filtered = subset_taxa(ps_filtered, !is.na(Kingdom)) #969 Species
ps_filtered = subset_taxa(ps_filtered, !is.na(Class)) #964 Species
ps_filtered
View(tax_table(ps_filtered))

ps_filtered_prop = transform_sample_counts(ps_filtered, function(x) x/sum(x)*100) #Transforming to proportions
ps_filtered_taxa = filter_taxa(physeq = ps_filtered_prop, function(x) mean(x) > 0.00005, prune =  TRUE) #Filtering anything below 0.005% to remove rare taxa
ps_filtered_taxa #803 Species

ps_filtered_taxa_genus = tax_glom(ps_filtered_taxa, taxrank="Genus")
#otu_table(ps_filtered_taxa_genus) Sanity check
ps_filtered_taxa_genus #543 taxa (genus level)

#SILVA second
colnames(taxa_table_silva)[ncol(taxa_table_silva)] = "Species 2"
ps_object_silva = merge_phyloseq(ps_object_silva, dna)
taxa_names(ps_object_silva) = paste0("ASV",seq(ntaxa(ps_object_silva)))
ps_object_silva = subset_taxa(ps_object_silva, !is.na(Kingdom))
ps_object_silva #Down to 3,244 ASVs
ps_object_silva = subset_taxa(ps_object_silva, !is.na(Class))
ps_object_silva #Down to 3,204 ASVs
ps_filtered_silva = tax_glom(ps_object_silva, taxrank = "Species", NArm = TRUE)
View(tax_table(ps_filtered_silva))

tax_table(ps)

tax_table(ps_filtered_silva)[,"Species"] = paste(tax_table(ps_filtered_silva)[, "Genus"],
                                               tax_table(ps_filtered_silva)[, "Species"],
                                               sep=" ") #This merges genus and species into the species column

ps_filtered_silva #867 species clustered when NAs are not removed, and 234 species when removed
ps_filtered_silva_genus = tax_glom(ps_object_silva, taxrank = "Genus", NArm = TRUE) #640 taxa when NA not removed and #504 if removed (genus level)
ps_filtered_silva_genus
ps_filtered_prop_silva = transform_sample_counts(ps_filtered_silva, function(x) x/sum(x)*100) #Transforming species to proportions
ps_filtered_taxa_silva = filter_taxa(physeq = ps_filtered_prop_silva, function(x) mean(x) > 0.00005, prune=TRUE)
ps_filtered_taxa_silva #Still 731 taxa after filtering
ps_filtered_prop_genus_silva = transform_sample_counts(ps_filtered_silva_genus, function(x) x/sum(x)*100)
ps_filtered_taxa_silva_genus = filter_taxa(physeq = ps_filtered_prop_genus_silva, function(x) mean(x) > 0.00005, prune=TRUE)
ps_filtered_taxa_silva_genus #Still 504 after filtering

#TO USE:
ps_filtered_taxa_silva_genus #proportion normalized
ps_filtered_taxa_silva #proportion normalized

painter_plot = function(ps_to_do){
  top = names(sort(taxa_sums(ps_to_do), decreasing=TRUE))[1:40]
  ps_relative_abundance_genus = prune_taxa(top, ps_to_do)
  df = psmelt(ps_relative_abundance_genus)
  #df Sanity check
  my_colors = c("#58643d","#6a42c3","#7ed05b","#c853bc","#cdb756",
              "#4b2e4b","#90c9b6","#c25d39","#8e92c2","#c35774",
              "red","blue","yellow","green","gray","pink",
              "purple","salmon","cyan","black")
  df$Species = factor(df$Species, levels = names(sort(tapply(df$Abundance, df$Species, sum), decreasing=TRUE))) #Change taxonomy rank if needed
  unique_genera = levels(df$Species) #Change here too
  repeated_colors = rep(my_colors, length.out = length(unique_genera))
  color_mapping = setNames(repeated_colors, unique_genera)

  plot_genus_top = ggplot(df, aes(x = Sample, y = Abundance, fill=Species)) + #here too
  geom_bar(stat = "identity", position="stack") +
  theme_minimal() + labs(y = "Relative abundance (%)") + 
  theme(axis.text.x = element_text(angle=90)) +
  scale_fill_manual(values=color_mapping)
  print(plot_genus_top)
}

painter_plot(ps_filtered_taxa_silva_genus)
painter_plot(ps_filtered_taxa_silva)
painter_plot(ps_filtered_taxa)
painter_plot(ps_filtered_taxa_genus)

length(unique(tax_table(ps_filtered_taxa_silva)[, "Species"]))
length(unique(tax_table(ps_filtered_silva_genus)[, "Genus"]))
ps_filtered_silva_genus


#Counting the proportion of ASVs in a table
ps_mycoplasmopsis_gtdb = subset_taxa(ps_object, Genus=="Mycoplasmopsis")
ps_mycoplasmopsis_bovis_gtdb = subset_taxa(ps_object, Species=="Mycoplasmopsis_A_bovirhinis(RS_GCF_900660515_1")
sample_sums(ps_object)
sample_sums(ps_mycoplasmopsis_gtdb)
sample_sums(ps_mycoplasmopsis_bovis_gtdb)
View(tax_table(ps_object))


tax_table(ps_object_silva)[,"Species"] = paste(tax_table(ps_object_silva)[, "Genus"],
                                                 tax_table(ps_object_silva)[, "Species"],
                                                 sep=" ") #This merges genus and species into the species column

ps_mycoplasmopsis_silva = subset_taxa(ps_object_silva, Genus=="Mycoplasmopsis")
ps_mycoplasmopsis_bovis_silva = subset_taxa(ps_object_silva, Species =="Mycoplasmopsis bovis")
sample_sums(ps_object_silva)
sample_sums(ps_mycoplasmopsis_silva)
sample_sums(ps_mycoplasmopsis_bovis_silva)

ps_bacteria_gtdb = subset_taxa(ps_object, Kingdom=="Bacteria")
sample_sums(ps_bacteria_gtdb)
sample_sums(ps_object)

ps_bacteria_silva = subset_taxa(ps_object_silva, Kingdom=="Archaea")
sample_sums(ps_bacteria_silva)
sample_sums(ps_object_silva)

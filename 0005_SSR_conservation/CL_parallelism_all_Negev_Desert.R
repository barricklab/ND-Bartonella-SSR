library(tidyverse)

## All output files will be placed here
output_path = "output"
input_path = "input"
dir.create(output_path, showWarnings = FALSE)

#########################################################################
## Load and tidy the orthogroup output from Orthofinder

orthogroups = read_tsv(file.path(input_path, "orthogroups.tsv"))
orthogroups = orthogroups %>% pivot_longer(!Orthogroup, names_to="species", values_to="locus_tag")
orthogroups = orthogroups %>% filter(!is.na(locus_tag))

orthogroups = orthogroups %>% separate_rows(locus_tag,sep=",\\s+")
orthogroups = orthogroups %>% rename(orthogroup = Orthogroup)
write_csv(orthogroups, file.path(output_path, "orthogroups_processed.csv"))

genes_per_orthogroup = orthogroups %>% group_by(orthogroup) %>% summarize(num_genes = n())


#########################################################################
## Load CL

contingency_loci = read.csv(file.path(input_path, "contingency_loci.csv"))

contingency_loci_merge = contingency_loci %>% select(genome_id, position, base, length, gene, location_type, gene_position, gene_length)
contingency_loci_merge = contingency_loci_merge %>% rename(locus_tag = gene)

orthogroups_merge = orthogroups %>% select(locus_tag, orthogroup)

X = contingency_loci_merge %>% left_join(orthogroups_merge, by="locus_tag") 

printable_list = function(in_list) {
  new_list = in_list
  new_list = new_list[new_list!=""]
  new_list = new_list[!is.na(new_list)]
  return (paste(new_list, collapse = "; "))
}

#########################################################################
## Add gene descriptions from Prokka annotation
gene_descriptions = read.csv(file.path(input_path, "gene_descriptions.csv"))
orthogroup_descriptions = orthogroups %>% left_join(gene_descriptions, by="locus_tag")
orthogroup_descriptions = orthogroup_descriptions %>% group_by(orthogroup) %>% summarize(gene = printable_list(unique(gene)), product = printable_list(unique(product)), avg_bp_length=mean(bp_length))
write_csv(orthogroup_descriptions, file.path(output_path, "orthogroup_descriptions.csv"))

#########################################################################
## Load species descriptions and categories

species_descriptions = read.csv(file.path(input_path, "strain_metadata.csv"))

species_descriptions = species_descriptions %>% mutate(category="other")
#to keep or not to keep recombinants?
#species_descriptions$category[(species_descriptions$genus_species=="Bartonella gerbillinarum") & !species_descriptions$recombinant] = "Bg"
#species_descriptions$category[(species_descriptions$genus_species=="Bartonella krasnovii") & !species_descriptions$recombinant] = "Bk"

species_descriptions$category[(species_descriptions$genus_species=="Bartonella gerbillinarum")] = "Bg"
species_descriptions$category[(species_descriptions$genus_species=="Bartonella krasnovii")] = "Bk"
species_descriptions$category[(species_descriptions$genus_species=="Bartonella khokhlovae")] = "Bh"
species_descriptions$category[(species_descriptions$genus_species=="Bartonella negeviensis")] = "Bn"

Bk_underscore_species = species_descriptions$full_name_underscores[species_descriptions$category=="Bk"]
Bg_underscore_species = species_descriptions$full_name_underscores[species_descriptions$category=="Bg"]
Bh_underscore_species = species_descriptions$full_name_underscores[species_descriptions$category=="Bh"]
Bn_underscore_species = species_descriptions$full_name_underscores[species_descriptions$category=="Bn"]

ND_underscore_species = c(Bk_underscore_species, Bg_underscore_species, Bh_underscore_species, Bn_underscore_species)

total_num_Bg = nrow(species_descriptions %>% filter(category == "Bg"))
total_num_Bk = nrow(species_descriptions %>% filter(category == "Bk"))
total_num_Bh = nrow(species_descriptions %>% filter(category == "Bh"))
total_num_Bn = nrow(species_descriptions %>% filter(category == "Bn"))
total_num_other = nrow(species_descriptions %>% filter(category == "other"))
write_csv(species_descriptions, "strain_metadata_processed.csv")

species_descriptions = species_descriptions %>% rename(genome_id = full_name_underscores)

# Add them to the contingency loci table
X = X %>% left_join(species_descriptions %>% select(genome_id, category, genus_species), by="genome_id")

#########################################################################
# Output tables of Negev Desert Orthogroups with >= 2 examples per species

orthogroup_avg_bp_lengths = orthogroups %>% filter(species %in% ND_underscore_species) %>% left_join(gene_descriptions, by="locus_tag") %>% group_by(orthogroup) %>% summarize(avg_bp_length=mean(bp_length))

Bg_represented_orthogroups = orthogroups %>% filter(species %in%  Bg_underscore_species) %>% group_by(orthogroup) %>% summarize(n=n()) %>% filter(n>=2) %>% left_join(orthogroup_avg_bp_lengths, by="orthogroup")

Bk_represented_orthogroups = orthogroups %>% filter(species %in%  Bk_underscore_species) %>% group_by(orthogroup) %>% summarize(n=n()) %>% filter(n>=2) %>% left_join(orthogroup_avg_bp_lengths, by="orthogroup")

Bh_represented_orthogroups = orthogroups %>% filter(species %in%  Bh_underscore_species) %>% group_by(orthogroup) %>% summarize(n=n()) %>% filter(n>=2) %>% left_join(orthogroup_avg_bp_lengths, by="orthogroup")

Bn_represented_orthogroups = orthogroups %>% filter(species %in%  Bn_underscore_species) %>% group_by(orthogroup) %>% summarize(n=n()) %>% filter(n>=2) %>% left_join(orthogroup_avg_bp_lengths, by="orthogroup")


write_tsv(Bg_represented_orthogroups, file=file.path(output_path, "Bg_orthogroups_represented.tsv"), col_names=T)

write_tsv(Bk_represented_orthogroups, file=file.path(output_path, "Bk_orthogroups_represented.tsv"), col_names=T)

write_tsv(Bh_represented_orthogroups, file=file.path(output_path, "Bh_orthogroups_represented.tsv"), col_names=T)

write_tsv(Bn_represented_orthogroups, file=file.path(output_path, "Bn_orthogroups_represented.tsv"), col_names=T)

###########################################################
# SSR analysis

##Filter types of SSRs

contingency_loci_minimum_AT_length = 9
contingency_loci_minimum_GC_length = 9

contingency_loci_filtered_on_length = X %>% 
  filter( 
     ( ((base=="A") | (base == "T")) & (length >= contingency_loci_minimum_AT_length) ) |
     ( ((base=="G") | (base == "C")) & (length >= contingency_loci_minimum_GC_length) )
    ) 

output_name = "undefined"


## This just requires that it be within the gene (Remember: 0-indexed positions)
if (TRUE) {
  output_name = "orthogroup_results_genic"
  #contingency_loci_filtered_on_length = contingency_loci_filtered_on_length %>% filter( !(gene_position < 0))
  contingency_loci_filtered_on_length = contingency_loci_filtered_on_length %>% filter(location_type == "gene")
}

#### Alternatives used for testing
if (FALSE){
#distance and not in promoter
  output_name = "orthogroup_results_first_150bp_of_gene"
  distance_within_gene = 150 
  contingency_loci_filtered_on_length = contingency_loci_filtered_on_length %>%  filter( (location_type == "gene") & (gene_position <= distance_within_gene) )
}

if (FALSE) {
  output_name = "orthogroup_results_all"
  contingency_loci_filtered_on_length = contingency_loci_filtered_on_length %>% filter( !(gene_position < -60))
}


contingency_loci_filtered_on_length = contingency_loci_filtered_on_length %>% mutate(base_length = paste0(base,length))

#################################################
sd = species_descriptions %>% select(genome_id,category)

### Count and composition for graphing
contingency_loci_representation_totals = contingency_loci_filtered_on_length %>% group_by(genome_id) %>% summarize(n=n())
contingency_loci_representation_totals = contingency_loci_representation_totals %>% left_join(sd, by="genome_id")
contingency_loci_representation_totals = contingency_loci_representation_totals %>% filter(category != "other")
contingency_loci_representation_totals_averages =  contingency_loci_representation_totals %>% group_by(category) %>% summarize(m=mean(n))

contingency_loci_representation_totals$category = factor(contingency_loci_representation_totals$category, levels = rev(c("Bh","Bk","Bg","Bn")))
contingency_loci_representation_totals_averages$category = factor(contingency_loci_representation_totals_averages$category, levels = rev(c("Bh","Bk","Bg","Bn")))

ggplot(contingency_loci_representation_totals, aes(y=category, x=n)) + geom_col(aes(x=m), data=contingency_loci_representation_totals_averages) + geom_jitter(height = 0.25) + theme_bw() +  scale_x_continuous(expand = c(0, 0))
ggsave(file.path(output_path, "CL_representation_totals_by_species.pdf"), height=5, width=8)

contingency_loci_representation = contingency_loci_filtered_on_length %>% group_by(genome_id, base) %>% summarize(n=n())
contingency_loci_representation = contingency_loci_representation %>% left_join(sd, by="genome_id")
contingency_loci_representation = contingency_loci_representation %>% filter(category != "other")
contingency_loci_representation_by_species = contingency_loci_representation %>% group_by(category, base) %>% summarize(n=sum(n))
contingency_loci_representation_by_species_totals = contingency_loci_representation_by_species %>% group_by(category) %>% summarize(total=sum(n))
contingency_loci_representation_by_species = contingency_loci_representation_by_species %>% left_join(contingency_loci_representation_by_species_totals, by="category")
contingency_loci_representation_by_species$fr = contingency_loci_representation_by_species$n / contingency_loci_representation_by_species$total

contingency_loci_representation_by_species$base = factor(contingency_loci_representation_by_species$base, levels=rev(c("G", "C", "A", "T")))

ggplot(contingency_loci_representation_by_species, aes(y=category, x=fr, fill=base)) + geom_col() + theme_bw() +  scale_x_continuous(expand = c(0, 0))
ggsave(file.path(output_path, "CL_representation_bases_by_species.pdf"), height=5, width=8)

contingency_loci_representation_by_species %>% group_by(base) %>% summarize(mean_fr=mean(fr))

# Make a final graph that is the distribution across the length of the gene

#gene position is 0 indexed
contingency_loci_length_fractions = contingency_loci_filtered_on_length %>% mutate(gene_fraction=(gene_position+1)/gene_length) %>% filter(category != "other")

# There are a few examples of an AAAAAAAAATG with -8 for the position and overlapping a gene that we need to filter out. Similarly filter out any that go past the end of the gene.
contingency_loci_length_fractions = contingency_loci_length_fractions %>% filter(gene_position+1>=1)
contingency_loci_length_fractions = contingency_loci_length_fractions %>% filter(gene_position+length<=gene_length)

for(this_category in unique(contingency_loci_length_fractions$category)) {

  ggplot(contingency_loci_length_fractions %>% filter(category==this_category), aes(x=gene_fraction)) + 
    geom_histogram(aes(y = after_stat(count / sum(count))), breaks = seq(0, 1, by = 0.05)) + 
    theme_bw() +theme_classic() +
    scale_y_continuous(limits=c(0, 0.35), expand=expansion(mult=c(0, 0.1)), breaks=seq(0, 0.35, by=0.05)) + 
    scale_x_continuous(limits=c(0,1), expand = expansion(add=c(0, 0)), breaks=seq(0, 1, by = 0.1))
  
  ggsave(file.path(output_path, paste0(this_category, "_CL_representation_by_relative_gene_position.pdf")), height=5, width=8)

}

# Do binomial tests on counts being different in first 5% than expected

first_5_counts = contingency_loci_length_fractions %>% group_by(category) %>% filter(gene_fraction<=0.05) %>% summarize(first_5 = n())
last_95_counts = contingency_loci_length_fractions %>% group_by(category) %>% filter(gene_fraction>0.05) %>% summarize(last_95 = n())
fraction_counts = first_5_counts %>% left_join(last_95_counts, by="category")

fraction_counts$first_5_fraction = fraction_counts$first_5 / (fraction_counts$first_5 + fraction_counts$last_95)

for(this_category in fraction_counts$category) {
  this_fraction_counts = fraction_counts %>% filter(category==this_category)
  p_value = binom.test(c(this_fraction_counts$first_5, this_fraction_counts$last_95), p = 0.05, alternative="greater")$p.value
  cat(this_category, " ", p_value, "\n" )
}


#################################################

contingency_loci_with_different_base_lengths = contingency_loci_filtered_on_length %>% group_by(orthogroup, base_length) %>% summarize(num_base_length=n())
contingency_loci_summary = contingency_loci_with_different_base_lengths %>% mutate(count_string = paste0(base_length, "x", num_base_length))
contingency_loci_summary = contingency_loci_summary %>% arrange(desc(num_base_length)) %>% group_by(orthogroup) %>% summarize (contingency_loci_counts = toString(count_string))

orthogroup_overall_averages = contingency_loci_filtered_on_length %>% group_by(orthogroup) %>% summarize(num_contingency_loci=n())


num_genes_with_contingency_loci = contingency_loci_filtered_on_length %>% group_by(orthogroup) %>% summarize(num_genes_with_contingency_loci=length(unique(locus_tag)))

orthogroup_overall_averages = orthogroup_overall_averages %>% left_join(num_genes_with_contingency_loci, by="orthogroup")


orthogroup_overall_averages = orthogroup_overall_averages %>% left_join(genes_per_orthogroup)
orthogroup_overall_averages = orthogroup_overall_averages %>% left_join(contingency_loci_summary, by="orthogroup")

orthogroup_overall_averages = orthogroup_overall_averages %>% mutate(contingency_loci_per_gene = num_contingency_loci/num_genes)

orthogroup_overall_averages = orthogroup_overall_averages %>% arrange(desc(contingency_loci_per_gene))

#add orthogroup descriptions
orthogroup_overall_averages = orthogroup_overall_averages %>% left_join(orthogroup_descriptions, by="orthogroup")

orthogroup_overall_averages = orthogroup_overall_averages %>% filter(!is.na(orthogroup))

Bg_with_CL = contingency_loci_filtered_on_length %>% filter(category=="Bg") %>% group_by(orthogroup) %>% summarize(num_Bg_CL = length(unique(genome_id)))
Bk_with_CL = contingency_loci_filtered_on_length %>% filter(category=="Bk") %>% group_by(orthogroup) %>% summarize(num_Bk_CL = length(unique(genome_id)))
Bh_with_CL = contingency_loci_filtered_on_length %>% filter(category=="Bh") %>% group_by(orthogroup) %>% summarize(num_Bh_CL = length(unique(genome_id)))
Bn_with_CL = contingency_loci_filtered_on_length %>% filter(category=="Bn") %>% group_by(orthogroup) %>% summarize(num_Bn_CL = length(unique(genome_id)))

### Output orthogroups with CL in each species
Bg_orthogroups_CL = Bg_with_CL %>% filter(num_Bg_CL>=2) %>% filter(!is.na(orthogroup))
write_tsv(Bg_orthogroups_CL %>% rename(n=num_Bg_CL), file=file.path(output_path, "Bg_orthogroups_CL.tsv"), col_names=T)

Bk_orthogroups_CL = Bk_with_CL %>% filter(num_Bk_CL>=2) %>% filter(!is.na(orthogroup))
write_tsv(Bk_orthogroups_CL %>% rename(n=num_Bk_CL), file=file.path(output_path, "Bk_orthogroups_CL.tsv"), col_names=T)

Bh_orthogroups_CL = Bh_with_CL %>% filter(num_Bh_CL>=2) %>% filter(!is.na(orthogroup))
write_tsv(Bh_orthogroups_CL %>% rename(n=num_Bh_CL), file=file.path(output_path, "Bh_orthogroups_CL.tsv"), col_names=T)

Bn_orthogroups_CL = Bn_with_CL %>% filter(num_Bn_CL>=2) %>% filter(!is.na(orthogroup))
write_tsv(Bn_orthogroups_CL %>% rename(n=num_Bn_CL), file=file.path(output_path, "Bn_orthogroups_CL.tsv"), col_names=T)


other_with_CL = contingency_loci_filtered_on_length %>% filter(category=="other") %>% group_by(orthogroup) %>% summarize(num_other_CL = length(unique(genome_id)))

Bg_with_gene = orthogroups %>% filter(species %in% Bg_underscore_species) %>% group_by(orthogroup) %>% summarize(num_Bg_gene = length(unique(species)))
Bk_with_gene = orthogroups %>% filter(species %in% Bk_underscore_species) %>% group_by(orthogroup) %>% summarize(num_Bk_gene = length(unique(species)))
Bh_with_gene = orthogroups %>% filter(species %in% Bh_underscore_species) %>% group_by(orthogroup) %>% summarize(num_Bh_gene = length(unique(species)))
Bn_with_gene = orthogroups %>% filter(species %in% Bn_underscore_species) %>% group_by(orthogroup) %>% summarize(num_Bn_gene = length(unique(species)))
other_with_gene = orthogroups %>% filter( !(species %in% Bk_underscore_species) & !(species %in% Bg_underscore_species) & !(species %in% Bh_underscore_species) & !(species %in% Bn_underscore_species)) %>% group_by(orthogroup) %>% summarize(num_other_gene = length(unique(species)))

orthogroup_overall_averages = orthogroup_overall_averages %>% left_join(Bg_with_CL, by="orthogroup")
orthogroup_overall_averages$num_Bg_CL[is.na(orthogroup_overall_averages$num_Bg_CL)]=0
orthogroup_overall_averages = orthogroup_overall_averages %>% left_join(Bg_with_gene, by="orthogroup")
orthogroup_overall_averages$num_Bg_gene[is.na(orthogroup_overall_averages$num_Bg_gene)]=0

orthogroup_overall_averages = orthogroup_overall_averages %>% left_join(Bk_with_CL, by="orthogroup")
orthogroup_overall_averages$num_Bk_CL[is.na(orthogroup_overall_averages$num_Bk_CL)]=0
orthogroup_overall_averages = orthogroup_overall_averages %>% left_join(Bk_with_gene, by="orthogroup")
orthogroup_overall_averages$num_Bk_gene[is.na(orthogroup_overall_averages$num_Bk_gene)]=0

orthogroup_overall_averages = orthogroup_overall_averages %>% left_join(Bh_with_CL, by="orthogroup")
orthogroup_overall_averages$num_Bh_CL[is.na(orthogroup_overall_averages$num_Bh_CL)]=0
orthogroup_overall_averages = orthogroup_overall_averages %>% left_join(Bh_with_gene, by="orthogroup")
orthogroup_overall_averages$num_Bh_gene[is.na(orthogroup_overall_averages$num_Bh_gene)]=0

orthogroup_overall_averages = orthogroup_overall_averages %>% left_join(Bn_with_CL, by="orthogroup")
orthogroup_overall_averages$num_Bn_CL[is.na(orthogroup_overall_averages$num_Bn_CL)]=0
orthogroup_overall_averages = orthogroup_overall_averages %>% left_join(Bn_with_gene, by="orthogroup")
orthogroup_overall_averages$num_Bn_gene[is.na(orthogroup_overall_averages$num_Bn_gene)]=0

orthogroup_overall_averages = orthogroup_overall_averages %>% left_join(other_with_CL, by="orthogroup")
orthogroup_overall_averages$num_other_CL[is.na(orthogroup_overall_averages$num_other_CL)]=0
orthogroup_overall_averages = orthogroup_overall_averages %>% left_join(other_with_gene, by="orthogroup")
orthogroup_overall_averages$num_other_gene[is.na(orthogroup_overall_averages$num_other_gene)]=0

orthogroup_overall_averages$fr_Bg_genes_with_CL = orthogroup_overall_averages$num_Bg_CL /orthogroup_overall_averages$num_Bg_gene
orthogroup_overall_averages$fr_Bg_genes_with_CL[is.na(orthogroup_overall_averages$fr_Bg_genes_with_CL)] = 0

orthogroup_overall_averages$fr_Bk_genes_with_CL = orthogroup_overall_averages$num_Bk_CL /orthogroup_overall_averages$num_Bk_gene
orthogroup_overall_averages$fr_Bk_genes_with_CL[is.na(orthogroup_overall_averages$fr_Bk_genes_with_CL)] = 0

orthogroup_overall_averages$fr_Bh_genes_with_CL = orthogroup_overall_averages$num_Bh_CL /orthogroup_overall_averages$num_Bh_gene
orthogroup_overall_averages$fr_Bh_genes_with_CL[is.na(orthogroup_overall_averages$fr_Bh_genes_with_CL)] = 0

orthogroup_overall_averages$fr_Bn_genes_with_CL = orthogroup_overall_averages$num_Bn_CL /orthogroup_overall_averages$num_Bn_gene
orthogroup_overall_averages$fr_Bn_genes_with_CL[is.na(orthogroup_overall_averages$fr_Bn_genes_with_CL)] = 0

orthogroup_overall_averages$fr_other_genes_with_CL = orthogroup_overall_averages$num_other_CL /orthogroup_overall_averages$num_other_gene
orthogroup_overall_averages$fr_other_genes_with_CL[is.na(orthogroup_overall_averages$fr_other_genes_with_CL)] = 0

orthogroup_overall_averages$sum_fr_three_category_genes_with_CL = 
  orthogroup_overall_averages$fr_Bg_genes_with_CL +
  orthogroup_overall_averages$fr_Bk_genes_with_CL +
  orthogroup_overall_averages$fr_Bh_genes_with_CL +
  orthogroup_overall_averages$fr_Bn_genes_with_CL +
  orthogroup_overall_averages$fr_other_genes_with_CL

# Flag ones with at least two examples
orthogroup_overall_averages$Bg_orthogroup_represented = orthogroup_overall_averages$num_Bg_gene >2
orthogroup_overall_averages$Bg_orthogroup_CL = orthogroup_overall_averages$num_Bg_CL >=2

orthogroup_overall_averages$Bk_orthogroup_represented = orthogroup_overall_averages$num_Bk_gene >2
orthogroup_overall_averages$Bk_orthogroup_CL = orthogroup_overall_averages$num_Bk_CL >=2

orthogroup_overall_averages$Bh_orthogroup_represented = orthogroup_overall_averages$num_Bh_gene >2
orthogroup_overall_averages$Bh_orthogroup_CL = orthogroup_overall_averages$num_Bh_CL >=2

orthogroup_overall_averages$Bn_orthogroup_represented = orthogroup_overall_averages$num_Bn_gene >2
orthogroup_overall_averages$Bn_orthogroup_CL = orthogroup_overall_averages$num_Bn_CL >=2

orthogroup_overall_averages$num_ND_genes_with_CL = orthogroup_overall_averages$num_Bk_CL + orthogroup_overall_averages$num_Bn_CL + orthogroup_overall_averages$num_Bh_CL + orthogroup_overall_averages$num_Bg_CL

orthogroup_overall_averages$num_ND_genes = orthogroup_overall_averages$num_Bk_gene + orthogroup_overall_averages$num_Bn_gene + orthogroup_overall_averages$num_Bh_gene + orthogroup_overall_averages$num_Bg_gene

#Require there be two examples in at least one species
orthogroup_overall_averages = orthogroup_overall_averages %>% filter( num_Bg_CL>=2 | num_Bk_CL>=2 | num_Bh_CL>=2 | num_Bn_CL>=2)

print(n=50, orthogroup_overall_averages)
write_csv(orthogroup_overall_averages, file.path(output_path, paste0(output_name, ".csv")))

#Output detailed info for each group
folder_name = output_name
dir.create(file.path(output_path,folder_name), showWarnings=F)
for (i in 1:nrow(orthogroup_overall_averages) ) {
  this_orthogroup = orthogroup_overall_averages$orthogroup[i]
  orthogroup_results = orthogroups %>% filter(orthogroup == this_orthogroup)
  
  orthogroup_results = orthogroup_results %>% left_join(contingency_loci_filtered_on_length %>% select(locus_tag, gene_position, base_length), by="locus_tag")
  
  orthogroup_results = orthogroup_results %>% left_join(gene_descriptions, by="locus_tag")
  
  # Add missing species so it's easier to see what didn't have a gene in the orthogroup...
  add_species = species_descriptions$genome_id[!(species_descriptions$genome_id %in% orthogroup_results$species)]
  orthogroup_results = orthogroup_results %>% add_row(orthogroup=this_orthogroup, species=add_species)
  orthogroup_results = orthogroup_results %>% arrange(species, locus_tag)
  
 # write_csv(orthogroup_results, file.path(folder_name, paste0(sprintf("%03d",i),"_", this_orthogroup, ".csv")))
  write_csv(orthogroup_results, file.path(output_path, folder_name, paste0(this_orthogroup, ".csv")))
  
  #Filter to just Bg and Bk subsets
  
  #write_csv(orthogroup_results %>% filter(species %in% Bk_underscore_species), file.path(folder_name, paste0(sprintf("%03d",i),"_", this_orthogroup, ".Bk.csv")))
  #write_csv(orthogroup_results %>% filter(species %in% Bg_underscore_species), file.path(folder_name, paste0(sprintf("%03d",i),"_", this_orthogroup, ".Bg.csv")))
  #write_csv(orthogroup_results %>% filter( !(species %in% Bk_underscore_species) & !(species %in% Bg_underscore_species)), file.path(folder_name, paste0(sprintf("%03d",i),"_", this_orthogroup, ".Bg.csv")))
  
#  write_csv(orthogroup_results %>% filter(species %in% Bk_underscore_species), file.path(folder_name, paste0(this_orthogroup, ".Bk.csv")))
#  write_csv(orthogroup_results %>% filter(species %in% Bg_underscore_species), file.path(folder_name, paste0(this_orthogroup, ".Bg.csv")))
#  write_csv(orthogroup_results %>% filter( !(species %in% Bk_underscore_species) & !(species %in% Bg_underscore_species)), file.path(folder_name, paste0(this_orthogroup, ".Bg.csv")))
  
}

#### Venn diagram of CL orthogroups
library(ggVennDiagram)
combos = orthogroup_overall_averages %>% select(Bg_orthogroup_CL, Bk_orthogroup_CL, Bh_orthogroup_CL, Bn_orthogroup_CL)

z<- list(Bh=Bh_orthogroups_CL$orthogroup,
         Bk=Bk_orthogroups_CL$orthogroup,
         Bg=Bg_orthogroups_CL$orthogroup,
          Bn=Bn_orthogroups_CL$orthogroup)

ggVennDiagram(z, label_alpha=0)
ggsave(file.path(output_path,"venn_diagram.pdf"))

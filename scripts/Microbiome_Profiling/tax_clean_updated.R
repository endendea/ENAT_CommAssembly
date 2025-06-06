## Clean taxonomy information from non-informative tags in phyloseq
## clean taxonomy tags with no information
library(phyloseq)
library(tidyverse)
library(glue)

tax_clean <- function(psdata){
  
  psdata_name <- deparse(substitute(psdata))
  # save uncleaned psdata
  if(!dir.exists("output_data")){dir.create("output_data")}
  saveRDS(psdata, file = paste0("output_data/",psdata_name,"_uncleaned.rds"))
  
  # ASV count start
  psdata_in = psdata
  
  # specify NA taxon name tags to last known taxon names
  tax.clean <- data.frame(phyloseq::tax_table(psdata))
  
  tax.clean2 = 
    tax.clean %>% 
    mutate_if(is.factor, as.character) %>% 
    mutate(across(everything(),
                  ~ str_replace_all(.x,  "Ambiguous_taxa|metagenome|uncultured archeaon|uncultured bacterium|uncultured prokaryote|uncultured soil bacterium|uncultured rumen bacterium|uncultured compost bacterium|uncultured organism|uncultured$",
                                    replacement = NA_character_
                  ))) %>%
    replace(is.na(.), NA_character_) %>% 
    mutate(Phylum = if_else(is.na(Phylum), paste0("Phylum of ", Kingdom), Phylum),
           Class = if_else(is.na(Class), paste0("Class of ", Phylum), Class),
           Order = if_else(is.na(Order), paste0("Order of ", Class), Order),
           Family = if_else(is.na(Family), paste0("Family of ", Order), Family),
           Genus = if_else(is.na(Genus), paste0("Genus of ", Family), Genus),
           Species = if_else(is.na(Species), paste0("Species of ", Genus), Species)
    ) %>% 
    mutate(across(.cols = Kingdom:Species, .fns = ~if_else(str_detect(.,'\\bof\\b.*\\bof\\b'), paste0(word(., 1)," ", word(., 2)," ", word(., -1)), .))) %>% 
    mutate(Kingdom = if_else(Kingdom == "d__Bacteria", "Bacteria", Kingdom))
  
  # put cleaned tax_table into phyloseq object
  phyloseq::tax_table(psdata) <- phyloseq::tax_table(as.matrix(tax.clean2))
  
  # remove unwanted taxa such as Mitochondria, Chloroplasts, Unclassified Kingdom, Eukaryota, etc.
  psdata <-
    psdata %>% subset_taxa(
        Class != "Chloroplast" & 
        Order != "Chloroplast" &
        Family != "Mitochondria" &
        Kingdom != "d__Eukaryota" &
        Kingdom != "Unassigned" &
        Phylum != "Phylum of d__Bacteria" &  
        Phylum != "Phylum of Bacteria"   
    )
  
  # count ASVs after cleaning
  psdata_out = psdata
  
  # difference
  removed_ASV_count = ntaxa(psdata_in) - ntaxa(psdata_out)
  
  #save psdata after cleaning as RDS object
  saveRDS(psdata, file = glue::glue("output_data/",{psdata_name},"_cleaned.rds", .sep = ""))
  print(
    paste("taxonomy table cleaned: ",  removed_ASV_count, " ASVs removed.  Phyloseq object saved as .rds object in output_data")
  )
  return(psdata)
}



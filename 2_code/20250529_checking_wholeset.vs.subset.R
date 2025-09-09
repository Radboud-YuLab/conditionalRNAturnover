

used_data_static_rev.run1 <- readRDS("data_KC_A/derived_data/compiled_data/revision/static/20250402_static.hl.basic.codonfreq_clean.RDS")
used_data_Static_subset.run1 <- readRDS("data_KC_A/derived_data/compiled_data/20240822_hl.mean.all_.basic.genetic.RDS")

used_data_dynamic_rev.run2 <- readRDS("data_KC_A/derived_data/compiled_data/revision/dynamic/20250515_dynamic.hl_unscaled.progeny.RDS")
used_data_dynamic_subset_run <- readRDS("data_KC_A/derived_data/compiled_data/20240910_hl.basic.genetic.unscaled.progeny.RDS")


# the different genes that are in the revised (with a cutoff) and subset df in dynamic model
length(setdiff(unique(used_data_dynamic_rev.run2$ensembl_gene_id),unique(used_data_dynamic_subset_run$ensembl_gene_id)))

# the different genes that are in the revised (with a cutoff) and subset df in static model
length(setdiff(unique(used_data_Static_subset.run1$ensembl_gene_id),unique(used_data_static_rev.run1$ensembl_gene_id)))

# there is some difference in which genes are incorporated in the df even though the sizes are similar



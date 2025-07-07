################################################################################
##### SUBSCRIPT_FIG3A
################################################################################
## sub script to run heatmap on data.  Need pre_processed set1 data (train1) and set1 Pheno

library(ComplexHeatmap)
setwd("~/Desktop/DESKTOP/ALZHEIMERS/EG_Pred/Heatmap")

# Step 2: load the centered and scaled transformed train SET1 data
gene_matrix_scaled <- train1
heatmap_pheno <- SET1_pheno
#subset SET1 by IDs
heatmap_pheno <-heatmap_pheno[index,]

# Prepare clinical annotations
heatmap_pheno$pTau.Aβ42<-heatmap_pheno$pTAU.ABETA
heatmap_pheno$pTau<-heatmap_pheno$pTAU
heatmap_pheno$tTau<-heatmap_pheno$tTAU
heatmap_pheno$Aβ42<-heatmap_pheno$ABETA

clinical_annot <- heatmap_pheno[, c("pTau", "tTau", "Aβ42", "pTau.Aβ42")]

# Sort clinical_annot based on pTAU.ABETA in ascending order
sorted_indices <- order(clinical_annot$pTau.Aβ42)
clinical_annot <- clinical_annot[sorted_indices, , drop = FALSE]  # Keep structure

# Reorder gene expression matrix rows to match sorted clinical_annot
gene_matrix_scaled <- gene_matrix_scaled[sorted_indices, , drop = FALSE]

# Ensure row names remain consistent
rownames(clinical_annot) <- rownames(gene_matrix_scaled)

summary(clinical_annot$pTau.Aβ42)

heatmap_pheno$Clinical_Condition<-heatmap_pheno$DX
unique(heatmap_pheno$Clinical_Condition)

# Load necessary libraries
# Convert clinical data into a heatmap annotation object
heatmap_pheno$Clinical_Condition[heatmap_pheno$Clinical_Condition=="MS"| heatmap_pheno$Clinical_Condition=="Imp"|heatmap_pheno$Clinical_Condition=="IMP"| heatmap_pheno$Clinical_Condition=="Other"]="Other"
# Define color mapping for clinical conditions (Clinical_Condition)
Clinical_Condition_colors <- c("Control" = "lightyellow", "MCI" = "grey","AD"="black","Other"="lightblue")
# Recreate heatmap annotation with sorted data
ha <- HeatmapAnnotation(
  Clinical_Condition = heatmap_pheno$Clinical_Condition[sorted_indices],  # Apply sorting to Clinical_Condition labels
  df = clinical_annot,  
  col = list(Clinical_Condition = Clinical_Condition_colors,
    pTau = colorRamp2(c(min(clinical_annot$pTau), median(clinical_annot$pTau), max(clinical_annot$pTau)), c("blue", "white", "red")),
    tTau = colorRamp2(c(min(clinical_annot$tTau), median(clinical_annot$tTau), max(clinical_annot$tTau)), c("blue", "white", "red")),
    Aβ42 = colorRamp2(c(min(clinical_annot$Aβ42), median(clinical_annot$Aβ42), max(clinical_annot$Aβ42)), c("blue", "white", "red")),
    pTau.Aβ42 = colorRamp2(c(min(clinical_annot$pTau.Aβ42), 0.025, max(clinical_annot$pTau.Aβ42)), c("darkblue", "white", "yellow"))
  ),
  annotation_name_gp = gpar(fontsize = 14, fontface = "bold"), # Make annotation labels bold and bigger
  annotation_legend_param = list(
    Clinical_Condition = list(title = "Clinical_Condition")
  ))


tiff("PLOTS/SET1_44CORGenes_pTAU_heatmap_sorted.tiff", width = 7, height = 4, units = "in", res = 500)

Heatmap(
  t(gene_matrix_scaled), 
  name = "Gene Expression", 
  top_annotation = ha, 
  column_names_gp = gpar(fontsize = 4),         # Increased column label font size
  row_names_gp = gpar(fontsize = 4),            # Increased row label font size
  cluster_rows = FALSE, 
  cluster_columns = FALSE,
  show_row_names = TRUE,
  heatmap_legend_param = list(
    title_gp = gpar(fontsize = 2),#, fontface = "bold"),  # Bigger, bold legend title
    labels_gp = gpar(fontsize = 2)                     # Bigger legend label text
  )
)
dev.off()


# Save as a TIFF file
tiff("PLOTS/SET1_sig_pTAU_heatmap_sorted.tiff", width = 12, height = 3, units = "in", res = 300)
# Generate the heatmap with the updated annotation
Heatmap(t(gene_matrix_scaled), 
        name = "Gene Expression", 
        top_annotation = ha, 
        column_names_gp = gpar(fontsize = 18),
        row_names_gp = gpar(fontsize = 14),
        cluster_rows = FALSE, 
        cluster_columns = FALSE,
        #column_title = "SET2 significant genes vs diagnostic biomarkers", 
        show_row_names = TRUE,
        heatmap_legend_param = list(
          title_gp = gpar(fontsize = 12, fontface = "bold")  # Make legend title bold and bigger
        )
)

# Close the TIFF device
dev.off()




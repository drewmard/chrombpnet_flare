library(data.table)
library(rhdf5)
library(reshape2)
library(ggplot2)
library(ggseqlogo)
library(patchwork)
library(dplyr)
library(tidyr)

########################################

# Output:
outdir = "/Users/amarderstein/Library/Mobile Documents/com~apple~CloudDocs/Documents/Research/chrombpnet_flare/output/visualize_variants"

# SHAP:
# f = "/data1/offitk/mardera1/chrombpnet_flare/output/variant_shap/HSCs/myb_motif_snps/Trisomy_Controls/fold_0/HSCs.fold_0.variant_shap.counts.h5"
# f = "/Users/amarderstein/Downloads/HSCs.fold_0.variant_shap.counts.h5"
# f = "/Users/amarderstein/Downloads/genetics/trevino_2021.c7.fold_0.variant_shap.counts.h5"
# f = "/Users/amarderstein/Downloads/ahi1_test/HSCs.fold_0.variant_shap.counts.h5"
f_shap = "/Users/amarderstein/Downloads/ahi1_mhe/HSCs.fold_0.variant_shap.counts.h5"

# PRED:
# f_pred = "/Users/amarderstein/Downloads/ahi1_test/HSCs.fold_0.scores.variant_predictions.h5"
f_pred = "/Users/amarderstein/Downloads/ahi1_mhe/HSCs.fold_0.scores.variant_predictions.h5"

# ISM:
f_scores = "~/Downloads/ahi1_mhe/HSCs.fold_0.scores.variant_scores.tsv"

# PHYLOP:
f_phylop = "/Users/amarderstein/Downloads/ahi1_test/phylop_region.bedGraph"

###########
# pos_of_interest = 135323454
pos_of_interest = 135323396; idx = 208
# pos_of_interest = 135323392; idx = 196
# pos_of_interest = 135323380; idx = 162
plot_title = "H-Me"
###################################################

# predictions = fread(f_scores,data.table = F,stringsAsFactors = F)
# subset(predictions,pos==135323396-16)
# as.data.frame(predictions)[idx <- which.max(predictions$abs_logfc),];idx
# pos_of_interest
# subset(as.data.frame(predictions),pos==135323392)[,c("pos","allele1","allele2","logfc")]

##############################
# f = "~/Downloads/HSCs.fold_0.scores.variant_scores.tsv"
# df=fread(f,data.table = F,stringsAsFactors = F);df
width_to_use=20
# dir.create("~/Downloads/myb_childhood_leukemia")
dir.create(outdir)

visualize_ism_heatmap=TRUE; ism_heatmap_prefix = "h_me"
visualize_pred=TRUE
visualize_phylop=TRUE
visualize_ism=TRUE
save_plot=TRUE
print_plot=TRUE
variant_line = TRUE
for (index in idx:idx) {
  # for (index in 1:1) {
  # for (index in 1:5) {
  flank      <- 75 # 50
  flank2      <- 200
  center_shift = 0
  
  alleles <- h5read(f_shap, "/alleles")
  seq_arr <- h5read(f_shap, "projected_shap/seq")
  # seq_arr <- h5read(f, "shap/seq")
  idx1 <- alleles == 0                     # TRUE/FALSE mask for allele 1
  idx2 <- alleles == 1                     # mask for allele 2
  allele1_shap <- seq_arr[,,idx1 , drop = FALSE]   # keeps all other dims
  allele2_shap <- seq_arr[,,idx2 , drop = FALSE]
  
  # Check new dimensions
  dim(allele1_shap)
  dim(allele2_shap)
  
  pos_center <- 1057 + center_shift
  
  # low_pos_idx = pos_center - (flank) - 1
  # hi_pos_idx = pos_center + (flank) + 1
  low_pos_idx = pos_center - (flank) + 1
  hi_pos_idx = pos_center + (flank) + 1
  
  shap_ref_window <- allele1_shap[low_pos_idx:hi_pos_idx,,index]     
  shap_alt_window <- allele2_shap[low_pos_idx:hi_pos_idx,,index]
  
  #––– 4) compute y‑limits (10% padding)
  all_vals <- c(shap_ref_window, shap_alt_window)
  # ylim <- c(min(all_vals), max(all_vals)) * 1.1
  # ylim <- c(-0.03,0.1)
  
  #––– 6) panels 2 & 3: sequence‑logo of SHAP scores
  #    ggseqlogo wants a matrix [letters × positions], so transpose
  mat_ref <- t(shap_ref_window)
  rownames(mat_ref) <- c("A","C","G","T")
  
  mat_alt <- t(shap_alt_window)
  rownames(mat_alt) <- c("A","C","G","T")
  
  # if plotting on reverse strand:
  # n <- ncol(mat_ref)
  # mat_ref <- rbind(
  #   A = mat_ref["T", ],  # now shows T→A contributions
  #   C = mat_ref["G", ],  # G→C
  #   G = mat_ref["C", ],  # C→G
  #   T = mat_ref["A", ]   # A→T
  # )[, n:1]              # then flip the columns
  # mat_alt <- rbind(
  #   A = mat_alt["T", ],  # now shows T→A contributions
  #   C = mat_alt["G", ],  # G→C
  #   G = mat_alt["C", ],  # C→G
  #   T = mat_alt["A", ]   # A→T
  # )[, n:1]              # then flip the columns
  # # x0 <- n - x0 + 1                       # mirrored x
  
  ylim=c(min(min(mat_ref),min(mat_alt)),max(max(mat_ref),max(mat_alt)))
  
  n = ncol(mat_ref)
  p2 <- ggseqlogo(mat_ref, method = "custom") +
    coord_cartesian(xlim = c(1, ncol(mat_ref)), ylim = ylim) +
    labs(x = NULL, y = "SHAP") +
    # scale_x_continuous(breaks = seq(1, ncol(mat_ref), by = 50),
    #                    labels = seq(-flank, flank, by = 50)) +
    theme_bw() +
    theme(
      panel.grid = element_blank(),
      axis.text.x = element_text(size = 14),
      axis.text.y = element_text(size = 14),
      plot.margin = margin(t = 5, b = 5),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 1)
    ) +
    scale_x_continuous(
      breaks = c(1,n),
      labels = c(-flank,flank)
    ) +
    geom_hline(yintercept = 0,col='lightgrey',lty='dotted');#p2
  
  if (variant_line) {
    p2 <- p2 + geom_vline(xintercept = flank+1 - center_shift,lty='dashed',col='grey')
  }
  
  
  
  p3 <- ggseqlogo(mat_alt, method = "custom") +
    coord_cartesian(xlim = c(1, ncol(mat_alt)), ylim = ylim) +
    labs(x = NULL, y = "SHAP") +
    theme_bw() +
    theme(
      panel.grid = element_blank(),
      axis.text.x = element_text(size = 14),
      axis.text.y = element_text(size = 14),
      plot.margin = margin(t = 5, b = 5),
      panel.border = element_rect(color = "red", fill = NA, linewidth = 1)
    ) +
    scale_x_continuous(
      breaks = c(1,n),
      labels = c(-flank,flank)
    ) +
    geom_hline(yintercept = 0,col='lightgrey',lty='dotted');#p2
  
  if (variant_line) {
    p3 <- p3 + geom_vline(xintercept = flank+1 - center_shift,lty='dashed',col='grey')
  }
  
  # final_plot
  
  if (visualize_pred) {
    # f_pred = "/Users/amarderstein/Downloads/ahi1_test/HSCs.fold_0.scores.variant_predictions.h5"
    # f_pred <- "~/Downloads/HSCs.fold_0.scores.variant_predictions.h5"
    # h5ls(f_pred)
    # c(h5read(f_pred, "/observed/allele1_pred_counts"),h5read(f_pred, "/observed/allele2_pred_counts"))
    # log2(h5read(f_pred, "/observed/allele2_pred_counts")/h5read(f_pred, "/observed/allele1_pred_counts"))
    
    # read 1D profile for this variant (index)
    a1 <- h5read(f_pred, "/observed/allele1_pred_profiles")[, index]
    a2 <- h5read(f_pred, "/observed/allele2_pred_profiles")[, index]
    
    pos_center <- 500
    
    low_pos_idx = pos_center - (flank2) + 1
    hi_pos_idx = pos_center + (flank2) + 1
    
    df_pred = data.frame(seq=1:length(seq(low_pos_idx:hi_pos_idx)),
                         a1=as.numeric(a1[low_pos_idx:hi_pos_idx]),
                         a2=a2[low_pos_idx:hi_pos_idx])
    
    n = nrow(df_pred)
    spanval=0.5
    p1 <- ggplot(df_pred, aes(x = seq)) + 
      geom_line(aes(y=a1),alpha=1,col='black') +
      geom_line(aes(y=a2),col='red',alpha=1) +
      # geom_smooth(aes(y=a1),col='black',se=F,span=spanval,alpha=0.1) +
      # geom_smooth(aes(y=a2),col='red',se=F,span=spanval,alpha=0.1) +
      # geom_line(aes(y=a1),alpha=0.1,col='black') +
      # geom_line(aes(y=a2),col='red',alpha=0.1) +
      # geom_smooth(aes(y=a1),col='black',se=F,span=spanval) +
      # geom_smooth(aes(y=a2),col='red',se=F,span=spanval) +
      # scale_x_continuous(
      #   breaks = c(1,flank+1,n),
      #   labels = c(-flank,0,flank)
      # ) +
      scale_x_continuous(
        breaks = c(1,n),
        labels = c(-flank2,flank2)
      ) +
      labs(x = NULL, y = "Pred. ATAC") +
      theme_bw() +
      theme(
        panel.grid = element_blank(),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        plot.margin = margin(t = 5, b = 5)
      ) + scale_y_continuous(breaks = scales::breaks_width(1))
  }
  
  if (visualize_phylop) {
    phylop = fread(f_phylop,data.table = F,stringsAsFactors = F)
    pos_center = which(phylop$V3==pos_of_interest)
    pos_center <- pos_center + center_shift
    low_pos_idx = pos_center - (flank) + 1
    hi_pos_idx = pos_center + (flank) + 1
    
    phylop_plot <- ggplot(phylop,aes(x=V3,y=V4)) +
      geom_hline(yintercept = 0,col='darkgrey',lty='dashed') +
      # geom_line() + 
      geom_smooth(span=0.05,se=F,method='loess',col="black") + 
      theme_bw() +
      theme(
        panel.grid = element_blank(),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        plot.margin = margin(t = 5, b = 5)
      ) + scale_y_continuous(breaks = scales::breaks_width(4)) +
      labs(x="Position",y="PhyloP") +
      scale_x_continuous(
        breaks = c(min(phylop$V3),max(phylop$V3)),
        labels = c(-1000,+1000)
      )
    
    
    # x-range for ribbon
    xmin <- phylop$V3[low_pos_idx]
    xmax <- phylop$V3[hi_pos_idx]
    
    # y-range for ribbon (full plot height)
    ymin <- min(phylop$V4, na.rm = TRUE)
    ymax <- max(phylop$V4, na.rm = TRUE)
    
    ribbon_df <- data.frame(
      xmin = xmin,
      xmax = xmax,
      ymin = ymin,
      ymax = ymax
    )
    
    phylop_plot = phylop_plot +
      geom_rect(
        data = ribbon_df,
        aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
        inherit.aes = FALSE,
        fill = "grey80",
        alpha = 0.4
      )
  }
  
  if (visualize_ism) {
    predictions = fread(f_scores,data.table = F,stringsAsFactors = F)
    predictions_max <- predictions %>%
      group_by(chr, pos) %>%
      slice_max(order_by = abs_logfc, n = 1, with_ties = FALSE) %>%
      ungroup()
    ism_plot = ggplot(predictions_max,aes(x=pos,y=abs_logfc)) + 
      geom_bar(stat='identity',fill='darkorange',col='white') +
      theme_bw() +
      labs(x = NULL, y = "Max |logFC| ") +
      theme(
        panel.grid = element_blank(),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        plot.margin = margin(t = 5, b = 5)#,
        # panel.border = element_rect(color = "black", fill = NA, linewidth = 1)
      ) + 
      scale_y_continuous(breaks = scales::breaks_width(1)) +
      xlim(pos_of_interest-flank,pos_of_interest+flank) +
      scale_x_continuous(
        breaks = c(pos_of_interest-flank,pos_of_interest+flank),
        labels = c(-flank,flank),
        limits = c(pos_of_interest-flank,pos_of_interest+flank)
      )
    if (variant_line) {
      ism_plot <- ism_plot + geom_vline(xintercept = pos_of_interest,lty='dashed',col='grey')
    }
    
  }
  
  # plot_title = index
  #––– 7) stitch them together & save
  if (visualize_ism) {
    final_plot <- p1 / p2 / p3 / ism_plot / phylop_plot +
      plot_layout(heights = c(1, 1, 1, 1, 1)) +
      plot_annotation(
        title=plot_title,
        theme = theme(
          plot.margin = margin(10,10,10,10),
          plot.title = element_text(hjust=0.5))
      )
  } else if (!visualize_phylop & !visualize_pred) {
    final_plot <- p2 / p3 +
      plot_layout(heights = c(1, 1, 1)) +
      plot_annotation(
        title=plot_title,
        theme = theme(
          plot.margin = margin(10,10,10,10),
          plot.title = element_text(hjust=0.5))
      )
  } else if (!visualize_phylop & visualize_pred) {
    final_plot <- p1 / p2 / p3 +
      plot_layout(heights = c(1, 1, 1)) +
      plot_annotation(
        title=gsub(":","_",df$variant_id[index]),
        theme = theme(
          plot.margin = margin(10,10,10,10),
          plot.title = element_text(hjust=0.5))
      )
  } else if (visualize_phylop & !visualize_pred) {
    final_plot <- p2 / p3 / phylop_plot +
      plot_layout(heights = c(1, 1, 1)) +
      plot_annotation(
        title=plot_title,
        theme = theme(
          plot.margin = margin(10,10,10,10),
          plot.title = element_text(hjust=0.5))
      )
  } else if (visualize_phylop & visualize_pred) {
    final_plot <- p1 / p2 / p3 / phylop_plot +
      plot_layout(heights = c(1, 1, 1, 1)) +
      plot_annotation(
        title=plot_title,
        theme = theme(
          plot.margin = margin(10,10,10,10),
          plot.title = element_text(hjust=0.5))
      )
  }
  
  if (print_plot) {
    print(final_plot)
  }
  
  # final_plot
  if (save_plot) {
    ggsave(
      filename = file.path(paste0(outdir,"/",gsub(":","_",predictions$variant_id[index]),".pdf")),
      plot     = final_plot,
      width    = 15, height = 8.2, units = "in"
    )
  }
}

print(final_plot)

##################################################

if (visualize_ism_heatmap) {
  
  heatmap_df <- predictions %>%
    mutate(pos_id = pos) %>%
    select(pos_id, allele = allele2, logfc) %>%
    complete(
      pos_id,
      allele = c("A", "C", "G", "T")
    )
  ism_heatmap <- ggplot(heatmap_df, aes(x = allele, y = pos_id, fill = logfc)) +
    geom_tile(color = "grey90") +
    scale_fill_gradient2(
      low = "blue",
      mid = "white",
      high = "red",
      midpoint = 0,
      na.value = "white",
      name = "logFC"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      # axis.text.y = element_text(size = 8),
      panel.grid = element_blank(),
      legend.position = 'top',
      axis.text.y = element_blank()
    ) +
    labs(
      x = "Allele",
      y = "Genomic position"
    ) 
  # + scale_y_continuous(breaks=c(min(heatmap_df$pos_id),max(heatmap_df$pos_id)))
  
  ggsave(
    filename = file.path(paste0(outdir,"/",ism_heatmap_prefix,".ism.pdf")),
    plot     = ism_heatmap,
    width    = 3*0.7, height = 10*0.7, units = "in"
  )
}
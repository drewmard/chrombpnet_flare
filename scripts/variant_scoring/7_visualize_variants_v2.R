library(data.table)
library(rhdf5)
library(reshape2)
library(ggplot2)
library(ggseqlogo)
library(patchwork)
library(dplyr)
library(tidyr)

# PARAMS
visualize_ism_heatmap=FALSE; ism_heatmap_prefix = ""
visualize_pred=TRUE
visualize_phylop=FALSE
visualize_ism=FALSE
save_plot=TRUE
print_plot=TRUE #TRUE
print_ism_plot=FALSE
variant_line = TRUE
reverse_complement = FALSE
num_repeated_runs = NULL # NULL
flank      <- 50 # for shap
flank2      <- 100 # 300 # for prediction
center_shift = 0

###########################

motif_disruption_score <- function(index, allele1_shap, allele2_shap, predictions, pos_of_interest) {
  
  ELF1  <- c(27100459, 27100468)
  GATA1 <- c(27100451, 27100454)
  RUNX1 <- c(27100440, 27100447)
  FLI1  <- c(27100432, 27100438)
  
  center_idx <- 1058
  
  motifs <- list(
    ELF1  = ELF1,
    GATA1 = GATA1,
    RUNX1 = RUNX1,
    FLI1  = FLI1
  )
  
  genomic_to_idx <- function(gr, pos_of_interest, center_idx) {
    (gr - pos_of_interest) + center_idx
  }
  
  get_motif_sums <- function(shap1, shap2, idx, motifs, pos_of_interest, center_idx) {
    out <- lapply(names(motifs), function(m) {
      gr <- motifs[[m]]
      rows <- genomic_to_idx(gr, pos_of_interest, center_idx)
      
      mat1 <- shap1[rows[1]:rows[2], , idx, drop = FALSE]
      mat2 <- shap2[rows[1]:rows[2], , idx, drop = FALSE]
      
      score <- abs(sum(mat2) / sum(mat1))
      
      data.frame(
        motif = m,
        disruption_score = score
      )
    })
    
    do.call(rbind, out)
  }
  
  res <- get_motif_sums(
    shap1 = allele1_shap,
    shap2 = allele2_shap,
    idx = index,
    motifs = motifs,
    pos_of_interest = pos_of_interest,
    center_idx = center_idx
  )
  
  out <- data.frame(
    variant_id = predictions$variant_id[index],
    t(res$disruption_score)
  )
  
  colnames(out)[-1] <- res$motif
  
  out
}
########################################

# Output:
outdir = "/Users/amarderstein/Library/Mobile Documents/com~apple~CloudDocs/Documents/Research/chrombpnet_flare/output/variant_shap/AoU_abbrev/Trisomy_Controls/DCs/fold_0/visualize_variants"

# SHAP:
f_shap = "/Users/amarderstein/Library/Mobile Documents/com~apple~CloudDocs/Documents/Research/chrombpnet_flare/output/variant_shap/AoU_abbrev/Trisomy_Controls/DCs/fold_0/DCs.fold_0.variant_shap.counts.h5"

# PRED:
f_pred = "/Users/amarderstein/Library/Mobile Documents/com~apple~CloudDocs/Documents/Research/chrombpnet_flare/output/variant_shap/AoU_abbrev/Trisomy_Controls/DCs/fold_0/DCs.fold_0.scores.variant_predictions.h5"

# ISM:
f_scores = "/Users/amarderstein/Library/Mobile Documents/com~apple~CloudDocs/Documents/Research/chrombpnet_flare/output/variant_shap/AoU_abbrev/Trisomy_Controls/DCs/fold_0/DCs.fold_0.scores.variant_scores.tsv"

###########
# pos_of_interest = 27100413; idx = 2655
# pos_of_interest = 27100444; idx = 3595
# pos_of_interest = 27100453; allele_of_interest = "A"
# pos_of_interest = 27100445; allele_of_interest = "C"
# idx=100
# pos_of_interest <- df$pos[idx]
plot_title = "AoU eQTL"
###################################################

# DATA LOAD:
dir.create(outdir)
predictions = fread(f_scores,data.table = F,stringsAsFactors = F)
alleles <- h5read(f_shap, "/alleles")
seq_arr <- h5read(f_shap, "projected_shap/seq")
allele1_pred_profiles <- h5read(f_pred, "/observed/allele1_pred_profiles")
allele2_pred_profiles <- h5read(f_pred, "/observed/allele2_pred_profiles")
df=fread(f_scores,data.table = F,stringsAsFactors = F);df[1,]

##############################################################

# INITIALIZE:
idx1 <- alleles == 0                     # TRUE/FALSE mask for allele1
idx2 <- alleles == 1                     # mask for allele2
allele1_shap <- seq_arr[,,idx1 , drop = FALSE]   # keeps all other dims
allele2_shap <- seq_arr[,,idx2 , drop = FALSE]
if (!is.null(num_repeated_runs)) {
  initial_instance = seq(1, nrow(predictions), by = num_repeated_runs)
  predictions <- predictions[initial_instance, ]
  df <- df[initial_instance,]
  allele1_pred_profiles <- allele1_pred_profiles[,initial_instance]
  allele2_pred_profiles <- allele2_pred_profiles[,initial_instance]
  n_rep <- 10
  n_var <- dim(allele1_shap)[3] / num_repeated_runs
  allele1_shap <- array(
    apply(
      array(allele1_shap, dim = c(dim(allele1_shap)[1],
                                  dim(allele1_shap)[2],
                                  n_rep,
                                  n_var)),
      c(1,2,4),
      mean
    ),
    dim = c(dim(allele1_shap)[1], dim(allele1_shap)[2], n_var)
  )
  allele2_shap <- array(
    apply(
      array(allele2_shap, dim = c(dim(allele2_shap)[1],
                                  dim(allele2_shap)[2],
                                  n_rep,
                                  n_var)),
      c(1,2,4),
      mean
    ),
    dim = c(dim(allele2_shap)[1], dim(allele2_shap)[2], n_var)
  )
  
}

idx <- which(predictions$pos == pos_of_interest)[1]
idx <- which(predictions$pos == pos_of_interest & predictions$allele2==allele_of_interest)[1]
# idx = ceiling((idx-1)/num_repeated_runs)

#########################
# RUN:

motif_disruption_score_results=list()
for (index in 1:nrow(predictions)) {
  print(index)
# for (index in 355:370) {
  pos_of_interest = predictions$pos[index]
  
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
  if (reverse_complement) {
    n <- ncol(mat_ref)
    mat_ref <- rbind(
      A = mat_ref["T", ],  # now shows T→A contributions
      C = mat_ref["G", ],  # G→C
      G = mat_ref["C", ],  # C→G
      T = mat_ref["A", ]   # A→T
    )[, n:1]              # then flip the columns
    mat_alt <- rbind(
      A = mat_alt["T", ],  # now shows T→A contributions
      C = mat_alt["G", ],  # G→C
      G = mat_alt["C", ],  # C→G
      T = mat_alt["A", ]   # A→T
    )[, n:1]              # then flip the columns
    # # x0 <- n - x0 + 1                       # mirrored x
  }
  
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
    
    # read 1D profile for this variant (index)
    # a1 <- h5read(f_pred, "/observed/allele1_pred_profiles")[, index]
    # a2 <- h5read(f_pred, "/observed/allele2_pred_profiles")[, index]
    a1 = allele1_pred_profiles[,index]
    a2 = allele2_pred_profiles[,index]
    
    pos_center <- 500
    
    low_pos_idx = pos_center - (flank2) + 1
    hi_pos_idx = pos_center + (flank2) + 1
    
    df_pred = data.frame(seq=1:length(seq(low_pos_idx:hi_pos_idx)),
                         a1=as.numeric(a1[low_pos_idx:hi_pos_idx]),
                         a2=a2[low_pos_idx:hi_pos_idx])
    
    n = nrow(df_pred)
    spanval=0.5
    
    # if (reverse_complement) {
    #   df_pred$seq = df_pred$seq[nrow(df_pred):1]
    # }
    if (reverse_complement) {
      df_pred <- df_pred[nrow(df_pred):1, ]
      df_pred$seq <- seq_len(nrow(df_pred))
    }
    
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
    
    # if (reverse_complement) {
    #   phylop$V3 = phylop$V3[nrow(phylop):1]
    # }
    if (reverse_complement) {
      phylop <- phylop[nrow(phylop):1, ]
    }
    
    # x-range for ribbon
    xmin <- phylop$V3[low_pos_idx]
    xmax <- phylop$V3[hi_pos_idx]
    
    # y-range for ribbon (full plot height)
    ymin <- min(phylop$V4, na.rm = TRUE)
    ymax <- max(phylop$V4, na.rm = TRUE)
    # lo <- loess(V4 ~ V3, data = phylop, span = 0.01)
    # smooth_pred <- predict(lo, newdata = phylop$V3)
    # ymin <- min(smooth_pred, na.rm = TRUE)
    # ymax <- max(smooth_pred, na.rm = TRUE)
    
    ribbon_df <- data.frame(
      xmin = xmin,
      xmax = xmax,
      ymin = ymin,
      ymax = ymax
    )
    
    height_tick <- ifelse(max(ribbon_df$ymax) > 2, floor(ribbon_df$ymax),
                          ifelse(max(ribbon_df$ymax) > 1, 1,
                                 ifelse(max(ribbon_df$ymax) > 0.5, 0.5,
                                        floor(max(ribbon_df$ymax) * 10) / 10)))
    
    phylop$x_axis_position <- (phylop$V3 - pos_of_interest) * ifelse(reverse_complement, -1, 1)
    ribbon_df$xmin_axis_position = (ribbon_df$xmin - pos_of_interest) * ifelse(reverse_complement, -1, 1)  
    ribbon_df$xmax_axis_position = (ribbon_df$xmax - pos_of_interest) * ifelse(reverse_complement, -1, 1)  
    
    phylop_plot <- ggplot(phylop,aes(x=x_axis_position,y=V4)) +
      geom_hline(yintercept = 0,col='darkgrey',lty='dashed') +
      geom_line() +
      # geom_smooth(span=0.01,se=F,method='loess',col="black") +
      theme_bw() +
      theme(
        panel.grid = element_blank(),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        plot.margin = margin(t = 5, b = 5)
      ) + scale_y_continuous(breaks = scales::breaks_width(height_tick)) +
      labs(x="Position",y="PhyloP") +
      scale_x_continuous(
        breaks = c(min(phylop$x_axis_position),max(phylop$x_axis_position)),
        labels = c(-1000,+1000)
      )
    
    phylop_plot = phylop_plot +
      geom_rect(
        data = ribbon_df,
        aes(xmin = xmin_axis_position, xmax = xmax_axis_position, ymin = ymin, ymax = ymax),
        inherit.aes = FALSE,
        fill = "grey80",
        alpha = 0.4
      )
    
    # phylop_plot <- ggplot(phylop,aes(x=V3,y=V4)) +
    #   geom_hline(yintercept = 0,col='darkgrey',lty='dashed') +
    #   # geom_line() + 
    #   geom_smooth(span=0.05,se=F,method='loess',col="black") + 
    #   theme_bw() +
    #   theme(
    #     panel.grid = element_blank(),
    #     axis.text.x = element_text(size = 14),
    #     axis.text.y = element_text(size = 14),
    #     plot.margin = margin(t = 5, b = 5)
    #   ) + scale_y_continuous(breaks = scales::breaks_width(height_tick)) +
    #   labs(x="Position",y="PhyloP") +
    #   scale_x_continuous(
    #     breaks = c(min(phylop$V3),max(phylop$V3)),
    #     labels = c(-1000,+1000)
    #   )
    # 
    # phylop_plot = phylop_plot +
    #   geom_rect(
    #     data = ribbon_df,
    #     aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
    #     inherit.aes = FALSE,
    #     fill = "grey80",
    #     alpha = 0.4
    #   )
  }
  
  if (visualize_ism) {
    # predictions = fread(f_scores,data.table = F,stringsAsFactors = F)
    predictions_max <- predictions %>%
      group_by(chr, pos) %>%
      slice_max(order_by = abs_logfc, n = 1, with_ties = FALSE) %>%
      ungroup()
    
    # if (reverse_complement) {
    #   predictions_max$pos = predictions_max$pos[nrow(predictions_max):1]
    # }
    
    if (reverse_complement) {
      predictions_max <- predictions_max[nrow(predictions_max):1, ]
    }
    
    predictions_max_sub = subset(predictions_max,pos >= (pos_of_interest-flank) & pos <= (pos_of_interest+flank))
    
    
    height_tick <- ifelse(max(predictions_max_sub$abs_logfc) > 2, 2,
                          ifelse(max(predictions_max_sub$abs_logfc) > 1, 1,
                                 ifelse(max(predictions_max_sub$abs_logfc) > 0.5, 0.5,
                                        floor(max(predictions_max_sub$abs_logfc) * 10) / 10)))
    
    # ism_plot = ggplot(predictions_max,aes(x=pos,y=abs_logfc)) + 
    #   geom_bar(stat='identity',fill='darkorange',col='white') +
    #   theme_bw() +
    #   labs(x = NULL, y = "Max |logFC| ") +
    #   theme(
    #     panel.grid = element_blank(),
    #     axis.text.x = element_text(size = 14),
    #     axis.text.y = element_text(size = 14),
    #     plot.margin = margin(t = 5, b = 5)#,
    #     # panel.border = element_rect(color = "black", fill = NA, linewidth = 1)
    #   ) + 
    #   # scale_y_continuous(breaks = scales::breaks_width(1)) +
    #   scale_y_continuous(breaks = scales::breaks_width(height_tick)) +
    #   xlim(pos_of_interest-flank,pos_of_interest+flank) +
    #   scale_x_continuous(
    #     breaks = c(pos_of_interest-flank,pos_of_interest+flank),
    #     labels = c(-flank,flank),
    #     limits = c(pos_of_interest-flank,pos_of_interest+flank)
    #   )
    # predictions_max_sub$x_axis_position <- predictions_max_sub$pos - pos_of_interest
    predictions_max_sub$x_axis_position <- (predictions_max_sub$pos - pos_of_interest) * ifelse(reverse_complement, -1, 1)
    
    ism_plot = ggplot(predictions_max_sub,aes(x=x_axis_position,y=abs_logfc)) + 
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
      # scale_y_continuous(breaks = scales::breaks_width(1)) +
      scale_y_continuous(breaks = scales::breaks_width(height_tick)) +
      # xlim(pos_of_interest-flank,pos_of_interest+flank) +
      scale_x_continuous(
        breaks = c(0-flank,0+flank),
        labels = c(-flank,flank),
        limits = c(0-flank,0+flank)
      )
    if (variant_line) {
      # ism_plot <- ism_plot + geom_vline(xintercept = pos_of_interest,lty='dashed',col='grey')
      ism_plot <- ism_plot + geom_vline(xintercept = 0,lty='dashed',col='grey')
    }
    
  }
  
  # plot_title = index
  #––– 7) stitch them together & save
  if (visualize_ism & visualize_phylop & visualize_pred) {
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
  
  if (print_ism_plot) {
    print(ism_heatmap)
  }
  
  if (save_plot) {
    ggsave(
      filename = file.path(paste0(outdir,"/",ism_heatmap_prefix,".ism.pdf")),
      plot     = ism_heatmap,
      width    = 3*0.7, height = 10*0.7, units = "in"
    )
  }
}

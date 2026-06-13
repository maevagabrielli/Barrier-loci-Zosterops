# ============================================================
# Pairwise circular bootstrap overlap
# ============================================================

library(data.table)
library(ggplot2)

# -----------------------------
# 1. DATA
# -----------------------------

fst   <- as.data.table(outfst)
sweep <- as.data.table(outsweep)
dils  <- as.data.table(outdils)

setnames(fst,   c("chr","start","end","method"))
setnames(sweep, c("chr","start","end","method"))
setnames(dils,  c("chr","start","end","method"))

fst[,   `:=`(start = as.numeric(start), end = as.numeric(end))]
sweep[, `:=`(start = as.numeric(start), end = as.numeric(end))]
dils[,  `:=`(start = as.numeric(start), end = as.numeric(end))]

# -----------------------------
# 2. CHROMOSOME LENGTHS
# -----------------------------

agp <- fread("DeCoSTAR_27avian_ADseq+scaff_Boltz_kT0.1_Lin0.1_M2_Zosterops_borbonicus_ZeFi-ZoBo+DeCoSTAR.agp",
             header = FALSE)

agp <- agp[, .(chr = V1,
               start = as.numeric(V2),
               end   = as.numeric(V3))]

chr_lengths <- agp[grepl("^chromosome_", chr),
                   .(length = max(end, na.rm = TRUE)),
                   by = chr]

chr_len_vec <- setNames(chr_lengths$length, chr_lengths$chr)

# -----------------------------
# 3. OVERLAP FUNCTION
# -----------------------------

overlap_two <- function(a, b) {
  
  if (nrow(a) == 0 || nrow(b) == 0) return(0)
  
  total <- 0
  
  for (i in seq_len(nrow(a))) {
    
    ai <- a[i]
    
    hits <- b[
      chr == ai$chr &
        end > ai$start &
        start < ai$end
    ]
    
    if (nrow(hits) == 0) next
    
    total <- total + sum(
      pmax(0,
           pmin(ai$end, hits$end) -
             pmax(ai$start, hits$start))
    )
  }
  
  total
}

# -----------------------------
# 4. CIRCULAR SHIFT
# -----------------------------

circular_shift <- function(df, shift, L) {
  
  if (nrow(df) == 0) return(df)
  
  df <- copy(df)
  
  df[, start := (start + shift) %% L]
  df[, end   := (end + shift) %% L]
  
  df[, `:=`(
    start = pmin(start, end),
    end   = pmax(start, end)
  )]
  
  df
}

# -----------------------------
# 5. PAIRWISE TEST FUNCTION
# -----------------------------

pairwise_test <- function(a, b, label_a, label_b, B = 1000) {
  
  obs <- overlap_two(a, b)
  boot <- numeric(B)
  
  for (i in seq_len(B)) {
    
    a_p <- copy(a)
    
    for (chr_name in names(chr_len_vec)) {
      
      L <- chr_len_vec[[chr_name]]
      shift <- sample(0:L, 1)
      
      idx <- a_p$chr == chr_name
      if (!any(idx)) next
      
      a_p[idx] <- circular_shift(a_p[idx], shift, L)
    }
    
    boot[i] <- overlap_two(a_p, b)
  }
  
  list(
    label = paste(label_a, "vs", label_b),
    obs   = obs,
    p     = mean(boot >= obs),
    boot  = boot
  )
}

# -----------------------------
# 6. RUN ANALYSES
# -----------------------------

res_fst_sweep <- pairwise_test(
  fst[method == "fst"],
  sweep[method == "sweep"],
  "FST", "Sweep"
)

res_fst_dils <- pairwise_test(
  fst[method == "fst"],
  dils[method == "dils"],
  "FST", "DILS"
)

res_sweep_dils <- pairwise_test(
  sweep[method == "sweep"],
  dils[method == "dils"],
  "Sweep", "DILS"
)

# -----------------------------
# 7. SUMMARY TABLE
# -----------------------------

results <- data.table(
  pair = c(res_fst_sweep$label,
           res_fst_dils$label,
           res_sweep_dils$label),
  
  obs  = c(res_fst_sweep$obs,
           res_fst_dils$obs,
           res_sweep_dils$obs),
  
  p    = c(res_fst_sweep$p,
           res_fst_dils$p,
           res_sweep_dils$p)
)

# -----------------------------
# EFFECT SIZE METRICS
# -----------------------------

results[, mean_null := c(
  mean(res_fst_sweep$boot),
  mean(res_fst_dils$boot),
  mean(res_sweep_dils$boot)
)]

results[, sd_null := c(
  sd(res_fst_sweep$boot),
  sd(res_fst_dils$boot),
  sd(res_sweep_dils$boot)
)]

# Fold enrichment
results[, fold_enrichment := obs / mean_null]

# Z-score (standardized effect size)
results[, z_score := (obs - mean_null) / sd_null]

print(results)

results[, label := paste0(
  "p = ", signif(p, 3),
  "\nZ = ", round(z_score, 2),
  "\nFold = ", round(fold_enrichment, 2)
)]

# -----------------------------
# 8. FIGURE
# -----------------------------

plot_df <- rbind(
  data.table(x = res_fst_sweep$boot,  pair = "FST vs Sweep"),
  data.table(x = res_fst_dils$boot,   pair = "FST vs DILS"),
  data.table(x = res_sweep_dils$boot, pair = "Sweep vs DILS")
)

obs_df <- data.frame(
  pair = c("FST vs Sweep",
           "FST vs DILS",
           "Sweep vs DILS"),
  
  obs = c(res_fst_sweep$obs,
          res_fst_dils$obs,
          res_sweep_dils$obs)
)

p <- ggplot(plot_df, aes(x = x)) +
  
  geom_histogram(aes(y = after_stat(density)),
                 bins = 40,
                 fill = "grey80",
                 color = "black") +
  
  geom_density(color = "black", linewidth = 0.6) +
  
  facet_wrap(~pair, scales = "free") +
  
  # observed values
  geom_vline(
    data = obs_df,
    aes(xintercept = obs),
    color = "red",
    linewidth = 1.2
  ) +
  
  # mean null
  geom_vline(
    data = results,
    aes(xintercept = mean_null),
    color = "blue",
    linetype = "dashed",
    linewidth = 1
  ) +
  
  theme_classic(base_size = 14) +
  
  labs(
    title = "Circular bootstrap overlap between genomic outlier sets",
    x = "Overlap (bp)",
    y = "Density (null)"
  )

p <- p +
  geom_text(
    data = results,
    aes(x = Inf, y = Inf, label = label),
    hjust = 1.1, vjust = 1.1,
    size = 3
  )

print(p)

ggsave(
  filename = "circular_bootstrap_overlap_TREE_pairwise.png",
  plot = p,
  width = 10,
  height = 6,
  dpi = 300
)


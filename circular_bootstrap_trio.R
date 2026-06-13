### read outlier files

setwd("/Users/mgabrielli/Documents/Doctorat/reviewZosteropsEvolution/")

library(ggplot2)
library(dplyr)

outNEPN=read.table("pos.island.highrec.NE_PN.merged.chr.bed")
outSNE=read.table("pos.island.highrec.S_NE.merged.chr.bed")
outWNE=read.table("pos.island.highrec.W_NE.merged.chr.bed")
outWPN=read.table("pos.island.highrec.W_PN.merged.chr.bed")
outWS=read.table("pos.island.highrec.W_S.merged.chr.bed")

outfst=bind_rows(outNEPN,outSNE,outWNE,outWPN,outWS)

outHN=read.table("Out.HN.SNP.top.100kb.bed")
outHS=read.table("Out.HS.SNP.top.100kb.bed")
outNE=read.table("Out.NE.SNP.top.100kb.bed")
outS=read.table("Out.S.SNP.top.100kb.bed")
outW=read.table("Out.W.SNP.top.100kb.bed")

outsweep=bind_rows(outHN,outHS,outNE,outS,outW)
outsweep=outsweep[,c(1:3)]
outdils=read.table("pos.barriers.p0.95.window100kb.chr.bed")

# Add method labels
outfst$method <- "fst"
outsweep$method <- "sweep"
outdils$method <- "dils"

# ============================================================
# Circularization bootstrap: overlap among DILS / sweep / FST
# (Heidbreder-style chromosome-level circular permutation)
# ============================================================

library(data.table)
library(ggplot2)

# ============================================================
# 1. INPUT DATA
# ============================================================

fst   <- as.data.table(outfst)
sweep <- as.data.table(outsweep)
dils  <- as.data.table(outdils)

setnames(fst,   c("chr", "start", "end", "method"))
setnames(sweep, c("chr", "start", "end", "method"))
setnames(dils,  c("chr", "start", "end", "method"))

# ensure numeric coordinates
fst[,   `:=`(start = as.numeric(start), end = as.numeric(end))]
sweep[, `:=`(start = as.numeric(start), end = as.numeric(end))]
dils[,  `:=`(start = as.numeric(start), end = as.numeric(end))]

# ============================================================
# 2. CHROMOSOME LENGTHS
# ============================================================

agp <- fread("DeCoSTAR_27avian_ADseq+scaff_Boltz_kT0.1_Lin0.1_M2_Zosterops_borbonicus_ZeFi-ZoBo+DeCoSTAR.agp",
             header = FALSE)

agp <- agp[, .(chr = V1, start = as.numeric(V2), end = as.numeric(V3))]

chr_lengths <- agp[grepl("^chromosome_", chr),
                   .(length = max(end, na.rm = TRUE)),
                   by = chr]

chr_len_vec <- setNames(chr_lengths$length, chr_lengths$chr)

# ============================================================
# 3. OVERLAP FUNCTION (BEDTOOLS-EQUIVALENT)
# ============================================================

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
    
    overlap_start <- pmax(ai$start, hits$start)
    overlap_end   <- pmin(ai$end,   hits$end)
    
    total <- total + sum(pmax(0, overlap_end - overlap_start))
  }
  
  total
}

obs_overlap <- function(d1, d2, d3) {
  overlap_two(d1, d2) +
    overlap_two(d1, d3) +
    overlap_two(d2, d3)
}

# ============================================================
# 4. OBSERVED VALUE
# ============================================================

obs <- obs_overlap(
  fst[method == "fst"],
  sweep[method == "sweep"],
  dils[method == "dils"]
)

cat("Observed overlap (bp):", obs, "\n")

# ============================================================
# 5. CIRCULAR SHIFT FUNCTION
# ============================================================

circular_shift <- function(df, shift, L) {
  
  if (nrow(df) == 0) return(df)
  
  df <- copy(df)
  
  df[, start := (start + shift) %% L]
  df[, end   := (end   + shift) %% L]
  
  # reorder after wrap-around
  df[, `:=`(
    start = pmin(start, end),
    end   = pmax(start, end)
  )]
  
  df
}

# ============================================================
# 6. CIRCULARIZATION BOOTSTRAP
# ============================================================

set.seed(1)

B <- 1000
boot <- numeric(B)

d1 <- fst[method == "fst"]
d2 <- sweep[method == "sweep"]
d3 <- dils[method == "dils"]

for (i in seq_len(B)) {
  
  boot_sum <- 0
  
  for (chr_name in names(chr_len_vec)) {
    
    L <- chr_len_vec[[chr_name]]
    shift <- sample(0:L, 1)
    
    a <- circular_shift(d1[chr == chr_name], shift, L)
    b <- circular_shift(d2[chr == chr_name], shift, L)
    c <- circular_shift(d3[chr == chr_name], shift, L)
    
    if (nrow(a) == 0 && nrow(b) == 0 && nrow(c) == 0) next
    
    boot_sum <- boot_sum +
      overlap_two(a, b) +
      overlap_two(a, c) +
      overlap_two(b, c)
  }
  
  boot[i] <- boot_sum
}

# ============================================================
# 7. STATISTICAL TEST
# ============================================================

p_enrich <- mean(boot >= obs)

cat("P-value (enrichment):", p_enrich, "\n")

# ============================================================
# 8. FIGURE
# ============================================================

df_plot <- data.frame(overlap = boot)

p <- ggplot(df_plot, aes(x = overlap)) +
  geom_histogram(bins = 40,
                 fill = "grey70",
                 color = "black") +
  geom_vline(xintercept = obs,
             color = "red",
             linewidth = 1.2) +
  annotate("text",
           x = obs,
           y = Inf,
           label = "Observed",
           vjust = 2,
           color = "red") +
  annotate("text",
    x = Inf,
    y = Inf,
    hjust = 1.1,
    vjust = 1.5,
    label = paste0(
    "P = ", signif(p_enrich, 3)
    ),
    size = 5) +
  theme_classic(base_size = 14) +
  labs(
    title = "Circularization bootstrap of genomic overlap",
    subtitle = "DILS barriers ∩ sweep regions ∩ FST outliers",
    x = "Total overlap (bp)",
    y = "Frequency (null distribution)"
  )

print(p)

ggsave("circular_bootstrap_overlap_TREE.png",
       p, width = 10, height = 6, dpi = 300)

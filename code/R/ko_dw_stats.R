# Statistical analysis to investigate differences on dry weight between wild type
# and knock-out lines growth at 17 °C

library("readxl")
library("car")
library("emmeans")
library("multcomp")
library("lme4")
library("nlme")

### ---------------- Data import and cleaning ----------------

# read data from Excel spreadsheet
top_dir <- "AraTModel/dry-weights-tdna-lines"
rgr_file_name <- paste(top_dir, "rgr_measurements_20231005.txt",
                       sep = "/")
data <- read.table(rgr_file_name, header = TRUE, sep = "\t")

# remove samples with comment (whole pot)
idx_comment <- data$comment != ""
# do not exclude samples that occur in more than one batch / experiment
# idx_comment[data$comment=="Philipp: not considered, repeated in dataRAW4th"] = F
data <- data[!idx_comment,]

# read additional table with update predictions (expected phenotypes)
si_file_name <- paste(top_dir, "si_table_ko_predictions.xlsx", sep = "/")
effect_data <- read_xlsx(si_file_name, skip = 3)

# update predicted effects on RGR
genes_uniq <- unique(data$Gene)
for (i in 1:length(genes_uniq)) {
  gene_idx_data <- data$Gene %in% genes_uniq[i]
  gene_idx_si <- effect_data$`Gene ID` %in% genes_uniq[i]
  gene_effect <- unique(effect_data$`expected phenotype`[gene_idx_si])
  if (length(gene_effect) > 1) {
    if ("Reduced growth at 17 °C" %in% gene_effect &&
        "negative control" %in% gene_effect) {
      # if knock-out lines has both effect and no effect predicted, assign 
      # "Reduced growth at 17 °C"
        gene_effect <- "Reduced growth at 17 °C"
    }
  } else if (genes_uniq[i] == "WT") {
    gene_effect <- "WT"
  }
  if (length(gene_effect) == 0) {
    # some KO lines have no effect because the updated model did no longer
    # predict an effect or no effect on RGR for them 
    gene_effect <- ""
  }
  data$effect[gene_idx_data] <- gene_effect
}

# reduce table to be unique with respect to Set, Gene, and DW_final
data <- data[!duplicated(data[,c(1, 2, 6)]),]

# remove genes with less than three replicates
n_reps <- vector(mode = "double", length = nrow(data))
for (i in 1:length(genes_uniq)) {
  gene_idx <- data$Gene %in% genes_uniq[i]
  n_reps[gene_idx] <- length(data$Gene[gene_idx])
}
data <- data[n_reps>2,]

# remove genes without predicted effect
has_effect <- data$effect != ""
data <- data[has_effect,]

### ---------------- type II ANOVA ----------------

# check if data set is unbalanced
with(data, table(Gene, Set))

# perform 2-way ANOVA
print("Gene * Set")
anova2.results.gs <- aov(DW_final ~ Gene * Set, data = data)
summary(anova2.results.gs)

print("Set * Gene")
anova2.results.sg <- aov(DW_final ~ Set * Gene, data = data)
summary(anova2.results.sg)

# post-hoc test
p_gs <- TukeyHSD(anova2.results.gs, which = "Gene")
p_sg <- TukeyHSD(anova2.results.sg, which = "Gene")

# significant differences
rownames(p_gs$Gene)[p_gs$Gene[,4] < 0.05]
rownames(p_sg$Gene)[p_sg$Gene[,4] < 0.05]

# type II ANOVA
anova22.results <- Anova(anova2.results.gs, type = 2)
anova22.results

### ---------------- Linear mixed-effects models ----------------

# change type of "Gene" to factor to enable use of glht on lme model
data$Gene <- factor(data$Gene,
                    levels = unique(data$Gene[order(data$effect,
                                                    decreasing = T)]))

# Model 1: random intercepts and fixed slopes with lmer function
m1 <- lmer(DW_final ~ Gene + (1|Set) - 1, data = data)
Anova(m1)
# Model 2: random intercepts and fixed slopes with lme function
m2 <- lme(DW_final ~ Gene - 1, random = ~1|Set, data = data)
Anova(m2)

# contrast definitions for glht for pairwise comparisons of WT against KOs
contr_glht <- unlist(lapply(levels(data$Gene)[-1], function(x) {paste0(x, " - WT == 0")}))

# contrast definitions for emmeans for pairwise comparisons of WT against KOs
contr_mat <- matrix(rep(0, length(levels(data$Gene))*(length(levels(data$Gene))-1)),
                    nrow = length(levels(data$Gene)))
contr_mat[1,] <- -1
diag(contr_mat[2:nrow(contr_mat),]) <- 1
contr_emm = list()
for (i in 1:ncol(contr_mat)) {
  contr_emm[i] = list(contr_mat[,i])
}

# define MHT correction method
mhtcorr_method <- "BH"

# --- Model 1 contrasts

# comparisons using glht function
comp1 <- glht(m1, mcp(Gene = contr_glht))

# without MHT correction
comp1.summary <- summary(comp1, test = adjusted("none"))
comp1.pvalues <- as.numeric(comp1.summary$test$pvalues)

# with MHT correction
comp1.adj.summary <- summary(comp1, test = adjusted(mhtcorr_method))
comp1.adj.pvalues <- as.numeric(comp1.adj.summary$test$pvalues)

# comparisons using emmeans function
comp2 <- emmeans(m1, pairwise ~ Gene)

# without MHT correction
comp2.pvalues <- as.data.frame(
  contrast(comp2, adjust = "none", method = contr_emm))$p.value

# with MHT correction
comp2.adj.pvalues <- as.data.frame(
  contrast(comp2, adjust = mhtcorr_method, method = contr_emm))$p.value

# --- Model 2 contrasts

# comparisons using glht function
comp3 <- glht(m2, mcp(Gene = contr_glht))

# without MHT correction
comp3.summary <- summary(comp3, test = adjusted("none"))
comp3.pvalues <- as.numeric(comp3.summary$test$pvalues)

# with MHT correction
comp3.adj.summary <- summary(comp3, test = adjusted(mhtcorr_method))
comp3.adj.pvalues <- as.numeric(comp3.adj.summary$test$pvalues)

# comparisons using emmeans function
comp4 <- emmeans(m2, pairwise ~ Gene)

# without MHT correction
comp4.pvalues <- as.data.frame(
  contrast(comp4, adjust = "none", method = contr_emm))$p.value

# with MHT correction
comp4.adj.pvalues <- as.data.frame(
  contrast(comp4, adjust = mhtcorr_method, method = contr_emm))$p.value

# combine p-values into data frame
p_values_combined <- cbind(comp1.pvalues,
                           comp1.adj.pvalues,
                           comp2.pvalues,
                           comp2.adj.pvalues,
                           comp3.pvalues,
                           comp3.adj.pvalues,
                           comp4.pvalues,
                           comp4.adj.pvalues)
# construct data frame with p-values and log-likelihoods
logl.reml <- unlist(lapply(list(m1, m1, m1, m1, m2, m2, m2, m2),
                     function(x) {logLik(x, REML = T)}))
logl <- unlist(lapply(list(m1, m1, m1, m1, m2, m2, m2, m2),
                     function(x) {logLik(x, REML = F)}))

result_df <- data.frame(rbind(p_values_combined, logl, logl.reml),
                          row.names = c(levels(data$Gene)[2:length(levels(data$Gene))],
                                        "log-li.", "log-li. (REML)"))
names(result_df) <- c("lmer_glht", "lmer_glht_adj", 
                        "lmer_emmeans", "lmer_emmeans_adj",
                        "lme_glht", "lme_glht_adj",
                        "lme_emmeans", "lme_emmeans_adj")

# write data frame to file
out_file = paste(top_dir, "p_values_lme.csv", sep = "/")
write.table(result_df, out_file, sep = ",")

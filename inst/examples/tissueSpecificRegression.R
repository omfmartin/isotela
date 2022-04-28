nsamples <- 500
ngenes <- 150
ntissueuniq <- 50
nercc <- 50
ntissues <- 2

sample_names <- paste0("Sample", seq_len(nsamples))
gene_names <- paste0("Gene", seq_len(ngenes))
ercc_names <- paste0("ERCC", seq_len(nercc))
tissue_names <- paste0("Tissue", seq_len(ntissues))

# simulate gene mean
mu_gene <- rlnorm(ngenes, meanlog = 5)
mu_ercc <- rlnorm(ngenes, meanlog = 5)

# simulate size and technical factors
size_factors <- sapply(seq_len(ntissues), function(i) {
  k <- rlnorm(nsamples, 0, 1)
  k/mean(k)
})
rownames(size_factors) <- sample_names
colnames(size_factors) <- tissue_names

tech_factors <- rlnorm(nsamples, 0, 1)
names(tech_factors) <- sample_names
tech_factors <- tech_factors/mean(tech_factors)

# genes are from which tissue
tissues <- sample(tissue_names, size = ntissueuniq, replace = TRUE)
names(tissues) <- sample(gene_names, ntissueuniq)

# simulate tissue proportions
prop <- sapply(seq_len(ntissues), function(i) rnorm(ngenes, 0, 1))
prop <- t(apply(prop, 1, function(x) x^2/sum(x^2)))
rownames(prop) <- gene_names
colnames(prop) <- tissue_names
for (i in tissue_names) {
  prop[names(which(tissues == i)), i] <- 1
  prop[names(which(tissues == i)), ! colnames(prop) %in% i] <- 0
}

# simulate over-dispersion
disp <- rlnorm(ngenes, -1)
names(disp) <- gene_names

# simulate counts and ERCC matrices
counts <- t(sapply(seq_len(ngenes), function(i) {
  eta <- exp(log(mu_gene[i]) + log(tech_factors) + log(size_factors %*% prop[i, ]))
  rnbinom(nsamples, mu = eta, size = 1/disp[i])
}))
colnames(counts) <- sample_names
rownames(counts) <- gene_names

ercc <- t(sapply(seq_len(nercc), function(i) {
  rpois(nsamples, exp(log(mu_gene[i]) + log(tech_factors)))
}))
colnames(ercc) <- sample_names
rownames(ercc) <- ercc_names

annots <- data.frame(Sample = sample_names)

# run standard normalization
nf <- normalizationFactor(counts)
norm <- normalizeCounts(counts, nf)

# run tissue-specific normalization
nf_tissues <- normalizationFactorTissueSpecific(counts, tissues)
norm_tissues <- normalizeCountsTissueSpecific(counts, nf_tissues, tissues)

# run regression
res <- tissueSpecificRegression(counts = counts,
                                ercc = ercc,
                                annots = annots,
                                tissues = tissues,
                                mean_design = ~1,
                                wald_test = FALSE,
                                cores = 1)

# check we recovered parameters
plot(exp(res$coefficients), mu_gene, log="xy")
abline(b=1, a=0)

plot(res$dispersions, disp, log="xy")
abline(b=1, a=0)

plot(res$tech_factors, tech_factors, log="xy")
abline(b=1, a=0)

plot(res$proportions[, 1], prop[, 1])
abline(b=1, a=0)


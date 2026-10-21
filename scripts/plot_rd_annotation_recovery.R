#!/usr/bin/env Rscript
# Plot the output of simulate_rd_annotations.py without additional R packages.
args <- commandArgs(trailingOnly=TRUE)
stopifnot(length(args) == 1)
folder <- args[1]
recovery <- read.table(file.path(folder, "annotation_recovery.tsv"), header=TRUE)
truth <- read.table(file.path(folder, "annotation_truth.tsv"), header=TRUE)
models <- unique(recovery$Model)
levels <- paste0("p", 2:5)
colors <- c("#828a94", "#9bc2e5", "#4989c3", "#e29b35", "#bf483b")
rc_generated <- any(truth$Component[truth$Annotation != "Intercept"] == "normal")
render <- function() {
  par(mfrow=c(length(models), 4), mar=c(4,4,2.5,1), oma=c(2,0,1,0))
  for (model in models) for (level in levels) {
    d <- subset(recovery, Model == model & Level == level)
    all <- subset(recovery, Level == level)
    limits <- range(c(all$TrueEffect, all$Mean))
    limits <- limits + c(-1,1)*diff(limits)*0.08
    t <- subset(truth, Level == level & Annotation != "Intercept")
    point_colors <- if (rc_generated)
      ifelse(t$Component[match(d$Annotation, t$Annotation)] == "normal", "#6655aa", "#828a94") else
      colors[as.integer(t$Component[match(d$Annotation, t$Annotation)])+1]
    plot(d$TrueEffect, d$Mean, xlim=limits, ylim=limits, asp=1,
         pch=16, col=point_colors,
         xlab="True annotation effect", ylab="Posterior mean",
         main=paste(model, level, if(max(d$SplitRhat) > 1.05) "(poor mixing)" else ""), cex=0.9)
    abline(0,1,lty=2,col="#444444")
    abline(h=0,v=0,lty=3,col="#dddddd")
    mtext(sprintf("RMSE %.3f", sqrt(mean((d$TrueEffect-d$Mean)^2))),
          side=3, line=0.2, cex=0.75, adj=0)
  }
  description <- if (rc_generated)
    paste("RC-generated effects; shared Gaussian variances p2-p5:",
          paste(sapply(levels, function(level) unique(truth$Variance[truth$Level == level & truth$Component == "normal"])), collapse=", "),
          "(grey: forced nulls)") else
    "Generating annotation variances: grey 0; light blue 0.0001; blue 0.001; orange 0.01; red 0.1"
  mtext(description,
        side=1, outer=TRUE, line=0.5, cex=0.8)
}
pdf(file.path(folder, "annotation_recovery.pdf"), width=14, height=3.7*length(models))
render()
dev.off()
png(file.path(folder, "annotation_recovery.png"), width=2100, height=560*length(models), res=150)
render()
dev.off()

#!/usr/bin/env Rscript
# Plot the any-level inclusion probabilities from check_annotation_pip_calibration.py.
args <- commandArgs(trailingOnly=TRUE)
stopifnot(length(args) >= 1)
for (folder in args) {
  d <- read.table(file.path(folder, "annotation_pip_calibration.tsv"), header=TRUE)
  render <- function() {
    par(mfrow=c(1,2), mar=c(4,4,3,1))
    boxplot(RDAnyPIP ~ factor(TrueAnyNonzero, levels=c(0,1), labels=c("Wholly null", "Causal")),
            data=d, ylim=c(0,1), ylab="RD posterior P(any effect nonzero)",
            col=c("#cccccc", "#acd5e9"), main="Inclusion probabilities by truth")
    abline(h=.95, lty=2, col="#c04c37")
    legend("bottomright", "PIP threshold 0.95", lty=2, col="#c04c37", bty="n")
    plot(c(0,1), c(0,1), type="n", xlim=c(0,1), ylim=c(0,1),
         xlab="Mean predicted inclusion probability", ylab="Observed causal fraction",
         main="Reliability by PIP bin")
    abline(0,1,lty=2,col="#666666")
    edges <- c(0,.5,.8,.9,.95,.99,1.00001)
    for (i in seq_len(length(edges)-1)) {
      selected <- d$RDAnyPIP >= edges[i] & d$RDAnyPIP < edges[i+1]
      if (any(selected)) {
        x <- mean(d$RDAnyPIP[selected]); y <- mean(d$TrueAnyNonzero[selected])
        points(x,y,pch=16,col="#2585aa",cex=1.2)
        text(x,y,paste0("n=",sum(selected)),pos=2,cex=.75,col="#2585aa")
      }
    }
    points(1,mean(d$TrueAnyNonzero),pch=17,col="#c04c37",cex=1.2)
    legend("topleft",c("RD any-level PIP", "RC implied PIP=1"),
           pch=c(16,17),col=c("#2585aa","#c04c37"),bty="n")
  }
  pdf(file.path(folder,"pip_calibration.pdf"),width=10,height=5)
  render(); dev.off()
  png(file.path(folder,"pip_calibration.png"),width=1500,height=750,res=150)
  render(); dev.off()
}

#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly=TRUE)
stopifnot(length(args)>0)
for (folder in args) {
  d <- read.table(file.path(folder,"rd_joint_calibration.tsv"),header=TRUE)
  render <- function() {
    par(mfrow=c(1,2),mar=c(4,4,3,1))
    b <- rbind(data.frame(PIP=d$IndependentAnyPIP,Model="Independent RD",Truth=d$TrueAnyNonzero),
               data.frame(PIP=d$JointPIP,Model="Joint RD",Truth=d$TrueAnyNonzero))
    b$Group <- factor(paste(b$Model,ifelse(b$Truth==0,"null","causal")),
                      levels=c("Independent RD null","Joint RD null","Independent RD causal","Joint RD causal"))
    boxplot(PIP~Group,data=b,ylim=c(0,1),names=c("Indep\nnull","Joint\nnull","Indep\ncausal","Joint\ncausal"),
            col=c("#aaaaaa","#2585aa","#aaaaaa","#2585aa"),ylab="Annotation inclusion probability",main="PIPs by annotation truth")
    abline(h=.95,lty=2,col="#c04c37")
    plot(c(0,1),c(0,1),type="n",xlim=c(0,1),ylim=c(0,1),xlab="Mean predicted inclusion probability",
         ylab="Observed causal fraction",main="Reliability by PIP bin")
    abline(0,1,lty=2,col="#666666")
    edges <- c(0,.5,.8,.9,.95,.99,1.00001)
    for (model in c("Independent RD","Joint RD")) {
      v <- b[b$Model==model,]
      color <- if(model=="Joint RD") "#2585aa" else "#888888"
      for (i in seq_len(length(edges)-1)) {
        chosen <- v$PIP>=edges[i] & v$PIP<edges[i+1]
        if(any(chosen)) {
          x <- mean(v$PIP[chosen]); y <- mean(v$Truth[chosen])
          points(x,y,pch=16,col=color)
          text(x,y,paste0("n=",sum(chosen)),pos=2,cex=.7,col=color)
        }
      }
    }
    legend("topleft",c("Independent RD","Joint RD"),col=c("#888888","#2585aa"),pch=16,bty="n")
  }
  pdf(file.path(folder,"rd_joint_calibration.pdf"),width=10,height=5);render();dev.off()
  png(file.path(folder,"rd_joint_calibration.png"),width=1500,height=750,res=150);render();dev.off()
}

# script to make the legend for manhattan plots for manuscript
# run this, then cut/paste the legend from here into the otehr plots
# trust me its easier this way

x<-1:10 
y1=x*x 
y2=2*y1

png("legend.png")
plot(x, y1, type="b", pch=19, col="red", xlab="x", ylab="y")
# Add a line
lines(x, y2, pch=18, col="blue", type="b", lty=2)
# Add a legend


legend(1, 95, legend=c("Novel", "Replicated",  "Not Replicated"),
      col =  c("black", "black", "black"), 
       pt.bg=c("dodgerblue","green2",  "red"), pch = c(22, 21, 23), cex=1.5)
dev.off()
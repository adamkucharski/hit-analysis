# Vaccine analysis


# Load functions -------------------------------------------------------------------

library(tidyverse)

# setwd("~/Documents/GitHub/hit-analysis/")

# Load data
data_vacc <- read_csv("data/vaccine_data.csv")
data_vacc$effectiveness <- data_vacc$effectiveness/100



# Plots -------------------------------------------------------------------

par(mfcol=c(2,1),mar=c(4,4,1,1),mgp=c(2,0.6,0),las=0)



# Plot herd immunity curve

plot(0,0,xlim=c(1,20),ylim=c(0,1),type="l",col="blue",xlab="R0",
     ylab="herd immunity threshold",yaxs="i",xaxs="i")

rr_seq <- seq(1,20,0.1)

lines(rr_seq,1-1/rr_seq)
points(data_vacc$r0, data_vacc$effectiveness,col="blue")
text(x=data_vacc$r0+0.5,y=data_vacc$effectiveness,labels=data_vacc$Pathogen,cex=0.8,adj=0,col="blue")

title(LETTERS[1],adj=0)

# Plot proportion required for SARS-CoV-2
plot(0,0,xlim=c(0,1),ylim=c(0,1),type="l",col="blue",xlab="proportion immune following natural SARS-CoV-2 infection",
     ylab="% vaccine coverage required for HI",yaxs="i",xaxs="i")

for(ii in 1:3){
  for(jj in 1){
  
  if(ii==1){v_eff <- 0.8; col_pick <- "dark green"}
  if(ii==2){v_eff <- 0.7; col_pick <- "blue"}
  if(ii==3){v_eff <- 0.6; col_pick <- "purple"}
  
  xx <- seq(0,1,0.01)

  if(jj==1){r0 <- 4;type_p <- 1}
  if(jj==2){r0 <- 4;type_p <- 2}
  
  yy1 <- pmin(1,((1-1/r0) - xx)/((1-xx)*v_eff))
  
  lines(xx,yy1,lty=type_p,col=col_pick)
  
  if(jj==1){  text(x=0.1,y=0.1+0.1*ii,labels=paste0("Veff = ",v_eff),adj=0,col=col_pick)}

  
  }
  
}

title(LETTERS[2],adj=0)


# Output plots
dev.copy(png,paste0("outputs/herd_immunity.png"),units="cm",width=15,height=25,res=200)
dev.off()


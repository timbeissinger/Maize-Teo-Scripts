library("plyr")
library("reshape2")
library("ggplot2")

#read all BKN bootstrap result files
bootstrap.BKN <- lapply(Sys.glob("BKN_bootstrap_*.result.txt"), read.delim)
length(bootstrap.BKN)
names(bootstrap.BKN) <- paste("BKN_bootstrap", seq(1, 20, 1), sep="")
names(bootstrap.BKN)
#convert list to data frame
bootstrap.BKN.df <- ldply(bootstrap.BKN, data.frame)
head(bootstrap.BKN.df)

#read all TIL bootstrap result files
bootstrap.TIL <- lapply(Sys.glob("TIL_bootstrap_*.result.txt"), read.delim)
length(bootstrap.TIL)
names(bootstrap.TIL) <- paste("TIL_bootstrap", seq(1, 20, 1), sep="")
bootstrap.TIL.df <- ldply(bootstrap.TIL, data.frame)
head(bootstrap.TIL.df)

# read the BKN direct estimate result file
direct.BKN <- read.delim("Hapmap2.BKN.6Hap.p2.result.txt")
#keep the header same as the previous dataframe for a convenience of rbind
direct.BKN[".id"] <- rep("direct.BKN", length(direct.BKN$time))
head(direct.BKN)
#reorder the columns
direct.BKN <- direct.BKN[, c(3, 1, 2)]
head(direct.BKN)

# read the TIL direct estimate result file
direct.TIL <- read.delim("Hapmap2.TIL.6Hap.p2.result.txt")
direct.TIL[".id"] <- rep("direct.TIL", length(direct.TIL$time))
direct.TIL <- direct.TIL[, c(3, 1, 2)]
head(direct.TIL)

#add a new column for all the data frames in order for the color assignment
bootstrap.TIL.df["group"] <- "teosinte bootstrapping estimate"
bootstrap.BKN.df["group"] <- "maize bootstrapping estimate"
direct.BKN["group"] <- "maize direct estimate"
direct.TIL["group"] <- "teosinte direct estimate"

#combine into a big dataframe
df <- rbind(bootstrap.TIL.df, bootstrap.BKN.df, direct.BKN, direct.TIL)

##try to find out ggplot2 default color assignment
#gg_color_hue <- function(n) {
#  hues = seq(15, 375, length=n+1)
#  hcl(h=hues, l=65, c=100)[1:n]
#}
#cols <- gg_color_hue(2)
#cols

#try to feed the same color scheme for maize and teosinte as those used in the Tim's figures
ggplot(data=df, aes(x=time, y=Ne, group=.id, color=group)) + geom_line() + scale_x_log10() + scale_y_log10() + xlab("years (u=3e-8, generation=1)") + ylab("effective population size") + theme_bw() + scale_color_manual(name="", values=c("maize direct estimate"="black", "maize bootstrapping estimate"="#F8766D", "teosinte direct estimate"="black", "teosinte bootstrapping estimate"="#3300FF")) + theme(legend.position=c(0.8, 0.9))
ggsave("TIL.BKN.boostrapping.msmc.pdf")
savehistory("msmc.bootstrapping.figures.R")

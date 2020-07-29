library(gridExtra)
library(grid)
library(ggplot2)
data<-read.csv("Family_GO_bar.csv")
data$family = factor(data$family, levels=c( "Caulobacteraceae", "Chitinophagaceae", "Chlamydiales", "Clostridiaceae", "Nocardiaceae", "Peptostreptococcaceae", 
                                            "Betaproteobacteria",  "Geodermatophilaceae", "Gp10", "Incertae_Sedis_XI", "Methylobacteriaceae", "Spartobacteria", 
                                            "Halomonadaceae","Orbaceae","Rubrobacteraceae"))
data$role = factor(data$role, levels=c("Positive", "Neutral", "Negative"))

ggplot(data, aes(x=reorder(name, pval),y=pval, fill=role)) + geom_bar(stat = "identity") +
  scale_fill_manual(values=c("#F21A00", "grey", "#3B9AB2"), name = "Immune Role") +
  facet_wrap(~family, ncol=6, scales="free") +
  coord_flip() +
  geom_col(color = "black", size = .2) +
  ylab("-log(padj)")+xlab(NULL) + theme_bw() + 
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(),
        axis.title.x = element_text(color="black", face="bold", size = 15),
        strip.text.x = element_text(size = 11, color = "black", face = "bold"),
        strip.background = element_rect(fill = "#F0F0F0", size=1, linetype="solid", 
                                        colour ="black"),
        panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
        panel.border = element_rect(fill=NA, colour = "black", size=1),
        panel.background = element_rect(fill = "#F0F0F0"),
        legend.position = c(.93, .15),
        legend.background = element_rect(fill="#F0F0F0",
                                         size=0.5, linetype="solid", 
                                         colour ="black"),
        legend.key = element_rect(fill="#F0F0F0"),
        legend.title = element_text(colour="black", size=16, 
                                    face="bold"),
        legend.text=element_text(size=14),
        legend.key.size = unit(.8,"cm"))

# g <- ggplotGrob(p1)
# # get the grobs that must be removed
# rm_grobs <- g$layout$name %in% c("panel-4-3", "strip-t-6-1")
# # remove grobs
# g$grobs[rm_grobs] <- NULL
# g$layout <- g$layout[!rm_grobs, ]
# 
# grid.newpage()
# grid.draw(g)



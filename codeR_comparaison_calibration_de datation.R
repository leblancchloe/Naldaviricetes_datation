library(tidyverse)

# Charge les 4 fichiers .dates
d1 <- read.table("chrono_both2_sample.dates", header=TRUE, 
                 col.names=c("node","meandate","stderr","inf95",
                             "sup95","instant_rate","stderr2",
                             "avg_rate","stderr3"))
d2 <- read.table("chrono_bd_both1_sample.dates", header=TRUE,
                 col.names=c("node","meandate","stderr","inf95",
                             "sup95","instant_rate","stderr2",
                             "avg_rate","stderr3"))
d3 <- read.table("chrono_euco2_sample.dates", header=TRUE,
                 col.names=c("node","meandate","stderr","inf95",
                             "sup95","instant_rate","stderr2",
                             "avg_rate","stderr3"))
d4 <- read.table("chrono_micro2_sample.dates", header=TRUE,
                 col.names=c("node","meandate","stderr","inf95",
                             "sup95","instant_rate","stderr2",
                             "avg_rate","stderr3"))

# Filtre uniquement les noeuds internes (ceux qui ont un age > 0)
d1_nodes <- d1 %>% filter(meandate > 0)
print(d1_nodes$node)  # affiche les numéros de noeuds

# noeud 59 = racine des Filamentoviridae 
noeud59 <- data.frame(
  modele = c("2 calibrations", "2 cal + birth-death", 
             "Eucoilini seul", "Bracovirus seul"),
  mean = c(d1$meandate[d1$node==59], d2$meandate[d2$node==59],
           d3$meandate[d3$node==59], d4$meandate[d4$node==59]),
  inf95 = c(d1$inf95[d1$node==59], d2$inf95[d2$node==59],
            d3$inf95[d3$node==59], d4$inf95[d4$node==59]),
  sup95 = c(d1$sup95[d1$node==59], d2$sup95[d2$node==59],
            d3$sup95[d3$node==59], d4$sup95[d4$node==59]))



#Figure age racine Filamentoviridae
noeud59$modele <- factor(noeud59$modele, 
                         levels=c("2 calibrations",
                                  "2 cal + birth-death",
                                  "Bracovirus seul",
                                  "Eucoilini seul"))


#Figure avec age médian 
ggplot(noeud59, aes(x=modele, y=mean)) +
  geom_errorbar(aes(ymin=inf95, ymax=sup95), 
                width=0.15, linewidth=0.7, color="grey40") +
  geom_point(size=4, color="black", fill="white", 
             shape=21, stroke=1.2) +
  geom_text(aes(label=paste0(round(mean, 1), " Ma")),
            vjust=-1.2, hjust=0.5, size=3.8, fontface="bold") +
  ylab("Âge estimé (Ma)") +
  xlab("") +
  ggtitle("Âge de la racine des Filamentoviridae") +
  scale_y_continuous(expand=expansion(mult=c(0.05, 0.15))) +
  theme_classic(base_size=13) +
  theme(
    axis.text.x = element_text(angle=20, hjust=1, size=11),
    axis.title.y = element_text(size=12),
    plot.title = element_text(hjust=0.5, face="bold", size=13),
    panel.grid.major.y = element_line(color="grey90", linewidth=0.4)
  ) +
  geom_hline(yintercept=noeud59$mean[noeud59$modele=="2 calibrations"],
             linetype="dashed", color="grey60", linewidth=0.5)



#noeud 39 = racine des arbres 

noeud39 <- data.frame(
  modele = c("2 calibrations", "2 cal + birth-death", 
             "Eucoilini seul","Bracovirus seul"),
  mean = c(d1$meandate[d1$node==39], d2$meandate[d2$node==39],
           d3$meandate[d3$node==39], d4$meandate[d4$node==39]),
  inf95 = c(d1$inf95[d1$node==39], d2$inf95[d2$node==39],
            d3$inf95[d3$node==39], d4$inf95[d4$node==39]),
  sup95 = c(d1$sup95[d1$node==39], d2$sup95[d2$node==39],
            d3$sup95[d3$node==39], d4$sup95[d4$node==39])
)


# Deviation par rapport aux 2 calibrations
ref39 <- noeud39$mean[noeud39$modele=="2 calibrations"]
noeud39$deviation <- round((noeud39$mean - ref39) / ref39 * 100, 1)
print(noeud39)

# Figure age racine Naldaviricetes 
noeud39$modele <- factor(noeud39$modele,
                         levels=c("2 calibrations",
                                  "2 cal + birth-death",
                                  "Bracovirus seul",
                                  "Eucoilini seul"))

ggplot(noeud39, aes(x=modele, y=mean)) +
  geom_errorbar(aes(ymin=inf95, ymax=sup95),
                width=0.15, linewidth=0.7, color="grey40") +
  geom_point(size=4, color="black", fill="white",
             shape=21, stroke=1.2) +
  geom_text(aes(label=paste0(round(mean, 1), " Ma")),
            vjust=-1.2, hjust=0.5, size=3.8, fontface="bold") +
  ylab("Âge estimé (Ma)") +
  xlab("") +
  ggtitle("Âge de la racine des Naldaviricetes") +
  scale_y_continuous(expand=expansion(mult=c(0.05, 0.15))) +
  theme_classic(base_size=13) +
  theme(
    axis.text.x = element_text(angle=20, hjust=1, size=11),
    axis.title.y = element_text(size=12),
    plot.title = element_text(hjust=0.5, face="bold", size=13),
    panel.grid.major.y = element_line(color="grey90", linewidth=0.4)
  ) +
  geom_hline(yintercept=ref39,
             linetype="dashed", color="grey60", linewidth=0.5)





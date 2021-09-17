#-------------------------------------------------------------------------------
#        generating all figures and supplementary tables for this study
#-------------------------------------------------------------------------------
#------------ false positives rate in the null setting (Figure 2) -------
fdr_null <- fread("./results/FPR_NULL.txt", data.table = F)
fdr_null  %>% mutate(cor2 = if_else(cor == TRUE, "0.5", "0.0")) %>%
  ggplot(aes(x = species, fill = cor2, y = sum)) +
  geom_bar(position = "dodge", stat="identity") +
  facet_grid( ~ size) + ylim(0, 100) +
  scale_fill_manual(values=c( "#006400", "#C1FFC1")) +
  xlab("") +  ylab("False Positives (%)") +  theme_bw(base_size = 20) +
  theme(
    panel.border = element_rect( fill = NA, size = 0.5, linetype = "solid", colour = "black"),
    strip.background =  element_rect(fill = "white", colour = "white"),
    axis.text.x = element_text( colour = "black", family = "Times"),
    axis.ticks = element_blank(),
    axis.title.y = element_text( vjust = 1.3, colour = "black", family = "Times", size = 20),
    axis.text.y = element_text(colour = "black", family = "Times"),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.spacing = unit(.5, "lines"),
    strip.text.x = element_text(size = 18, hjust = 0.5, vjust = 0.5),
    strip.text.y = element_text(size = 18, hjust = 0.5, vjust = 0.5),
    legend.title = element_text(hjust = 1),
    legend.text = element_text(size = 20),
    legend.key.width = unit(0.5, 'cm'),
    legend.box = "vertical",
    legend.direction = "vertical",
    legend.position = "bottom"
  ) + guides(
    fill = guide_legend(
      title = "Pearson correlation between traits:",
      nrow = 1,
      byrow = TRUE,
      order = 2,
      title.hjust = 0
    )
  )
# save figure 2 as PDF
ggsave(file="./results/Figure_2.pdf", dpi=300,family="Times")

#------------ false positives rate in the non-null setting  (Figure 3) -------
fp <- fread("./results/FPR.txt", data.table = F)
fp  %>%  
  pivot_longer(cols = c("gemma", "mt", "t1", "t2", "t3"), names_to = "model") %>%
  mutate(
    model = recode(
      model,
      "mt" = "MT",
      "gemma" = "mvLMM",
      "t1" = "T1",
      "t2" = "T2",
      "t3" = "T3"
    ),
    architecture = as.factor(architecture),
    architecture = recode(
      architecture,
      "LD" = "Spurious Pleiotropy",
      "partially" = "Partial Pleiotropy",
      "pleiotropic" = "Pleiotropy"
    )) %>% 
  ggplot(aes(x = species, fill = model, y = value)) +
  geom_bar(position = "dodge", stat="identity") +
  facet_grid(size ~ architecture, scales = "free_x", space = "free_x")+
  scale_fill_manual(values=c( "#006400","#66CDAA", "#A2CD5A", "#ADFF2F", "#8B8B00")) +
  xlab("") +  ylab("False Positives (%)") +  theme_bw(base_size = 20) + ylim(0, 100) +
  theme(
    panel.border = element_rect( fill = NA, size = 0.5, linetype = "solid", colour = "black"),
    strip.background =  element_rect(fill = "white", colour = "white"),
    axis.text.x = element_text( colour = "black", family = "Times"),
    axis.ticks = element_blank(),
    axis.title.y = element_text( vjust = 1.3, colour = "black", family = "Times", size = 20),
    axis.text.y = element_text(colour = "black", family = "Times"),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.spacing = unit(.5, "lines"),
    strip.text.x = element_text(size = 18, hjust = 0.5, vjust = 0.5),
    strip.text.y = element_text(size = 18, hjust = 0.5, vjust = 0.5),
    legend.title = element_text(hjust = 0.5),
    legend.text = element_text(size = 20),
    legend.key.width = unit(0.5, 'cm'),
    legend.box = "vertical",
    legend.direction = "vertical",
    legend.position = "bottom"
  ) + guides(fill = guide_legend(title ="GWAS Model:", nrow = 1, byrow = TRUE, order = 2)  )
# save figure 2 as PDF
ggsave(file="./results/Figure_3.pdf", dpi=300,family="Times")

#------------ true positives rate (Figure 4) -------
tp <- fread("./results/TP.txt", data.table = F)
tp  <- tp %>%
  mutate(
    model = recode(
      model,
      "detected_mt_1" = "MT",
      "detected_mt_2" = "MT",
      "detected_mt_3" = "MT",
      "detected_gemma_1" = "mvLMM",
      "detected_gemma_2" = "mvLMM",
      "detected_gemma_3" = "mvLMM",
      "detected_t1" = "T1",
      "detected_t2" = "T2",
      "detected_t3" = "T3"
    ),
    architecture = as.factor(architecture),
    species = as.factor(species),
    size = as.factor(size),
    trait = as.factor(trait),
    model = as.factor(model),
    snp = as.factor(snp),
    trait = recode(
      trait,
      "1" = "Trait 1",
      "2" = "Trait 2",
      "3" = "Trait 3"
    ),
    snp = recode(snp, "1" = "QTN 1", "2" = "QTN 2", "3" = "QTN 3"),
    architecture = recode(
      architecture,
      "LD" = "Spurious Pleiotropy",
      "partially" = "Partial Pleiotropy",
      "pleiotropic" = "Pleiotropy"
    ),
    alpha = ifelse(
      architecture == "Pleiotropy" |
        (architecture == "Partial Pleiotropy" &
           snp == "QTN 2" & trait != "Trait 3"),
      "Pleiotropic QTN",
      "Trait-Specific QTN"
    )
  )

tp  %>%  filter(architecture == "Spurious Pleiotropy" & snp !=  "QTN 1")%>% 
  ggplot(aes(x = trait, fill = model, y = sum)) +
  geom_col(position = "dodge") +
  facet_grid(snp + species ~ size, scales = "free_x", space = "free_x") +
  scale_fill_manual(values=c( "#006400", "#66CDAA","#A2CD5A", "#ADFF2F", "#8B8B00")) +
  xlab("") +  ylab("QTN Detection (%)") +  theme_bw(base_size = 20) + ylim(0, 100) +
  theme(
    panel.border = element_rect( fill = NA, size = 0.5, linetype = "solid", colour = "black"),
    strip.background =  element_rect(fill = "white", colour = "white"),
    axis.text.x = element_text( colour = "black", family = "Times"),
    axis.ticks = element_blank(),
    axis.title.y = element_text( vjust = 1.3, colour = "black", family = "Times", size = 20),
    axis.text.y = element_text(colour = "black", family = "Times"),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.spacing = unit(.5, "lines"),
    strip.text.x = element_text(size = 18, hjust = 0.5, vjust = 0.5),
    strip.text.y = element_text(size = 18, hjust = 0.5, vjust = 0.5),
    legend.title = element_text(hjust = 0.5),
    legend.text = element_text(size = 20),
    legend.key.width = unit(0.5, 'cm'),
    legend.box = "vertical",
    legend.direction = "vertical",
    legend.position = "bottom"
  ) + guides(fill = guide_legend(title ="Models:", nrow = 1, byrow = TRUE, order = 2)  )

ggsave(file="./results/Figure_4.pdf", width = 15, height = 20, units = "cm", dpi=300,family="Times")

#------------ true positives rate (Figure 5) -------
tp  %>%  filter(architecture == "Pleiotropy" & snp !=  "QTN 1")%>% 
  ggplot(aes(x = trait, fill = model, y = sum)) +
  geom_col(position = "dodge") +
  facet_grid(snp + species ~ size, scales = "free_x", space = "free_x") +
  scale_fill_manual(values=c( "#006400", "#66CDAA","#A2CD5A", "#ADFF2F", "#8B8B00")) +
  xlab("") +  ylab("QTN Detection (%)") +  theme_bw(base_size = 20) + ylim(0, 100) +
  theme(
    panel.border = element_rect( fill = NA, size = 0.5, linetype = "solid", colour = "black"),
    strip.background =  element_rect(fill = "white", colour = "white"),
    axis.text.x = element_text( colour = "black", family = "Times"),
    axis.ticks = element_blank(),
    axis.title.y = element_text( vjust = 1.3, colour = "black", family = "Times", size = 20),
    axis.text.y = element_text(colour = "black", family = "Times"),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.spacing = unit(.5, "lines"),
    strip.text.x = element_text(size = 18, hjust = 0.5, vjust = 0.5),
    strip.text.y = element_text(size = 18, hjust = 0.5, vjust = 0.5),
    legend.title = element_text(hjust = 0.5),
    legend.text = element_text(size = 20),
    legend.key.width = unit(0.5, 'cm'),
    legend.box = "vertical",
    legend.direction = "vertical",
    legend.position = "bottom"
  ) + guides(fill = guide_legend(title ="Models:", nrow = 1, byrow = TRUE, order = 2)  )

ggsave(file="./results/Figure_5.pdf", width = 15, height = 20, units = "cm", dpi=300,family="Times")

#------------ true positives rate (Figure 6) -------
tp  %>%  filter(architecture == "Partial Pleiotropy" & snp !=  "QTN 1")%>% 
  ggplot(aes(x = trait, fill = model, y = sum)) +
  geom_col(position = "dodge") +
  facet_grid(snp + species ~ size, scales = "free_x", space = "free_x") +
  scale_fill_manual(values=c( "#006400", "#66CDAA","#A2CD5A", "#ADFF2F", "#8B8B00")) +
  xlab("") +  ylab("QTN Detection (%)") +  theme_bw(base_size = 20) + ylim(0, 100) +
  theme(
    panel.border = element_rect( fill = NA, size = 0.5, linetype = "solid", colour = "black"),
    strip.background =  element_rect(fill = "white", colour = "white"),
    axis.text.x = element_text( colour = "black", family = "Times"),
    axis.ticks = element_blank(),
    axis.title.y = element_text( vjust = 1.3, colour = "black", family = "Times", size = 20),
    axis.text.y = element_text(colour = "black", family = "Times"),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.spacing = unit(.5, "lines"),
    strip.text.x = element_text(size = 18, hjust = 0.5, vjust = 0.5),
    strip.text.y = element_text(size = 18, hjust = 0.5, vjust = 0.5),
    legend.title = element_text(hjust = 0.5),
    legend.text = element_text(size = 20),
    legend.key.width = unit(0.5, 'cm'),
    legend.box = "vertical",
    legend.direction = "vertical",
    legend.position = "bottom"
  ) + guides(fill = guide_legend(title ="Models:", nrow = 1, byrow = TRUE, order = 2)  )

ggsave(file="./results/Figure_6.pdf", width = 15, height = 20, units = "cm", dpi=300,family="Times")


#------------ true positives rate (Figure S1) -------
tp  %>%  filter(architecture == "Spurious Pleiotropy")%>% 
  ggplot(aes(x = trait, fill = model, y = sum)) +
  geom_col(position = "dodge") +
  facet_grid(snp + species ~ size, scales = "free_x", space = "free_x") +
  scale_fill_manual(values=c( "#006400", "#66CDAA","#A2CD5A", "#ADFF2F", "#8B8B00")) +
  xlab("") +  ylab("QTN Detection (%)") +  theme_bw(base_size = 20) + ylim(0, 100) +
  theme(
    panel.border = element_rect( fill = NA, size = 0.5, linetype = "solid", colour = "black"),
    strip.background =  element_rect(fill = "white", colour = "white"),
    axis.text.x = element_text( colour = "black", family = "Times"),
    axis.ticks = element_blank(),
    axis.title.y = element_text( vjust = 1.3, colour = "black", family = "Times", size = 20),
    axis.text.y = element_text(colour = "black", family = "Times"),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.spacing = unit(.5, "lines"),
    strip.text.x = element_text(size = 18, hjust = 0.5, vjust = 0.5),
    strip.text.y = element_text(size = 18, hjust = 0.5, vjust = 0.5),
    legend.title = element_text(hjust = 0.5),
    legend.text = element_text(size = 20),
    legend.key.width = unit(0.5, 'cm'),
    legend.box = "vertical",
    legend.direction = "vertical",
    legend.position = "bottom"
  ) + guides(fill = guide_legend(title ="Models:", nrow = 1, byrow = TRUE, order = 2)  )

ggsave(file="./results/Figure_S1.pdf", width = 15, height = 20, units = "cm", dpi=300,family="Times")

#------------ true positives rate (Figure S2) -------
tp  %>%  filter(architecture == "Pleiotropy")%>% 
  ggplot(aes(x = trait, fill = model, y = sum)) +
  geom_col(position = "dodge") +
  facet_grid(snp + species ~ size, scales = "free_x", space = "free_x") +
  scale_fill_manual(values=c( "#006400", "#66CDAA","#A2CD5A", "#ADFF2F", "#8B8B00")) +
  xlab("") +  ylab("QTN Detection (%)") +  theme_bw(base_size = 20) + ylim(0, 100) +
  theme(
    panel.border = element_rect( fill = NA, size = 0.5, linetype = "solid", colour = "black"),
    strip.background =  element_rect(fill = "white", colour = "white"),
    axis.text.x = element_text( colour = "black", family = "Times"),
    axis.ticks = element_blank(),
    axis.title.y = element_text( vjust = 1.3, colour = "black", family = "Times", size = 20),
    axis.text.y = element_text(colour = "black", family = "Times"),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.spacing = unit(.5, "lines"),
    strip.text.x = element_text(size = 18, hjust = 0.5, vjust = 0.5),
    strip.text.y = element_text(size = 18, hjust = 0.5, vjust = 0.5),
    legend.title = element_text(hjust = 0.5),
    legend.text = element_text(size = 20),
    legend.key.width = unit(0.5, 'cm'),
    legend.box = "vertical",
    legend.direction = "vertical",
    legend.position = "bottom"
  ) + guides(fill = guide_legend(title ="Models:", nrow = 1, byrow = TRUE, order = 2)  )

ggsave(file="./results/Figure_S2.pdf", width = 15, height = 20, units = "cm", dpi=300,family="Times")

#------------ true positives rate (Figure 6) -------
tp  %>%  filter(architecture == "Partial Pleiotropy")%>% 
  ggplot(aes(x = trait, fill = model, y = sum)) +
  geom_col(position = "dodge") +
  facet_grid(snp + species ~ size, scales = "free_x", space = "free_x") +
  scale_fill_manual(values=c( "#006400", "#66CDAA","#A2CD5A", "#ADFF2F", "#8B8B00")) +
  xlab("") +  ylab("QTN Detection (%)") +  theme_bw(base_size = 20) + ylim(0, 100) +
  theme(
    panel.border = element_rect( fill = NA, size = 0.5, linetype = "solid", colour = "black"),
    strip.background =  element_rect(fill = "white", colour = "white"),
    axis.text.x = element_text( colour = "black", family = "Times"),
    axis.ticks = element_blank(),
    axis.title.y = element_text( vjust = 1.3, colour = "black", family = "Times", size = 20),
    axis.text.y = element_text(colour = "black", family = "Times"),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.spacing = unit(.5, "lines"),
    strip.text.x = element_text(size = 18, hjust = 0.5, vjust = 0.5),
    strip.text.y = element_text(size = 18, hjust = 0.5, vjust = 0.5),
    legend.title = element_text(hjust = 0.5),
    legend.text = element_text(size = 20),
    legend.key.width = unit(0.5, 'cm'),
    legend.box = "vertical",
    legend.direction = "vertical",
    legend.position = "bottom"
  ) + guides(fill = guide_legend(title ="Models:", nrow = 1, byrow = TRUE, order = 2)  )

ggsave(file="./results/Figure_S3.pdf", width = 15, height = 20, units = "cm", dpi=300,family="Times")


#--- LD heatmaps (Figures S4 and S5)------
# The phisical distance was removed manually by editing the PDF as it is not
# accurate in these figures nor it is necessary in the context used here
gdsfmt::showfile.gds(closeall = TRUE, verbose = F)
genofile <- snpgdsOpen("./data/maize/ames_2815_mac5_ld09_imputed.gds")
snps <- read.gdsn(index.gdsn(genofile, "snp.id"))
rgb.palette <- colorRampPalette(rev(c("blue", "yellow", "red")), space = "rgb")
# selecting QTN2 for trait 1 and trait 2
set <- which( snps %in% "S7_106651946" )

pdf("./results/Figure_S4.pdf")
# calculating the LD with the 20 SNPs upstream and downstream
ld1 <- snpgdsLDMat(genofile, snp.id = snps[(set-20) : (set + 20)],  slide = -1, method = "r", verbose=F)$LD  ^ 2
ld1 <- LDheatmap(ld1, color=rgb.palette(50))
LDheatmap.marks(ld1,  21,  21,  gp=grid::gpar(cex=2),  pch = "*")
dev.off()
# selecting QTN2 for trait 3
set <- which( snps %in% "S1_116015247" )
pdf("./results/Figure_S5.pdf")
ld <- snpgdsLDMat(genofile, snp.id = snps[(set-20) : (set + 20)],  slide = -1, method = "r", verbose=F)$LD  ^ 2
ld2 <- LDheatmap(ld, color=rgb.palette(50))
LDheatmap.marks(ld2,  21,  21,  gp=grid::gpar(cex=2),  pch = "*")
dev.off()
gdsfmt::showfile.gds(closeall = TRUE, verbose = F)
#------------ Supplementary Table S1------------
tp <- tp[,-10]
colnames(tp) <- c("Scenario", "Species", "Sample Size",  "Trait", "QTN", "Model", "MAF", "Allelic Effect",  "QTN Detection rate (%)" )
tp <- tp %>% mutate(Scenario = recode(Scenario, 
                                      "Spurious Pleiotropy" = "SP", 
                                      "Pleiotropy" = "P",
                                      "Partial Pleiotropy" = "PP")) %>%
  arrange(Scenario, Species, `Sample Size`, Trait, QTN, Model)
tp$QTN <- gsub("QTN", "", tp$QTN)
tp$Trait <- gsub("Trait", "", tp$Trait)
tp[,7] <- round(tp[,7],2)
ft <- flextable(tp)
ft <- line_spacing(ft, space = 1, part = "header")
ft <- line_spacing(ft, space = 0.3, part = "body")
ft <- add_header_lines(ft,
                       values = c("TABLE S1: Quantitative trait nucleotide (QTN) detection rate of all the 20 settings simulated."))
ft <- add_footer_lines(ft, values = "SP: Spurious Pleiotropy; P: Pleiotropy; PP: Partial Pleiotropy;
T1, T2, T3: Univariate stepwise model selection in traits 1, 2, and 3, repectively.;
MT: Multivariate Multi-locus Stepwise model selection; 
mvLMM: Multivariate linear mixed model.", top = FALSE)
ft <- align(ft, align = "right", part = "body")
ft <- align(ft, align = "center", part = "header")
ft <- font(ft, fontname = "Times New Roman", part = "all")
ft <- fontsize(ft, size = 12, part = "all")
# page properties
sect_properties <- prop_section(
  page_size = page_size(orient = 'portrait'),
  type = "continuous",
  page_margins = page_mar(  bottom = 0.5,
                            top = 0.5,
                            right = 0.5,
                            left = 0.5)
)
save_as_docx(ft, path = "./results/Table_S1.docx", pr_section = sect_properties)


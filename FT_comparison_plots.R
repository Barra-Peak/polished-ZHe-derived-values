### Code to generate synthetic polished zircon dataset
### Compare derived data using [insert name of r script] and other methods to calculate FT
### Generate plots published in Peak (2025, Geochronology). Plot titles, legends, and annotations modified in Adobe Illustrator for publication.
# Barra Peak June 10, 2025

#### Load necessary R packages ----
library(openxlsx) 
library(tidyverse)
library(ggplot2)
library(ggh4x)


##### Synthetic Data
####  Generate synthetic dataset ----
initial <- read.xlsx("/Users/barrapeak/Mirror/Moved_Files/grain_polishing/Geochron_submission/Reviews/Resubmit_Text/Repository_code/initial_synthetic.xlsx",
                            startRow = 1, colNames = TRUE, check.names = TRUE, sep.names = "_") # load formatted data file with correct column names using full file path
# set up synthetic grain parameters:
eU <- unique(initial[, 14:21]) # take eU from formatted initial file
eU <- eU[1:3, ] %>%
  mutate(X4He = 1000,
         X._1σ = 2,
         U = 1,
         X._1σ.1 = 0.01,
         Th = U,
         X._1σ.2 = X._1σ.1,
         X147Sm = 0,
         X._1σ.3 = 0)
eU <- unique(eU)
L <- data.frame("Smallest" = 60, 
                "Small" = 100,
                "Medium" = 150,
                "Large" = 200) # c-axis parallel length, one of 4 options
A.Ratio <- data.frame("0.3" = 0.3,"0.4" = 0.4,"0.5" = 0.5,"0.6" = 0.6, 
                      "0.7" = 0.7, "0.8" = 0.8, "0.9" = 0.9, "1" = 1) # aspect ratio between c-axis parallel and perpendicular lengths
W.Ratio <- data.frame("W.0.5" = 0.5, "W.0.6" = 0.6, "W.0.7" = 0.7,
                      "W.0.8" = 0.8, "W.0.9" = 0.9, "W.1.0" = 1) # width ratio between two c-axis perpendicular crystal width measurements
Orientation <- data.frame("Perpendicular" = 1,
                          "Parallel" = 2) # grinding and polished orientation relative to the crystal c-axis
g.percent <- data.frame("None" = 0,
                        "P25" = 0.25, "P30" = 0.3,"P35" = 0.35,"P40" = 0.4,
                        "P45" = 0.45, "P50" = 0.5, "P55" = 0.55, "P60" = 0.6,
                        "P65" = 0.65, "P70" = 0.7,"P75" = 0.75) # percent of original grain length or width removed by grinding and polishing
Shape <- data.frame("Ellipsoid" = 1, 
                    "Cylinder" = 2,
                    "Tetragon" = 3) # base geometry
Np <- data.frame("Zero" = 0,
                 "One" = 1,
                 "Two" = 2) # number of possible terminations for "tetragonal" grains

### create full dataset of all possible parameter combinations
Base <- cbind(eU, L, A.Ratio, W.Ratio, Orientation, g.percent, Shape, Np)
Base.2 <- pivot_longer(Base, cols = c("Smallest", "Small", "Medium", "Large"), names_to = "Size",
                       values_to = "Length_2")
Base.3 <- pivot_longer(Base.2, cols = names(A.Ratio),
                       names_to = "X1", values_to = "Aspect.Ratio") %>%
  select(!c("X1"))
Base.4 <- pivot_longer(Base.3, cols = c("Perpendicular", "Parallel"),
                       names_to = "Orientation.Name",
                       values_to = "Orientation")
Base.5 <- pivot_longer(Base.4, cols = names(g.percent),
                       names_to = "Grind.Name", values_to = "Grind.Percent") 
Base.6 <- pivot_longer(Base.5, cols = c("Ellipsoid", "Cylinder", "Tetragon"),
                       names_to = "Shape", values_to = "Geometry")
Base.7 <- pivot_longer(Base.6, cols = c("Zero", "One", "Two"),
                       names_to = "XNp", values_to = "Np") %>%
  select(!c("XNp"))
Base.8 <- pivot_longer(Base.7, cols = names(W.Ratio),
                       names_to = "W", values_to = "W.Ratio") %>%
  select(!c(W))

Clean.1 <- Base.8 %>%
  filter((Shape == "Ellipsoid" & Np == 0) |
           (Shape == "Cylinder" & Np == 0) |
           (Shape == "Tetragon" & Np == 0) |
           (Shape == "Tetragon" & Np == 1) |
           (Shape == "Tetragon" & Np == 2))

Clean.2 <- Clean.1 %>%
  mutate(Length_1 = Length_2) %>%
  mutate(Width_1 = Length_2*Aspect.Ratio) %>%
  mutate(Width_2 = Width_1*W.Ratio) %>%
  mutate(Grind_Depth = case_when(Orientation == 1 ~ Length_1*Grind.Percent,
                                 .default = Width_2*Grind.Percent))

All.Data <- Clean.2 %>%
  mutate(Length_2 = case_when((Orientation == 1 & Grind_Depth != 0) ~ Length_2 - Grind_Depth,
                              .default = Length_2)) %>%
  mutate(Length_1 = case_when((Orientation == 1 & Shape == "Ellipsoid" & Grind_Depth != 0) ~
                                Width_2,
                              .default = Length_2)) %>%
  mutate(Width_2 = case_when((Orientation == 2 & Grind_Depth != 0) ~ Width_2 - Grind_Depth,
         .default = Width_2)) %>%
  mutate(Length_P = case_when((Orientation == 2 & Grind_Depth != 0) ~ 0.95*Length_1,
                              .default = NA)) %>%
  mutate(Width_P = case_when((Orientation == 2 & Grind_Depth != 0) ~ 0.95*Width_1,
                             .default = NA)) %>%
  mutate("Crystal_Fragment." = "No",
         Analysis_Session = Geometry,
         Grain = paste(Shape, as.character(Np), Size, as.character(Aspect.Ratio), as.character(W.Ratio),
                       Grind.Name, Orientation.Name, sep = "_")) %>%
  select(names(initial))

# save synthetic dataset, replace pathname with local pathname
write.xlsx(All.Data, "/Users/barrapeak/Desktop/synthetic_zircon.xlsx")

# save synthetic dataset with grinding depth set to 0 for all analyses (used to evaluate Ketcham et al., 2011 protocol)
All.Data.g0 <- All.Data %>%
  mutate(Grind_Depth = 0)
write.xlsx(All.Data.g0, "/Users/barrapeak/Desktop/synthetic_zircon_g0.xlsx")


#### Synthetic Data Results Comparison
### Load synthetic data results - new protocol ----
synthetic_data <- read.xlsx("/Users/barrapeak/Mirror/Moved_Files/grain_polishing/Geochron_submission/Reviews/Resubmit_Text/Repository_code/synthetic_zircon_new_protocol.xlsx", check.names = TRUE) %>%
  mutate(across(-c("Grain", "Crystal.Fragment."), as.numeric)) %>%
  separate_wider_delim(Grain, "_", names = c("Shape", "XNp", "Size", "Aspect.Ratio", "W.Ratio", "Ground", "Orientation.Name")) %>%
  mutate(Shape = case_when(Geometry == 3 ~ paste(Shape, XNp, "Np"),
                           .default = Shape)) %>%
  select(!c("XNp")) %>%
  mutate(Aspect.Ratio = as.numeric(Aspect.Ratio)) %>%
  mutate(W.Ratio = as.numeric(W.Ratio)) %>%
  mutate(Percent.unc = X..1σ.16/Combined.FT)

### Load synthetic data results - Ketcham et al. (2011) protocol ----
synthetic_ketcham <- read.xlsx("/Users/barrapeak/Mirror/Moved_Files/grain_polishing/Geochron_submission/Reviews/Resubmit_Text/Repository_code/synthetic_zircon_ketcham_protocol.xlsx", check.names = TRUE) %>%
  mutate(across(-c("Grain", "Crystal.Fragment."), as.numeric)) %>%
  separate_wider_delim(Grain, "_", names = c("Shape", "XNp", "Size", "Aspect.Ratio", "W.Ratio", "Ground", "Orientation.Name")) %>%
  mutate(Shape = case_when(Geometry == 3 ~ paste(Shape, XNp, "Np"),
                           .default = Shape)) %>%
  select(!c("XNp")) %>%
  mutate(Aspect.Ratio = as.numeric(Aspect.Ratio)) %>%
  mutate(W.Ratio = as.numeric(W.Ratio)) %>%
  mutate(Percent.unc = X..1σ.16/Combined.FT)

synthetic_data_merge <- merge(synthetic_data, synthetic_ketcham, by = c(
  "Analysis.Session", "Shape", "Size", "Aspect.Ratio", "W.Ratio", "Orientation.Name",
  "Length.1", "Width.1", "Length.2", "Width.2", "Geometry", "Np", "Length.P", "Width.P",
  "Orientation", "Ground", "Crystal.Fragment.", "X4He", "X..1σ", "U", "X..1σ.1","Th", "X..1σ.2", "X147Sm", "X..1σ.3")
) %>%
select(!c("Crystal.Fragment.")) %>%
  rename("Syn.FT" = "Combined.FT.x", "Syn.FT.Unc" = "X..1σ.16.x",
         "Syn.RFT" = "RFT.x", "Syn.RFT.Unc" = "X..1σ.17.x",
         "Syn.Percent.Unc" = "Percent.unc.x",
         "NC.FT" = "Combined.FT.y", "NC.FT.Unc" = "X..1σ.16.y",
         "NC.RFT" = "RFT.y", "NC.RFT.Unc" = "X..1σ.17.y",
         "NC.Percent.Unc" = "Percent.unc.y")

### Calculate Reiners et al. (2007) FT for synthetic data ----
a1.U238 <- -4.31 # fit parameters from Farley 2002
a2.U238 <- 4.92
a1.Th232 <- -5
a2.Th232 <- 6.8

intermediate_vals <- synthetic_data_merge %>%
  mutate(r_c = (Width.1 + (Width.2 + Grind.Depth.x))/4) %>%
  mutate(l = (Length.1 + Length.2)/2) %>%
  mutate(r_t1 = Width.1/2) %>%
  mutate(r_t2 = Width.2/2) %>%
  mutate(w = case_when(Grind.Depth.x == r_c ~ 0,
                       .default = acos(1-(Grind.Depth.x/r_c)))) %>%
  mutate(Beta_c = case_when(Grind.Depth.x == 0 ~ (2/l) + (2/(r_c)),
                            .default = (2/l) + 
                              ((2*r_c*(pi-w))/(((pi-w)*(r_c^2)) + 
                                                 ((r_c - Grind.Depth.x)*sqrt((r_c^2) - ((r_c-Grind.Depth.x)^2))))))) %>%
  mutate(Beta_t = case_when(Grind.Depth.x == 0 ~ (2/l) + (1/r_t1) + (1/r_t2),
                            .default = (2/l) + (1/(2*r_t2)) + (1/r_t1)))

all_synthetic_data <- cbind(synthetic_data_merge, select(intermediate_vals, c("Beta_c", "Beta_t"))) %>%
  mutate(Beta = case_when((Orientation == 2 & 
                             Geometry == 2 & 
                             between(Grind.Depth.x, 0.1, (Grind.Depth.x + Width.2)/2)) ~ Beta_c,
                          (Orientation == 2 & 
                             Geometry == 2 & 
                             Grind.Depth.x == 0) ~ Beta_c,
                          (Orientation == 2 & 
                             Geometry == 3 & 
                             between(Grind.Depth.x, 0.1, (Grind.Depth.x + Width.2)/2) &
                             Np == 0) ~ Beta_t,
                          (Orientation == 2 & 
                             Geometry == 3 & 
                             Grind.Depth.x == 0 &
                             Np == 0) ~ Beta_t,
                          .default = NA)) %>%
  select(!c("Beta_c", "Beta_t")) %>%
  mutate(FT_238U = 1 + (a1.U238*Beta) + (a2.U238*(Beta^2)),
         FT_232Th = 1 + (a1.Th232*Beta) + (a2.Th232*(Beta^2))) %>%
  mutate(A = 1/(1.04 + (0.245*(Th/U)))) %>%
  mutate(A = case_when(is.na(Beta) == F ~ A,
                       .default = NA)) %>%
  mutate(FTc = (A*FT_238U) + ((1-A)*FT_232Th))

# save all synthetic data method comparisons as excel file
write.xlsx(all_synthetic_data, "/Users/barrapeak/Desktop/all_synthetic_data_comparisons.xlsx")

### Plot synthetic data and comparisons ----
# color, shape, and size scales to use across plots
geom_color <- c("1" = "#E58606", "2" = "#5D69B1", "3" = "#52BCA3")
size_color <- c("Smallest" = "#72190E", "Small" = "#332288", "Medium" = "#225555", "Large" = "#997700")
size_scale <- c("Smallest" = 1, "Small" = 2, "Medium" = 3, "Large" = 4)
shape_scale <- c("Ellipsoid" = 1, "Cylinder" = 10, "Tetragon 0 Np" = 0,
                 "Tetragon 1 Np" = 2, "Tetragon 2 Np" = 5)

# prepare dataframe for plotting
all_data_comp <- all_synthetic_data %>%
  rename("Reiners.FT" = "FTc") %>%
  group_by(Shape, Size, Aspect.Ratio, W.Ratio, Orientation.Name) %>%
  arrange(Shape, Size, Aspect.Ratio, W.Ratio, Orientation.Name) %>%
  mutate(Syn.Diff.From.Whole = ((Syn.FT/Syn.FT[which.min(Grind.Depth.x)])-1)*100) %>%
  mutate(Syn.Diff.From.Ketcham = ((Syn.FT/NC.FT) -1)*100) %>%
  mutate(Syn.Diff.From.Reiners = ((Syn.FT/Reiners.FT) -1)*100) %>%
  mutate(Reiners.Diff.From.Whole = ((Reiners.FT/Reiners.FT[which.min(Grind.Depth.x)])-1)*100) %>%
  mutate(Shape.fac = factor(Shape, levels = c("Ellipsoid", "Cylinder",
                                              "Tetragon 0 Np", "Tetragon 1 Np",
                                              "Tetragon 2 Np"))) %>%
  mutate(Size = factor(Size, levels = c("Smallest",
                                        "Small",
                                        "Medium",
                                        "Large"))) %>%
  mutate(Refs = case_when(Ground == "P50" ~ "same",
                          Ground == "None" ~ "same",
                          .default = "diff")) %>%
  mutate(Refs = as.factor(Refs)) %>%
  mutate(Orientation.fac = case_when(Orientation == 1 ~ "Perpendicular",
                                     .default = "Parallel"))

## Plot percent difference combined FT vs grind depth for synthetic data ----
# Maximum symmetry case
syn.diffFt.g1 <- all_data_comp %>%
  filter(W.Ratio == 1) %>%
  filter(Aspect.Ratio == 1) %>%
  filter(Volume.x > 0) %>%
  ggplot(aes(x = Grind.Depth.x, y = Syn.Diff.From.Whole, group = interaction(Size, as.factor(Aspect.Ratio)), shape = as.factor(Shape), color = Size)) +
  geom_hline(yintercept = 0, linetype = "dotted", linewidth = 0.5, color = "gray") +
  geom_line() +
  scale_color_manual(values = size_color) +
  coord_cartesian(
    xlim = NULL,
    ylim = c(-50, 50),
    expand = TRUE,
    default = FALSE,
    clip = "on"
  ) +
  scale_y_continuous(breaks = seq(-50, 50, 10)) +
  scale_x_continuous(limits = c(0, 150), breaks = seq(0, 150, 50)) +
  xlab("Grinding Depth (um)") +
  ylab("% Difference FT from whole grain") + 
  ggtitle(paste("FT Difference from Whole Grain ((Syn.FT/Syn.FT0)-1)*100)", paste("Width Ratio = 1", sep = " "), sep = "\n")) +
  theme(panel.background = element_rect(fill = "transparent", colour = "transparent", linewidth = 1),
        panel.border = element_rect(fill = NA, color = "black", linewidth = 1),
        axis.ticks = element_line(colour = "black", linewidth = 0.5),
        axis.ticks.length = unit(-0.1,"cm"),
        axis.text = element_text(size = 7),
        axis.text.x = element_text(size = 7),
        axis.title = element_text(size = 9),
        plot.margin = grid::unit(c(0,0,0,0), "mm"),
        legend.background = element_rect(fill = NA, colour = NA),
        legend.key = element_blank(),
        legend.position = "bottom",
        legend.text = element_text(size = 7),
        legend.title = element_text(size = 9),
        plot.title = element_text(hjust = 0.5)) +
  facet_grid2(as.factor(Shape) ~ Orientation.fac, scales = "free_y", axes ="margins") 
syn.diffFt.g1

# save plot as pdf, replace path name with local path name
ggsave("/Users/barrapeak/Desktop/Synthetic.FT.Max.Symmetry.pdf", syn.diffFt.g1, dpi = 300, width = 3.25, height = 7)

# Minimum symmetry
syn.diffFt.g2 <- all_data_comp %>%
  filter(W.Ratio == 0.5) %>%
  filter(Aspect.Ratio == 0.3) %>%
  filter(Volume.x > 0) %>%
  ggplot(aes(x = Grind.Depth.x, y = Syn.Diff.From.Whole, group = interaction(Size, as.factor(Aspect.Ratio)), shape = as.factor(Shape), color = Size)) +
  geom_hline(yintercept = 0, linetype = "dotted", linewidth = 0.5, color = "gray") +
  geom_line(linetype = "dashed") +
  #geom_point(size = 1) +
  scale_color_manual(values = size_color) +
  scale_shape_manual(values = shape_scale) +
  coord_cartesian(
    xlim = NULL,
    ylim = c(-50, 50),
    expand = TRUE,
    default = FALSE,
    clip = "on"
  ) +
  scale_y_continuous(breaks = seq(-50, 50, 10)) +
  scale_x_continuous(limits = c(0, 150), breaks = seq(0, 150, 20))+
  xlab("Grinding Depth (um)") +
  ylab("% Difference FT from whole grain") + 
  ggtitle(paste("FT Difference from Whole Grain ((Syn.FT/Syn.FT0)-1)*100)", paste("Width Ratio = 0.5", sep = " "), sep = "\n")) +
  theme(panel.background = element_rect(fill = "transparent", colour = "transparent", linewidth = 1),
        panel.border = element_rect(fill = NA, color = "black", linewidth = 1),
        axis.ticks = element_line(colour = "black", linewidth = 0.5),
        axis.ticks.length = unit(-0.1,"cm"),
        axis.text = element_text(size = 7),
        axis.text.x = element_text(size = 7),
        axis.title = element_text(size = 9),
        plot.margin = grid::unit(c(0,0,0,0), "mm"),
        legend.background = element_rect(fill = NA, colour = NA),
        legend.key = element_blank(),
        legend.position = "bottom",
        legend.text = element_text(size = 7),
        legend.title = element_text(size = 9),
        plot.title = element_text(hjust = 0.5)) +
  facet_grid2(as.factor(Shape) ~ Orientation.fac, scales = "free", axes ="margins") 
syn.diffFt.g2

# save plot as pdf, replace path name with local path name
ggsave("/Users/barrapeak/Desktop/Synthetic.FT.Min.Symmetry.pdf", syn.diffFt.g2, dpi = 300, width = 3.25, height = 7)

## Compare combined FT between different methods for synthetic data ----
# Comparison between new protocol and Ketcham et al. (2011) protocol
# Maximum symmetry case
K.comparison.plot.1 <- all_data_comp %>%
  filter(Ground != "None") %>%
  filter(W.Ratio == 1) %>%
  filter(Volume.x > 0) %>%
  filter(Volume.y > 0) %>%
  filter(Aspect.Ratio == 1) %>%
  ggplot(aes(x = NC.FT, y = Syn.FT, color = Size, xmin = NC.FT - NC.FT.Unc,
             xmax = NC.FT + NC.FT.Unc, ymin = Syn.FT - Syn.FT.Unc,
             ymax = Syn.FT + Syn.FT.Unc, label = Ground)) +
  geom_abline(intercept = 0, slope = 1) +
  geom_hline(yintercept = 0.5) +
  geom_vline(xintercept = 0.5) +
  geom_errorbar(width = 0) +
  geom_errorbarh(width = 0) +
  geom_point(size = 1) +
  geom_text(aes(label = Ground, y = Syn.FT + 0.1), size = 5, size.unit = "pt") +
  scale_color_manual(values = size_color) +
  coord_cartesian(
    xlim = c(0,1),
    ylim = c(0, 1),
    expand = TRUE,
    default = FALSE,
    clip = "on"
  ) +
  scale_x_continuous(name = "Combined FT, Ketcham", breaks = seq(0, 1, 0.1)) +
  scale_y_continuous(name = "Combined FT, This Study", breaks = seq(0, 1, 0.1)) +
  ggtitle(paste("This Study vs. Ketcham et al. (2011)", "Width Ratio = 1, Aspect Ratio = 1", sep = "\n")) +
  theme(panel.background = element_rect(fill = "transparent", colour = "transparent", linewidth = 1),
        panel.border = element_rect(fill = NA, color = "black", linewidth = 1),
        axis.ticks = element_line(colour = "black", linewidth = 0.5),
        axis.ticks.length = unit(-0.1,"cm"),
        axis.text = element_text(size = 7),
        axis.text.x = element_text(size = 7),
        axis.title = element_text(size = 9),
        plot.margin = grid::unit(c(0,0,0,0), "mm"),
        legend.background = element_rect(fill = NA, colour = NA),
        legend.key = element_blank(),
        legend.position = "bottom",
        legend.text = element_text(size = 7),
        legend.title = element_text(size = 9),
        plot.title = element_text(hjust = 0.5)) +
  facet_grid2(as.factor(Shape) ~ Orientation.fac, scales = "fixed")
K.comparison.plot.1

# save plot as pdf, replace path name with local path name
ggsave("/Users/barrapeak/Desktop/Ketcham.1to1.Comp.Max.Symmetry.pdf", K.comparison.plot.1, dpi = 300, width = 3.25, height = 7)

# Comparison between new protocol and Ketcham et al. (2011) protocol
# Minimum symmetry case
K.comparison.plot.2 <- all_data_comp %>%
  filter(Ground != "None") %>%
  filter(W.Ratio == 0.5) %>%
  filter(Volume.x > 0) %>%
  filter(Volume.y > 0) %>%
  filter(Aspect.Ratio == 0.3) %>%
  ggplot(aes(x = NC.FT, y = Syn.FT, color = Size, xmin = NC.FT - NC.FT.Unc,
             xmax = NC.FT + NC.FT.Unc, ymin = Syn.FT - Syn.FT.Unc,
             ymax = Syn.FT + Syn.FT.Unc, label = Ground)) +
  geom_abline(intercept = 0, slope = 1) +
  geom_hline(yintercept = 0.5) +
  geom_vline(xintercept = 0.5) +
  geom_errorbar(width = 0) +
  geom_errorbarh(width = 0) +
  geom_point(size = 1) +
  geom_text(aes(label = Ground, y = Syn.FT + 0.1), size = 5, size.unit = "pt") +
  scale_color_manual(values = size_color) +
  coord_cartesian(
    xlim = c(0,1),
    ylim = c(0, 1),
    expand = TRUE,
    default = FALSE,
    clip = "on"
  ) +
  scale_x_continuous(name = "Combined FT, Ketcham", breaks = seq(0, 1, 0.1)) +
  scale_y_continuous(name = "Combined FT, This Study", breaks = seq(0, 1, 0.1)) +
  ggtitle(paste("This Study vs. Ketcham et al. (2011)", "Width Ratio = 0.5, Aspect Ratio = 0.3", sep = "\n")) +
  theme(panel.background = element_rect(fill = "transparent", colour = "transparent", linewidth = 1),
        panel.border = element_rect(fill = NA, color = "black", linewidth = 1),
        axis.ticks = element_line(colour = "black", linewidth = 0.5),
        axis.ticks.length = unit(-0.1,"cm"),
        axis.text = element_text(size = 7),
        axis.text.x = element_text(size = 7),
        axis.title = element_text(size = 9),
        plot.margin = grid::unit(c(0,0,0,0), "mm"),
        legend.background = element_rect(fill = NA, colour = NA),
        legend.key = element_blank(),
        legend.position = "bottom",
        legend.text = element_text(size = 7),
        legend.title = element_text(size = 9),
        plot.title = element_text(hjust = 0.5)) +
  facet_grid2(as.factor(Shape) ~ Orientation.fac, scales = "fixed")
K.comparison.plot.2

# save plot as pdf, replace path name with local path name
ggsave("/Users/barrapeak/Desktop/Ketcham.1to1.Comp.Min.Symmetry.pdf", K.comparison.plot.2, dpi = 300, width = 3.25, height = 7)

# Comparison between new protocol and Reiners et al. (2007) protocol
# Maximum symmetry case
R.comparison.plot.1 <- all_data_comp %>%
  drop_na(Reiners.FT) %>%
  filter(Volume.x > 0) %>%
  filter(Ground != "None") %>%
  filter(W.Ratio == 0.5) %>%
  filter(Aspect.Ratio == 0.3) %>%
  ggplot(aes(x = Reiners.FT, y = Syn.FT, color = Size,
             ymin = Syn.FT - Syn.FT.Unc,
             ymax = Syn.FT + Syn.FT.Unc, label = Ground)) +
  geom_abline(intercept = 0, slope = 1) +
  geom_hline(yintercept = 0.5) +
  geom_vline(xintercept = 0.5) +
  geom_errorbar(width = 0) +
  geom_point(size = 1) +
  geom_text(aes(label = Ground, y = Syn.FT + 0.1), size = 5, size.unit = "pt") +
  scale_color_manual(values = size_color) +
  coord_cartesian(
    xlim = c(0,1),
    ylim = c(0, 1),
    expand = TRUE,
    default = FALSE,
    clip = "on"
  ) +
  scale_x_continuous(name = "Combined FT, Reiners", breaks = seq(0, 1, 0.1)) +
  scale_y_continuous(name = "Combined FT, This Study", breaks = seq(0, 1, 0.1)) +
  ggtitle(paste("This Study vs. No Correction", "Width Ratio = 1, Aspect Ratio = 1", sep = "\n")) +
  theme(panel.background = element_rect(fill = "transparent", colour = "transparent", linewidth = 1),
        panel.border = element_rect(fill = NA, color = "black", linewidth = 1),
        axis.ticks = element_line(colour = "black", linewidth = 0.5),
        axis.ticks.length = unit(-0.1,"cm"),
        axis.text = element_text(size = 7),
        axis.text.x = element_text(size = 7),
        axis.title = element_text(size = 9),
        plot.margin = grid::unit(c(0,0,0,0), "mm"),
        legend.background = element_rect(fill = NA, colour = NA),
        legend.key = element_blank(),
        legend.position = "bottom",
        legend.text = element_text(size = 7),
        legend.title = element_text(size = 9),
        plot.title = element_text(hjust = 0.5)) +
  facet_grid2(as.factor(Shape) ~ Orientation.fac, scales = "fixed")
R.comparison.plot.1

# save plot as pdf, replace path name with local path name
ggsave("/Users/barrapeak/Desktop/Reiners.1to1.Max.Symmetry.pdf", R.comparison.plot.1, dpi = 300, width = 2, height = 3.5)

# Comparison between new protocol and Reiners et al. (2007) protocol
# Minimum symmetry case
R.comparison.plot.2 <- all_data_comp %>%
  drop_na(Reiners.FT) %>%
  filter(Volume.x > 0) %>%
  filter(Ground != "None") %>%
  filter(W.Ratio == 0.5) %>%
  filter(Aspect.Ratio == 0.3) %>%
  ggplot(aes(x = Reiners.FT, y = Syn.FT, color = Size,
             ymin = Syn.FT - Syn.FT.Unc,
             ymax = Syn.FT + Syn.FT.Unc, label = Ground)) +
  geom_abline(intercept = 0, slope = 1) +
  geom_hline(yintercept = 0.5) +
  geom_vline(xintercept = 0.5) +
  geom_errorbar(width = 0) +
  geom_point(size = 1) +
  geom_text(aes(label = Ground, y = Syn.FT + 0.1), size = 5, size.unit = "pt") +
  scale_color_manual(values = size_color) +
  coord_cartesian(
    xlim = c(0,1),
    ylim = c(0, 1),
    expand = TRUE,
    default = FALSE,
    clip = "on"
  ) +
  scale_x_continuous(name = "Combined FT, Reiners", breaks = seq(0, 1, 0.1)) +
  scale_y_continuous(name = "Combined FT, This Study", breaks = seq(0, 1, 0.1)) +
  ggtitle(paste("This Study vs. No Correction", "Width Ratio = 0.5, Aspect Ratio = 0.3", sep = "\n")) +
  theme(panel.background = element_rect(fill = "transparent", colour = "transparent", linewidth = 1),
        panel.border = element_rect(fill = NA, color = "black", linewidth = 1),
        axis.ticks = element_line(colour = "black", linewidth = 0.5),
        axis.ticks.length = unit(-0.1,"cm"),
        axis.text = element_text(size = 7),
        axis.text.x = element_text(size = 7),
        axis.title = element_text(size = 9),
        plot.margin = grid::unit(c(0,0,0,0), "mm"),
        legend.background = element_rect(fill = NA, colour = NA),
        legend.key = element_blank(),
        legend.position = "bottom",
        legend.text = element_text(size = 7),
        legend.title = element_text(size = 9),
        plot.title = element_text(hjust = 0.5)) +
  facet_grid2(as.factor(Shape) ~ Orientation.fac, scales = "fixed")
R.comparison.plot.2

# save plot as pdf, replace path name with local path name
ggsave("/Users/barrapeak/Desktop/Reiners.1to1.Min.Symmetry.pdf", R.comparison.plot.2, dpi = 300, width = 2, height = 3.5)


#### Real Data Results Comparison

##### Real Data
### Load real data results - new protocol ----
real_data <- read.xlsx("/Users/barrapeak/Mirror/Moved_Files/grain_polishing/Geochron_submission/Reviews/Resubmit_Text/Repository_code/real_data_new_protocol.xlsx", check.names = TRUE) %>%
  mutate(across(-c("Grain", "Crystal.Fragment."), as.numeric)) %>%
  mutate(Shape = case_when(Geometry == 3 ~ paste("Tetragon", Np, "Np"),
                           Geometry == 2 ~ "Cylinder",
                           .default = "Ellipsoid")) %>%
  mutate(Size = case_when(Orientation == 1 ~ Length.2 + Grind.Depth,
                          .default = (Length.1 + Length.2)/2)) %>%
  mutate(o.W1 = Width.1,
         o.W2 = case_when(Orientation == 1 ~ (Width.2 + Length.1)/2,
                          .default = Width.2 + Grind.Depth),
         o.L = case_when(Orientation == 1 ~ Length.2 + Grind.Depth,
                         .default = (Length.1 + Length.2)/2)) %>%
  rowwise() %>%
  mutate(W.Ratio = min(o.W1, o.W2)/max(o.W1, o.W2)) %>%
  mutate(Aspect.Ratio = max(o.W1, o.W2)/o.L) %>%
  mutate(Ground = case_when(Orientation == 2 ~ signif(Grind.Depth/o.W2, 2),
                            .default = signif(Grind.Depth/o.L, 2))) %>%
  mutate(Percent.unc = X..1σ.16/Combined.FT)

# summarize dataset size, aspect ratio, and width ratio
data_summary <- list("Size" = summary(real_data$Size), "Aspect.Ratio" = summary(real_data$Aspect.Ratio), "Width.Ratio" = summary(real_data$W.Ratio))

# summarize dataset geometries
Shape_summary <- real_data %>%
  group_by(Shape, Orientation) %>%
  tally()

### Load real data results - Ketcham et al. (2011) protocol ----
ketcham_real_data <- read.xlsx("/Users/barrapeak/Mirror/Moved_Files/grain_polishing/Geochron_submission/Reviews/Resubmit_Text/Repository_code/real_data_ketcham_protocol.xlsx", check.names = TRUE) %>%
  mutate(across(-c("Grain", "Crystal.Fragment."), as.numeric)) %>%
  mutate(Percent.unc = X..1σ.16/Combined.FT)

real_data_merge <- merge(real_data, ketcham_real_data, by = c(
  "Analysis.Session","Grain", "Length.1", "Width.1", "Length.2", "Width.2", "Geometry", "Np", "Length.P", "Width.P",
  "Orientation", "Crystal.Fragment.", "X4He", "X..1σ", "U", "X..1σ.1","Th", "X..1σ.2", "X147Sm", "X..1σ.3")
) %>%
  select(!c("Crystal.Fragment.")) %>%
  rename("NP.FT" = "Combined.FT.x", "NP.FT.Unc" = "X..1σ.16.x",
         "NP.RFT" = "RFT.x", "NP.RFT.Unc" = "X..1σ.17.x",
         "NP.Percent.Unc" = "Percent.unc.x",
         "NC.FT" = "Combined.FT.y", "NC.FT.Unc" = "X..1σ.16.y",
         "NC.RFT" = "RFT.y", "NC.RFT.Unc" = "X..1σ.17.y",
         "NC.Percent.Unc" = "Percent.unc.y")

### Calculate Reiners et al. (2007) FT for real data ----
a1.U238 <- -4.31 # fit parameters from Farley 2002
a2.U238 <- 4.92
a1.Th232 <- -5
a2.Th232 <- 6.8

intermediate_vals <- real_data_merge %>%
  mutate(r_c = (Width.1 + (Width.2 + Grind.Depth.x))/4) %>%
  mutate(l = (Length.1 + Length.2)/2) %>%
  mutate(r_t1 = Width.1/2) %>%
  mutate(r_t2 = Width.2/2) %>%
  mutate(w = case_when(Grind.Depth.x == r_c ~ 0,
                       .default = acos(1-(Grind.Depth.x/r_c)))) %>%
  mutate(Beta_c = case_when(Grind.Depth.x == 0 ~ (2/l) + (2/(r_c)),
                            .default = (2/l) + 
                              ((2*r_c*(pi-w))/(((pi-w)*(r_c^2)) + 
                                                 ((r_c - Grind.Depth.x)*sqrt((r_c^2) - ((r_c-Grind.Depth.x)^2))))))) %>% # Reiners 2007 cylinder eq (modified for no terminations)
  mutate(Beta_t = case_when(Grind.Depth.x == 0 ~ (2/l) + (1/r_t1) + (1/r_t2),
                            .default = (2/l) + (1/(2*r_t2)) + (1/r_t1)))

all_real_data <- cbind(real_data_merge, select(intermediate_vals, c("Beta_c", "Beta_t"))) %>%
  mutate(Beta = case_when((Orientation == 2 & 
                             Geometry == 2 & 
                             between(Grind.Depth.x, 0.1, (Grind.Depth.x + Width.2)/2)) ~ Beta_c,
                          (Orientation == 2 & 
                             Geometry == 2 & 
                             Grind.Depth.x == 0) ~ Beta_c,
                          (Orientation == 2 & 
                             Geometry == 3 & 
                             between(Grind.Depth.x, 0.1, (Grind.Depth.x + Width.2)/2) &
                             Np == 0) ~ Beta_t,
                          (Orientation == 2 & 
                             Geometry == 3 & 
                             Grind.Depth.x == 0 &
                             Np == 0) ~ Beta_t,
                          .default = NA)) %>%
  select(!c("Beta_c", "Beta_t")) %>%
  mutate(FT_238U = 1 + (a1.U238*Beta) + (a2.U238*(Beta^2)),
         FT_232Th = 1 + (a1.Th232*Beta) + (a2.Th232*(Beta^2))) %>%
  mutate(A = 1/(1.04 + (0.245*(Th/U)))) %>%
  mutate(FTc = (A*FT_238U) + ((1-A)*FT_232Th))

# save all real data method comparisons as excel file, replace path name with local path name
write.xlsx(all_real_data, "/Users/barrapeak/Desktop/all_real_data_comparisons.xlsx")

### Plot real data comparisons ----
real_size_color <- c("Small" = "#332288", "Medium" = "#225555", "Large" = "#997700", "Extra Large" = "black") # color scale across plots

# prepare data frome for plotting
all_real_data_comp <- all_real_data %>%
  rename("Reiners.FT" = "FTc") %>%
  group_by(Shape, Size, Aspect.Ratio, W.Ratio, Orientation) %>%
  arrange(Shape, Size, Aspect.Ratio, W.Ratio, Orientation) %>%
  mutate(Syn.Diff.From.Ketcham = ((NP.FT/NC.FT) -1)*100) %>%
  mutate(Syn.Diff.From.Reiners = ((NP.FT/Reiners.FT) -1)*100) %>%
  mutate(Shape.fac = factor(Shape, levels = c("Ellipsoid", "Cylinder",
                                              "Tetragon 0 Np", "Tetragon 1 Np"))) %>%
  mutate(Size.bins = case_when(between(Size, 100, 150) ~ "Small",
                               between(Size, 150, 200) ~ "Medium",
                               between(Size, 200, 250) ~ "Large",
                               .default = "Extra Large")) %>%
  mutate(Size.fac = factor(Size.bins, levels = c("Small",
                                        "Medium",
                                        "Large", "Extra Large"))) %>%
  mutate(Orientation.fac = case_when(Orientation == 1 ~ "Perpendicular",
                                     .default = "Parallel"))

# summary of number of each geometry present in real Reiners FT dataset
Shape_summary_Reiners <- all_real_data_comp %>%
  drop_na(Reiners.FT) %>%
  group_by(Shape, Orientation) %>%
  tally()

## Compared combined FT between different methods for real data ----
# Combined FT comparison between new protocol and Ketcham et al. (2011)
K.real.comparison.plot <- all_real_data_comp %>%
  filter(Volume.x > 0) %>%
  filter(Volume.y > 0) %>%
  ggplot(aes(x = NC.FT, y = NP.FT, color = Size.fac, xmin = NC.FT - NC.FT.Unc,
             xmax = NC.FT + NC.FT.Unc, ymin = NP.FT - NP.FT.Unc,
             ymax = NP.FT + NP.FT.Unc, label = as.character(Ground))) +
  geom_abline(intercept = 0, slope = 1) +
  geom_hline(yintercept = 0.5) +
  geom_vline(xintercept = 0.5) +
  geom_errorbar(width = 0) +
  geom_errorbarh(width = 0) +
  geom_point(size = 1) +
  geom_text(aes(label = as.character(Ground), y = NP.FT + 0.1), size = 5, size.unit = "pt") +
  scale_color_manual(values = real_size_color) +
  coord_cartesian(
    xlim = c(0,1),
    ylim = c(0, 1),
    expand = TRUE,
    default = FALSE,
    clip = "on"
  ) +
  scale_x_continuous(name = "Combined FT, Ketcham", breaks = seq(0, 1, 0.1)) +
  scale_y_continuous(name = "Combined FT, This Study", breaks = seq(0, 1, 0.1)) +
  ggtitle(paste("This Study vs. Ketcham et al. (2011)", "All Sizes", sep = "\n")) +
  theme(panel.background = element_rect(fill = "transparent", colour = "transparent", linewidth = 1),
        panel.border = element_rect(fill = NA, color = "black", linewidth = 1),
        axis.ticks = element_line(colour = "black", linewidth = 0.5),
        axis.ticks.length = unit(-0.1,"cm"),
        axis.text = element_text(size = 7),
        axis.text.x = element_text(size = 7),
        axis.title = element_text(size = 9),
        plot.margin = grid::unit(c(0,0,0,0), "mm"),
        legend.background = element_rect(fill = NA, colour = NA),
        legend.key = element_blank(),
        legend.position = "bottom",
        legend.text = element_text(size = 7),
        legend.title = element_text(size = 9),
        plot.title = element_text(hjust = 0.5)) +
  facet_grid2(Shape.fac ~ Orientation.fac, scales = "fixed")
K.real.comparison.plot

# save plot as pdf, replace path name with local path name
ggsave("/Users/barrapeak/Desktop/real_data_Ketcham.1to1.pdf", K.real.comparison.plot, dpi = 300, width = 3.25, height = 5.8)

# Combined FT comparison between new protocol and Reiners et al. (2007)
R.real.comparison.plot <- all_real_data_comp %>%
  drop_na(Reiners.FT) %>%
  filter(Volume.x > 0) %>%
  ggplot(aes(x = Reiners.FT, y = NP.FT, color = Size.fac,
             ymin = NP.FT - NP.FT.Unc,
             ymax = NP.FT + NP.FT.Unc, label = as.character(Ground))) +
  geom_abline(intercept = 0, slope = 1) +
  geom_hline(yintercept = 0.5) +
  geom_vline(xintercept = 0.5) +
  geom_errorbar(width = 0) +
  geom_point(size = 1) +
  geom_text(aes(label = as.character(Ground), y = NP.FT + 0.1), size = 5, size.unit = "pt") +
  scale_color_manual(values = real_size_color) +
  coord_cartesian(
    xlim = c(0,1),
    ylim = c(0, 1),
    expand = TRUE,
    default = FALSE,
    clip = "on"
  ) +
  scale_x_continuous(name = "Combined FT, Reiners", breaks = seq(0, 1, 0.1)) +
  scale_y_continuous(name = "Combined FT, This Study", breaks = seq(0, 1, 0.1)) +
  ggtitle(paste("This Study vs. No Correction", "Width Ratio = 0.5, Aspect Ratio = 0.3", sep = "\n")) +
  theme(panel.background = element_rect(fill = "transparent", colour = "transparent", linewidth = 1),
        panel.border = element_rect(fill = NA, color = "black", linewidth = 1),
        axis.ticks = element_line(colour = "black", linewidth = 0.5),
        axis.ticks.length = unit(-0.1,"cm"),
        axis.text = element_text(size = 7),
        axis.text.x = element_text(size = 7),
        axis.title = element_text(size = 9),
        plot.margin = grid::unit(c(0,0,0,0), "mm"),
        legend.background = element_rect(fill = NA, colour = NA),
        legend.key = element_blank(),
        legend.position = "bottom",
        legend.text = element_text(size = 7),
        legend.title = element_text(size = 9),
        plot.title = element_text(hjust = 0.5)) +
  facet_grid2(Shape.fac ~ Orientation.fac, scales = "fixed")
R.real.comparison.plot

# save plot as pdf, replace path name with local path name
ggsave("/Users/barrapeak/Desktop/real_data_Reiners.1to1.pdf", R.real.comparison.plot, dpi = 300, width = 2, height = 3.5)

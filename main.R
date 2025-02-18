library(yaml)

# Load config vars
config <- yaml.load_file("./config.yaml", eval.expr = TRUE)

NZSCC <- raster(config$data$nzscc.groups)
BI.richness <- raster(config$data$richness$BI)
DF.richness <- raster(config$data$richness$DF)
RF.richness <- raster(config$data$richness$RF)
MA.richness <- raster(config$data$richness$MA)
hoki.boffff <- raster(config$data$boffff)
VME.index <- raster(config$data$vme)

# Hoki Boffff layer cells to 1 where gt than 0
hoki.boffff[hoki.boffff > 0] <- 1
# reproject VME layer
VME.index <- projectRaster(from = VME.index,
                           to = hoki.boffff,
                           crs = crs(hoki.boffff))

# NZSCC Groups dataframe (extent and percentage)
DF.Groups <- setNames(
    data.frame(matrix(ncol = 4, nrow = (75))),
    c(
        "GF.grp",
        "Group.extent",
        "Group.percentage",
        "Between-Group.Sim.percentage"
    )
)
DF.Groups$GF.grp <- seq(1:75)

DF.Overlap <- setNames(
    data.frame(matrix(ncol = 8, nrow = 1)),
    c(
        "BI.rich.percentage",
        "DF.rich.percentage",
        "RF.rich.percentage",
        "MA.rich.percentage",
        "VME.0.percentage",
        "VME.1.percentage",
        "VME.2.percentage",
        "VME.3.percentage"
    )
)

# For all NZSCC groups
for (i in 1:nrow(DF.Groups)) {
    # load group within similarity layer
    within.sim <- raster(paste(
        config$data$nzscc.intra.sim,
        config$sim.intra.prefix,
        i,
        ".tif",
        sep = ""
    ))
    # Isolate NZSCC group
    Grp.R <- NZSCC
    Grp.R[Grp.R != i] <- 0
    Grp.R[Grp.R == i] <- 1
    
    # Overlap with Hoki BOFFFF
    P <- Grp.R * hoki.boffff
    DF.Groups[i, 2] <- sum(na.omit(values(P)))
    DF.Groups[i, 3] <- round((sum(na.omit(values(
        P
    ))) / sum(na.omit(values(
        Grp.R
    )))) * 100, 1)
    print(paste(
        "Finished calculating hoki BOFFFF overlap with NZSCC Group: ",
        i
    ))
    
    # Overlap with group within similarity
    within.sim[within.sim == 0] <- NA
    Intra.P <- within.sim * hoki.boffff
    
    # percentage
    DF.Groups[i, 4] <- round((sum(na.omit(
        values(Intra.P)
    )) / sum(na.omit(values(
        within.sim
    )))) * 100, 1)
    
}

# Overlap with Richness
j <- 1
for (richness in c(BI.richness, DF.richness, RF.richness, MA.richness)) {
    P <- richness * hoki.boffff
    
    DF.Overlap[1, j] <- round((sum(na.omit(values(
        P
    ))) / sum(na.omit(
        values(richness)
    ))) * 100, 1)
    j <- j + 1
}

# Overlap with VME index (0 to 3)

for (i in 0:3) {
    VME.subset <- VME.index
    VME.subset[VME.subset != i] <- 0
    VME.subset[VME.subset == i] <- 1
    
    P <- VME.subset * hoki.boffff
    DF.Overlap[1, j] <- round((sum(na.omit(values(
        P
    ))) / sum(na.omit(
        values(VME.subset)
    ))) * 100, 1)
    j <- j + 1
}

write.csv(DF.Groups, sep = ";", file = "./output/results_NZSCC_overlap.csv")
write.csv(DF.Overlap, sep = ";", file = "./output/results_richness_VME_overlap.csv")
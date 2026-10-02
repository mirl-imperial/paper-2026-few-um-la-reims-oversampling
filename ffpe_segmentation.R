library(rawsewage)
library(ggplot2)
library(imager)
library(dplyr)
# recent versions of genieclust have bugs with certain kinds of data
library(lumbermark)
scale_values <- function(x){(x-min(x))/(max(x)-min(x))}

normalise_raster <- function(input_df, image_dims, method = c("mean", "median"), offset_to_zero = FALSE) {
  method <- match.arg(method)
  # Preallocate output matrix
  normalised_df <- matrix(data = NA, nrow = nrow(input_df), ncol = ncol(input_df))
  for (i in 1:ncol(input_df)) {
    #cat("\r", "Processing column", i, "of", ncol(input_df))
    # Convert column to matrix representing an image
    image <- matrix(data = input_df[, i], nrow = image_dims[1], ncol = image_dims[2])
    # Preallocate vectors to store means/medians for each raster (column)
    raster_values <- numeric(image_dims[2])
    for (z in 1:image_dims[2]) {
      # Compute the mean or median for all pixels in the current raster
      if (method == "mean") {
        raster_values[z] <- mean(image[, z])
      } else {
        raster_values[z] <- median(image[, z])
      }
    }
    # Normalise each raster by subtracting its mean/median
    for (z in 1:image_dims[2]) {
      image[, z] <- image[, z] - raster_values[z]
    }
    # Optionally offset values so they are all non-negative
    if (offset_to_zero) {
      min_value <- min(image)
      image <- image + abs(min_value)
    }
    # Flatten the normalised image back into a column and store it
    normalised_df[, i] <- as.vector(image)
  }
  return(normalised_df)
}

# Paths: set RAW_DIR, PEAKS_TXT and OUT_PDF before source()-ing this file, or run
#   Rscript ffpe_segmentation.R <MS19 .raw dir> <HDI 1000-peak .txt> [cluster_plot.pdf]
if (!exists("RAW_DIR")) {
  args <- commandArgs(trailingOnly = TRUE)
  RAW_DIR <- args[1]
  PEAKS_TXT <- args[2]
  OUT_PDF <- if (length(args) >= 3) args[3] else "cluster_plot.pdf"
}

dims <- estDims(RAW_DIR)

# peak picked via hdi with 1000 peaks
dat <- read.csv(PEAKS_TXT, sep = "\t", skip = 3, header = FALSE)

mz <- dat[1,4:ncol(dat)]
dat <- dat[2:nrow(dat),4:(ncol(dat)-2)]

# "TIC" (but not really) of sum of peak picked data
cimgplot <- matrix(data = rowSums(dat), nrow = dims[1], ncol = dims[2])
cimgplot <- imager::as.cimg(cimgplot, x = dims[1], y = dims[2])
p <- ggplot(data = as.data.frame(cimgplot)) +
  geom_raster(aes(x = x, y = rev(y), fill = (value))) +
  scale_fill_viridis_c(option = "magma", name = "Intensity (A.U)") + coord_fixed() +
  scale_x_continuous(expand=c(0,0))+scale_y_continuous(expand=c(0,0)) + xlab(NULL) + ylab(NULL) + ggtitle("TIC Image") +
  theme(plot.title = element_text(hjust = 0.5))
p

library(irlba)
set.seed(12345)
pca_obj <- prcomp_irlba(log1p(dat), n = 30)
pca_img <- matrix(data = pca_obj$x, nrow = prod(dims), ncol = ncol(pca_obj$x))
colnames(pca_img) <- colnames(pca_obj$x)

pca_img_df <- imager::as.cimg(array(1,c(dims[1],dims[2])))
pca_img_df <- as.data.frame(pca_img_df)
# and bind
pca_img_df <- cbind(pca_img_df[1:nrow(pca_img),1:2], pca_img)

p <- ggplot(data = pca_img_df) +
  geom_raster(aes(x = x, y = y, fill = PC1)) +
  scale_fill_viridis_c(option = "magma") + coord_fixed() + scale_x_continuous(expand=c(0,0))+scale_y_continuous(expand=c(0,0), trans=scales::reverse_trans()) + xlab(NULL) + ylab(NULL) + theme(plot.title = element_text(hjust = 0.5))
p
p <- ggplot(data = pca_img_df) +
  geom_raster(aes(x = x, y = y, fill = PC2)) +
  scale_fill_viridis_c(option = "magma") + coord_fixed() + scale_x_continuous(expand=c(0,0))+scale_y_continuous(expand=c(0,0), trans=scales::reverse_trans()) + xlab(NULL) + ylab(NULL) + theme(plot.title = element_text(hjust = 0.5))
p
p <- ggplot(data = pca_img_df) +
  geom_raster(aes(x = x, y = y, fill = PC3)) +
  scale_fill_viridis_c(option = "magma") + coord_fixed() + scale_x_continuous(expand=c(0,0))+scale_y_continuous(expand=c(0,0), trans=scales::reverse_trans()) + xlab(NULL) + ylab(NULL) + theme(plot.title = element_text(hjust = 0.5))
p
df_log_pca_cac <- pca_img_df %>% mutate(rgb.val=rgb(scale_values(PC1),scale_values(PC2),scale_values(PC3)))
p <- ggplot(df_log_pca_cac,aes(x,y))+geom_raster(aes(fill=rgb.val))+scale_fill_identity()
p +coord_fixed()+scale_x_continuous(expand=c(0,0))+scale_y_continuous(expand=c(0,0), trans=scales::reverse_trans())+
  xlab(NULL) + ylab(NULL) + ggtitle("PCA image (PCs 1-3)") +
  theme(plot.title = element_text(hjust = 0.5))

# estimate tissue mask so as not to incorporate variance from the background
set.seed(12345)
tree <- genieclust::gclust(pca_img)
cut <- cutree(tree, k = 2)

cimgplot <- matrix(data = cut, nrow = dims[1], ncol = dims[2])
cimgplot <- imager::as.cimg(cimgplot)
cimgplot <- imager::mirror(cimgplot, "y")
# tinker with rotations
p <- ggplot(data = as.data.frame(cimgplot)) +
  geom_raster(aes(x = x, y = y, fill = factor(value))) +
  scale_fill_manual(values = c("black", scales::hue_pal()(nlevels(as.factor(cut))))) + coord_fixed()
p

# select relevant cluster
tissuepix <- which(cut == 2)

rastnorm <- normalise_raster(log1p(dat), dims, method = "mean")

rastnorm <- expm1(rastnorm)

clr_df <- vegan:::.calc_clr(rastnorm+1, na.rm = TRUE)

pca_obj <- prcomp_irlba((clr_df[tissuepix,]), n = 10)
pca_img <- matrix(data = 0, nrow = prod(dims), ncol = ncol(pca_obj$x))
colnames(pca_img) <- colnames(pca_obj$x)
pca_img[tissuepix,] <- pca_obj$x

pca_img_df <- imager::as.cimg(array(1,c(dims[1],dims[2])))
pca_img_df <- as.data.frame(pca_img_df)
# and bind
pca_img_df <- cbind(pca_img_df[1:nrow(pca_img),1:2], pca_img)

p <- ggplot(data = pca_img_df) +
  geom_raster(aes(x = x, y = y, fill = PC1)) +
  scale_fill_viridis_c(option = "magma") + coord_fixed() + scale_x_continuous(expand=c(0,0))+scale_y_continuous(expand=c(0,0), trans=scales::reverse_trans()) + xlab(NULL) + ylab(NULL) + theme(plot.title = element_text(hjust = 0.5))
p
p <- ggplot(data = pca_img_df) +
  geom_raster(aes(x = x, y = y, fill = PC2)) +
  scale_fill_viridis_c(option = "magma") + coord_fixed() + scale_x_continuous(expand=c(0,0))+scale_y_continuous(expand=c(0,0), trans=scales::reverse_trans()) + xlab(NULL) + ylab(NULL) + theme(plot.title = element_text(hjust = 0.5))
p
p <- ggplot(data = pca_img_df) +
  geom_raster(aes(x = x, y = y, fill = PC3)) +
  scale_fill_viridis_c(option = "magma") + coord_fixed() + scale_x_continuous(expand=c(0,0))+scale_y_continuous(expand=c(0,0), trans=scales::reverse_trans()) + xlab(NULL) + ylab(NULL) + theme(plot.title = element_text(hjust = 0.5))
p
df_log_pca_cac <- pca_img_df %>% mutate(rgb.val=rgb(scale_values(PC1),scale_values(PC2),scale_values(PC3)))
p <- ggplot(df_log_pca_cac,aes(x,y))+geom_raster(aes(fill=rgb.val))+scale_fill_identity()
p +coord_fixed()+scale_x_continuous(expand=c(0,0))+scale_y_continuous(expand=c(0,0), trans=scales::reverse_trans())+
  xlab(NULL) + ylab(NULL) + ggtitle("PCA image (PCs 1-3)") +
  theme(plot.title = element_text(hjust = 0.5))

screeplot(pca_obj)
library(uwot)

temp_umap <- umap(pca_obj$x[,1:3], n_components = 2, metric = "cosine", seed = 12345)

cut <- lumbermark(temp_umap, k = 3, distance = "cosine")

tempalette <- c("#5CB047", "#E41A1C","#377EB8","#FC8D62","#E6C948","#A65628","#984EA3",
                "#E78AC3","#4DAF4A","#6A3D9A","#FF7F00")

cimgplot <- matrix(data = 0, nrow = dims[1], ncol = dims[2])
cimgplot[tissuepix] <- cut
cimgplot <- imager::as.cimg(cimgplot)
cimgplot <- imager::mirror(cimgplot, "y")
p <- ggplot(data = as.data.frame(cimgplot)) +
  geom_raster(aes(x = x, y = y, fill = factor(value))) +
  scale_fill_manual(values = c("black", tempalette)) + coord_fixed() +
  theme_void() + guides(fill = guide_legend(title = "Cluster", override.aes = list(size = 3)))

pdf(OUT_PDF, height = 8, width = 8)
p
dev.off()
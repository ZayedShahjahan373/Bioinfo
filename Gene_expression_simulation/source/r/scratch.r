# Install and load required packages
install.packages("glmnet")
library(glmnet)


# Example data (replace this with your actual dataset)
set.seed(123)
n <- 100
p <- 10
strat <- rep(1:5, each = n/5)
x <- matrix(rnorm(n * p), ncol = p)
y <- rnorm(n)

# Combine data into a data frame
data <- data.frame(strat, x, y)


# Function to fit elastic net model with k-fold cross-validation
fit_enet_cv <- function(data, indices, k = 5) {
  # Subset data based on bootstrap indices
  boot_data <- data[indices, ]
  
  # Extract predictors and response
  x_boot <- boot_data[, -1]  # Assuming the first column is the stratifying variable
  y_boot <- boot_data[, 1]   # Assuming the first column is the response variable
  
  
  
  # Check if the response variable has variability
  if (length(unique(y_boot)) <= 1) {
    # If y is constant, return NA or an appropriate placeholder
    return(rep(NA, length(y_boot)))
  }
  
  # Create a fold index vector for cross-validation
  fold_indices <- sample(1:k, length(y_boot), replace = TRUE)
  
  # Initialize an empty vector to store the cross-validated predictions
  cv_preds <- numeric(length = length(y_boot))
  
  # Perform k-fold cross-validation
  for (fold in 1:k) {
    # Split data into training and validation sets
    x_train <- x_boot[fold_indices != fold, ]
    y_train <- y_boot[fold_indices != fold]
    x_valid <- x_boot[fold_indices == fold, ]
    
    # Fit elastic net model on training data
    enet_model <- cv.glmnet(x_train, y_train, alpha = 0.5)
    
    # Make predictions on the validation set
    cv_preds[fold_indices == fold] <- predict(enet_model, newx = x_valid, s = "lambda.min")
  }
  
  # Return cross-validated predictions
  return(cv_preds)
}

# Set up bootstrap
n_bootstrap <- 1000  # Number of bootstrap samples
bootstrap_indices <- unlist(lapply(1:n_bootstrap, 
                            function(i) sample(1:n, replace = TRUE)))

# Perform bootstrap and fit elastic net models with k-fold CV
boot_results <- lapply(bootstrap_indices, fit_enet_cv, data = data)

# Display results (for example, predictions for the first observation in each bootstrap sample)
head(sapply(boot_results, function(cv_preds) cv_preds[1]))

#retrieve the first element of the simulated gene expression dataset 
gnex.dir = file.path(anaDir,"01")

gnex.df.list = readRDS(file.path(gnex.dir,"bile_sim_gnexp.RDS"))
sim.gene.exp = gnex.df.list[["x"]]

#Scaled data on which pca is to be performed
dta.pca = apply(sim.gene.exp,2,scale)

pca.res = prcomp(dta.pca)
pca.res$pvar=round(pca.res$sdev)
#Store only first 3 

pca.x = as.data.frame(pca.res$x)
pca.rotation = as.data.frame(pca.res$rotation)
pca.rotation$GENE.SYM = rownames(pca.rotation)


rownames(pca.x) = rownames(sim.gene.exp)
pca.x$SAMID = rownames(pca.x)

#Keep first 2 PCs and the sample IDs
pca.x = pca.x[,c("SAMID","PC1","PC2")]
#merge with 
pca.mta = merge(pca.x,mta.df,by="SAMID")

ggplot(pca.mta,aes(PC1,PC2,color=treatType))+
  geom_point(alpha=0.5)+
  stat_ellipse(type = "norm",level = 0.99)+
  theme_bw()
#Arrange the loadings from largest to smallest in the pca rotation object
ggplot(pca.rotation,aes(PC1,GENE.SYM))+geom_bar(stat = "identity")

pca.var.expl = data.frame(Var_exp=round(pca.res$sdev^2/sum(pca.res$sdev^2)*100,1),
                             components=colnames(pca.res$rotation))
pca.var.expl$components = factor(pca.var.expl$components,
                                 levels = pca.var.expl[order(pca.var.expl$Var_exp),"components"])

ggplot(pca.var.expl,aes(components,Var_exp))+
  geom_segment(aes(x=components,xend=components,y=0,yend=Var_exp),color="dodgerblue")+
  geom_point()+
  coord_flip()+
  labs(x="Principal Components",
       y="Variance Explained",
       title = "Gene Expression Principal Components")+theme_bw()

summary(pca.res)

gen.dist=dist(dta.pca,method =)
plot(as.dendrogram(hclust(gen.dist)),)
 
###GGPLOT2

#Why use ggplot2:
#Without a grammar, you can only do what you know.
#For each plot, learn a different function.
#But you can't create novel plots for new visualizations

#BAsic structure of a plot. Axes+mapping+shape
#Axes+points+shapes=easthetics (things we can perceive on plot)
#going from data to visual is the mapping 
#kind of plot -> geometry of the plot
#Apply a scale transformation to print into the window of the plot
#Statistical transformation.
#group observations into ranges=binning
#Coordinate system.
#Polar coordinates for CCA plot
#Facet to groups according to comparison variables

#Mapping from data to aesthetics

#Components of the Grammar of Graphics:
# A plot needs:
#   A default dataset+A mapping from data to aesthetics+
#   At least one layer:
#     -One geometry
#     -One statistical transformations
#     -One position adjustment
#     -Secondary mapping for datapoints
#   One scale for each aesthetic mapping
#   One coordinate system
#   One facet specification

library(GenomicRanges)


library(readxl)

kohlman <- read_excel("C:/Users/momom/Downloads/PICharterRecord_62000600_2012-06-21.xlsx", sheet = 'Kohlman data refined')
View(kohlman)
colnames(kohlman)[1] = "dates"
colnames(kohlman)
par(mfrow=c(1,1))


# p.crispus x m. spicatum
plot(kohlman$potamogeton_crispus, kohlman$myriophyllum_spicatum, xlab = 'P. crispus', ylab = 'M. spicatum Population')
myr.crisp.mod <- lm(kohlman$myriophyllum_spicatum~kohlman$potamogeton_crispus)
abline(myr.crisp.mod, col = 'red')
plot(myr.crisp.mod)

plot(kohlman$potamogeton_crispus, kohlman$myriophyllum_spicatum, xlab = 'P. crispus', ylab = 'M. spicatum Population')
myr.crisp.mod <- lm(kohlman$myriophyllum_spicatum~kohlman$potamogeton_crispus)
abline(myr.crisp.mod, col = 'red')
plot(myr.crisp.mod)


# shows competition
plot(kohlman$elodea_canadensis, kohlman$myriophyllum_spicatum, xlab = 'E. canadensis', ylab = 'M. spicatum Population')
myr.elo.can.mod <- lm(kohlman$myriophyllum_spicatum~kohlman$elodea_canadensis)
abline(myr.elo.can.mod, col = 'red')
summary(myr.elo.can.mod)
# y-intercept -> 49.9580
# x-intercept -> 115.6971

# shows competition
plot(kohlman$potamogeton_pusillus, kohlman$myriophyllum_spicatum, xlab = 'P. fusillus', ylab = 'M. spicatum Population')
myr.pot.fus.mod <- lm(kohlman$myriophyllum_spicatum~kohlman$potamogeton_pusillus)
abline(myr.pot.fus.mod, col = 'red')

plot(kohlman$potamogeton_pusillus, kohlman$myriophyllum_spicatum, xlab = 'P. fusillus', ylab = 'M. spicatum Population')
myr.pot.fus.mod <- lm(kohlman$myriophyllum_spicatum~kohlman$potamogeton_pusillus)
abline(myr.pot.fus.mod, col = 'red')




silver <- read.csv("C:/Users/momom/Downloads/silver data.csv")
# Amy's code
#install.packages("glmnet")
library(glmnet)

setwd("~/Desktop/UROP_Momo")
plant_data <- read.csv("UROP_lake_data_refined.csv")
plant_data <- silver[,-c(1,42,41)]
# Load data (assuming `plant_data` is your abundance matrix)
head(plant_data)
plant_matrix <- as.matrix(plant_data)
View(plant_matrix) # -> exclude dates column for the surveys
# for .csv -> values stored in [[x,y,z], [a,b,c]] ? to loop through indices

# Fit LASSO regression for each taxon
network_edges <- list()
for (j in 1:ncol(plant_matrix)) {
  for (i in 41:nrow(plant_matrix)) {
    target <- plant_matrix[, j]  
    predictors <- plant_matrix[, -j]  
    
    # Fit LASSO
    lasso_fit <- cv.glmnet(predictors, target, alpha = .5)# how does alpha impact the outputs?
    
    # Extract coefficients
    coef_values <- coef(lasso_fit, s = "lambda.min")
    
    # Store nonzero coefficients as a named vector
    nonzero_coefs <- coef_values[coef_values[,1] != 0, , drop = FALSE]
    network_edges[[colnames(plant_matrix)[j]]] <- 
      setNames(as.list(nonzero_coefs[-1]), colnames(plant_matrix)[-j][coef_values[-1] != 0])
  }
}

network_edges

library(igraph)

# Convert the network_edges list to a data frame

edges_df <- do.call(rbind, lapply(names(network_edges), function(target) {
  if (length(network_edges[[target]]) > 0) {
    # Create the data frame with all edges
    all_edges <- data.frame(
      from = target,
      to = names(network_edges[[target]]),
      weight = unlist(network_edges[[target]])
    )
    
    # Filter to only include edges below the cutoff
    cutoff_value <- 0.5  # Set your desired cutoff here
    filtered_edges <- all_edges[all_edges$weight < cutoff_value, ]
    
    # Return the filtered data frame
    if (nrow(filtered_edges) > 0) {
      return(filtered_edges)
    }
  }
}))

edges_df <- do.call(rbind, lapply(names(network_edges), function(target) {
  if (length(network_edges[[target]]) > 0) {
    data.frame(
      from = target,
      to = names(network_edges[[target]]),
      weight = unlist(network_edges[[target]])
    )
  }
}))

# Print to verify
print(edges_df)

# Create an undirected graph from the data frame
g <- graph_from_data_frame(edges_df, directed = FALSE)

# Assign edge weights
E(g)$width <- E(g)$weight * 0.5  # Scale edge thickness for visibility

# Label network
label_map <- c('potamogeton_pusillus' = 'P. crispus', 'elodea_canadensis' = 'E. canadensis', 
               'potamogeton_foliosus' = 'P. foliosus', 'nuphar_advena' = 'N. advena',
               'myriophyllum_spicatum' = 'M. spicatum', 'elodea_sp' = 'Elodea species',
               'nitella_sp' = 'Nitella species', 'chara_sp' = 'Chara species')

# Plot the network
plot(g,
     edge.width = E(g)$width,  # Edge thickness represents weight
     vertex.size = 10,         # Node size
     vertex.label = label_map[V(g)$name], 
     vertex.label.cex = 1.8,   # Label size
     vertex.color = "lightgreen",
     edge.color = "gray50",
     vertex.label.family = 'sans',
     main = "Weighted Co-occurrence Network")

plot(g,
     layout = layout_with_fr,  # Force-directed layout
     edge.width = E(g)$weight * 10,  
     vertex.size = 15,
     vertex.label.cex = 1,
     vertex.color = "lightgreen",
     edge.color = "black",
     vertex.label.dist = 4, 
     main = "Weighted Network of Plant Co-occurrence")




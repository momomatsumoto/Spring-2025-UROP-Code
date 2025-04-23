# Load dataset where each row is an observation and each column is a taxon
kohlman <- read.csv("C:/Users/momom/Downloads/kohlman data.csv")

plant_data <- kohlman[,-1]
# assign lake data to variable plant_data 
# remove column 1 (dates will give NA error)

install.packages("glmnet")
library(glmnet)

# Load data (assuming 'plant_data' is your abundance matrix)
plant_matrix <- as.matrix(plant_data)
View(plant_matrix)

# Fit LASSO regression for each taxon
network_edges <- list()
for (j in 1:ncol(plant_matrix)) {
  for (i in 41:nrow(plant_matrix)) {
    target <- plant_matrix[, j]  
    predictors <- plant_matrix[, -j]  
    
    # Fit LASSO
    lasso_fit <- cv.glmnet(predictors, target, alpha = .5)
    
    # Extract coefficients
    coef_values <- coef(lasso_fit, s = "lambda.min")
    
    # Store nonzero coefficients as a named vector
    nonzero_coefs <- coef_values[coef_values[,1] != 0, , drop = FALSE]
    network_edges[[colnames(plant_matrix)[j]]] <- 
      setNames(as.list(nonzero_coefs[-1]), colnames(plant_matrix)[-j][coef_values[-1] != 0])
  }
}


# Convert the network_edges list to a data frame
edges_df <- do.call(rbind, lapply(names(network_edges), function(target) {
  if (length(network_edges[[target]]) > 0) {
    data.frame(
      from = target,
      to = names(network_edges[[target]]),
      weight = unlist(network_edges[[target]])
    )
  }
}))

install.packages("igraph")
library(igraph)

# Assign positive or negative relationships
edges_df$relationship<-ifelse(edges_df$weight>0, "positive", "negative")
edges_df$weight_orig<-edges_df$weight
edges_df$weight<-abs(edges_df$weight) # assign the absolute variables of the weights

# Print to verify it is still sender, receiver, weight
print(edges_df)

# Assign a 'sign' attribute

# Create an undirected graph from the data frame
g <- graph_from_data_frame(edges_df, directed = FALSE)

# Create a directed graph from the data frame
g2 <- graph_from_data_frame(edges_df, directed = TRUE)

# Set edge color based on 'sign'
E(g2)$color <- ifelse(E(g2)$relationship == "positive", "lightgreen", "red")

# Assign edge weights
E(g)$width <- E(g)$weight * 10 # Scale edge thickness for visibility

# Label map
label_map <- c('potamogeton_pusillus' = 'P. pusillus', 'elodea_canadensis' = 'E. canadensis', 
               'potamogeton_foliosus' = 'P. foliosus', 'nuphar_advena' = 'N. advena',
               'myriophyllum_spicatum' = 'M. spicatum', 'elodea_sp' = 'Elodea species',
               'nitella_sp' = 'Nitella species', 'chara_sp' = 'Chara species', 'algae' = 'algae',
               'myriophyllum_sibiricum' = 'M. sibiricum', 'nymphaea_odorata' = 'N. odorata',
               'eleocharis_acicularis' = 'E. acicularis', 'elodea_ca0densis' = 'E. canadensis',
               'ceratophyllum_demersum' = 'C. demersum', 'characeae_taxa' = 'Characeae',
               'potamogeton_crispus' = 'P. crispus')


# Plot the network

# Specify layout
my_layout_g <- layout_in_circle(g)
my_layout_g2 <- layout_in_circle(g2)

# Undirected, weighted network
plot(g,
     layout = my_layout_g,
     edge.width = E(g)$width, 
     vertex.size = 10,
     vertex.label = label_map[V(g)$name], 
     vertex.label.cex = 1,  
     vertex.color = "darkseagreen1",
     vertex.label.family = 'serif',
     edge.color = "gray50",
     vertex.label.color = "black",
     vertex.frame.color = 'darkseagreen4',
     edge.curved = 0.3,
     main = "Weighted Plant Co-occurrence Network")

# Directed, unweighted network
plot(g2,
     layout = my_layout_g2, 
     edge.width = 3,  
     vertex.size = 10,
     vertex.label = label_map[V(g)$name],
     vertex.label.family = 'serif',
     vertex.label.cex = 1,
     vertex.color = "pink",
     vertex.frame.color ='pink4',
     edge.color = E(g)$color, 
     vertex.label.color = 'black',
     edge.arrow.size = 0.5,
     edge.curved = 0.3,
     main = "Directed Unweighted Plant Co-occurrence Network")

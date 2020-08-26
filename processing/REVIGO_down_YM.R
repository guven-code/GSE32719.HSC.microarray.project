

# A plotting R script produced by the REVIGO server at http://revigo.irb.hr/
# If you found REVIGO useful in your work, please cite the following reference:
# Supek F et al. "REVIGO summarizes and visualizes long lists of Gene Ontology
# terms" PLoS ONE 2011. doi:10.1371/journal.pone.0021800


# --------------------------------------------------------------------------
# If you don't have the ggplot2 package installed, uncomment the following line:
# install.packages( "ggplot2" );
library( ggplot2 );
# --------------------------------------------------------------------------
# If you don't have the scales package installed, uncomment the following line:
# install.packages( "scales" );
library( scales );


# --------------------------------------------------------------------------
# Here is your data from REVIGO. Scroll down for plot configuration options.

revigo.names <- c("term_ID","description","frequency_%","plot_X","plot_Y","plot_size","log10_p_value","uniqueness","dispensability");
revigo.data <- rbind(c("GO:0000122"," ", 4.380,-5.125, 1.530, 2.881,-4.4473,0.392,0.000),
c("GO:0048511","rhythmic process", 1.812, 0.925,-5.429, 2.498,-3.6882,0.955,0.000),
c("GO:0045444","fat cell differentiation", 1.183,-1.483,-5.761, 2.314,-3.2168,0.952,0.003),
c("GO:0071407","cellular response to organic cyclic compound", 3.104, 4.393, 5.410, 2.732,-3.5376,0.794,0.003),
c("GO:0042221","response to chemical",23.993, 4.637,-0.537, 3.619,-4.2358,0.873,0.146),
c("GO:0009719","response to endogenous stimulus", 9.175, 3.918,-2.202, 3.202,-3.9281,0.884,0.188),
c("GO:0032922","circadian regulation of gene expression", 0.317,-4.771,-1.296, 1.748,-3.2168,0.723,0.233),
c("GO:0007584","response to nutrient", 0.998, 3.759, 6.866, 2.241,-3.2168,0.841,0.269),
c("GO:0010604","positive regulation of macromolecule metabolic process",16.613,-5.567, 4.112, 3.459,-3.9872,0.481,0.294),
c("GO:0051716","cellular response to stimulus",40.358, 4.848, 0.352, 3.845,-3.5331,0.865,0.307),
c("GO:0009892","negative regulation of metabolic process",14.605,-4.177, 2.698, 3.403,-4.6536,0.518,0.381),
c("GO:0010033","response to organic substance",16.515, 4.816, 2.925, 3.457,-4.5735,0.785,0.411),
c("GO:0009889","regulation of biosynthetic process",25.216,-5.271, 2.953, 3.641,-3.8182,0.539,0.417),
c("GO:0009893","positive regulation of metabolic process",17.715,-4.652, 2.298, 3.487,-3.5575,0.550,0.459),
c("GO:2000112","regulation of cellular macromolecule biosynthetic process",23.220,-5.651, 3.346, 3.605,-4.2104,0.443,0.476),
c("GO:0051252","regulation of RNA metabolic process",21.864,-5.011, 3.509, 3.579,-3.2628,0.474,0.485),
c("GO:0010468","regulation of gene expression",25.084,-5.636, 2.450, 3.638,-3.7959,0.513,0.509),
c("GO:0032870","cellular response to hormone stimulus", 3.612, 4.523, 5.041, 2.797,-3.0701,0.766,0.523),
c("GO:0071310","cellular response to organic substance",12.977, 4.762, 3.731, 3.352,-4.7620,0.760,0.582),
c("GO:0051172","negative regulation of nitrogen compound metabolic process", 8.632,-4.158, 3.884, 3.175,-5.1085,0.426,0.589),
c("GO:0070887","cellular response to chemical stimulus",15.701, 4.911, 3.181, 3.435,-4.3363,0.784,0.594),
c("GO:0033993","response to lipid", 5.055, 4.843, 4.529, 2.943,-3.0762,0.789,0.598),
c("GO:0014070","response to organic cyclic compound", 5.338, 4.601, 4.288, 2.967,-3.6517,0.788,0.605),
c("GO:0010629","negative regulation of gene expression", 8.477,-4.880, 1.289, 3.167,-4.4353,0.420,0.668),
c("GO:0009890","negative regulation of biosynthetic process", 8.546,-4.255, 4.374, 3.171,-4.5918,0.410,0.669),
c("GO:0010557","positive regulation of macromolecule biosynthetic process", 9.406,-5.361, 4.537, 3.212,-3.6696,0.461,0.685),
c("GO:0045935","positive regulation of nucleobase-containing compound metabolic process", 9.763,-5.048, 4.851, 3.229,-3.3979,0.485,0.691),
c("GO:0010628","positive regulation of gene expression", 9.833,-5.845, 1.959, 3.232,-4.3809,0.506,0.692));

one.data <- data.frame(revigo.data);
names(one.data) <- revigo.names;
one.data <- one.data [(one.data$plot_X != "null" & one.data$plot_Y != "null"), ];
one.data$plot_X <- as.numeric( as.character(one.data$plot_X) );
one.data$plot_Y <- as.numeric( as.character(one.data$plot_Y) );
one.data$plot_size <- as.numeric( as.character(one.data$plot_size) );
one.data$log10_p_value <- as.numeric( as.character(one.data$log10_p_value) );
one.data$frequency <- as.numeric( as.character(one.data$frequency) );
one.data$uniqueness <- as.numeric( as.character(one.data$uniqueness) );
one.data$dispensability <- as.numeric( as.character(one.data$dispensability) );
#head(one.data);


# --------------------------------------------------------------------------
# Names of the axes, sizes of the numbers and letters, names of the columns,
# etc. can be changed below

p1 <- ggplot( data = one.data );
p1 <- p1 + geom_point( aes( plot_X, plot_Y, colour = log10_p_value, size = plot_size), alpha = I(0.6) ) + scale_size_area();
p1 <- p1 + scale_colour_gradientn( colours = c("blue", "green", "yellow", "red"), limits = c( min(one.data$log10_p_value), 0) );
p1 <- p1 + geom_point( aes(plot_X, plot_Y, size = plot_size), shape = 21, fill = "transparent", colour = I (alpha ("black", 0.6) )) + scale_size_area();
p1 <- p1 + scale_size( range=c(5, 30)) + theme_bw(); # + scale_fill_gradientn(colours = heat_hcl(7), limits = c(-300, 0) );
ex <- one.data [ one.data$dispensability < 0.15, ]; 
p1 <- p1 + geom_text( data = ex, aes(plot_X, plot_Y, label = description), colour = I(alpha("black", 0.85)), size = 3 );
p1 <- p1 + labs (y = "semantic space x", x = "semantic space y");
p1 <- p1 + theme(legend.key = element_blank()) ;
one.x_range = max(one.data$plot_X) - min(one.data$plot_X);
one.y_range = max(one.data$plot_Y) - min(one.data$plot_Y);
p1 <- p1 + xlim(min(one.data$plot_X)-one.x_range/10,max(one.data$plot_X)+one.x_range/10);
p1 <- p1 + ylim(min(one.data$plot_Y)-one.y_range/10,max(one.data$plot_Y)+one.y_range/10);



# --------------------------------------------------------------------------
# Output the plot to screen

p1;

pdf("plots/revigo-down_YM.pdf")
print(p1)     # Plot 1 --> in the first page of PDF
#print(myplot2)     # Plot 2 ---> in the second page of the PDF
dev.off() 

# Uncomment the line below to also save the plot to a file.
# The file type depends on the extension (default=pdf).

# ggsave("C:/Users/path_to_your_file/revigo-plot.pdf");

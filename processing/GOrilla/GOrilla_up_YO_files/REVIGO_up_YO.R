

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
revigo.data <- rbind(c("GO:0007010","cytoskeleton organization", 0.786,-0.166,-5.404, 5.004,-3.2464,0.921,0.000),
c("GO:0014070","response to organic cyclic compound", 0.227, 5.142,-2.853, 4.464,-4.1838,0.873,0.000),
c("GO:0032502","developmental process", 2.812,-4.714,-3.918, 5.557,-3.2487,0.944,0.000),
c("GO:0051241","negative regulation of multicellular organismal process", 0.211,-5.329, 3.031, 4.432,-5.6655,0.414,0.000),
c("GO:0048869","cellular developmental process", 1.896, 3.615, 4.764, 5.386,-3.9706,0.388,0.081),
c("GO:1903530","regulation of secretion by cell", 0.134,-1.555, 3.004, 4.234,-3.5114,0.509,0.207),
c("GO:0048519","negative regulation of biological process", 1.984,-4.464, 5.319, 5.406,-3.4214,0.670,0.263),
c("GO:0060326","cell chemotaxis", 0.060, 4.664,-0.903, 3.888,-3.0448,0.744,0.500),
c("GO:0055002","striated muscle cell development", 0.037, 3.721, 6.090, 3.676,-3.9586,0.341,0.592),
c("GO:0048638","regulation of developmental growth", 0.085, 0.517, 5.702, 4.039,-3.0711,0.336,0.637),
c("GO:0051239","regulation of multicellular organismal process", 0.628,-5.082, 3.674, 4.906,-4.1221,0.525,0.650),
c("GO:0002064","epithelial cell development", 0.065, 4.052, 5.532, 3.923,-3.0186,0.390,0.674));

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
pdf("plots/revigo-up_YO.pdf")
print(p1)     # Plot 1 --> in the first page of PDF
#print(myplot2)     # Plot 2 ---> in the second page of the PDF
dev.off() 
# Uncomment the line below to also save the plot to a file.
# The file type depends on the extension (default=pdf).

# ggsave("C:/Users/path_to_your_file/revigo-plot.pdf");

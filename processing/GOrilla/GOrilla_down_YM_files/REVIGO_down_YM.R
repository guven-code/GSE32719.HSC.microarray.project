

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
revigo.data <- rbind(c("GO:0008152","metabolic process",75.387, 0.017, 5.450, 6.986,-4.3307,0.990,0.000),
c("GO:0009987","cellular process",63.780, 0.943, 2.916, 6.913,-3.4609,0.986,0.000),
c("GO:0071840","cellular component organization or biogenesis", 8.568, 3.112, 4.905, 6.041,-3.2518,0.964,0.000),
c("GO:2000113","negative regulation of cellular macromolecule biosynthetic process", 0.715, 3.065,-5.245, 4.963,-4.8894,0.328,0.000),
c("GO:0001570","vasculogenesis", 0.016, 6.575, 0.771, 3.300,-3.8153,0.856,0.024),
c("GO:0044237","cellular metabolic process",53.061,-4.479,-0.648, 6.833,-6.0376,0.861,0.058),
c("GO:1901564","organonitrogen compound metabolic process",17.886,-5.971,-3.799, 6.361,-3.6126,0.779,0.067),
c("GO:0006807","nitrogen compound metabolic process",38.744,-2.180, 4.220, 6.696,-5.2692,0.911,0.088),
c("GO:0071704","organic substance metabolic process",58.357,-4.288, 2.950, 6.874,-4.8268,0.908,0.119),
c("GO:0044238","primary metabolic process",53.743,-6.016, 1.669, 6.839,-4.5670,0.909,0.120),
c("GO:0044242","cellular lipid catabolic process", 0.180, 5.543, 3.563, 4.363,-3.1024,0.868,0.124),
c("GO:0044260","cellular macromolecule metabolic process",34.276,-2.601,-5.343, 6.643,-5.0362,0.690,0.191),
c("GO:0010810","regulation of cell-substrate adhesion", 0.036, 6.652,-4.275, 3.668,-3.2725,0.649,0.206),
c("GO:0043170","macromolecule metabolic process",39.491,-4.752,-6.698, 6.705,-4.4034,0.834,0.224),
c("GO:0010559","regulation of glycoprotein biosynthetic process", 0.008, 4.502,-6.478, 2.993,-3.0472,0.604,0.277),
c("GO:0043412","macromolecule modification", 9.785,-1.377,-8.199, 6.099,-4.9666,0.770,0.331),
c("GO:0044267","cellular protein metabolic process",14.293,-1.514,-6.149, 6.263,-4.8297,0.685,0.373),
c("GO:0034641","cellular nitrogen compound metabolic process",34.137,-4.422,-3.148, 6.641,-3.3468,0.748,0.415),
c("GO:0030500","regulation of bone mineralization", 0.015, 5.875,-1.508, 3.294,-3.0472,0.645,0.492),
c("GO:0050793","regulation of developmental process", 1.205, 5.554,-2.459, 5.189,-3.1931,0.549,0.526),
c("GO:2000112","regulation of cellular macromolecule biosynthetic process",10.683, 1.950,-5.232, 6.137,-4.5317,0.287,0.555),
c("GO:1903018","regulation of glycoprotein metabolic process", 0.009, 4.673,-6.157, 3.042,-3.0472,0.607,0.597),
c("GO:0019222","regulation of metabolic process",11.942, 5.183,-4.119, 6.185,-3.5287,0.459,0.607),
c("GO:0006464","cellular protein modification process", 7.726,-1.056,-6.854, 5.996,-4.5622,0.702,0.611));

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

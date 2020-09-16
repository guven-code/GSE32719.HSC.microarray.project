

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
revigo.data <- rbind(c("GO:0008152","metabolic process",75.387,-3.843,-0.480, 6.986,-3.6253,0.982,0.000),
c("GO:0009987","cellular process",63.780,-5.041, 2.567, 6.913,-4.1107,0.973,0.000),
c("GO:0016192","vesicle-mediated transport", 1.085,-6.396,-1.408, 5.144,-3.5528,0.928,0.000),
c("GO:0071840","cellular component organization or biogenesis", 8.568, 3.459, 5.710, 6.041,-4.8356,0.933,0.000),
c("GO:2000113","negative regulation of cellular macromolecule biosynthetic process", 0.715, 4.086,-5.776, 4.963,-4.1427,0.402,0.000),
c("GO:0001570","vasculogenesis", 0.016,-0.386,-8.065, 3.300,-4.1124,0.826,0.024),
c("GO:0044237","cellular metabolic process",53.061, 4.636, 3.269, 6.833,-5.0315,0.856,0.058),
c("GO:0006807","nitrogen compound metabolic process",38.744,-1.107, 1.710, 6.696,-3.7328,0.896,0.088),
c("GO:0016043","cellular component organization", 7.239,-3.278,-6.640, 5.968,-4.4306,0.859,0.100),
c("GO:0071704","organic substance metabolic process",58.357, 0.503, 5.213, 6.874,-4.0565,0.899,0.119),
c("GO:0044238","primary metabolic process",53.743,-2.038, 4.896, 6.839,-3.5901,0.898,0.120),
c("GO:0044242","cellular lipid catabolic process", 0.180,-4.888,-4.514, 4.363,-3.3507,0.843,0.124),
c("GO:0044260","cellular macromolecule metabolic process",34.276, 6.022,-1.298, 6.643,-4.1296,0.697,0.191),
c("GO:0043170","macromolecule metabolic process",39.491, 7.396, 0.803, 6.705,-3.3335,0.828,0.224),
c("GO:0043412","macromolecule modification", 9.785, 5.505,-2.276, 6.099,-3.7878,0.755,0.331),
c("GO:0044267","cellular protein metabolic process",14.293, 6.631,-2.721, 6.263,-3.7055,0.671,0.373),
c("GO:0050793","regulation of developmental process", 1.205, 1.430,-7.419, 5.189,-3.3233,0.679,0.526),
c("GO:0060255","regulation of macromolecule metabolic process",11.716, 4.291,-4.637, 6.177,-3.0297,0.520,0.528),
c("GO:0010769","regulation of cell morphogenesis involved in differentiation", 0.055,-0.265,-6.872, 3.852,-3.3215,0.673,0.583),
c("GO:0006464","cellular protein modification process", 7.726, 6.824,-3.547, 5.996,-3.5452,0.675,0.611));

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
pdf("plots/revigo-up_YM.pdf")
print(p1)     # Plot 1 --> in the first page of PDF
#print(myplot2)     # Plot 2 ---> in the second page of the PDF
dev.off() 

# Uncomment the line below to also save the plot to a file.
# The file type depends on the extension (default=pdf).

# ggsave("C:/Users/path_to_your_file/revigo-plot.pdf");

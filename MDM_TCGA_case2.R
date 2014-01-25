## need help with a function?
## just type '?' followed by the function name: e.g. ?plot

########### BASICS
	L <- 5
	i <- 1:10
	s <- matrix(1:10, nrow=2)
	a <- data.frame(first=1:5, second=letters[1:5])
		#Matrices and data.frames must be square
	M <- list(L1 = 1:20, L2 = letters[1:3])
	c <- list(v1 = 20:1, L1 = M, m1 = s)

	L;i;s;a;M;c;
		# semi colons ; separate statements for execution, # designates comment for everything after it
		# colon : acts as a "to" function 1:5 generates numbers 1 to 5
	
	## Create number sequences
		seq(1,3); seq(5,2)		# create a sequence 
		seq(1,20, by=2)			
		rep(5, 3)				# rep(x,y) repeats x y times
		c(rep(c(2,4), 5), seq(1,10, by=2), K =letters[15:19])
			# create a vector by joining elements, can be mixed types and can be named
	
################
## is.installed
##		checks to see if library or package is installed already
################
	is.installed <- function(mypkg) is.element(mypkg, installed.packages()[,1]) 

	is.installed("rgl")

################
## BioConductor as source for obtaining packages
##		biocLite command makes for easy install of dependencies
################

	source("http://bioconductor.org/biocLite.R")
    mypkg <- "rgl"; # mypkg <- c("rgl", "gplots", "VennDiagram", "RColorBrewer")
    biocLite(mypkg)


#######################
## Easy plots in 2D ###

	x <- rnorm(1000)		#generate 1000 random normal values with mean 0 sd 1
	y <- rnorm(1000)
	z <- 0.2*x - 0.3*y + rnorm(1000, sd=0.3)
	fit <- lm(z ~ x + y)	# evaluate the linear regression
	
## Scatterplot
	plot(x,y, pch=19, col="blue", xlab="foo", ylab="bar", main="Easy as that")
		#pch defines the point type, see points() for the variety of options
	abline(fit$coefficients[1:2], lwd=4, lty=2, col="grey")
		# draw a line from with intercept a and slope b; lwd defines the line width, and lty the line type
	points(y,z, pch=21, col="red", cex=0.6, bg="black")
		# add smaller points (0.6% size) to existing plot with red border and black center

## Barplot	
	FooBar_Matrix <- rbind(foo=x, bar=y)[,1:20]
		# join two vectors to by row (cbind forms columns), and only keep first 20 rows
	dimnames(FooBar_Matrix) <- list(c("foo", "bar"),letters[1:20])
		# name the rows and columns
	fb_plot <- barplot(FooBar_Matrix, beside=T, col=c("red", "blue"), main="This is a barplot")
	legend(x=25, y=-1.2, legend=c("foo", "bar"), fill=c("red", "blue"), bty="n", ncol=2)
	
	lines(fb_plot, z[1:40], lwd=4, lty=3, col="green")
		# use the x coordinates from the barplot and add a line to the image
	legend(x=25, y=-1.5, legend=c("random"), col=c("green"), bty="n",lwd=4, lty=3)
	
## Both in side by side panes

	par(mfcol=c(2,1))	# par contains all your plotting information, mfcol defines grid of rows x cols
	plot(x,y, pch=19, col="blue", xlab="foo", ylab="bar", main="Scatterplot")
	fb_plot <- barplot(FooBar_Matrix, beside=T, col=c("red", "blue"), main="Barplot")

################
## Copy image to File
##		save the current image to a specified file
################

	 dev.copy2pdf(file="MDM_2014-1-25_Scatter_Barplot.pdf")
		# saves image on screen to file
		# can also do things like png(), pdf(), postscript(), dev.print(), dev.copy2eps


################
## error.bar
##		add error bars to barplot
################
	error.bar <- function(x, y, upper, lower=upper, length=0.1,...){
		# the "..." allows for other parameters to be passed to other functions
		# in this case, they are pertinent to arrows. 
	
    	if(length(x) != length(y) | length(y) !=length(lower) | length(lower) != length(upper))
	    stop("vectors must be same length")

    	arrows(x,y+upper, x, y-lower, angle=90, code=3, length=length, ...)
	}
	
	require(RColorBrewer)
		# Package useful for creating pallettes of color shades
	
	# make with function calls like
	y <- sapply(1:5,function(i) rnorm(1000, mean=i)); 
	y.means<-colMeans(y); 
	y.sd <- apply(y, 2, sd)
	barx <- barplot(y.means, ylim=c(0,6), col=brewer.pal(5, "Blues"))
	error.bar(barx,y.means, 1.96*y.sd/10)

################
## Zoomable 3D plots
##		library rgl
################
	library(rgl)

	x <- rnorm(100)		#generate 100 random normal values with mean 0 sd 1
	y <- rnorm(100)
	z <- 0.2*x - 0.3*y + rnorm(1000, sd=0.3)
	fit <- lm(z ~ x + y)	# evaluate the linear regression
	plot3d(x,y,z, type="s", col="red", size=1)
		# plot a zoomable and rotating 3D visual of the data, type=spheres

	names(fit)		# get the names of variables returned with the linear fit
	coefs <- coef(fit)	# extract the linear coefficients; fit$coefficients synonymous
	
	a <- coefs["x"];	b <- coefs["y"]
	c <- -1;			d <- coefs["(Intercept)"]
	planes3d(a, b, c, d, alpha=0.5)


################
## ls, rm
##		view and remove variables from environment
################

	ls()
	rm(x, y, z, a,b,c,d)
	rm(list=ls())	#removes pretty much everything from the working environment

##########################
######## TCGA DATA #######
##########################

################
## setwd()
##		set your working directory to a certain folder
################

	setwd("/Volumes/homes/Lisa/MDM/TCGA")

### Define Parameters ###

	Date="2013-12-18"
	TumorNames <- c("LUNG","BRCA","GBM", "LGG") 

################
## read.table, read.delim, read.csv
##		read in structured data
################

	foo <- read.delim(file="2013-12-18/TCGA_BRCA_exp_HiSeqV2-2013-12-18/genomicMatrix",  header=T,sep="\t", row.names=1, stringsAsFactors=F, check.names=F)
		# reads in the genomicMatrix file from the subfolders
		# stores the first row and column as identifiers
		# keeps strings as text instead of factors
		# keeps column names as-is (e.g. doesn't convert dashes to dots)
	# also read.table(), read.csv(); use scan() if you have very large matrices and readLines() if you want single line feeds
	
	### Find out information about an object
	dim(foo)
	class(foo)
	names(foo)[1:10]
	head(foo)[,1:10]	# grabs first couple rows of an object

	foo[1:5,1:5]
	foo[1:5, -c(10:dim(foo)[2])]
		# the minus removes those columns from the display
	# NOTE: nothing in foo is changed since it has not been assigned back to foo
			
	RNF_foo <- foo[grep( "RNF", rownames(foo)),]
		# find all rows that have RNF in the name and save it separately
		
	dim(RNF_foo)
	RNF_foo[, 1:3]
	
	rm(foo, RNF_foo)

################
## apply, sapply, lapply, vapply
##		iteratively traverse data
##		apply for matrix, lapply for lists
##		sapply, vapply wrappers for lapply
##
##    	This is similar to for loops in other functions, which also exist in R
##		But this class of functions is much faster as the data type and size is 
##		defined prior to processing the function, so there aren't multiple allocation steps
##
################

#### Iteratively read in data and store into a list		

	TumorList<- list()
	TumorList<-sapply(TumorNames, function(Name){
		# iterate through TumorNames by assigning them to Name and performing a function

		cat("Getting full ", Name, "File\n")	# cat prints to terminal for progress report, can also be exported to file
   		folder<- paste(Date, "/TCGA_", Name, "_exp_HiSeqV2-", Date, sep="")
   				# paste joins variable values and strings, separated by the sep= token: \t is a tab, \n new line, \r carriage return
		Expr <- read.delim(file=paste(folder, "/genomicMatrix", sep=""), header=T,sep="\t", row.names=1, stringsAsFactors=F, check.names=F)
		Expr		#returns the Expression data read from file
	})
	names(TumorList) <- TumorNames
	# TumorList now contains a list of matrices 
	
	length(TumorList)
	names(TumorList)
	
	class(TumorList[[1]])	# For lists: [[""]]	accesses single element, while [""] is a general subscripting operator
	class(TumorList[1])

	dim(TumorList[["LUNG"]]) # is the same as
	dim(TumorList[[1]])		 #  but not the same as
	dim(TumorList[1])		 # this returns the first top level element, basically the pointer to the object and not the object itself
	length(TumorList[1])

	# Accessing nested data can be done using $ as a delimiter
	TumorList[["LUNG"]][1:5,1:10]
	TumorList$LUNG[1:5,1:3]
	
################
## save/load
##		save: writes data, with variable names, to file
##		load: retrieves variables and data
################

	save(TumorList, file="TCGA_TumorList_2013-12-18.RData")

#	load("TCGA_TumorList_2013-12-18.RData")
	#		TumorList

#	save(TumorNames, TumorList, GeneExpr, Gene_PCA, Tumor_PCA, Tumor_GeneSymbols, TumorVenn, file="MDM_TCGA_AllVariables_2014-1-23.RData")

####################################################
########### PRINCIPAL COMPONENT ANALYSIS ###########
####################################################

# GOAL: Use dimension reduction analysis to identify gene groups (variables) that stratify patients (observations)
# princomp and prcomp use similar but use different methods of eigenvector analysis

### For deeper understanding of PCA, Factor Analysis, or Discriminant Analysis see
## 	www4.ncsu.edu/~lgmcferr

	GeneExpr <- t(TumorList$LUNG)
		# PCA uses observations as rows and columns for variables 
		#	- trying to reduce the number of variables yet still explain variability in observations
		
	Gene_PCA <- prcomp(x=GeneExpr)
		# PCA reports a transformation matrix known as rotation (here) or loadings (princomp) that weight the variables
		#				scores (labeled 'x' here) which are the transformed observations
		# Together the loadings describe which variables are important to the sources of variability (axes of priniciple components)
		# Scores reorient the data according to these similarity measures
		
	screeplot(Gene_PCA, npcs = dim(GeneExpr)[2],type = "lines")		#shows the ordered eigenvector contribution
	screeplot(Gene_PCA, npcs = 20,type = "lines")
		# around 3-6 PCs explain a decent proportion of the variability among patients (20501 variables [genes], 1081 observations [patients])
		# This can be estimated from the "elbow" of the curve or when the standard deviation drops below 1

	summary(Gene_PCA)$importance[,1:20]
	# The first PC (principal Component) explains 19.2% of the variability, the second adds 7%, so 26% of variability can be explained by 2 components
	# Gene_PCA$loadings; Gene_PCA$scores;		# in case you wanted to see all the actual values (keep in mind dimensions)

	min <- min(Gene_PCA$x[,1:2]); max <- max(Gene_PCA$x[,1:2])		# get min and max values to define viewing area
	plot(min:max, min:max, type="n", xlab="PC1", ylab="PC2")
	text(Gene_PCA$x[,1], Gene_PCA$x[,2], labels=rownames(Gene_PCA$x), cex=0.4)
		# Plotting the scores shows how the genes relate to each other according to the first 2 components
	
	min <- min(Gene_PCA$rotation[,1:2]); max <- max(Gene_PCA$rotation[,1:2])
	plot(seq(min, max, by=0.005), seq(min, max, by=0.005), type="n", xlab="PC1", ylab="PC2")
	text(Gene_PCA$rotation[,1], Gene_PCA$rotation[,2], labels=rownames(Gene_PCA$rotation), cex=0.8)
		# Plotting the loadings shows how the samples are separated/related using the first 2 components
	
	####
	##  This plots all component comparisons in different screens for the loadings to observe how the variation is explained among samples

	split.screen( figs = c( 2,3 )); scr <- 1
	for(i in 1:3){	for(j in (i+1):4){
		screen(scr); scr <- scr+1
		min <- min(Gene_PCA$x[,1:2]); max <- max(Gene_PCA$x[,1:2])
		plot(seq(min, max, by=0.005), seq(min, max, by=0.005), type="n", xlab=paste("PC",i,sep=""), ylab=paste("PC",j,sep=""))
		text(Gene_PCA$x[,i], Gene_PCA$x[,j], labels=rownames(Gene_PCA$x), cex=0.5)
	} }
	close.screen(all=TRUE)
	

################
############ Plot in Zoomable 3D!!! ######
################

	text3d(Gene_PCA$x[,1], Gene_PCA$x[,2], Gene_PCA$x[,3], texts=rownames(Gene_PCA$x), cex=0.5)
		# Visualize the relationship among patients

	# NOTE: text3d doesn't form new plotting window, must refresh/close screen to redraw
#	text3d(Gene_PCA$rotation[,1], Gene_PCA$rotation[,2], Gene_PCA$rotation[,3], texts=rownames(Gene_PCA$rotation), cex=0.5)
		# Visualize the contribution of each gene to describing the variability among patients
		# Has Many genes with zero contribution
		
	plot(density(Gene_PCA$rotation[,1:3]), main="Density plot of First 3 loadings")
		# removing genes with rotation coefficient less than 0.02 seems reasonable
	Gene_IDpos <- apply(Gene_PCA$rotation[,1:3], 2, function(loading) { which(abs(loading)>0.02)})
	Gene_IDpos <- unique(unlist(Gene_IDpos))
	
	min <- min(Gene_PCA$rotation[,1:3]); max <- max(Gene_PCA$rotation[,1:3])
	plot3d(seq(min, max, by=0.005), seq(min, max, by=0.005),seq(min, max, by=0.005), type="n", xlab="PC1", ylab="PC2", zlab="PC3", box=F)
	text3d(Gene_PCA$rotation[Gene_IDpos,1], Gene_PCA$rotation[Gene_IDpos,2], Gene_PCA$rotation[Gene_IDpos,3], texts=rownames(Gene_PCA$rotation[Gene_IDpos,]), cex=0.5)


################ Iterate PCA for Tumor Types
# Run PCA for each tumor type and save to a list

Threshold <- 0.02

	Tumor_PCA<- lapply(TumorNames, function(Name){
	
		GeneExpr <- t(TumorList[[Name]])
		Gene_PCA <- prcomp(x=GeneExpr)
		Gene_PCA
	})
	names(Tumor_PCA) <- TumorNames
	
	par(mfcol=c(1,4))
		#another way of separating the plotting window into different frames
		# par() contains all the parameters for defining plots
	Tumor_GeneSymbols<- sapply(TumorNames, function(Name){
		PCA_type <- Tumor_PCA[[Name]]
		plot(density(PCA_type$rotation[,1:3]), main= paste(Name, " Density plot of First 3 loadings", sep=""))
				
		Gene_row <- apply(PCA_type$rotation[,1:3], 2, function(loading) { which(abs(loading)>Threshold)})
		Gene_row <- unique(unlist(Gene_row))
		unique(rownames(Tumor_PCA[[Name]]$rotation[Gene_row,]))
	})
	# line up and plot the 4 loading distributions for comparison in a 1x4 grid
	# and return the row number for which genes contribute to the axes transformation

	Name <- "GBM"; 	#change the name to other TumorNames to visualize each analysis
	text3d(Tumor_PCA[[Name]]$x[,1], Tumor_PCA[[Name]]$x[,2], Tumor_PCA[[Name]]$x[,3], texts=rownames(Tumor_PCA[[Name]]$x), cex=0.5)
		# Visualize the relationship among patients

	text3d(Tumor_PCA[[Name]]$rotation[Tumor_GeneSymbols[[Name]],1], Tumor_PCA[[Name]]$rotation[Tumor_GeneSymbols[[Name]],2], Tumor_PCA[[Name]]$rotation[Tumor_GeneSymbols[[Name]],3], texts=rownames(Tumor_PCA[[Name]]$rotation[Tumor_GeneSymbols[[Name]],]), cex=0.5)
		# Visualize the contribution of each gene to describing the variability among patients

	
################
## Venn Diagram
##		visually compare overlap among sets
################

	if(!is.installed("gplots")) biocLite("gplots")
	require(gplots)
	
	TumorVenn<- venn(Tumor_GeneSymbols)
	# Draw Venn Diagram of Overlapping Gene Symbols from Tumor Types

	universe <- unique(unlist(Tumor_GeneSymbols))
	Universe_In <- sapply(TumorNames, function(Name) universe %in% Tumor_GeneSymbols[[Name]])
	Universe_In[1:10,]
	# Define Universal table for Tumor types and membership of symbols

	universe[Universe_In[,"LUNG"] & Universe_In[,"BRCA"]& Universe_In[,"GBM"]]
	# compare symbol sets and return corresponding gene symbols

	if(!is.installed("VennDiagram")) biocLite("VennDiagram")
	require(VennDiagram)
	
	venn.diagram(Tumor_GeneSymbols, file="TCGA_PCAgeneSym_T02_VennDiagram.tiff", fill=c("blue", "yellow", "green", "red"), cat.fontface=2, cat.fontfamily="serif")
		# make high quality, publication worthy venn diagrams
	


##########################
####### TRICKS ###########
## and additional info

#### Viewing/Extracting Sructured Call info
## S3 classes uses $ for getting values
## S4 classes uses @


#str 
	library(GenomicRanges)
	gr <- GRanges()
	getClass(class(gr))
	str(gr)  # compactly display structure of an object (instance in S3)

	gr@strand
	strand(gr)	# access strand information from gr object

library(Biostrings)
	showMethods(reverseComplement)
	showMethods(class="DNAString", where=search())
  	#	 gives what you can do with that object, only returns functions for what packages are loaded

#xtabs
  # do cross tabulation for multiple parameter queries of datatable
#with
  # Conveniently access columns of a data frame or other element without having to
  # repeat the name of the data frame

# Restructure List
	l <- list(a=c(1,2,3), b=c(2,3,4,5))
	logl <- log(unlist(l, use.names=F))
	relist(logl, l)
	  # uses the skeleton of list l to reorganize the "flesh" of logl

### Parallel Computing
	cl = makeCluster(2L)
	  # for parallel computing
	parLapply(cl, l, log)
	  # run each of the list elements on different cluster
	  # platform independent
	mclapply(l, log)
	  # little more clever, don't have to manage clusters


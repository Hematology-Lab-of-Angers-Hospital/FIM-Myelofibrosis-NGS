# Generation of circos for figure S3 and 1-E
# Author:RENARD Maxime
# Date: 16/12/20
# laboratory of hematology

# ***********************
# Library
library(circlize)
library(dplyr)
library("RColorBrewer")

# ***********************
# Function
# ***********************

# ***********************
# Circos by type of category of gene
# ***********************
Circos_figure <- function(file_matrix,matrix_col,color,category,order_circos,gene_type,seuil) {
	# link treshold
	color[file_matrix < seuil ] = NA
	if (category == "GRINFELD SIGNATURE"){
		svg(filename ="../Figure/Circos_Figure_1-E.svg")
		circos.par(start.degree = 90, clock.wise = TRUE)

		cmd = chordDiagram(file_matrix, annotationTrack = "grid",grid.border="black", preAllocateTracks = 1, annotationTrackHeight = c(0.03, 0.01),small.gap = 1.5,
		,transparency = 0.5,order = order_circos,directional=FALSE,grid.col= matrix_col,col=color)
		title(paste0('\n\n',"4-tier classification",'\n\n'))
		circos.trackPlotRegion(track.index = 1, panel.fun = function(x, y) {
			xlim = get.cell.meta.data("xlim")
			ylim = get.cell.meta.data("ylim")
			sector.name = get.cell.meta.data("sector.index")
		circos.text(mean(xlim), ylim[1] + .1, sector.name,facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5), col = "black",cex=0.4)
		}, bg.border = NA)
		dev.off()
	}
	else{
		
		cmd = chordDiagram(file_matrix, annotationTrack = "grid",grid.border="black", preAllocateTracks = 1, annotationTrackHeight = c(0.03, 0.01),small.gap = 1.5,
		,transparency = 0.5,order = order_circos,directional=FALSE,grid.col= matrix_col,col=color)
		title(paste0('\n\n',category,'\n\n'))
		circos.trackPlotRegion(track.index = 1, panel.fun = function(x, y) {
			xlim = get.cell.meta.data("xlim")
			ylim = get.cell.meta.data("ylim")
			sector.name = get.cell.meta.data("sector.index")
		circos.text(mean(xlim), ylim[1] + .1, sector.name,facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5), col = "black",cex=0.4)
		}, bg.border = NA) 
	}
}

# Main 
main_function_circos <- function(name_matrix,list_category,seuil) {

	# Read file
	file_table = read.table(name_matrix,header=TRUE, sep="\t",row.names = "X")
	# Convert to matrix
	file_matrice = as.matrix(file_table)
	# Replace NA value by 0
	file_matrix = replace(file_matrice, file_matrice == 0.0,NA)
	# Gene show in circos 
	gene_seuil = colnames(file_matrix)

	# Function order circos by category
	order_circos = order_by_category(list_category,gene_seuil)
	# Figure
	svg(filename = paste("../Figure/Circos_Figure_S3.svg",sep=""))
	par(mfrow = c(3, 2))
	circos.par(start.degree = 90, clock.wise = TRUE)
	# Gene of Grinfeld study
	Gene_mauvais_progno = c("EZH2", "IDH1", "IDH2", "ASXL1", "PHF6", "CUX1", "ZRSR2","SRSF2", "U2AF1", "KRAS", "NRAS", "GNAS", "CBL", "RUNX1", "STAG2","BCOR")

	# For each category
	for (category in list_category){

		if (category == "GRINFELD SIGNATURE"){
			
			circos.clear()
			dev.off()
			
			list_gene_cate = Gene_mauvais_progno
		}
		else{
			list_gene_cate = gene_category(category)
		}
		# Intersection with list of gene
		list_gene = intersect(gene_seuil,list_gene_cate)
		# Create matrix with color by category
		matrice = function_Implementation_matrix_color(file_matrix,list_gene) 
		# Recover function result of function_Implementation_matrix_color
		grid_color = matrice[[1]]
		matrix_color = matrice[[2]]

		if(category == "OTHER" | category == "COHESINE"){
			#print(paste0("Category of gene ",category," are no shown in figure S3"))
		}
		else{
			# Call Circos figure
			Circos_figure(file_matrix,grid_color,matrix_color,category,order_circos,list_gene_cate,seuil)
		}
	}
}

# ***************
# Sub-Fonction 
# ***************
#Recover gene by category
gene_category <- function(categorie) {

	# Read file
	file_group_gene = "../Data/Groupe_gene.csv"
	file_cate_gene = read.csv(file_group_gene,header=TRUE, sep=",")
	group = file_cate_gene %>% filter( GROUPE == categorie)%>%select(Gene)
	# List of gene
	list_group = as.vector(group[, "Gene", drop=TRUE])
	
	return(list_group)
}

# Order circos by gene category
order_by_category <- function(liste_category,gene_seuil) {
	# File of category gene
	file_group_gene = "../Data/Groupe_gene.csv"
	file_cate_gene = read.csv(file_group_gene,header=TRUE, sep=",")
	order_circos = c()
	# Course
	for (category in liste_category){
		group = file_cate_gene %>% filter( GROUPE == category)%>%select(Gene)
		list_group = as.vector(group[, "Gene", drop=TRUE])
		for (gene_cate in list_group){
			# Write gene in order for representation
			if (is.element(gene_cate,gene_seuil)){
				order_circos = c(order_circos,gene_cate)
			}
		}
	}
	return(order_circos)
}

# Implementation matrix
function_Implementation_matrix_color<- function(matrix_file,gene_liste) {
	# Color vector
	col = c("#E41A1C","#377EB8","#4DAF4A","#984EA3","#FF7F00","#F781BF","#A65628","#F781BF")
	# Matrix initialisation
	# Grid.color
	grid_color = matrix(nrow=1,ncol=length(colnames(matrix_file)))
	colnames(grid_color) = colnames(matrix_file)

	# Matrix color
	matrix_color = matrix(nrow = length(rownames(matrix_file)),ncol = length(colnames(matrix_file)))
	colnames(matrix_color) = colnames(matrix_file)
	rownames(matrix_color) = rownames(matrix_file)

	# list of gene by category
	list_ONCO = gene_category("TUMOR SUPPRESSOR GENES")
	list_Transcri = gene_category("TRANSCRIPTION FACTORS")
	list_METHYL = gene_category("DNA METHYLATION")
	list_HISTONE = gene_category("HISTONE MODIFICATION")
	list_SIGNAL = gene_category("SIGNALLING")
	list_SPLICE = gene_category("SPLICEOSOME")

	# Color implementation
	for (col_gene in gene_liste){
		if (is.element(col_gene,list_ONCO)){
			color = col[1]
		}
		else if (is.element(col_gene,list_Transcri)){
			color = col[2]
		}
		else if (is.element(col_gene,list_METHYL)){
			color = col[3]
		}
		else if (is.element(col_gene,list_HISTONE)){
			color = col[4]
		}
		else if (is.element(col_gene,list_SIGNAL)){
			color = col[5]
		}
		else if (is.element(col_gene,list_SPLICE)){
			color = col[7]
		}
		else{
			color = NA
		}
		# add value
		grid_color[1,col_gene] = color
		matrix_color[,col_gene] = color
		matrix_color[col_gene,] = color
	}
	# Converse matrix to character 
	converse_color = setNames(grid_color[1,], colnames(grid_color))
	# Return result
	result = list(converse_color,matrix_color)
	return(result)
}

# ********************************************
# Main
# ********************************************

# path directory
path = getwd()
setwd("Figure")

# *******************************************
# Function
# Parameter
type_category =c("TUMOR SUPPRESSOR GENES","TRANSCRIPTION FACTORS","DNA METHYLATION","HISTONE MODIFICATION", "SIGNALLING","SPLICEOSOME","COHESINE","OTHER","GRINFELD SIGNATURE")

# Remove link with no interest
treshold_recurrence = 3
# file 
file = "Matrix_statistic_A_B_5-All.csv"
# Circos
main_function_circos(file,type_category,treshold_recurrence)

setwd("..")

# Statistic on figure
# Author:RENARD Maxime
# Date: 16/12/20
# laboratory of hematology

# ***********************
# Statistic Function
# ***********************
# Statistic test Fischer
file_fischer<- function(file,type) {

	# Initialisation vector
	statistique = c()
	liste_gene = c()
	# recover line repartition
	if(type != "figure_s4_B" && type != "figure_s4_A"){
		All = t(file["All",])
	}
	else if(type == "figure_s4_A" || type == "figure_s4_B"){
		
		All = t(file["Sum",])
	}
	else{
		print("Error of type")
	}
	# Course Data
	for (line in (2:dim(file)[1]-1)){
		data = file[line,]
		Result_gene = t(data)
		tab_contingence = cbind(Result_gene,All)
		# CImplementation of matrix
		if (type == "figure_1_B"){
			tab_contingence["Sum_MFP","All"] = tab_contingence["Sum_MFP","All"] - tab_contingence["Sum_MFP",1]
			tab_contingence["Sum_MFS","All"] = tab_contingence["Sum_MFS","All"] - tab_contingence["Sum_MFS",1]
		}
		else if (type == "figure_1_C"){
			tab_contingence["Sum_JAK2","All"] = tab_contingence["Sum_JAK2","All"] - tab_contingence["Sum_JAK2",1]
			tab_contingence["Sum_CALR","All"] = tab_contingence["Sum_CALR","All"] - tab_contingence["Sum_CALR",1]
		}
		else if (type == "figure_1_A"){
			tab_contingence["Sum_ABC","All"] = tab_contingence["Sum_ABC","All"] - tab_contingence["Sum_ABC",1]
			tab_contingence["Sum_AB","All"] = tab_contingence["Sum_AB","All"] - tab_contingence["Sum_AB",1]
		}
		else if (type == "figure_s4_A"){
			tab_contingence["F","Sum"] = tab_contingence["F","Sum"] - tab_contingence["F",1]
			tab_contingence["M","Sum"] = tab_contingence["M","Sum"] - tab_contingence["M",1]
			
		}
		else if (type == "figure_s4_B"){
			for ( i in 1:dim(tab_contingence)[1]){
				tab_contingence[i,"Sum"] = tab_contingence[i,"Sum"] - tab_contingence[i,1]
			}
			print("Plusieurs Catégories")
		}
		else{
			print("Error")
		}
		tab = t(tab_contingence)
		name_gene = row.names(file[line,])
		# statistic
		stat = fisher.test(tab,alternative = "two.sided")
		statistique = c(statistique,stat$p.value) 
		liste_gene = c(liste_gene,name_gene)
	}
	statistic = setNames(statistique, liste_gene)
	hochberg_ajustement = p.adjust(statistic,method = "hochberg",n = length(statistic))
	matrix_statistic_result = rbind(statistic,hochberg_ajustement)

	write.table(matrix_statistic_result, file = paste("../Figure/Statistic_Fischer_",type,".csv",sep = ""), sep = ",",, row.names=TRUE, col.names=TRUE)
}

# ***********************
# Main
# ***********************

path = getwd()
setwd("File_statistic")
# File 
Files = c("Data_statistic_Figure_1_A.csv","Data_statistic_Figure_1_B.csv","Data_statistic_Figure_1_C.csv")
# Type of file
type = c("figure_1_A","figure_1_B","figure_1_C")

i = 1
for (name_file in Files){
	if (type[i]== "figure_s4_A" || type[i] == "figure_s4_B"){
		file_table = read.table(name_file,header=TRUE, sep=",",row.names = "Classif_count")
		file_table[is.na(file_table)] <- 0
	}
	else{
		file_table = read.table(name_file,header=TRUE, sep=",",row.names = "Gene.refGene")
		
	}
	file_fischer(file_table,type[i])
	i= i + 1
}
setwd("../")



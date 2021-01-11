# -*- coding: utf-8 -*-
#!/usr/bin/python3
# Author Maxime RENARD, Bioinformatician
# Contact: maxime.renard
# Date: 28/08/20
# Laboratory of d'Hematology, CHU
# Generation of figure of article

# ***************************************************************
# Import Library
# ***************************************************************

import os
# Treat Data
import numpy as np
import pandas as pd
# Module for figure
import matplotlib.pyplot as plt
import seaborn as sns

# ***************************************************************
# Function
# ***************************************************************

# Fonction to read data 
def read_files(files,car):
    """
        IN : Name of files
        OUT : Read of files
    """
    data = pd.read_csv(files, sep = car,na_values = ['.'])
    return(data.reset_index())


# 2 - Main Function to realise all figure
# Main Fonction to create figure 
def figure_genomic(data,size_label,quality_dpi,size_figure,seuil_gene,name_repertory):
    """
        IN: Data and parameters
        OUT: figure on csv format
    """

    # Figure 1:
    # Figure 1-A
    # Proportion of patients with at least one mutation per gene according to the classification of mutations
    figure_1_A(data,size_label,seuil_gene,quality_dpi,size_figure,name_repertory)
    # Figure 1-B
    # Proportion of patients with at least one mutation per gene according to  the type of myelofibrosis 
    figure_1_B(data,size_label,seuil_gene,quality_dpi,size_figure,name_repertory)
    # Figure 1-C
    # Proportion of patients with at least one mutation per gene according to the driver mutation (JAK2 or CALR)
    figure_1_C(data,size_label,seuil_gene,quality_dpi,size_figure,name_repertory)

    # Figure 1-D
    # Distribution of allelic burden of additionnal mutations 
    # Treshold of number variant per gene is 25 and  gene IDH1/2
    seuil_rep = 25
    figure_1_D(data,size_label,seuil_rep,quality_dpi,size_figure,name_repertory)
    
    # Figure 1-E
    # association between mutations on the genes of the 4-tier classification were shown in a circos plot
    #Figure S3
    # Association of pairwise mutations represented by circus plot for each functional category of mutated genes. 
    treshold = 5
    # Creation of matrix use by R code to generate new file
    matrix_stat_circos(data,int(treshold),name_repertory)
    os.system("Rscript Code/Generation_circos.R")
    # **********************
    # Figure supplementary :
    # Figure S2 
    # Allele Burden of driver mutation according to the type of myelofibrosis
    # *********
    figure_s2(data,size_label,quality_dpi,size_figure,name_repertory)

    # Figure S4A 
    # Distribution of the number of additionnal mutations according to gender
    figure_s4A(data,size_label,quality_dpi,size_figure,name_repertory)
    # Figure S4B 
    # Distribution of the number of additionnal mutations according to age at diagnosis myelofibrosis
    figure_s4B(data,quality_dpi,size_figure,name_repertory)

    # Figure S7
    # Figure with ASXL1 
    # *********
    # Preparation of data
    data_group = Prepare_Data_ASXL1(data)
    # Distribution of allele burden ASXL1 according to 4-tiers classification
    figure_s7(data_group,size_label,seuil_rep,quality_dpi,size_figure,name_repertory)

# # 3- Function Figure     
def figure_1_A(data,size_label,seuil,quality_dpi,size_figure,name_repertory):
    """
       IN : Data
       OUT : Figure_1A Proportion of patients with at least one mutation per gene according to the classification of mutations
    """
    # Preparation Data
    index_gene = ['NOTCH1', 'RUNX1', 'CBL', 'MPL', 'KMT2D', 'CUX1', 'KMT2C', 'SH2B3',
       'TP53', 'ATM', 'DNMT3A', 'NFE2', 'SF3B1', 'EZH2', 'U2AF1', 'SRSF2',
       'TET2', 'ASXL1']
    
    # All Class A B C 
    classe=["A","B","C"]
    #  Information Total of 433 patient with at least one additionnal mutation
    data_ABC = data[data['Classification_mutation'].isin(classe)]
    df_plot = data_ABC.groupby(['PATIENT','Gene.refGene']).size().reset_index().pivot(columns='PATIENT', index='Gene.refGene')
    # Replace NA to 0
    df_plot.fillna(0,inplace= True)
    # If patient have more 1 mutation in same gene we count only once number of variant in this gene
    df_plot[df_plot != 0] = 1
    df_plot["Sum_ABC"]= df_plot.apply(np.sum, axis=1)
    df_plot_ABC = df_plot[(df_plot["Sum_ABC"] >= seuil) ].sort_values(by='Sum_ABC', ascending=True)

     
    # Class A and B only
    classe=["A","B"]
    data_AB = data[data['Classification_mutation'].isin(classe)]
    df_plot_AB = data_AB.groupby(['PATIENT','Gene.refGene']).size().reset_index().pivot(columns='PATIENT', index='Gene.refGene')
    # Replace NA to 0
    df_plot_AB.fillna(0,inplace= True)
    # If patient have more 1 mutation in same gene we count only once number of variant in this gene
    df_plot_AB[df_plot_AB != 0] = 1
    index_ABC = list(df_plot_ABC.index)
    
    df_plot_AB["Sum_AB"]= df_plot_AB.apply(np.sum, axis=1)
    df_plot_sum_AB = df_plot_AB[(df_plot_AB.index.isin(index_ABC))].sort_values(by='Sum_AB', ascending=True)
    
    # Merge file ABC clasification and AB
    combine_file = pd.merge(left=df_plot_ABC["Sum_ABC"],right=df_plot_sum_AB["Sum_AB"], left_index =True, right_index=True,how="inner")
    # Convert count to proportion
    combine_file_stat = combine_file
    combine_file_stat.loc["All","Sum_ABC"] = 479
    combine_file_stat.loc["All","Sum_AB"] = 479
   
    percentage_class = combine_file.apply(lambda x: (100. * x / 479),axis=1).sort_values(by='Sum_ABC', ascending=True)
    
     # Reindexation
    percentage_class_reindex = percentage_class.reindex(index_gene)
    #Figure 
    # Parameter
    sns.set_style("ticks")
    sns.set_context('paper')
    # Main
    percentage_class_reindex[["Sum_AB","Sum_ABC"]].plot.barh(color=["#e74c3c","#3498db"],width=0.8,xlim=(0,40),xticks=[2,5,10,20,30,40],grid = [2,5,10,20,30,40],align="center",figsize=(7,10),logx=False)
    # Labels
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    plt.ylabel('')
    plt.xlabel("Proportion of patients with mutation (%)", fontsize=18)
    plt.legend(loc='lower right',prop={'size':12},title="",title_fontsize = 16, borderpad=0.5, labelspacing=0.5,labels= ["Pathogenic and likely pathogenic mutations",'All mutations'])

    sns.despine(top = True, bottom = False, left = False, right = True)
    plt.savefig(name_repertory + 'Figure_1_A.svg' ,dpi = quality_dpi,bbox_inches='tight')

    plt.close()
    percentage_class.to_csv(name_repertory + "Data_Describe_Figure_1_A.csv",sep = ",")

def figure_1_B(data,size_label,seuil,quality_dpi,size_figure,name_repertory):
    """
        IN : Data
        OUT : Figure 1_B Proportion of patients with at least one mutation per gene according to  the type of myelofibrosis
    """
    # Preparation Data
    
    # Primary Myelofibrosis : 305 patients
    data_MFP = data[data['TYPE_MYELOFIBROSIS'] == "Primary myelofibrosis"]
    total_MFP = data_MFP["PATIENT"].nunique()
    index_gene = ['NOTCH1', 'RUNX1', 'CBL', 'MPL', 'KMT2D', 'CUX1', 'KMT2C', 'SH2B3',
       'TP53', 'ATM', 'DNMT3A', 'NFE2', 'SF3B1', 'EZH2', 'U2AF1', 'SRSF2',
       'TET2', 'ASXL1']
    df_plot = data_MFP.groupby(['PATIENT','Gene.refGene']).size().reset_index().pivot(columns='PATIENT', index='Gene.refGene')
    # Replace NA to 0
    df_plot.fillna(0,inplace= True)
    # If patient have more 1 mutation in same gene we count only once number of variant in this gene
    df_plot[df_plot != 0] = 1
    df_plot["Sum_MFP"]= df_plot.apply(np.sum, axis=1)
    df_plot_MFP = df_plot[(df_plot.index.isin(index_gene))]
    index_MFP = list(df_plot_MFP.index)
     
    # Secondary Mielofibrosis : 174 patients
    data_MFS = data[data['TYPE_MYELOFIBROSIS'] == "Secondary myelofibrosis"]
    total_MFS = data_MFS["PATIENT"].nunique()
    df_plot_MFS = data_MFS.groupby(['PATIENT','Gene.refGene']).size().reset_index().pivot(columns='PATIENT', index='Gene.refGene')
     # Replace NA to 0
    df_plot_MFS.fillna(0,inplace= True)
    # If patient have more 1 mutation in same gene we count only once number of variant in this gene
    df_plot_MFS[df_plot_MFS != 0] = 1
    df_plot_MFS["Sum_MFS"]= df_plot_MFS.apply(np.sum, axis=1)
    df_plot_sum_MFS = df_plot_MFS[df_plot_MFS.index.isin(index_gene)]

    # Merge Data of MFP ans MFS
    combine_file = pd.merge(left=df_plot_MFP["Sum_MFP"],right=df_plot_sum_MFS["Sum_MFS"], left_index =True, right_index=True,how="inner")
    # Calcul
    combine_file_stat = combine_file
    combine_file_stat.loc["All","Sum_MFP"] = total_MFP
    combine_file_stat.loc["All","Sum_MFS"] = total_MFS
  
    # Percentage
    combine_file["Sum_MFP"] = (combine_file["Sum_MFP"]/total_MFP)*100
    combine_file["Sum_MFS"] = (combine_file["Sum_MFS"]/total_MFS)*100
 
    percentage_class = combine_file.reindex(index_gene)
    colors = [(0.9882352941176471, 0.5529411764705883, 0.3843137254901961),(0.4, 0.7607843137254902, 0.6470588235294118)]
    # Figure 
    # Parameters
    sns.set_style("ticks")
    sns.set_context('paper')
    #Main
    percentage_class[["Sum_MFS","Sum_MFP"]].plot.barh(color=colors,width=0.8,xlim=(0,42),xticks=[2,5,10,20,30,40],grid = [2,5,10,20,30,40],align="center",figsize=(7,10),logx=False)
    #Labels
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    plt.ylabel('')
    plt.xlabel("Proportion of patients with mutation (%)", fontsize=18)
    label = ['Secondary myelofibrosis',"Primary myelofibrosis"] 
    plt.legend(loc='lower right',prop={'size':12},title="",title_fontsize = 16, borderpad=0.5, labelspacing=0.5,labels= label)
    sns.despine(top = True, bottom = False, left = False, right = True)
    # Save figure
    plt.savefig(name_repertory + "Figure_1_B.svg" ,dpi = quality_dpi,bbox_inches='tight')
    # Close Figure
    plt.close()
    # Save statistic
    percentage_class.to_csv(name_repertory + "Data_Describe_Figure_1_B.csv",sep = ",")

def figure_1_C(data,size_label,seuil,quality_dpi,size_figure,name_repertory):
    """
        IN: Data
        OUT: Figure 1_C Proportion of patients with at least one mutation per gene according to the driver mutation (JAK2 or CALR)
    """
    # Preparation Data

    # JAK2 V617F: 312 patient
    data_JAK2 = data[data['DRIVER_USE_Figure'] == "JAK2"]
    total_JAK2 = data_JAK2["PATIENT"].nunique()
    index_gene = ['NOTCH1', 'RUNX1', 'CBL', 'MPL', 'KMT2D', 'CUX1', 'KMT2C', 'SH2B3',
       'TP53', 'ATM', 'DNMT3A', 'NFE2', 'SF3B1', 'EZH2', 'U2AF1', 'SRSF2',
       'TET2', 'ASXL1']
    df_plot = data_JAK2.groupby(['PATIENT','Gene.refGene']).size().reset_index().pivot(columns='PATIENT', index='Gene.refGene')
    # Replace NA to 0
    df_plot.fillna(0,inplace= True)
    # If patient have more 1 mutation in same gene we count only once number of variant in this gene
    df_plot[df_plot != 0] = 1
    df_plot["Sum_JAK2"]= df_plot.apply(np.sum, axis=1)
    df_plot_JAK2 = df_plot[(df_plot.index.isin(index_gene))]

     # CALR: 113 patients
    data_CALR = data[data['DRIVER_USE_Figure'] == "CALR"]
    total_CALR = data_CALR["PATIENT"].nunique()
    df_plot_CALR = data_CALR.groupby(['PATIENT','Gene.refGene']).size().reset_index().pivot(columns='PATIENT', index='Gene.refGene')
     # Replace NA to 0
    df_plot_CALR.fillna(0,inplace= True)
    # If patient have more 1 mutation in same gene we count only once number of variant in this gene
    df_plot_CALR[df_plot_CALR != 0] = 1
    
    df_plot_CALR["Sum_CALR"]= df_plot_CALR.apply(np.sum, axis=1)
    df_plot_sum_CALR = df_plot_CALR[df_plot_CALR.index.isin(index_gene)]

    # Combine column to JAK2 and CALR
    combine_file = pd.merge(left=df_plot_JAK2["Sum_JAK2"],right=df_plot_sum_CALR["Sum_CALR"], left_index =True, right_index=True,how="inner")
    combine_file_stat = combine_file
    combine_file_stat.loc["All","Sum_JAK2"] = total_JAK2
    combine_file_stat.loc["All","Sum_CALR"] = total_CALR
    

    # Convert count to proportion
    combine_file["Sum_JAK2"] = (combine_file["Sum_JAK2"]/total_JAK2)*100
    combine_file["Sum_CALR"] = (combine_file["Sum_CALR"]/total_CALR)*100
    # Reindexation
    percentage_class = combine_file.reindex(index_gene)
    
    #Figure 
    # Parameter
    sns.set_style("ticks")
    sns.set_context('paper')
    colors = ["#FFD700","#9b59b6"]
    # Main
    percentage_class[["Sum_CALR","Sum_JAK2"]].plot.barh(color=colors,width=0.8,xlim=(0,42),xticks=[2,5,10,20,30,40],grid = [2,5,10,20,30,40],align="center",figsize=(7,10),logx=False)
    # Labels
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    plt.ylabel('')
    plt.xlabel("Proportion of patients with mutation (%)", fontsize=18)
    label = ["CALR",'JAK2 V617F'] 
    plt.legend(loc='lower right',prop={'size':12},title="",title_fontsize = 16, borderpad=0.5, labelspacing=0.5,labels= label)
    sns.despine(top = True, bottom = False, left = False, right = True)
    # Save fig
    plt.savefig(name_repertory + 'Figure_1_C.svg' ,dpi = quality_dpi,bbox_inches='tight')
    # Close Figure
    plt.close()
    # Write Statistic
    percentage_class.to_csv(name_repertory + "Data_Describe_Figure_1_C.csv",sep = ",")

# Repartition de la VAF des genes avec seuil de 25 en ajoutant IDH1 et IDH2 
def figure_1_D(data,size_label,seuil,quality_dpi,size_figure,name_repertory):
    """
        IN: Data
        OUT: Figure 1_D Distribution of allelic burden of additionnal mutations
    """
    df_plot = data.groupby(['Gene.refGene']).size().sort_values(ascending=False)
    
    # Filter accordind to treshold
    df_plot_filter = df_plot[(df_plot > seuil)].sort_values(ascending=False)

    # Gene which pass treshold
    gene_filtre = list(df_plot_filter.index)
    # add gene IDH1/2
    gene_filtre.append("IDH1")
    gene_filtre.append("IDH2")
    # Data use for figure 
    data_filter = data[(data["Gene.refGene"].isin(gene_filtre))]
    
    # Parameter
    plt.figure(figsize=size_figure)
    sns.set(style="whitegrid")
    sns.set_context('paper')
    
    # Call figure
    ax = sns.violinplot(x="Gene.refGene", y="Mean_VAF",order=gene_filtre, data=data_filter,palette= "Paired",cut=0,width=0.8,linewidth=1.5)
      
    # Subtitle Effectif
    nobs = data_filter['Gene.refGene'].value_counts()
    nobs = [str(x) for x in nobs.tolist()]
    nobs = ["n= " + i for i in nobs]
 
    # Add it to the plot
    pos = range(len(nobs))
    for tick,label in zip(pos,ax.get_xticklabels()):
        ax.text(pos[tick], -2, nobs[tick], horizontalalignment='center', size='12', color='black', weight='semibold')
    
    # Labels
    plt.xlabel(' ', fontsize=size_label)
    plt.ylabel("Distribution of VAF (%)", fontsize=size_label)
    plt.xticks(fontsize=14)
    plt.yticks(fontsize=size_label)
    sns.despine(top = True, bottom = False, left = False, right = True)
    # Save Figure 
    plt.savefig(name_repertory + 'Figure_1_D.svg',dpi=quality_dpi,bbox_inches='tight')
    # Close Figure
    plt.close()
    # Write statistics figure
    describe_figure =  data.groupby('Gene.refGene')['Mean_VAF'].agg(["count","median","min","max"]).sort_values(by="count",ascending=False)
    describe_figure.to_csv(name_repertory + "Data_Describe_Fig_1_D.csv",sep = ",")

def figure_s2(data,size_label,quality_dpi,size_figure,name_repertory):
    """
        IN : Data
        OUT : Figure_S2 Allele burden of driver mutations according to the type of myelofibrosis
    """
    # Extraction of data : one line per patient
    data_patient = data.groupby(['PATIENT'])[["PATIENT","DRIVER_USE_Figure","TYPE_DRIVER_USE_Figure","VAF_DRIVER_USE_Figure","TYPE_MYELOFIBROSIS"]].first()
    # Parameter figure
    plt.figure(figsize=size_figure)
    sns.set(style="whitegrid")
    sns.set_context('paper')
    
    ax = sns.violinplot(y = "VAF_DRIVER_USE_Figure", x = "TYPE_DRIVER_USE_Figure" , data=data_patient, order = ["V617F","TYPE1","TYPE2","Others","W515"],hue ="TYPE_MYELOFIBROSIS",split=True,palette="Set2",scale="width",cut=0,scale_hue="True", width=0.7,linewidth=4)
 
   # Subtitle effectif
    nobs_size = data_patient.groupby(['TYPE_DRIVER_USE_Figure','TYPE_MYELOFIBROSIS'])['VAF_DRIVER_USE_Figure'].count()
    
    # reindexation in order graphics                     
    nobs_reindex = nobs_size.reindex([("V617F","Primary myelofibrosis"),("V617F","Secondary myelofibrosis"),
        ("TYPE1","Primary myelofibrosis"),("TYPE1","Secondary myelofibrosis"),
        ("TYPE2","Primary myelofibrosis"),("TYPE2","Secondary myelofibrosis"),
        ("Others","Primary myelofibrosis"),("Others","Secondary myelofibrosis"),
        ("W515","Primary myelofibrosis"),("W515","Secondary myelofibrosis")])
    nobs = [str(x) for x in nobs_reindex.tolist()]
    nobs_2 = []
    # Notation of value by type of myelofibrosis
    i = 0
    while i <= len(nobs)-2:
        nobs_2.append("n= " + str(nobs[i]) + " n= " + str(nobs[i+1]))
        i+=2
    # Add it to the plot
    pos = range(len(nobs_2))
   
    for tick,label in zip(pos,ax.get_xticklabels()):

        ax.text(pos[tick], -2.5, nobs_2[tick], horizontalalignment='center', size='16', color='black', weight='semibold')
       
    # Labels
    plt.xlabel('', fontsize=size_label)
    plt.ylabel("Distribution of VAF driver (%)", fontsize=size_label)
    plt.legend(fontsize=size_label,title = "", title_fontsize=size_label,bbox_to_anchor=(0, 0, 0.5, 1))
    axes = plt.gca()
    axes.xaxis.set_ticklabels(['\nJAK2 V617F', '\nCALR Type 1', '\nCALR Type 2', '\nCALR Others', '\nMPL W515'], color = 'black', fontsize = 20, verticalalignment = 'center')
    plt.xticks(fontsize=size_label)
    plt.yticks(fontsize=size_label)
    sns.despine(top = True, bottom = False, left = False, right = True)
    # Save figure
    plt.savefig(name_repertory + 'Fig_supplementary_S2.svg',dpi = quality_dpi,bbox_inches='tight')
    # Close figure
    plt.close()
    # Write file for statistic
    describe_figure = data_patient.groupby(["TYPE_DRIVER_USE_Figure","TYPE_MYELOFIBROSIS"])["VAF_DRIVER_USE_Figure"].agg(["count","median","min","max"])
    describe_figure["sum"]= sum(describe_figure["count"])
    describe_figure.to_csv(name_repertory + "Data_Describe_Fig_supplementary_S2.csv",sep = ",")


def figure_s4A(data,size_label,quality_dpi,size_figure,name_repertory):
    """
        IN : Data and figure parameters
        OUT : Figure of supplementary 4-A Distribution of the number of additionnal mutations according to gender
    """
    # Data preparation
    dataAB = data [data["Classification_mutation"] != "C"]
    # keep all info patient : patient without additionnal mutation classification A or B
    data_info = data.groupby("PATIENT")[["FILTER","Age_category","SEX","TYPE_MYELOFIBROSIS"]].first().reset_index()
    colonne = ["PATIENT","FILTER","Age_category","SEX","TYPE_MYELOFIBROSIS"]
    # Merge data 
    dataAB_all = data_info.merge(dataAB,how="outer",on = colonne)
    # ***************
    # Calcul result to show 
    count_patient_sexe = dataAB_all.groupby(["PATIENT","SEX"])["ID"].count().reset_index()
    # Cut information by range
    ranges = [0,1,3,6,15]
    labels = ["0","1-2","3-5",">5"] 
    count_patient_sexe["Classif_count"] =pd.cut(count_patient_sexe["ID"], bins = ranges,right=False,labels=labels)
    
    # Repartition of classification  
    class_count = count_patient_sexe.groupby(['SEX', 'Classif_count'])[['PATIENT']].count().reset_index().pivot(columns='Classif_count', index='SEX',values= 'PATIENT')
    class_count.columns = class_count.columns.add_categories(['Sum'])
    class_count['Sum'] = class_count.sum(axis = 1)
    
    # Subtitle of population repartition
    nobs = [str(x) for x in class_count["Sum"].tolist()]
    nobs = ["n= " + i for i in nobs]
 
    # Save data to use for Statistique
    #class_count.T.to_csv(name_repertory_stat + "Data_Statistic_Fig_supplementary4-A.csv",sep = ",")
    
    # Deletion of Sum column : Columns use only for calcul
    del class_count["Sum"]
    
    # Convert count result in percentage
    class_count_reindex = class_count
    percentage_class_count = class_count_reindex.apply(lambda x: (100. * x / x.sum()),axis=1)
    percentage_class_count = percentage_class_count.fillna(0)
    
    # ***************
    # Figure
    # Parameters
    plt.figure(figsize=size_figure)
    sns.set_style("ticks")
    sns.set_context('paper')
    percentage_class_count[['0', '1-2', '3-5', '>5']].plot.bar(stacked=True,width=1,color=sns.color_palette("Greys",4))
    # Labels
    plt.ylabel("Proportion of patient with additionnal mutation (%)",fontsize=14 ,rotation="vertical")
    plt.xlabel("")
    axes = plt.gca()
    axes.xaxis.set_ticklabels(["Women\n"+ str(nobs[0]),"Men\n" + str(nobs[1])], color = 'black', fontsize = 14,horizontalalignment = 'center')
    plt.xticks(rotation = 'horizontal')
    plt.legend(bbox_to_anchor=(1.05, 1), title= "Number Mutations",title_fontsize = 18,loc='upper left',prop={'size':14}, borderpad=0.25, labelspacing=0.25)
    sns.despine(top = True, bottom = False, left = False, right = True)
    
    # save figure in svg format
    plt.savefig(name_repertory + 'Fig_supplementary_4-A.svg',dpi = quality_dpi,bbox_inches='tight')
    
    # Close graphics
    plt.close()
    percentage_class_count.to_csv(name_repertory + "Data_Describe_Fig_supplementary4-A.csv",sep = ",")


def figure_s4B(data,quality_dpi,size_figure,name_repertory):
    """
        IN: Data
        OUT : Fig s4B Distribution of the number of additionnal mutations according to age at diagnosis myelofibrosis
    """
    # Prepare Data
    dataC = data[data["Classification_mutation"]!= "C"]
    # keep all patient 
    data_info = data.groupby("PATIENT")[["FILTER","Age_category","SEX","TYPE_MYELOFIBROSIS"]].first().reset_index()
    colonne = ["PATIENT","FILTER","Age_category","SEX","TYPE_MYELOFIBROSIS"]
    dataAB = data_info.merge(dataC,how="outer",on = colonne)

    # ***************************
    # Calcul Data
    count_patient_age = dataAB.groupby(["PATIENT","Age_category"])["ID"].count().reset_index()
    # Cut Data 
    ranges = [0,1,3,6,15]
    labels = ["0","1-2","3-5",">5"] 
    # Creation des coupures
    count_patient_age["Classif_count"] =pd.cut(count_patient_age["ID"], bins = ranges,right=False,labels=labels)
    
    # Repartition de la classification  
    class_count = count_patient_age.groupby(['Age_category', 'Classif_count'])[['PATIENT']].count().reset_index().pivot(columns='Classif_count', index='Age_category',values= 'PATIENT')
    class_count.columns = class_count.columns.add_categories(['Sum'])
    class_count['Sum'] = class_count.sum(axis = 1)
    # Sub title of reparition in effectif
    nobs = [str(x) for x in class_count["Sum"].tolist()]
    nobs = ["n= " + i for i in nobs]
    
    # Save_data
    #class_count.T.to_csv(name_repertory_stat + "Data_Statistic_Fig_supplementary4-B.csv",sep = ",")
    # Suppresion column Sum no use in figure
    del class_count["Sum"]
    
    class_count_reindex = class_count.reindex(["[80,90[","[70,80[","[60,70[","[50,60[","[40,50[","< 40"])
    # Convert count in Proportion
    percentage_class_count = class_count_reindex.apply(lambda x: (100. * x / x.sum()),axis=1)
    percentage_class_count = percentage_class_count.fillna(0)
    # *****************
    # Figure
    # Parameters
    plt.figure(figsize=size_figure)
    sns.set_style("ticks")
    sns.set_context('paper')
    percentage_class_count[['0', '1-2', '3-5', '>5']].plot.barh(stacked=True,width=1,color=sns.color_palette("Greys",4))
    # Legend
    plt.xlabel("Proportion of patient with additionnal mutation (%)",fontsize=12)
    plt.ylabel("")
    axes = plt.gca()
    axes.yaxis.set_ticklabels(["[80,90[\n" + str(nobs[5]),"[70,80[\n"+ str(nobs[4]),"[60,70[\n" + str(nobs[3]),"[50,60[\n" + str(nobs[2]),"[40,50[\n" + str(nobs[1]),"< 40\n" + str(nobs[0])], color = 'black', fontsize = 10, verticalalignment = 'center')
    
    plt.legend(bbox_to_anchor=(1.05, 1), title= "Number Mutations",title_fontsize = 18,loc='upper left',prop={'size':14}, borderpad=0.25, labelspacing=0.25)
    sns.despine(top = True, bottom = False, left = False, right = True)
    plt.savefig(name_repertory + "Fig_supplementary_4-B.svg",dpi=quality_dpi,bbox_inches='tight')
    
    # Close figure
    plt.close()
    percentage_class_count.to_csv(name_repertory + "Data_Describe_Fig_supplementary_S4B.csv",sep = ",")
   

def Prepare_Data_ASXL1(data):
    """
        IN: All Data
        OUT : Creation of data for figure_7s
    """ 
    # All patient with no ariant are no representatated in this figure
    # a patient must have at least ASXL1 mutation
    # Keep only pathogenic variant
    data_AB = data[data['Classification_mutation'] != "C"]
    data_AB.rename(columns = {'Gene.refGene':'Gene_refGene'}, inplace = True)    
    # Group data by patient
    grouped = data_AB.groupby(['PATIENT'])
    # Variable to create group
    gene_interet = "ASXL1"
    groupe_figure_ASXL1 = ["TP53","SRSF2|RUNX1|IDH1/2|NRAS|KRAS","U2AF1","EZH2|CBL","Others"]
    # Read Data by group patient
    for name,group in grouped:

        # Repartition of gene by patient
        class_group = group['Gene_refGene'].value_counts()
       
        # If ASXL1 gene
        if gene_interet in class_group.index:
            # ASXL1  == unique groupe Class in Others group
            if len(class_group) == 1:
                groupe = groupe_figure_ASXL1[4]
                
            # Many gene
            elif len(class_group) >= 2:
                # Order in group by priority
                # group TP53 
                if "TP53" in class_group.index:
                    groupe = groupe_figure_ASXL1[0]
                     
                # Second Group
                elif "SRSF2" in class_group.index or "RUNX1" in class_group.index or "IDH1" in class_group.index or "IDH2" in class_group.index or "NRAS" in class_group.index or "KRAS" in class_group.index:
                    groupe = groupe_figure_ASXL1[1]
                    
                # Third group U2AF1
                elif "U2AF1" in class_group.index:
                    groupe = groupe_figure_ASXL1[2]
                    
                # Fourth group EZH2|CBL
                elif "EZH2" in class_group.index or "CBL" in class_group.index:
                    groupe = groupe_figure_ASXL1[3]
                    
                # Other mutations gene : groupe Others   
                else:
                    groupe = groupe_figure_ASXL1[4]
                    
            else:
                print("Error")
            # Assign value of group ASXL1
            # Warning copy
            data_AB.loc[(data_AB.PATIENT == name) & (data_AB.Gene_refGene  == gene_interet),"FILTER"] = groupe 
          
    # Keep only data use in graphics
    data_keep = data_AB[data_AB["FILTER"].isin(groupe_figure_ASXL1)]
    data_keep[["PATIENT","DRIVER_USE_Figure","TYPE_DRIVER_USE_Figure","VAF_DRIVER_USE_Figure","ID","Gene_refGene","Classification_mutation","Mean_VAF","FILTER"]].to_csv("Figure/Data_preparation_figs7_final.csv",sep=",")
    return(data_keep)

def figure_s7(data_group_ASXL1,size_label,seuil_rep,quality_dpi,size_figure,name_repertory):
    """
        IN : New Data
        OUT : Distribution of allele burden ASXL1 according to 4-tiers classification
    """
    # Parameter
    plt.figure(figsize=size_figure)
    order_graph = ["TP53","SRSF2|RUNX1|IDH1/2|NRAS|KRAS","U2AF1","EZH2|CBL","Others"]
    sns.set(style="whitegrid")
    sns.set_context('paper')
    
    ax = sns.violinplot(x= "FILTER", y="Mean_VAF",order=order_graph, data=data_group_ASXL1,palette= "Paired",cut=0,width=0.8,linewidth=3)
    
    # Labels
    plt.xlabel(' ', fontsize=size_label)
    plt.ylabel("Distribution of ASXL1 VAF (%)", fontsize=size_label)
    
    # Subtitle Effectif 
    nobs = data_group_ASXL1['FILTER'].value_counts()
    nobs_reindex = nobs.reindex(order_graph)
    nobs = [str(x) for x in nobs_reindex.tolist()]
    nobs = ["n= " + i for i in nobs]
    pos = range(len(nobs))
    for tick,label in zip(pos,ax.get_xticklabels()):
        ax.text(pos[tick], 0.3, nobs[tick], horizontalalignment='center', size='14', color='black', weight='semibold')
    plt.xticks(fontsize=size_label)
    plt.yticks(fontsize=size_label)
    sns.despine(top = True, bottom = False, left = False, right = True)
    
    # Save figure S7
    plt.savefig(name_repertory + 'Fig_supplementary_7.svg',dpi=quality_dpi,bbox_inches='tight')
    
    # Close figure
    plt.close()
    
    # Description of figure S7
    describe_figure =  data_group_ASXL1.groupby('FILTER')['Mean_VAF'].describe()
    describe_figure.to_csv( name_repertory + "Data_Describe_FigureS7.csv",sep = ",")


# Creation table comptage et du fichier de statistique
def creation_file_statistic(data,name_repertory):
    """ 
        IN: DATA
        OUT: Statistic file 
             New Data
    """
    # Extraction of patient data
    data_comptage = data.groupby(['PATIENT'])[["PATIENT","SEX","TYPE_MYELOFIBROSIS"]].first()
    # Extraction of classification data
    data_class = data[data['Classification_mutation'] != "C"]

    return(data_class,data_comptage) 

def initialisation_matrice(index_gene):
    """
        OUT :Initialisation of matrix
    """
    return(pd.DataFrame(0.0, index = sorted(index_gene), columns = sorted(index_gene)))


# Matrix implementation Circos
def implementation_matrix(gene_comptage,matrix_compte,data,liste_patient,driver):
    """
         IN : Data and Count table 
         OUT : implementation of  matrix 
    """
    # If patient have more 1 variant in same gene we count an occurence of 1
    gene_comptage[gene_comptage != 0] = 1
    
    # Group data by patient
    grouped = gene_comptage.groupby('PATIENT')
    # Fo reach patient
    for name, group in grouped:
        # Check patient Suppression attente review
        if name not in liste_patient:
            continue    
        # Take only columns which have variant in gene 
        group_value = group.loc[:, (group != 0).any(axis=0)]
        # recover name of gene
        col_name = group_value.columns
        # First course
        for i in range(len(group_value.columns)):
            # All link gene1-gene1 are no show in circos figure
            matrix_compte.loc[col_name[i],col_name[i]] = 0
            
            # Second course
            for j in range(i + 1,len(group_value.columns)):
                # Implementation
                minimum = np.minimum(group_value.loc[name,col_name[i]],group_value.loc[name,col_name[j]])
                matrix_compte.loc[col_name[j],col_name[i]] += minimum 
    
    return(matrix_compte)

def liste_patient_driver(file_stat):
    """
        IN: DATA
        OUT: List of patients
    """
    
    list_patient_driver = list(file_stat.index)

    return(list_patient_driver)

def Preparation_data_matrix(data_seuil,index_seuil):
    """
        IN: DATA
        OUT: Count Table of number of gene variant per patient
    """
    # Matrix creation
    gene_data_seuil = data_seuil.groupby(['PATIENT','Gene.refGene'],sort=True).size().reset_index().pivot(columns='Gene.refGene', index='PATIENT',values=0)    
    # Sort Data by name of column
    gene_data_seuil = gene_data_seuil.reindex(sorted(gene_data_seuil.columns), axis=1)
    # Replace null value by 0
    gene_data_seuil.fillna(0,inplace=True)

    return(gene_data_seuil)

def index_seuil_mutation(data,seuil):
    """
        IN: DATA
        OUT: FILTER DATA with only gene (Treshold >5 variants) 
    """
    df_plot = data.groupby('Gene.refGene').size().sort_values(ascending=False) 
    # filter data by treshold
    df_plot_filter = df_plot[(df_plot >= seuil)].sort_values(ascending=False)

    # List of genes which respect treshold
    index_gene = list(df_plot_filter.index)

    # data with filtered genes
    data_filter = data[(data["Gene.refGene"].isin(index_gene))]

    return(data_filter,index_gene)

def generation_matrix_seuil(data_stat,seuil_gene,type_driver,file_statistic,type_file,name_repertory): 
    """
       IN: DATA mody to create matrix 
       OUT: Matrix generation of link gene variant per patient
    """
    # Filter data with only gene more 5 variant ad index with all gene present show in matrix
    data_seuil,index_seuil = index_seuil_mutation(data_stat,seuil_gene)
    
    # Calcul Data of number of variant par gene per patient
    gene_data_seuil = Preparation_data_matrix(data_seuil,index_seuil)
    
    # Circos with driver Attente suppression
    liste_patient = liste_patient_driver(file_statistic)
        
    # Initialisation of matrix
    matrix_comptage_seuil = initialisation_matrice(index_seuil)
    
    # Main function function implementation of matrix
    matrix_driver = implementation_matrix(gene_data_seuil,matrix_comptage_seuil,data_seuil,liste_patient,type_driver)
    
    # Save matrix in repertory
    file_name_repertory = name_repertory + "Matrix_statistic_"  + type_file + "_" + str(seuil_gene) + "-" + type_driver  + ".csv"
    matrix_driver.to_csv(file_name_repertory, sep = '\t')


# ## Main function circos
def matrix_stat_circos(data,treshold_matrix,name_repertory_matrix):
    """
        IN: DATA
        OUT: Generation of statistic file and matrix file 
    """
    # Choice of category of classication mutation A or B
    type_cate = "A_B"
    driver_list = "All"
    # ************************* 
    # Generation of statistic file
    class_data,file_statistic = creation_file_statistic(data,name_repertory_matrix)      
    # ***********
    # Generation of matrix file
    generation_matrix_seuil(class_data,treshold_matrix,driver_list,file_statistic,type_cate,name_repertory_matrix)

# ******************************************************************
# ***************************** Main *******************************
# ******************************************************************
if __name__ == '__main__':
    # Path 
    path = os.getcwd()
    os.chdir(path)
    
    # Name of directory to save figure
    name_repertory_figure = "Figure/"
    # File Database
    name_files = "Data/Data_Variant_result_479MF.csv"
    # setting figure
    size_figure = (20,10)
    quality_dpi = 400
    size_label = 18
    # Treshold for gene figure 1-A, 1-B,1-C
    seuil_gene = 20
    
    # load files
    data = read_files(name_files,',')
    # *********
    # Call of main function Figure Genomic
    figure_genomic(data,size_label,quality_dpi,size_figure,seuil_gene,name_repertory_figure)

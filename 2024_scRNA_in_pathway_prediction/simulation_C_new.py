# coding=utf-8
import pandas as pd
import networkx as nx
import multiprocessing
import os
import sys
import itertools

Path = os.getcwd()
os.chdir(Path)

def clustering_coefficient(times, certain_int, edge_list_filename, clustering_coefficient_output):
    pathway_class_geneExp = pd.read_csv("all_gene_classification.csv")
    pathway_class_geneExp_unique = pathway_class_geneExp.drop(["Classification"], axis=1).drop_duplicates().reset_index(drop=True)
    wgcna = pd.read_csv(edge_list_filename, sep="\t") 
    wgcna = wgcna.loc[:, ["fromNode", "toNode", "weight"]]
    with open(clustering_coefficient_output, "a", encoding="UTF-8") as w:
        for n, m in pathway_class_geneExp.groupby("Classification"):
            num = certain_int
            m = m.drop(["Classification"], axis=1).drop_duplicates()
            rows = m.shape[0]
            for i in range(times):
                num = num + i
                m = pathway_class_geneExp_unique.sample(n=rows,random_state=num).reset_index(drop=True) 
                cro_b_df = wgcna[wgcna.fromNode.isin(m.Gene.to_list())].reset_index(drop=True) 
                cro_a_df = cro_b_df[cro_b_df.toNode.isin(m.Gene.to_list())].reset_index(drop=True)
                G = nx.from_pandas_edgelist(cro_a_df, 'fromNode', 'toNode', 'weight', create_using=nx.Graph())
                T = nx.average_clustering(G, weight='weight')
                S = nx.transitivity(G)
                w.write(f"{n}\t{T}\t{S}\t{num}\n")

def GeneExpr_density_plot_data(clustering_coefficient_output, density_plot_filename):
    cl_df = pd.read_csv(clustering_coefficient_output, sep='\t', header=None)
    cl_df.columns = ["classification", "value", "global_value", "random_state"]
    out = {}
    for each_class in cl_df.classification.value_counts().index.to_list():
        df = cl_df[cl_df.classification.isin([each_class])].reset_index(drop=True)
        value = df.value.to_list()
        colname = df.classification.unique()[0]
        out[colname] = value
    pd.DataFrame(out).to_csv(density_plot_filename, index=False)

def GeneExpr_real_C_value(GeneExpr_Refer_value, edge_list_filename):
    pathway_class_geneExp = pd.read_csv("all_gene_classification.csv")
    out = open(GeneExpr_Refer_value, "w", encoding="UTF-8")
    out.write("class, Refervalue, CheckPoint\n")
    wgcna = pd.read_csv(edge_list_filename, sep="\t")
    for n, m in pathway_class_geneExp.groupby("Classification"):
        m = m.drop(["Classification"], axis=1).drop_duplicates()
        cro_b_df = wgcna[wgcna.fromNode.isin(m.Gene.to_list())].reset_index(drop=True)
        cro_a_df = cro_b_df[cro_b_df.toNode.isin(m.Gene.to_list())].reset_index(drop=True)
        G = nx.from_pandas_edgelist(cro_a_df, 'fromNode', 'toNode', 'weight', create_using=nx.Graph())
        T = nx.average_clustering(G, weight='weight')
        S = nx.transitivity(G)
        out.write(f"{n},{T},{S}\n")
    out.close()

def Complete_network(edge_list_filename):
    all_gene = pd.read_csv("all_gene_classification.csv")
    wgcna = pd.read_csv(edge_list_filename, sep="\t")
    wgcna = wgcna.loc[:, ['fromNode', 'toNode', 'weight']]
    wgcna.loc[:, 'Node_pair'] = wgcna.loc[:, ["fromNode", "toNode"]].min(axis=1).map(str) + "," + wgcna.loc[:,["fromNode","toNode"]].max(axis=1)
    out = []
    for i in list(itertools.permutations(all_gene.Gene.unique().tolist(), 2)):
        i = min(i) + "," + max(i)
        out.append(i)
    all_gene_df = pd.DataFrame(out).drop_duplicates().reset_index(drop=True)
    all_gene_df.columns = ['Node_pair']
    complete_network = pd.merge(all_gene_df, wgcna, left_on="Node_pair", right_on="Node_pair", how="left").drop(['fromNode', 'toNode'], axis=1).fillna(0)
    test_df = complete_network.Node_pair.str.split(',', expand=True)
    test_df.columns = ["fromNode","toNode"]
    test_df['weight'] = complete_network['weight']
    test_df.to_csv(edge_list_filename, sep = '\t', index=False)



if __name__ == '__main__':
    path = os.getcwd()
    os.chdir(path)
    guest_geneExp_filename = sys.argv[1]
    ab_name = guest_geneExp_filename.split(".")[0]
    R_WGCNA_script = "Rscript WGCNA.R"+ " " + guest_geneExp_filename + " " + ab_name
    os.system(R_WGCNA_script)

    edge_list_filename = "CytoscapeInput_edges_" + ab_name + ".txt"
    clustering_coefficient_output = "clustering_coefficient_" + ab_name + ".txt"
    density_plot_filename = "plot_data_" + ab_name + ".csv"
    GeneExpr_Refer_value = "Refervalue_" + ab_name + ".csv"

    Complete_network(edge_list_filename)
    GeneExpr_real_C_value(GeneExpr_Refer_value, edge_list_filename)

    mission = multiprocessing.Pool(50)
    for i in range(50):
        certain_int = 200 * i
        mission.apply_async(clustering_coefficient, args=(200, certain_int, edge_list_filename, clustering_coefficient_output))
    mission.close()
    mission.join()

    GeneExpr_density_plot_data(clustering_coefficient_output, density_plot_filename)
    R_density_plot_script = "Rscript Simulation_C_density_plot_230328.R" + " " + density_plot_filename + " " + GeneExpr_Refer_value + " " + ab_name
    os.system(R_density_plot_script)
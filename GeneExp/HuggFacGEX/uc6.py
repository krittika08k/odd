import streamlit as st
import matplotlib.pyplot as plt
import scanpy as sc

# Add a dropdown to select the dataset
datasets = {
    "pbmc3k": sc.datasets.pbmc3k_processed,
    
    "Subsampled and processed 68k PBMCs": sc.datasets.pbmc68k_reduced,
    "Development of Myeloid Progenitors": sc.datasets.paul15
}
dataset_name = st.sidebar.selectbox("Select a dataset", list(datasets.keys()))

# Load the data
data_func = datasets[dataset_name]
adata = data_func()

# Compute nearest neighbors and UMAP coordinates
sc.pp.neighbors(adata)
sc.tl.umap(adata)

# Define widgets for user input
cluster_method = st.sidebar.selectbox('Clustering method', ['leiden', 'louvain'])
dim_reduction_method = st.sidebar.selectbox('Dimension reduction method', ['umap', 'tsne'])
gene = st.sidebar.text_input('Gene')
gene_set = st.sidebar.text_input('Gene set (comma separated)')

# Define plotting functions
def plot_cluster(adata, method):
    sc.pl.umap(adata, color=method, show=False)
    plt.title(f'Cluster method: {method}')
    st.set_option('deprecation.showPyplotGlobalUse', False)
    st.pyplot()

def plot_gene(adata, gene):
    sc.pl.umap(adata, color=gene, show=False)
    plt.title(f'Gene expression for {gene}')
    st.set_option('deprecation.showPyplotGlobalUse', False)
    st.pyplot()

def plot_gene_set(adata, gene_set):
    genes = gene_set.split(',')
    sc.pl.umap(adata, color=genes, show=False)
    plt.title(f'Gene set: {gene_set}')
    st.set_option('deprecation.showPyplotGlobalUse', False)
    st.pyplot()

# Define main function
def main(adata, cluster_method, dim_reduction_method, gene, gene_set):
    if cluster_method:
        if cluster_method == 'leiden':
            sc.tl.leiden(adata)
        elif cluster_method == 'louvain':
            sc.tl.louvain(adata)
        plot_cluster(adata, cluster_method)
    if dim_reduction_method:
        if dim_reduction_method == 'umap':
            plot_cluster(adata, cluster_method)
        elif dim_reduction_method == 'tsne':
            sc.tl.tsne(adata)
            plot_cluster(adata, cluster_method)
    if gene:
        try:
            if gene not in adata.var_names:
                st.error(f'Gene {gene} not found in dataset.')
            else:
                plot_gene(adata, gene)
        except KeyError:
            st.error(f'Gene {gene} not found in dataset.')
    if gene_set:
        plot_gene_set(adata, gene_set)

    

# Run the main function based on user input
if st.button('Update'):
    main(adata, cluster_method, dim_reduction_method, gene, gene_set)

if st.button('Refresh'):
    adata = sc.datasets.pbmc3k()
    sc.pp.neighbors(adata)
    sc.tl.umap(adata)
    st.success('Data refreshed.')
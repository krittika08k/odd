import scanpy as sc
import streamlit as st
from sklearnex import patch_sklearn
from sklearn.ensemble import RandomForestClassifier
from xgboost import XGBClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.tree import DecisionTreeClassifier
from sklearn.neighbors import KNeighborsClassifier
from sklearn.metrics import f1_score, confusion_matrix, accuracy_score
import plotly.express as px
import plotly.graph_objects as go

# Patch scikit-learn
patch_sklearn()

# Define the Streamlit app
st.title("Single-cell RNA sequencing analysis")
st.title("Powered by Intel oneAPI AI Analytics toolkit")

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

# Preprocess the data
sc.pp.filter_cells(adata, min_genes=200)
sc.pp.filter_genes(adata, min_cells=3)
sc.tl.pca(adata)
sc.pp.neighbors(adata)
sc.tl.umap(adata)
sc.tl.louvain(adata)
sc.tl.paga(adata)

# Select classification algorithm
clf = st.sidebar.selectbox(
    "Select a classification algorithm",
    ("Random Forest", "XGBoost", "Logistic Regression", "Decision Tree", "K-Nearest Neighbors"),
)

# Train the model
X = adata.X
y = adata.obs['louvain'].astype('category').cat.codes

if clf == "Random Forest":
    model = RandomForestClassifier(n_estimators=100, random_state=0)
elif clf == "XGBoost":
    model = XGBClassifier(use_label_encoder=False)
elif clf == "Logistic Regression":
    model = LogisticRegression(random_state=0)
elif clf == "Decision Tree":
    model = DecisionTreeClassifier(random_state=0)
else:
    model = KNeighborsClassifier()

model.fit(X, y)

cell_idx = st.sidebar.selectbox("Select a cell", adata.obs.index)

if st.button("Predict cancer"):
    cell_data = adata[adata.obs.index == cell_idx].X
    prediction = model.predict(cell_data)

    if prediction[0] == 0:
        st.write("This cell is not cancerous")
    else:
        st.write("This cell is cancerous")

    # Create a scatter plot with the predicted cancer status highlighted in red and blue
    adata.obs['predicted'] = model.predict(adata.X)
    fig1 = px.scatter(adata.obsm['X_umap'], x=0, y=1, color=adata.obs['louvain'])
    fig1.update_layout(title='Clusters (Louvain)')
    fig2 = px.scatter(adata.obsm['X_umap'], x=0, y=1, color=adata.obs['predicted'])
    fig2.update_layout(title='Predicted cancer status')
    st.plotly_chart(fig1)
    st.plotly_chart(fig2)
    
    # Add PAGA plot with DPT values
    # Add pseudotime as a color dimension to the scatter plot
    # Create PAGA plot
    
else:
   
    # Create a scatter plot with no highlighting
    fig = px.scatter(adata.obsm['X_umap'], x=0, y=1, color=adata.obs['louvain'])
    fig.update_layout(title='Clusters (Louvain)')
    st.plotly_chart(fig)

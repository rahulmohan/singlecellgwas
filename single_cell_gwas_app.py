import pandas as pd
import numpy as np
import plotly.graph_objects as go
import streamlit as st 
from streamlit_plotly_events import plotly_events

st.title("Single Cell GWAS Explorer")
st.write("Size of circles indicate enrichment score, opacity of circles indicate -log(p-value).")
pairs = pd.read_csv("trait_celltype_pairs.csv")
app_startup = True

def create_matrix(trait_categories, tissue_categories):
	traits = pairs[pairs["trait category"].isin(trait_categories)]["traits"].values
	pairs["celltypes"] = pairs["celltypes"].map(lambda x: str(x).split("(")[0].strip())
	celltypes = pairs[pairs["tissue"].isin(tissue_categories)]["celltypes"].values

	traits = np.unique(traits)
	celltypes = np.unique(celltypes)
	
	trait_celltype_escores = np.zeros((len(traits), len(celltypes)))
	trait_celltype_pvals = np.zeros((len(traits), len(celltypes)))
	trait_celltype_pvals_orig = np.zeros((len(traits), len(celltypes)))

	for i in range(len(traits)):
		trait_pairs = pairs.loc[pairs["traits"] == traits[i]]
		for j in range(len(celltypes)):
			if celltypes[j] in trait_pairs["celltypes"].values:
				trait_celltype_escores[i, j] = trait_pairs[trait_pairs["celltypes"] == celltypes[j]]["Escore"].values[0]
				trait_celltype_pvals[i, j] = -np.log(trait_pairs[trait_pairs["celltypes"] == celltypes[j]]["pEscore"].values[0])
				trait_celltype_pvals_orig[i, j] = trait_pairs[trait_pairs["celltypes"] == celltypes[j]]["pEscore"].values[0]
			else:
				trait_celltype_escores[i, j] = 0
				trait_celltype_pvals[i, j] = 0

	return trait_celltype_escores, trait_celltype_pvals, trait_celltype_pvals_orig, traits, celltypes

def generate_heatmap(traits_selected, tissues_selected):
	trait_celltype_escores, trait_celltype_pvals, original_pvals, y_axis, x_axis = create_matrix(traits_selected, tissues_selected)
	if trait_celltype_pvals.max() == trait_celltype_pvals.min():
		trait_celltype_pvals = (trait_celltype_pvals - trait_celltype_pvals.min()) / (1 - trait_celltype_pvals.min())
	else:
		trait_celltype_pvals = (trait_celltype_pvals - trait_celltype_pvals.min()) / (trait_celltype_pvals.max() - trait_celltype_pvals.min())

	x = []
	y = []
	escores = []
	pvals = []
	texts = []
	for i in range(len(y_axis)):
		for j in range(len(x_axis)):
			x.append(j)
			y.append(i)
			escores.append(trait_celltype_escores[i, j] * 5)
			pvals.append(trait_celltype_pvals[i, j])
			texts.append("P-Value: " + str(original_pvals[i, j]) + "<br>" + "E-Score: " + str(trait_celltype_escores[i, j]))

	fig = go.Figure(data=[go.Scatter(
	    x=x, y=y,
	    mode='markers',
	    text=texts,
	    marker=dict(
			color = ["red"] * len(x) * len(y),
	    	opacity = pvals,
	    	size = escores
	    )
	)])

	fig.update_layout(
	    xaxis = dict(
	        tickmode = 'array',
			tickvals=[i for i in range(len(x_axis))], 
			ticktext=x_axis
	    )
	)

	fig.update_layout(
	    yaxis = dict(
	        tickmode = 'array',
			tickvals=[i for i in range(len(y_axis))], 
			ticktext=y_axis
	    )
	)

	fig.update_layout(
	    autosize=True,
	 )

	def on_click(trace, points, selector):
	    inds = points.point_inds
	    st.write(str(inds[0]))


	fig.on_click(update_point)


	return fig

trait_categories = list(np.unique(pairs["trait category"].values))
tissue_categories = list(np.unique(pairs["tissue"].values))

traits_selected = st.sidebar.multiselect("Select Trait Categories: ", trait_categories, default=["blood test"])
tissues_selected = st.sidebar.multiselect("Select Tissue: ", tissue_categories, default=["blood"])

if len(traits_selected) > 0 and len(tissues_selected) > 0:
	fig = generate_heatmap(traits_selected, tissues_selected)
	selected_points = plotly_events(fig, click_event=True, hover_event=False)
	st.plotly_chart(fig)

st.markdown(
    """
    <style>
    [data-testid="stSidebar"][aria-expanded="true"] > div:first-child {
        width: 500px;
    }
    [data-testid="stSidebar"][aria-expanded="false"] > div:first-child {
        width: 500px;
        margin-left: -500px;
    }
    </style>
    """,
    unsafe_allow_html=True,
)

st.markdown(
        f"""
<style>
    .reportview-container .main .block-container{{
        margin-left: -275px;
    }}

</style>
""",
        unsafe_allow_html=True,
    )

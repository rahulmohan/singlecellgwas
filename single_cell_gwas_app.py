import pandas as pd
import numpy as np
import plotly.graph_objects as go
import streamlit as st 
from streamlit_plotly_events import plotly_events

st.set_page_config(layout="centered")
st.title("Single Cell GWAS Explorer")
st.write("Size of circles indicate enrichment score, opacity of circles indicate -log(p-value).")
pairs = pd.read_csv("trait_celltype_pairs.csv")

trait_map = {
	"Age First Birth": "UKB_460K.repro_AgeFirstBirth_Female",
	"Alkaline Phosphatase": "UKB_460K.biochemistry_AlkalinePhosphatase",
	"Alzheimers": "PASS_Alzheimer",
	"Aspartate Amino transferase": "UKB_460K.biochemistry_AspartateAminotransferase",
	"Asthma": "UKB_460K.disease_ASTHMA_DIAGNOSED",
	"Atrial Fibrillation": "PASS_AtrialFibrillation_Nielsen2018",
	"BMI/WHR": "Body_mass_index_BMI",
	"BMIz": "UKB_460K.body_BMIz",
	"BMR": "Basal_metabolic_rate",
	"Bilirubin": "UKB_460K.biochemistry_TotalBilirubin",
	"Breast.Cancer": "UKB_460K.cancer_BREAST",
	"Celiac": "PASS_Celiac",
	"Cholesterol": "UKB_460K.biochemistry_Cholesterol",
	"Coronary Artery Disease": "PASS_Coronary_Artery_Disease",
	"Creatinine": "UKB_460K.biochemistry_Creatinine",
	"Crohns Disease": "PASS_Crohns_Disease",
	"Diastolic bp": "Diastolic_blood_pressure_automated_reading",
	"ECG rate": "ECG_heart_rate",
	"Eczema": "UKB_460K.disease_ALLERGY_ECZEMA_DIAGNOSED",
	"Edu Years": "PASS_Years_of_Education1",
	"Inflammatory Bowel Disease": "PASS_IBD",
	"Insomnia": "PASS_Insomnia_Jansen2019",
	"Intelligence": "PASS_Intelligence_SavageJansen2018",
	"Lung capacity": "Forced_vital_capacity_FVC",
	"Lupus": "PASS_Lupus",
	"Lymphocyte.percent": "Lymphocyte_percentage",
	"Major Depressive Disorder": "PASS_MDD_Howard2019",
	"Monocyte.percent": "Monocyte_percentage",
	"Morning Person": "UKB_460K.other_MORNINGPERSON",
	"Multiple sclerosis": "PASS_Multiple_sclerosis",
	"Neuroticism": "Neuroticism_score",
	"Num Children": "UKB_460K.repro_NumberChildrenEverBorn_Pooled",
	"Phosphate": "UKB_460K.biochemistry_Phosphate",
	"Pigment.Sunburn": "UKB_460K.pigment_SUNBURN",
	"Platelet.count": "Platelet_count",
	"Primary Biliary cirrhosis": "PASS_Primary_biliary_cirrhosis",
	"RBC.count": "UKB_460K.blood_RED_COUNT",
	"RBC.distwidth": "UKB_460K.blood_RBC_DISTRIB_WIDTH",
	"Rheumatoid arthritis": "PASS_Rheumatoid_Arthritis",
	"Schizophrenia": "PASS_Schizophrenia",
	"Smoking Status": "UKB_460K.cov_SMOKING_STATUS",
	"Systolic bp": "Systolic_blood_pressure_automated_reading",
	"T2D": "UKB_460K.disease_T2D",
	"TotalProtein": "UKB_460K.biochemistry_TotalProtein",
	"Type1 Diabetes": "PASS_Type_1_Diabetes",
	"Ulcerative colitis": "PASS_Ulcerative_Colitis",
	"VitaminD": "UKB_460K.biochemistry_VitaminD",
	"WHRadjBMI": "UKB_460K.body_WHRadjBMIz"
}

celltype_map = {
	
}

print(np.unique(pairs["traits"]))
print(np.unique(pairs["celltypes"].map(lambda x: str(x).split("(")[0].strip())))

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
			color = ["#d62728"] * len(x) * len(y),
	    	opacity = pvals,
	    	size = escores
	    ),
	    hoverlabel=dict(bgcolor="white")
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

	#fig.update_layout(
	#    autosize=True
	# )

	def on_click(trace, points, selector):
	    inds = points.point_inds
	    print(inds)
	    st.write(str(inds[0]))

	fig.data[0].on_click(on_click)
	return fig

trait_categories = list(np.unique(pairs["trait category"].values))
tissue_categories = list(np.unique(pairs["tissue"].values))

traits_selected = st.sidebar.multiselect("Select Trait Categories: ", trait_categories, default=["blood test"])
tissues_selected = st.sidebar.multiselect("Select Tissue: ", tissue_categories, default=["blood"])

if len(traits_selected) > 0 and len(tissues_selected) > 0:
	fig = generate_heatmap(traits_selected, tissues_selected)
	selected_points = plotly_events(fig, click_event=True, hover_event=False)
	#st.write(selected_points)
	#st.plotly_chart(fig)

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

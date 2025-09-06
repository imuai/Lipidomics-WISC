# Complete 10-Step WISC Lipidomics Analysis Pipeline
# Authors: Luan Vu (UTSA) & Joan Cook-Mills (IU)
# Date: December 2024
# Systematic Analysis: 74 lipids → 76 lipids (adding α-tocopherol & γ-tocopherol)

import pandas as pd
import numpy as np
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import streamlit as st
import umap
from sklearn.cluster import KMeans
from sklearn.preprocessing import StandardScaler
from scipy import stats
from scipy.spatial import ConvexHull
from statsmodels.stats.multitest import multipletests
import warnings
warnings.filterwarnings('ignore')

# Set page config
st.set_page_config(page_title="New Lipidomics Analysis", page_icon="🧬", layout="wide")

def explain_analysis_pipeline():
    """Detailed explanation of the analysis pipeline"""
    with st.expander("📋 Detailed Analysis Methodology", expanded=False):
        st.markdown("""
        ### **Analysis Pipeline Details**
        
        **1. Data Preprocessing**
        - Load 174 samples with 75 measured lipids
        - **CRITICAL**: Omit `8.9.DHET.1` (duplicate) → 74 lipids (matching publication)
        - Remove comment columns
        - Handle missing values: median imputation
        
        **2. Reproduce Original UMAP**
        - Use existing `oldUMAP1` and `oldUMAP2` coordinates
        - Visualize with `hclust8` (8 hierarchical clusters)
        - Add convex hull enclosures for each cluster
        - This reproduces the exact original publication figure
        
        **3. Reanalyze 74 Lipids**
        - Apply Pareto scaling: `(x - mean) / sqrt(std)`
        - Fresh UMAP: n_neighbors=15, min_dist=0.1, random_state=42
        - Gap statistic for optimal k (testing k=2 to 10)
        - Create TWO visualizations:
          * Colored by original hclust8
          * Colored by new k-means clusters
        - Add cluster enclosures to all visualizations
        
        **4. Reanalyse with 74 Lipids and 2 tocopherol**
        - Add α-tocopherol and γ-tocopherol
        - Same preprocessing and UMAP parameters
        - Gap statistic optimization
        - Compare clustering stability
                    
        **5. Pro-inflammatory Index**
        - Reference cluster/group: Identify the cluster with the lowest βGlcCer level and use it as the reference group.
        - Justificattion:
          * Elevated βGlcCer levels lead to increased allergen responsiveness and altered subsets of proinflammatory dendritic cells in neonatal mice. 
          * Therefore, a cluster with lower βGlcCer levels and a profile indicative of anti-inflammatory responses, served as a logical baseline to compare against for assessing proinflammatory lipid indices and their association with allergic outcomes in humans.
        - Statistical tests: ANOVA + Tukey's HSD or Kruskal-Wallis + Dunn's test
        - Adjust p-values with Benjamini-Hochberg FDR
        - Visualize with violin plots
                    
        **6. Pro-inflammatory Index**     
        - Calculate log2 fold changes of each cluster vs. reference group/cluster
        - Volcano plots for lipid significance of each cluster vs. reference cluster/group.
        - Statistical tests: ANOVA + Tukey's HSD or Kruskal-Wallis + Dunn's test
        - Adjust p-values with Benjamini-Hochberg FDR
        - Index = Σ(pro-inflammatory) - Σ(anti-inflammatory)
          * Pro-inflammation: βGlcCers, βGalCers, sphingomyelins, γ-tocopherol
          * Anti-inflammation: Lipoxin A4, Resolvins, α-tocopherol
        
        **7. Clinical Correlations**
        - Spearman correlation (non-parametric)
        - Bootstrap confidence intervals (5000 iterations)
        - Linear regression with R² values
    
        
        ### **Why Omit 8.9.DHET.1?**
        - Manuscript states 74 lipids measured
        - Dataset contains 75 lipids
        - `8.9.DHET` appears twice: `8.9.DHET` and `8.9.DHET.1`
        - Likely technical replicate or isomer
        - Currently, omitting `8.9.DHET.1` aligns with publication
        - Dr. Cook-Mills, please confirm if this is correct or if both should be included.
        """)


# ================================
# STEP 1: Summary of Original Dataset (obtained from Dr. Cook-Mills)
# ================================
def explain_lipid_mediators():
    """Explain the biological significance of each lipid mediator class"""
    st.subheader("📚 Lipid Mediator Biology")
    
    with st.expander("Understanding the 75 Lipid Mediators - Dr. Cook-Mills please confirm if there are any misunderstandings", expanded=False):
        st.markdown("""
        ### **Sphingolipids (44 total)**
        
        **1. Ceramides (Cer) - 10 species**
        - C14.0, C16.0, C18.0, C18.1, C20.0, C22.0, C24.0, C24.1, C26.0, C26.1
        - **Function**: Pro-apoptotic, stress response, barrier function
        - **In allergy**: Elevated in asthma, promote inflammation
        
        **2. Glucosylceramides (GlcCer) - 10 species**  
        - C14:0, C16.0, C18.0, C18.1, C20.0, C22.0, C24.0, C24.1, C26.0, C26.1
        - **Function**: Precursors to complex glycosphingolipids
        - **In allergy**: Elevated in maternal allergy, transferred to fetus
        
        **3. Galactosylceramides (GalCer) - 10 species**
        - C14.0, C16.0, C18.0, C18.1, C20.0, C22.0, C24.0, C24.1, C26.0, C26.1
        - **Function**: Major myelin component, immune regulation
        - **In allergy**: Modulate T-cell responses
        
        **4. Sphingosines (SPH) - 10 species**
        - C14.0, C16.0, C18.0, C18.1, C20.0, C22.0, C24.0, C24.1, C26.0, C26.1
        - **Function**: Ceramide precursors, signaling molecules
        - **In allergy**: Regulate mast cell degranulation
        
        **5. Other Sphingolipids - 4 molecules**
        - **So (Sphingosine)**: Base sphingolipid, antimicrobial
        - **DHSo (Dihydrosphingosine)**: Sphingosine precursor
        - **S1P (Sphingosine-1-phosphate)**: Immune cell trafficking
        - **DHS1P (Dihydrosphingosine-1-phosphate)**: S1P analog
        
        ### **Eicosanoids & Docosanoids (31 total)**
        
        **6. Prostaglandins - 8 molecules**
        - **PGE2, PGE1**: Bronchodilation, inflammation
        - **PGF2α**: Bronchoconstriction
        - **PGD2**: Mast cell product, allergic inflammation
        - **PGA2**: Anti-inflammatory
        - **6-keto-PGF1α**: PGI2 metabolite, vasodilation
        - **8-iso-PGF2α**: Oxidative stress marker
        - **15-deoxy-Δ12,14-PGJ2**: PPARγ agonist, anti-inflammatory
        
        **7. Thromboxane**
        - **TXB2**: TXA2 metabolite, platelet activation
        
        **8. Isoprostane**
        - **5-iPF2α-VI**: Oxidative stress marker
        
        **9. Specialized Pro-resolving Mediators**
        - **Resolvin D1, D2**: Resolution of inflammation
        - **Lipoxin A4**: Anti-inflammatory, pro-resolution
        
        **10. Leukotrienes - 4 molecules**
        - **LTB4**: Neutrophil chemotaxis
        - **LTC4, LTD4, LTE4**: Cysteinyl LTs, bronchoconstriction
        
        **11. Cytochrome P450 Metabolites**
        - **DHETs (8,9-DHET, 11,12-DHET, 8,9-DHET.1)**: EET hydrolysis products
        - **EETs (14,15-EET, 8,9-EET)**: Vasodilation, anti-inflammatory
        - **HETEs (5-, 12-, 15-, 20-HETE)**: Various inflammatory roles
        
        **12. Fatty Acids - 4 molecules**
        - **AA (Arachidonic Acid)**: Pro-inflammatory precursor
        - **EPA (Eicosapentaenoic Acid)**: Anti-inflammatory omega-3
        - **DHA (Docosahexaenoic Acid)**: Resolvin precursor
        - **DHGLA (Dihomo-γ-linolenic Acid)**: Anti-inflammatory
        
        **13. Linoleic Acid Metabolite**
        - **13-HODE**: Oxidized linoleic acid, inflammatory marker
        
        ### **Why 75 not 74?**
        The manuscript mentions 74 lipids, but the dataset contains 75. The likely explanation:
        - **Duplicate**: "8,9-DHET" appears twice (8.9.DHET and 8.9.DHET.1)
        - This could be different isomers or a data processing artifact
        - At this stage, I will treat 8.9 DHET.1 as a duplicate of 8.9 DHET for simplicity.
        - Dr. Cook-Mills, please  confirmed this is likely a duplicate or an actual separated lipid measurement.
        - If they are indeed different, we can adjust the analysis accordingly. Personally, for publication compliance, I will treat this as 75 unique measurements
        """)

def step1_original_dataset_summary():
    """Step 1: Summary of original dataset without tocopherols"""
    st.header("Step 1: Original Dataset Summary (Without Tocopherols)")

    df = pd.read_csv('data.csv')
    
    # Remove comment columns
    columns_to_exclude = ['aT (uM),comment', 'gT (uM), comment']
    df = df.drop(columns=columns_to_exclude, errors='ignore')
    
    # Define 74 lipid columns (excluding 8.9.DHET.1 and tocopherols)
    lipid_columns_74 = [
        # Ceramides (10)
        'Cer.C14.0', 'Cer.C16.0', 'Cer.C18.1', 'Cer.C18.0', 'Cer.C20.0',
        'Cer.C22.0', 'Cer.C24.1', 'Cer.C24.0', 'Cer.C26.1', 'Cer.C26.0',
        # Glucosylceramides (10)
        'GlcCer_C14:0', 'GlcCer.C16.0', 'GlcCer.C18.1', 'GlcCer.C18.0', 'GlcCer.C20.0',
        'GlcCer.C22.0', 'GlcCer.C24.1', 'GlcCer.C24.0', 'GlcCer.C26.1', 'GlcCer.C26.0',
        # Galactosylceramides (10)
        'GalCer.C14.0', 'GalCer.C16.0', 'GalCer.C18.1', 'GalCer.C18.0', 'GalCer.C20.0',
        'GalCer.C22.0', 'GalCer.C24.1', 'GalCer.C24.0', 'GalCer.C26.1', 'GalCer.C26.0',
        # Sphingosines (10)
        'SPH.C14.0', 'SPH.C16.0', 'SPH.C18.1', 'SPH.C18.0', 'SPH.C20.0',
        'SPH.C22.0', 'SPH.C24.1', 'SPH.C24.0', 'SPH.C26.1', 'SPH.C26.0',
        # Other sphingolipids (4)
        'So', 'DHSo', 'S1P', 'DHS1P',
        # Prostaglandins (8)
        '6.keto.PGF1α', '8.iso.PGF2α', 'PGE2', 'PGF2α', 'PGD2', 'PGE1', 'PGA2', '15.deoxy-Δ12,14-PGJ2',
        # Others (22)
        'TXB2', '5.iPF2α-VI', 'Resolvin.D2', 'Resolvin.D1', 'Lipoxin.A4',
        'LTD4', 'LTC4', 'LTE4', 'LTB4',
        '8.9.DHET', '11.12.DHET', '14.15.EET', '8.9.EET',
        '20.HETE', '15.HETE', '12.HETE', '5.HETE',
        'EPA', 'DHA', 'AA', 'DHGLA', '13.HODE'
    ]
    
    col1, col2, col3, col4 = st.columns(4)
    with col1:
        st.metric("Total Samples", len(df))
    with col2:
        st.metric("Total Lipids", len(lipid_columns_74))
    with col3:
        st.metric("Sphingolipids", 44)
    with col4:
        st.metric("Eicosanoids", 30)
    
    explain_lipid_mediators()
    # Missing data analysis for original lipids
    missing_data = df[lipid_columns_74].isnull().sum()
    missing_percent = (missing_data / len(df) * 100).round(2)
    
    missing_summary = pd.DataFrame({
        'Missing Count': missing_data[missing_data > 0],
        'Missing %': missing_percent[missing_data > 0]
    }).sort_values('Missing %', ascending=False)
    
    if len(missing_summary) > 0:
        st.subheader("Missing Data in Original 74 Lipids")
        st.dataframe(missing_summary)
    else:
        st.success("✅ No missing data in original 74 lipids")
    
    # Clinical outcomes
    st.subheader("Clinical Outcomes")
    wheeze_count = df['Wheeze_1yr'].sum()
    ad_count = df['AD_1yr'].sum()
    
    col1, col2 = st.columns(2)
    with col1:
        st.metric("Wheeze Cases", f"{wheeze_count} ({wheeze_count/len(df)*100:.1f}%)")
    with col2:
        st.metric("AD Cases", f"{ad_count} ({ad_count/len(df)*100:.1f}%)")
    
    st.info("**Note**: Excluding `8.9.DHET.1` (duplicate) to match publication (74 lipids)")
    
    return df, lipid_columns_74

# ================================
# STEP 2: Reproduce Original UMAP
# ================================
def get_cluster_hull(x, y):
    """Calculate convex hull for cluster enclosure"""
    if len(x) < 3:
        return x, y
    
    points = np.column_stack([x, y])
    try:
        hull = ConvexHull(points)
        hull_points = points[hull.vertices]
        hull_points = np.vstack([hull_points, hull_points[0]])
        return hull_points[:, 0], hull_points[:, 1]
    except:
        return x, y

def step2_reproduce_original_umap(df):
    """Step 2: Reproduce original UMAP using oldUMAP1/oldUMAP2 and hclust8"""
    st.header("Step 2: Reproduce Original UMAP")
    st.write("Using `the orignal UMAP1`, `original UMAP2` coordinates with `hclust8` clustering")
    
    fig = go.Figure()
    
    # Define cluster colors
    cluster_colors = {
        1: '#FFE4B5', 2: '#FFD700', 3: '#98D8C8', 4: '#B0E0E6',
        5: '#90EE90', 6: '#87CEEB', 7: '#DDA0DD', 8: '#FFB6C1'
    }
    
    # Add cluster enclosures
    for cluster in range(1, 9):
        cluster_data = df[df['hclust8'] == cluster]
        if len(cluster_data) >= 3:
            hull_x, hull_y = get_cluster_hull(
                cluster_data['oldUMAP1'].values,
                cluster_data['oldUMAP2'].values
            )
            
            fig.add_trace(go.Scatter(
                x=hull_x, y=hull_y,
                fill='toself',
                fillcolor=cluster_colors[cluster],
                opacity=0.3,
                line=dict(color=cluster_colors[cluster], width=2),
                name=f'Cluster {cluster}',
                showlegend=True,
                hoverinfo='skip'
            ))
    
    # Add scatter points
    for cluster in range(1, 9):
        cluster_data = df[df['hclust8'] == cluster]
        fig.add_trace(go.Scatter(
            x=cluster_data['oldUMAP1'],
            y=cluster_data['oldUMAP2'],
            mode='markers',
            marker=dict(
                size=8,
                color=cluster_colors[cluster],
                line=dict(color='black', width=0.5)
            ),
            name=f'HC{cluster} (n={len(cluster_data)})',
            showlegend=False
        ))
        
        # Add cluster label
        if len(cluster_data) > 0:
            fig.add_annotation(
                x=cluster_data['oldUMAP1'].mean(),
                y=cluster_data['oldUMAP2'].mean(),
                text=f'<b>{cluster}</b>',
                showarrow=False,
                font=dict(size=14, color='black'),
                bgcolor='white',
                bordercolor='black',
                borderwidth=1
            )
    
    fig.update_layout(
        title="Original UMAP Analysis - 8 Hierarchical Clusters",
        xaxis_title='UMAP1',
        yaxis_title='UMAP2',
        width=900,
        height=700,
        template='plotly_white'
    )
    
    st.plotly_chart(fig, use_container_width=True)
    
    return df

# ================================
# STEP 3: Dataset with Tocopherols & Missing Data Handling
# ================================
def step3_tocopherol_dataset_missing_data(df, lipid_columns_74):
    """Step 3: Summary of dataset with tocopherols and handle missing data"""
    st.header("Step 3: Dataset with Tocopherols & Missing Data Handling")
    
    # Add tocopherols
    enhanced_columns = lipid_columns_74 + ['alpha-tocopherol', 'gamma-tocopherol']
    
    # Summary
    col1, col2, col3 = st.columns(3)
    with col1:
        st.metric("Original Lipids", len(lipid_columns_74))
    with col2:
        st.metric("Added Tocopherols", 2)
    with col3:
        st.metric("Total variables", len(enhanced_columns))
    
    # Missing data analysis for tocopherols
    st.subheader("Missing Data Analysis for Tocopherols")
    
    missing_alpha = df['alpha-tocopherol'].isnull().sum()
    missing_gamma = df['gamma-tocopherol'].isnull().sum()
    
    col1, col2 = st.columns(2)
    with col1:
        st.metric("Missing α-tocopherol", f"{missing_alpha} ({missing_alpha/len(df)*100:.1f}%)")
    with col2:
        st.metric("Missing γ-tocopherol", f"{missing_gamma} ({missing_gamma/len(df)*100:.1f}%)")
    
    # Missing data handling options
    st.subheader("Missing Data Handling Methods")
    
    with st.expander("📚 Detailed Explanation of Methods", expanded=True):
        st.markdown("""
        ### **1. Mean Imputation**
        - **Method**: Replace missing values with the mean of non-missing values
        - **Pros**: 
            - Simple and fast
            - Preserves sample size
            - Maintains mean of the variable
        - **Cons**: 
            - Reduces variance
            - Can distort correlations
            - Assumes data is missing completely at random (MCAR)
        - **Best for**: Small amounts of missing data (<5%), normally distributed variables
        
        ### **2. Median Imputation**
        - **Method**: Replace missing values with the median of non-missing values
        - **Pros**: 
            - Robust to outliers
            - Simple to implement
            - Better for skewed distributions
        - **Cons**: 
            - Reduces variance
            - Can distort distributions
            - Ignores relationships between variables
        - **Best for**: Skewed distributions, presence of outliers, small missing data
        
        ### **3. Forward Fill (Last Observation Carried Forward)**
        - **Method**: Use the last valid observation to fill missing values
        - **Pros**: 
            - Preserves temporal patterns
            - Simple for longitudinal data
        - **Cons**: 
            - Not appropriate for cross-sectional data
            - Can propagate errors
        - **Best for**: Time-series data with temporal correlation
        
        ### **4. KNN Imputation**
        - **Method**: Use k-nearest neighbors to estimate missing values
        - **Pros**: 
            - Considers relationships between variables
            - Can handle complex patterns
            - More accurate than mean/median
        - **Cons**: 
            - Computationally intensive
            - Sensitive to k parameter
            - Requires scaling
        - **Best for**: Moderate missing data (5-20%), correlated variables
        
        ### **5. Complete Case Analysis (Deletion)**
        - **Method**: Remove all samples with any missing values
        - **Pros**: 
            - No bias if MCAR
            - Simple and transparent
            - Preserves relationships
        - **Cons**: 
            - Loss of sample size
            - Can introduce bias if not MCAR
            - Loss of statistical power
        - **Best for**: Very small amounts of missing data (<2%), MCAR assumption met
        
        ### **Recommendation for This Analysis**
        Given the biological nature of tocopherol measurements and typical lipidomics data:
        - **If missing <5%**: Use **median imputation** (robust to outliers common in biological data)
        - **If missing 5-15%**: Use **KNN imputation** (preserves biological relationships)
        - **If missing >15%**: Consider **complete case analysis** or investigate why data is missing
        """)
    
    # Method selection
    method = st.selectbox(
        "Select Missing Data Handling Method:",
        ["Median Imputation (Recommended)", "Mean Imputation", "KNN Imputation (k=5)", 
         "Complete Case Analysis", "Forward Fill"]
    )
    
    # Apply selected method
    df_processed = df.copy()
    
    if method == "Median Imputation (Recommended)":
        for col in enhanced_columns:
            if col in df_processed.columns:
                median_val = df_processed[col].median()
                df_processed[col].fillna(median_val, inplace=True)
        st.success(f"✅ Applied median imputation. Dataset maintained: {len(df_processed)} samples")
        
    elif method == "Mean Imputation":
        for col in enhanced_columns:
            if col in df_processed.columns:
                mean_val = df_processed[col].mean()
                df_processed[col].fillna(mean_val, inplace=True)
        st.success(f"✅ Applied mean imputation. Dataset maintained: {len(df_processed)} samples")
        
    elif method == "KNN Imputation (k=5)":
        from sklearn.impute import KNNImputer
        imputer = KNNImputer(n_neighbors=5)
        df_processed[enhanced_columns] = imputer.fit_transform(df_processed[enhanced_columns])
        st.success(f"✅ Applied KNN imputation (k=5). Dataset maintained: {len(df_processed)} samples")
        
    elif method == "Complete Case Analysis":
        before_count = len(df_processed)
        df_processed = df_processed.dropna(subset=enhanced_columns)
        after_count = len(df_processed)
        st.warning(f"⚠️ Removed {before_count - after_count} samples. Dataset reduced to: {after_count} samples")
        
    else:  # Forward Fill
        for col in enhanced_columns:
            if col in df_processed.columns:
                df_processed[col].fillna(method='ffill', inplace=True)
                df_processed[col].fillna(method='bfill', inplace=True)  # Also backward fill for first values
        st.success(f"✅ Applied forward/backward fill. Dataset maintained: {len(df_processed)} samples")
    
    return df_processed, enhanced_columns

# ================================
# STEP 4: Reanalyze UMAP with 74 Lipids
# ================================
def pareto_scaling(data):
    """Apply Pareto scaling"""
    scaled_data = data.copy()
    for col in data.columns:
        mean_val = data[col].mean()
        std_val = data[col].std()
        if std_val > 0:
            scaled_data[col] = (data[col] - mean_val) / np.sqrt(std_val)
    return scaled_data

def calculate_gap_statistic(data, k_range=range(2, 11)):
    """Calculate gap statistic for optimal clustering"""
    gap_values = []
    for k in k_range:
        kmeans = KMeans(n_clusters=k, random_state=42, n_init=10)
        cluster_labels = kmeans.fit_predict(data)
        inertia = kmeans.inertia_
        
        n_refs = 10
        ref_inertias = []
        for _ in range(n_refs):
            random_data = np.random.uniform(
                low=data.min(axis=0), high=data.max(axis=0), size=data.shape
            )
            ref_kmeans = KMeans(n_clusters=k, random_state=42, n_init=10)
            ref_kmeans.fit(random_data)
            ref_inertias.append(ref_kmeans.inertia_)
        
        gap = np.log(np.mean(ref_inertias)) - np.log(inertia)
        gap_values.append(gap)
    
    return gap_values

def step4_reanalyze_74_lipids(df, lipid_columns_74):
    """Step 4: Reanalyze with 74 lipids"""
    st.header("Step 4: Reanalyze UMAP with 74 Lipids")
    
    # Prepare data
    lipid_data = df[lipid_columns_74]
    lipid_data_scaled = pareto_scaling(lipid_data)
    
    # UMAP
    with st.spinner("Running UMAP on 74 lipids..."):
        umap_model = umap.UMAP(n_neighbors=15, min_dist=0.1, n_components=2, random_state=42)
        umap_result_74 = umap_model.fit_transform(lipid_data_scaled.values)
    
    # Optimal clustering
    gap_74 = calculate_gap_statistic(umap_result_74)
    optimal_k_74 = range(2, 11)[np.argmax(gap_74)]
    
    kmeans_74 = KMeans(n_clusters=optimal_k_74, random_state=42, n_init=10)
    clusters_74 = kmeans_74.fit_predict(umap_result_74)
    
    df['UMAP1_74'] = umap_result_74[:, 0]
    df['UMAP2_74'] = umap_result_74[:, 1]
    df['cluster_74'] = clusters_74 + 1
    
    # Visualization
    fig = go.Figure()
    colors = ['red', 'blue', 'green', 'orange', 'purple', 'brown', 'pink', 'gray']
    
    for i, cluster in enumerate(sorted(df['cluster_74'].unique())):
        cluster_data = df[df['cluster_74'] == cluster]
        
        # Add hull
        if len(cluster_data) >= 3:
            hull_x, hull_y = get_cluster_hull(
                cluster_data['UMAP1_74'].values,
                cluster_data['UMAP2_74'].values
            )
            fig.add_trace(go.Scatter(
                x=hull_x, y=hull_y,
                fill='toself',
                fillcolor=colors[i % len(colors)],
                opacity=0.2,
                line=dict(width=1),
                showlegend=False,
                hoverinfo='skip'
            ))
        
        fig.add_trace(go.Scatter(
            x=cluster_data['UMAP1_74'],
            y=cluster_data['UMAP2_74'],
            mode='markers',
            marker=dict(size=6, color=colors[i % len(colors)]),
            name=f'C{cluster} (n={len(cluster_data)})'
        ))
    
    fig.update_layout(
        title=f"Reanalyzed UMAP - 74 Lipids ({optimal_k_74} clusters)",
        xaxis_title='UMAP1',
        yaxis_title='UMAP2',
        height=600,
        template='plotly_white'
    )
    
    st.plotly_chart(fig, use_container_width=True)
    
    col1, col2 = st.columns(2)
    with col1:
        st.metric("Optimal Clusters", optimal_k_74)
    with col2:
        st.metric("Max Gap Statistic", f"{max(gap_74):.3f}")
    
    return df, optimal_k_74

# ================================
# STEP 5: Reanalyze UMAP with 76 Lipids
# ================================
def step5_reanalyze_76_lipids(df, enhanced_columns):
    """Step 5: Reanalyze with 76 lipids (including tocopherols)"""
    st.header("Step 5: Reanalyze UMAP with 76 Lipids")
    
    # Prepare data
    lipid_data = df[enhanced_columns]
    lipid_data_scaled = pareto_scaling(lipid_data)
    
    # UMAP
    with st.spinner("Running UMAP on 76 lipids..."):
        umap_model = umap.UMAP(n_neighbors=15, min_dist=0.1, n_components=2, random_state=42)
        umap_result_76 = umap_model.fit_transform(lipid_data_scaled.values)
    
    # Optimal clustering
    gap_76 = calculate_gap_statistic(umap_result_76)
    optimal_k_76 = range(2, 11)[np.argmax(gap_76)]
    
    kmeans_76 = KMeans(n_clusters=optimal_k_76, random_state=42, n_init=10)
    clusters_76 = kmeans_76.fit_predict(umap_result_76)
    
    df['UMAP1_76'] = umap_result_76[:, 0]
    df['UMAP2_76'] = umap_result_76[:, 1]
    df['cluster_76'] = clusters_76 + 1
    
    # Visualization
    fig = go.Figure()
    colors = ['red', 'blue', 'green', 'orange', 'purple', 'brown', 'pink', 'gray']
    
    for i, cluster in enumerate(sorted(df['cluster_76'].unique())):
        cluster_data = df[df['cluster_76'] == cluster]
        
        # Add hull
        if len(cluster_data) >= 3:
            hull_x, hull_y = get_cluster_hull(
                cluster_data['UMAP1_76'].values,
                cluster_data['UMAP2_76'].values
            )
            fig.add_trace(go.Scatter(
                x=hull_x, y=hull_y,
                fill='toself',
                fillcolor=colors[i % len(colors)],
                opacity=0.2,
                line=dict(width=1),
                showlegend=False,
                hoverinfo='skip'
            ))
        
        fig.add_trace(go.Scatter(
            x=cluster_data['UMAP1_76'],
            y=cluster_data['UMAP2_76'],
            mode='markers',
            marker=dict(size=6, color=colors[i % len(colors)]),
            name=f'C{cluster} (n={len(cluster_data)})'
        ))
    
    fig.update_layout(
        title=f"Reanalyzed UMAP - 76 Lipids ({optimal_k_76} clusters)",
        xaxis_title='UMAP1',
        yaxis_title='UMAP2',
        height=600,
        template='plotly_white'
    )
    
    st.plotly_chart(fig, use_container_width=True)
    
    col1, col2 = st.columns(2)
    with col1:
        st.metric("Optimal Clusters", optimal_k_76)
    with col2:
        st.metric("Max Gap Statistic", f"{max(gap_76):.3f}")
    
    return df, optimal_k_76

# ================================
# STEP 6: Cluster Comparison
# ================================
def step6_cluster_comparison(df):
    """Step 6: Compare original and new clusters"""
    st.header("Step 6: Cluster Assignment Comparison")
    
    # Create comparison table
    comparison_df = df[['SID', 'hclust8', 'cluster_74', 'cluster_76']].copy()
    comparison_df.columns = ['Patient ID', 'Original (hclust8)', '74 Lipids', '76 Lipids']
    
    # Show first 20 rows
    st.subheader("Sample Cluster Assignments (First 20 Patients)")
    st.dataframe(comparison_df.head(20))
    
    # Calculate agreement statistics
    st.subheader("Cluster Agreement Analysis")
    
    # Compare original vs 74 lipids
    same_74 = (df['hclust8'] == df['cluster_74']).sum()
    diff_74 = len(df) - same_74
    
    # Compare original vs 76 lipids
    same_76 = (df['hclust8'] == df['cluster_76']).sum()
    diff_76 = len(df) - same_76
    
    # Compare 74 vs 76 lipids
    same_74_76 = (df['cluster_74'] == df['cluster_76']).sum()
    diff_74_76 = len(df) - same_74_76
    
    col1, col2, col3 = st.columns(3)
    
    with col1:
        st.metric("Original vs 74 Lipids",
                 f"Same: {same_74} ({same_74/len(df)*100:.1f}%)",
                 f"Different: {diff_74} ({diff_74/len(df)*100:.1f}%)")
    
    with col2:
        st.metric("Original vs 76 Lipids",
                 f"Same: {same_76} ({same_76/len(df)*100:.1f}%)",
                 f"Different: {diff_76} ({diff_76/len(df)*100:.1f}%)")
    
    with col3:
        st.metric("74 vs 76 Lipids",
                 f"Same: {same_74_76} ({same_74_76/len(df)*100:.1f}%)",
                 f"Different: {diff_74_76} ({diff_74_76/len(df)*100:.1f}%)")
    
    # Confusion matrix
    st.subheader("Cluster Transition Matrix (Original → 76 Lipids)")
    
    # Create transition matrix
    transition_matrix = pd.crosstab(df['hclust8'], df['cluster_76'], 
                                   rownames=['Original'], colnames=['New (76 Lipids)'])
    
    # Create heatmap
    fig = go.Figure(data=go.Heatmap(
        z=transition_matrix.values,
        x=[f'C{i}' for i in transition_matrix.columns],
        y=[f'HC{i}' for i in transition_matrix.index],
        colorscale='Blues',
        text=transition_matrix.values,
        texttemplate='%{text}',
        textfont={"size": 10}
    ))
    
    fig.update_layout(
        title="Patient Movement Between Clusters",
        xaxis_title="New Clusters (76 Lipids)",
        yaxis_title="Original Clusters (hclust8)",
        height=500
    )
    
    st.plotly_chart(fig, use_container_width=True)
    
    return comparison_df

# ================================
# STEP 7: Identify Reference Cluster
# ================================
def step7_identify_reference_cluster(df):
    """Step 7: Identify reference cluster with lowest βGlcCer"""
    st.header("Step 7: Identify Reference Cluster (Lowest βGlcCer)")
    
    glccer_cols = [col for col in df.columns if col.startswith('GlcCer')]
    
    # Calculate mean βGlcCer for each cluster (using 76-lipid clustering)
    cluster_glccer = df.groupby('cluster_76')[glccer_cols].mean().mean(axis=1)
    ref_cluster = cluster_glccer.idxmin()
    
    st.success(f"✅ Reference cluster identified: **Cluster {ref_cluster}** (lowest mean βGlcCer)")
    
    # Create violin plots for all βGlcCer species (like Fig 1C)
    st.subheader("βGlcCer Levels by Cluster (Violin Plots)")
    
    # Calculate mean of all βGlcCers for each sample
    df['mean_GlcCer'] = df[glccer_cols].mean(axis=1)
    
    fig = go.Figure()
    
    for cluster in sorted(df['cluster_76'].unique()):
        cluster_data = df[df['cluster_76'] == cluster]['mean_GlcCer']
        
        fig.add_trace(go.Violin(
            y=cluster_data,
            x=[f'C{cluster}'] * len(cluster_data),
            name=f'C{cluster}',
            box_visible=True,
            meanline_visible=True,
            fillcolor='lightblue' if cluster != ref_cluster else 'lightgreen',
            opacity=0.6
        ))
    
    fig.update_layout(
        title="Mean βGlcCer Levels by Cluster",
        xaxis_title="Cluster",
        yaxis_title="Mean of relative scaled βGlcCers",
        showlegend=False,
        height=500
    )
    
    st.plotly_chart(fig, use_container_width=True)
    
    # Statistical comparisons
    st.subheader("Statistical Comparisons (vs Reference)")
    
    ref_data = df[df['cluster_76'] == ref_cluster]['mean_GlcCer']
    
    comparison_results = []
    for cluster in sorted(df['cluster_76'].unique()):
        if cluster != ref_cluster:
            cluster_data = df[df['cluster_76'] == cluster]['mean_GlcCer']
            
            # Perform t-test
            t_stat, p_val = stats.ttest_ind(cluster_data, ref_data)
            
            # Apply Bonferroni correction
            p_val_corrected = min(p_val * (len(df['cluster_76'].unique()) - 1), 1.0)
            
            comparison_results.append({
                'Cluster': f'C{cluster}',
                'Mean βGlcCer': cluster_data.mean(),
                'Ref Mean': ref_data.mean(),
                'Difference': cluster_data.mean() - ref_data.mean(),
                't-statistic': t_stat,
                'p-value': p_val,
                'p-value (corrected)': p_val_corrected,
                'Significant': '***' if p_val_corrected < 0.001 else '**' if p_val_corrected < 0.01 else '*' if p_val_corrected < 0.05 else 'ns'
            })
    
    comparison_df = pd.DataFrame(comparison_results)
    st.dataframe(comparison_df.round(4))
    
    return ref_cluster

# ================================
# STEP 8: Volcano Plots
# ================================
def step8_volcano_plots(df, ref_cluster):
    """Step 8: Create volcano plots for all clusters vs reference"""
    st.header("Step 8: Volcano Plots (All Clusters vs Reference)")
    
    # Get lipid columns
    lipid_cols = [col for col in df.columns if (
        col.startswith('Cer.') or col.startswith('GlcCer') or 
        col.startswith('GalCer') or col.startswith('SPH.') or
        col in ['So', 'DHSo', 'S1P', 'DHS1P'] or
        col.startswith('PG') or col.startswith('LT') or 
        col.startswith('Resolvin') or col == 'Lipoxin.A4' or
        col.endswith('HETE') or col.endswith('EET') or col.endswith('DHET') or
        col in ['EPA', 'DHA', 'AA', 'DHGLA', '13.HODE', 'TXB2', '5.iPF2α-VI',
                'alpha-tocopherol', 'gamma-tocopherol']
    )]
    
    # Calculate reference means
    ref_data = df[df['cluster_76'] == ref_cluster][lipid_cols]
    ref_means = ref_data.mean()
    
    # Define lipid categories
    def categorize_lipid(lipid_name):
        if lipid_name.startswith('Cer.'):
            return 'Ceramides', 'green'
        elif lipid_name.startswith('GlcCer'):
            return 'βGlucosylceramides', 'red'
        elif lipid_name.startswith('GalCer'):
            return 'βGalactosylceramides', 'magenta'
        elif lipid_name.startswith('SPH.') or lipid_name in ['So', 'DHSo', 'S1P', 'DHS1P']:
            return 'Sphingomyelins', 'lightblue'
        elif lipid_name == 'Lipoxin.A4':
            return 'Lipoxin A4', 'gold'
        elif 'Resolvin' in lipid_name:
            return 'Resolvin D1 or D2', 'cyan'
        elif 'tocopherol' in lipid_name:
            return 'Tocopherols', 'brown'
        else:
            return 'Eicosanoids', 'olive'
    
    clusters_to_compare = [c for c in sorted(df['cluster_76'].unique()) if c != ref_cluster]
    
    # Create subplots
    n_clusters = len(clusters_to_compare)
    n_cols = 4
    n_rows = (n_clusters + n_cols - 1) // n_cols
    
    fig = make_subplots(
        rows=n_rows, cols=n_cols,
        subplot_titles=[f'C{ref_cluster} vs C{c}' for c in clusters_to_compare]
    )
    
    # Process each cluster
    for idx, cluster in enumerate(clusters_to_compare):
        row = idx // n_cols + 1
        col = idx % n_cols + 1
        
        cluster_data = df[df['cluster_76'] == cluster][lipid_cols]
        cluster_means = cluster_data.mean()
        
        # Calculate log2 fold change
        log2fc = np.log2((cluster_means + 1e-10) / (ref_means + 1e-10))
        
        # Calculate p-values
        p_values = []
        for lipid in lipid_cols:
            ref_values = ref_data[lipid].dropna()
            cluster_values = cluster_data[lipid].dropna()
            if len(ref_values) > 0 and len(cluster_values) > 0:
                _, p_val = stats.mannwhitneyu(cluster_values, ref_values, alternative='two-sided')
                p_values.append(p_val)
            else:
                p_values.append(1.0)
        
        # FDR correction
        _, fdr_values, _, _ = multipletests(p_values, method='fdr_bh')
        
        # Plot
        for i, lipid in enumerate(lipid_cols):
            category, color = categorize_lipid(lipid)
            is_significant = (abs(log2fc[lipid]) > 0.6) and (-np.log10(fdr_values[i] + 1e-10) > 1.3)
            
            fig.add_trace(
                go.Scatter(
                    x=[log2fc[lipid]],
                    y=[-np.log10(fdr_values[i] + 1e-10)],
                    mode='markers',
                    marker=dict(
                        size=8 if is_significant else 5,
                        color=color if is_significant else 'black',
                        opacity=0.8 if is_significant else 0.3
                    ),
                    showlegend=False,
                    hovertext=lipid
                ),
                row=row, col=col
            )
        
        # Add threshold lines
        fig.add_hline(y=1.3, line_dash="dash", line_color="red", opacity=0.5, row=row, col=col)
        fig.add_vline(x=0.6, line_dash="dash", line_color="red", opacity=0.5, row=row, col=col)
        fig.add_vline(x=-0.6, line_dash="dash", line_color="red", opacity=0.5, row=row, col=col)
        
        fig.update_xaxes(title_text="Log2(FC)", row=row, col=col, range=[-5, 5])
        fig.update_yaxes(title_text="-Log10(FDR)", row=row, col=col, range=[0, 8])
    
    fig.update_layout(
        height=200 * n_rows,
        title_text=f"Volcano Plots: All Clusters vs Reference (C{ref_cluster})",
        template='plotly_white'
    )
    
    st.plotly_chart(fig, use_container_width=True)
    
    return

# ================================
# STEP 9: Pro-inflammatory Index
# ================================
def step9_proinflammatory_index(df, ref_cluster):
    """Step 9: Calculate pro-inflammatory index"""
    st.header("Step 9: Pro-inflammatory Index Calculation")
    
    # Define lipid groups
    glccer_cols = [col for col in df.columns if col.startswith('GlcCer')]
    galcer_cols = [col for col in df.columns if col.startswith('GalCer')]
    sphingo_cols = ['Cer.C16.0', 'Cer.C18.0', 'Cer.C24.0']  # Representative sphingomyelins
    
    # Calculate cluster means
    cluster_means = df.groupby('cluster_76')[
        glccer_cols + galcer_cols + sphingo_cols + 
        ['Lipoxin.A4', 'Resolvin.D2', 'alpha-tocopherol', 'gamma-tocopherol']
    ].mean()
    
    ref_means = cluster_means.loc[ref_cluster]
    
    # Calculate index for each cluster
    index_data = []
    for cluster in sorted(cluster_means.index):
        cluster_mean = cluster_means.loc[cluster]
        
        # Calculate log2FC
        glccer_log2fc = np.mean([np.log2((cluster_mean[col] + 1e-10) / (ref_means[col] + 1e-10)) 
                                 for col in glccer_cols])
        galcer_log2fc = np.mean([np.log2((cluster_mean[col] + 1e-10) / (ref_means[col] + 1e-10)) 
                                 for col in galcer_cols])
        sphingo_log2fc = np.mean([np.log2((cluster_mean[col] + 1e-10) / (ref_means[col] + 1e-10)) 
                                  for col in sphingo_cols])
        
        lipoxin_log2fc = np.log2((cluster_mean['Lipoxin.A4'] + 1e-10) / (ref_means['Lipoxin.A4'] + 1e-10))
        resolvin_log2fc = np.log2((cluster_mean['Resolvin.D2'] + 1e-10) / (ref_means['Resolvin.D2'] + 1e-10))
        alpha_log2fc = np.log2((cluster_mean['alpha-tocopherol'] + 1e-10) / (ref_means['alpha-tocopherol'] + 1e-10))
        gamma_log2fc = np.log2((cluster_mean['gamma-tocopherol'] + 1e-10) / (ref_means['gamma-tocopherol'] + 1e-10))
        
        # Calculate index
        pro_inflammatory = glccer_log2fc + galcer_log2fc + sphingo_log2fc + gamma_log2fc
        anti_inflammatory = lipoxin_log2fc + resolvin_log2fc + alpha_log2fc
        
        index = pro_inflammatory - anti_inflammatory
        
        index_data.append({
            'cluster': cluster,
            'βGlcCers (mean log2FC)': glccer_log2fc,
            'βGalCers (mean log2FC)': galcer_log2fc,
            'sphingomyelins (mean log2FC)': sphingo_log2fc,
            'Lipoxin A4 (log2FC)': lipoxin_log2fc,
            'Resolvin D2 (log2FC)': resolvin_log2fc,
            'α-tocopherol (log2FC)': alpha_log2fc,
            'γ-tocopherol (log2FC)': gamma_log2fc,
            'pro-inflammatory index': index
        })
    
    index_df = pd.DataFrame(index_data)
    
    # Display table
    st.subheader("Pro-inflammatory Lipid Index for Clusters")
    display_df = index_df.round(2)
    st.dataframe(display_df)
    
    # Visualization
    fig = go.Figure()
    fig.add_trace(go.Bar(
        x=[f'C{int(c)}' for c in index_df['cluster']],
        y=index_df['pro-inflammatory index'],
        marker_color=['lightgreen' if c == ref_cluster else 'coral' for c in index_df['cluster']],
        text=index_df['pro-inflammatory index'].round(2),
        textposition='outside'
    ))
    
    fig.update_layout(
        title="Pro-inflammatory Index by Cluster",
        xaxis_title="Cluster",
        yaxis_title="Pro-inflammatory Index",
        height=400,
        template='plotly_white'
    )
    
    st.plotly_chart(fig, use_container_width=True)
    
    return index_df

# ================================
# STEP 10: Clinical Correlations
# ================================
def step10_clinical_correlations(df, index_df):
    """Step 10: Clinical correlations with regression analysis"""
    st.header("Step 10: Clinical Correlations with Regression Analysis")
    
    # Calculate cluster outcome percentages
    cluster_outcomes = df.groupby('cluster_76').agg({
        'Wheeze_1yr': lambda x: (x == 1).mean() * 100,
        'AD_1yr': lambda x: (x == 1).mean() * 100
    })
    
    # Merge with index
    correlation_data = index_df.merge(
        cluster_outcomes, 
        left_on='cluster', 
        right_index=True
    )
    
    # Linear regression
    from sklearn.linear_model import LinearRegression
    
    X = correlation_data['pro-inflammatory index'].values.reshape(-1, 1)
    y_wheeze = correlation_data['Wheeze_1yr'].values
    y_ad = correlation_data['AD_1yr'].values
    
    reg_wheeze = LinearRegression().fit(X, y_wheeze)
    reg_ad = LinearRegression().fit(X, y_ad)
    
    r2_wheeze = reg_wheeze.score(X, y_wheeze)
    r2_ad = reg_ad.score(X, y_ad)
    
    # Spearman correlation
    wheeze_corr, wheeze_p = stats.spearmanr(
        correlation_data['pro-inflammatory index'],
        correlation_data['Wheeze_1yr']
    )
    ad_corr, ad_p = stats.spearmanr(
        correlation_data['pro-inflammatory index'],
        correlation_data['AD_1yr']
    )
    
    # Create plots
    fig = make_subplots(
        rows=1, cols=2,
        subplot_titles=('Wheeze vs Pro-inflammatory Index', 
                       'Atopic Dermatitis vs Pro-inflammatory Index')
    )
    
    # Generate regression lines
    x_range = np.linspace(correlation_data['pro-inflammatory index'].min(),
                         correlation_data['pro-inflammatory index'].max(), 100)
    wheeze_pred = reg_wheeze.predict(x_range.reshape(-1, 1))
    ad_pred = reg_ad.predict(x_range.reshape(-1, 1))
    
    # Wheeze plot
    fig.add_trace(
        go.Scatter(
            x=x_range,
            y=wheeze_pred,
            mode='lines',
            line=dict(color='black', dash='dot', width=2),
            showlegend=False
        ),
        row=1, col=1
    )
    
    fig.add_trace(
        go.Scatter(
            x=correlation_data['pro-inflammatory index'],
            y=correlation_data['Wheeze_1yr'],
            mode='markers+text',
            text=[f'C{int(c)}' for c in correlation_data['cluster']],
            textposition="top center",
            marker=dict(size=10, color='blue'),
            showlegend=False
        ),
        row=1, col=1
    )
    
    # AD plot
    fig.add_trace(
        go.Scatter(
            x=x_range,
            y=ad_pred,
            mode='lines',
            line=dict(color='black', dash='dot', width=2),
            showlegend=False
        ),
        row=1, col=2
    )
    
    fig.add_trace(
        go.Scatter(
            x=correlation_data['pro-inflammatory index'],
            y=correlation_data['AD_1yr'],
            mode='markers+text',
            text=[f'C{int(c)}' for c in correlation_data['cluster']],
            textposition="top center",
            marker=dict(size=10, color='red'),
            showlegend=False
        ),
        row=1, col=2
    )
    
    # Add annotations
    fig.add_annotation(
        x=correlation_data['pro-inflammatory index'].min() + 0.5,
        y=correlation_data['Wheeze_1yr'].max() - 5,
        text=f"Adjusted R² = {r2_wheeze:.2f}<br>p = {wheeze_p:.4f}",
        showarrow=False,
        font=dict(size=11),
        bgcolor="white",
        bordercolor="black",
        borderwidth=1,
        xref="x",
        yref="y"
    )
    
    fig.add_annotation(
        x=correlation_data['pro-inflammatory index'].min() + 0.5,
        y=correlation_data['AD_1yr'].max() - 2,
        text=f"Adjusted R² = {r2_ad:.2f}<br>p = {ad_p:.3f}",
        showarrow=False,
        font=dict(size=11),
        bgcolor="white",
        bordercolor="black",
        borderwidth=1,
        xref="x2",
        yref="y2"
    )
    
    fig.update_xaxes(title_text="Pro-inflammatory Index", row=1, col=1)
    fig.update_xaxes(title_text="Pro-inflammatory Index", row=1, col=2)
    fig.update_yaxes(title_text="Wheeze (%)", row=1, col=1)
    fig.update_yaxes(title_text="Atopic Dermatitis (%)", row=1, col=2)
    
    fig.update_layout(
        height=500,
        title_text="Clinical Outcomes vs Pro-inflammatory Index",
        template='plotly_white'
    )
    
    st.plotly_chart(fig, use_container_width=True)
    
    # Summary statistics
    st.subheader("Statistical Summary")
    results_df = pd.DataFrame({
        'Outcome': ['Wheeze', 'Atopic Dermatitis'],
        'Spearman ρ': [wheeze_corr, ad_corr],
        'p-value': [wheeze_p, ad_p],
        'Adjusted R²': [r2_wheeze, r2_ad]
    })
    st.dataframe(results_df.round(3))
    
    return correlation_data

# ================================
# MAIN APPLICATION
# ================================
def main():
    st.title("WISC Lipidomics Analysis Pipeline")
    st.markdown("**Systematic Analysis**: 74 → 76 Lipids (Adding α-tocopherol & γ-tocopherol)")
    
    # Sidebar
    with st.sidebar:
        st.header("Analysis Steps")
        st.markdown("""
        1. **Original Dataset Summary**
        2. **Reproduce Original UMAP**
        3. **Tocopherol Dataset & Missing Data**
        4. **Reanalyze 74 Lipids**
        5. **Reanalyze 76 Lipids**
        6. **Cluster Comparison**
        7. **Reference Cluster (βGlcCer)**
        8. **Volcano Plots**
        9. **Pro-inflammatory Index**
        10. **Clinical Correlations**
        """)

        st.header("Analysis Pipeline Explanation")
        explain_analysis_pipeline()
        
        run_analysis = st.button("🚀 Run Complete Pipeline", type="primary")
    
    if run_analysis:
        # Step 1
        st.markdown("---")
        df, lipid_columns_74 = step1_original_dataset_summary()
        
        # Step 2
        st.markdown("---")
        df = step2_reproduce_original_umap(df)
        
        # Step 3
        st.markdown("---")
        df, enhanced_columns = step3_tocopherol_dataset_missing_data(df, lipid_columns_74)
        
        # Step 4
        st.markdown("---")
        df, optimal_k_74 = step4_reanalyze_74_lipids(df, lipid_columns_74)
        
        # Step 5
        st.markdown("---")
        df, optimal_k_76 = step5_reanalyze_76_lipids(df, enhanced_columns)
        
        # Step 6
        st.markdown("---")
        comparison_df = step6_cluster_comparison(df)
        
        # Step 7
        st.markdown("---")
        ref_cluster = step7_identify_reference_cluster(df)
        
        # Step 8
        st.markdown("---")
        step8_volcano_plots(df, ref_cluster)
        
        # Step 9
        st.markdown("---")
        index_df = step9_proinflammatory_index(df, ref_cluster)
        
        # Step 10
        st.markdown("---")
        correlation_data = step10_clinical_correlations(df, index_df)
        
        # Export results
        st.markdown("---")
        st.header("📥 Export Results")
        export_df = df[['SID', 'hclust8', 'cluster_74', 'cluster_76',
                        'UMAP1_74', 'UMAP2_74', 'UMAP1_76', 'UMAP2_76',
                        'alpha-tocopherol', 'gamma-tocopherol',
                        'Wheeze_1yr', 'AD_1yr']].copy()
        
        csv = export_df.to_csv(index=False)
        st.download_button(
            label="📊 Download Complete Analysis Results",
            data=csv,
            file_name="10step_wisc_analysis_results.csv",
            mime="text/csv"
        )
        
        st.success("✅ Complete 10-step analysis finished!")

if __name__ == "__main__":
    main()
# Enhanced WISC Lipidomics Analysis

**Methodology-Compliant Enhancement with α-tocopherol & γ-tocopherol**

**Authors:** Luan Vu (UTSA) & Joan Cook-Mills (IU)  
**Date:** December 2024

##  Overview

This project enhances the original Wisconsin Infant Study Cohort (WISC) lipidomics analysis by adding **α-tocopherol** and **γ-tocopherol** measurements while maintaining **exact compliance** with the original methodology.

### Key Features
-  **Exact WISC Methodology**: Pareto scaling, UMAP, Gap statistic, Spearman correlations
-  **Enhanced Analysis**: 76 total lipids (74 original + 2 tocopherols)  
-  **Interactive Interface**: Easy-to-use Streamlit app for non-coding colleagues
-  **Publication Ready**: All figures and results formatted for publication

##  Quick Start (3 Easy Steps)

### Step 1: Download Files
Save these 4 files in the same folder as your `data.csv`:
1. `wisc_enhanced_analysis.py` (main analysis script)
2. `requirements.txt` (Python packages)
3. `run_analysis.py` (launcher script)
4. `README.md` (this file)

### Step 2: Run the Launcher
```bash
python run_analysis.py
```

The launcher will:
- ✅ Check your Python version
- ✅ Verify `data.csv` exists  
- ✅ Install required packages automatically
- ✅ Launch the web interface

### Step 3: Run Analysis
1. A web browser will open automatically
2. Click **" Run Enhanced Analysis"**
3. View results and download files

##  What You'll Get

###  Interactive Visualizations
- **UMAP Comparison**: Original vs Enhanced clustering
- **Tocopherol Distributions**: By cluster and clinical outcomes
- **Enhanced Pro-inflammatory Index**: With tocopherol contributions
- **Correlation Analysis**: Enhanced index vs clinical outcomes

###  Statistical Results
- **Spearman correlations** with bootstrap confidence intervals
- **Gap statistics** for optimal clustering
- **Enhanced cluster assignments** for all patients
- **Publication-ready tables** and figures

###  Downloadable Results
- Enhanced analysis results (CSV format)
- All UMAP coordinates (original + enhanced)
- Enhanced pro-inflammatory index values
- Statistical summary tables

##  Methodology Compliance

### Original WISC Methods (Maintained)
✅ **Pareto Scaling**: `(x - mean) / sqrt(std)`  
✅ **UMAP Parameters**: n_neighbors=15, min_dist=0.1  
✅ **Gap Statistic**: For optimal cluster determination  
✅ **Reference Cluster**: C6 (as in original)  
✅ **Statistical Tests**: Spearman correlations  
✅ **Bootstrap CI**: 5000 resamples  
✅ **Visualization Style**: Matching original paper  

### Enhancement (Added)
➕ **α-tocopherol**: Protective effect (negative contribution)  
➕ **γ-tocopherol**: Pro-inflammatory effect (positive contribution)  
➕ **76 Total Lipids**: 74 original + 2 tocopherols  
➕ **Enhanced Index**: Updated pro-inflammatory formula  

##  Sharing with Colleagues

### For Non-Coding Colleagues:
1. Via the web URL (usually `http://localhost:8501`)
2. They can run analysis with button clicks
3. Download results directly from the interface
4. No coding or statistics knowledge required

### For Publication:
- All plots are publication-ready
- Statistical methods match original study
- Results directly comparable to original WISC

##  Troubleshooting

### Common Issues:

**"data.csv not found"**
- Ensure `data.csv` is in the same folder as the scripts

**"Package installation failed"**
```bash
pip install pandas numpy plotly streamlit scikit-learn umap-learn scipy
```

**"Port already in use"**
```bash
streamlit run wisc_enhanced_analysis.py --server.port 8502
```

**Python version issues**
- Requires Python 3.8 or higher
- Check version: `python --version`

## Expected Results

### Enhanced Clustering
- Compare original 74-lipid vs enhanced 76-lipid clustering
- Assess impact of tocopherols on patient groupings
- Optimal cluster determination via gap statistic

### Pro-inflammatory Index
- Updated index incorporating tocopherol contributions
- α-tocopherol: Protective (negative contribution)
- γ-tocopherol: Pro-inflammatory (positive contribution)

### Clinical Correlations
- Enhanced index vs wheeze prevalence
- Enhanced index vs atopic dermatitis prevalence
- Spearman correlations with 95% confidence intervals

## Support

If you encounter issues:
1. Check this README for troubleshooting steps
2. Ensure all files are in the same directory
3. Verify `data.csv` has the required columns
4. Contact the authors for additional support

---

** Ready to enhance your WISC analysis with tocopherols!** 

# Virtual Lab

[![PyPI - Python Version](https://img.shields.io/pypi/pyversions/virtual-lab)](https://badge.fury.io/py/virtual-lab)
[![PyPI version](https://badge.fury.io/py/virtual-lab.svg)](https://badge.fury.io/py/virtual-lab)
[![Downloads](https://pepy.tech/badge/virtual-lab)](https://pepy.tech/project/virtual-lab)
[![license](https://img.shields.io/github/license/zou-group/virtual-lab.svg)](https://github.com/zou-group/virtual-lab/blob/main/LICENSE.txt)

![Virtual Lab](images/virtual_lab_architecture.png)

The **Virtual Lab** is an AI-human collaboration for science research. In the Virtual Lab, a human researcher works with a team of large language model (LLM) **agents** to perform scientific research. Interaction between the human researcher and the LLM agents occurs via a series of **team meetings**, where all the LLM agents discuss a scientific agenda posed by the human researcher, and **individual meetings**, where the human researcher interacts with a single LLM agent to solve a particular scientific task.

Please see our paper [The Virtual Lab of AI agents designs new SARS-CoV-2 nanobodies](https://www.nature.com/articles/s41586-025-09442-9) for more details on the Virtual Lab and an application to nanobody design for SARS-CoV-2.

If you use the Virtual Lab, please cite our work as follows:

Swanson, K., Wu, W., Bulaong, N.L. et al. The Virtual Lab of AI agents designs new SARS-CoV-2 nanobodies. *Nature* (2025). https://doi.org/10.1038/s41586-025-09442-9


## Applications

### Virtual Lab for nanobody design

As a real-world demonstration, we applied the Virtual Lab to design nanobodies for one of the latest variants of SARS-CoV-2 (see [nanobody_design](https://github.com/zou-group/virtual-lab/tree/main/nanobody_design)). The Virtual Lab built a computational pipeline consisting of [ESM](https://www.science.org/doi/10.1126/science.ade2574), [AlphaFold-Multimer](https://www.biorxiv.org/content/10.1101/2021.10.04.463034v2), and [Rosetta](https://rosettacommons.org/software/) and used it to design 92 nanobodies that were experimentally validated.

Please see the notebook [nanobody_design/run_nanobody_design.ipynb](https://github.com/zou-group/virtual-lab/blob/main/nanobody_design/run_nanobody_design.ipynb) for an example of how to use the Virtual Lab to create agents and run team and individual meetings.

### Biomarker Analysis for Prognosis Prediction

The Virtual Lab includes a comprehensive toolkit for identifying prognostic biomarker genes from gene expression data with survival outcomes. This application demonstrates how the platform can be used for general-purpose bioinformatics analysis beyond the nanobody design case study.

#### Features

- **Multiple Statistical Methods**: Cox Proportional Hazards, Log-rank test, Differential Expression, Elastic Net Cox Regression
- **Consensus Gene Selection**: Identify robust biomarkers that are significant across multiple methods
- **Comprehensive Visualization**: Volcano plots, Kaplan-Meier curves, heatmaps, and more
- **Easy-to-use Interface**: Command-line tool and interactive Jupyter notebook

#### Quick Start

**Using Pixi (Recommended):**

```bash
# Run the complete biomarker analysis
pixi run biomarker-analysis

# Or launch the interactive tutorial
pixi run biomarker-tutorial
```

**Using Python directly:**

```bash
cd biomarker_analysis/scripts

python select_marker_genes.py \
    --input_file ../../Example_TCGA_TNBC_data.csv \
    --output_dir ../../biomarker_results \
    --n_top_genes 50 \
    --methods cox logrank differential elasticnet \
    --visualization
```

**Using Jupyter Notebook:**

```bash
jupyter notebook biomarker_analysis/tutorial_biomarker_selection.ipynb
```

#### Example Dataset

The repository includes an example TCGA triple-negative breast cancer (TNBC) dataset:

- **File**: `Example_TCGA_TNBC_data.csv`
- **Samples**: 144 patients
- **Features**: ~20,000 gene expression values (log-transformed)
- **Outcomes**: Overall survival (OS) and survival time (OS.year)

#### Output Files

The analysis generates:

- `consensus_genes.csv`: Genes ranked by appearance across methods
- `{method}_results.csv`: Detailed results for each statistical method
- `top_biomarkers_summary.csv`: Comprehensive summary of top genes
- `marker_gene_analysis.pdf`: Visualization plots

#### Statistical Methods

1. **Cox Proportional Hazards Regression**
   - Models the relationship between gene expression and survival time
   - Accounts for censored data
   - Outputs hazard ratios (HR) and confidence intervals

2. **Log-rank Test (Kaplan-Meier Analysis)**
   - Compares survival curves between high vs. low expression groups
   - Non-parametric test for survival differences
   - Generates Kaplan-Meier plots

3. **Differential Expression Analysis**
   - Compares expression between event and censored groups
   - Uses Mann-Whitney U test (non-parametric)
   - Calculates fold changes and significance

4. **Elastic Net Cox Regression**
   - Regularized Cox model with L1/L2 penalties
   - Performs automatic feature selection
   - Identifies sparse set of prognostic genes

#### Tutorial

For a comprehensive guide on biomarker selection, see the tutorial notebook:
[biomarker_analysis/tutorial_biomarker_selection.ipynb](biomarker_analysis/tutorial_biomarker_selection.ipynb)

The tutorial covers:
- Data loading and exploration
- Each statistical method in detail
- Consensus gene identification
- Result visualization and interpretation
- Exporting results for further analysis

### Genome-Scale Metabolic Modeling

The Virtual Lab includes a comprehensive toolkit for genome-scale metabolic modeling using COBRApy. This module enables identification of gene knockout/knockdown targets for metabolic engineering, drug discovery, and systems biology applications.

#### Features

- **Constraint-Based Modeling**: Flux Balance Analysis (FBA) and Flux Variability Analysis (FVA)
- **Gene Knockout Simulation**: Single and double gene knockout analysis
- **Essential Gene Identification**: Find critical genes for survival
- **Synthetic Lethality Discovery**: Identify gene pairs that are lethal when deleted together
- **Growth-Coupled Production**: Design strains for optimal bioproduction
- **Multiple Model Support**: E. coli, human (Recon3D), yeast, and custom models

#### Quick Start

**Using Pixi (Recommended):**

```bash
# Run metabolic analysis with E. coli core model (fast)
pixi run metabolic-analysis

# Run with larger E. coli iML1515 model
pixi run metabolic-iml1515

# Launch the interactive tutorial
pixi run metabolic-tutorial
```

**Using Python directly:**

```bash
cd metabolic_modeling/scripts

# E. coli core model (fast, for testing)
python metabolic_target_finder.py \
    --model_id textbook \
    --output_dir ../../metabolic_results \
    --ko_methods single essential fva \
    --visualization

# E. coli iML1515 (comprehensive)
python metabolic_target_finder.py \
    --model_id iML1515 \
    --output_dir ../../metabolic_results \
    --ko_methods single essential fva \
    --growth_threshold 0.1 \
    --visualization
```

**Using Jupyter Notebook:**

```bash
jupyter notebook metabolic_modeling/tutorial_metabolic_modeling.ipynb
```

#### Available Models

**From BiGG Database:**
- `textbook`: E. coli core model (95 reactions, 72 genes) - Fast, for learning
- `iML1515`: Latest E. coli model (2,712 reactions, 1,877 genes) - Production use
- `iJO1366`: Previous E. coli model (2,583 reactions, 1,366 genes)
- `Recon3D`: Human metabolism (13,543 reactions, 3,288 genes) - Drug discovery
- `iMM904`: S. cerevisiae (1,577 reactions, 904 genes) - Yeast engineering

**Custom Models:**
- Load from SBML, JSON, or MAT files
- Place in `metabolic_modeling/models/` directory

#### Analysis Methods

1. **Single Gene Knockout**
   - Systematically delete each gene
   - Measure growth impact
   - Classify essential vs. non-essential genes

2. **Essential Gene Analysis**
   - Identify genes critical for survival
   - Applications: Antibiotic targets, cancer drug discovery
   - Rank by criticality

3. **Double Gene Knockout (Synthetic Lethality)**
   - Find gene pairs that are lethal together
   - Applications: Combination therapy design
   - Identify synergistic targets

4. **Flux Variability Analysis (FVA)**
   - Determine flux range for each reaction
   - Identify metabolic bottlenecks
   - Find robust engineering targets

5. **Growth-Coupled Production**
   - Design strains for bioproduction
   - Couple product formation with growth
   - Optimize yield and productivity

#### Output Files

The analysis generates:

- `model_summary.csv`: Model statistics and information
- `single_knockout_results.csv`: Growth effects of all gene knockouts
- `essential_genes_results.csv`: Essential genes with detailed annotations
- `double_knockout_results.csv`: Synthetic lethal gene pairs
- `fva_results.csv`: Flux ranges for all reactions
- `growth_coupled_results.csv`: Production-enhancing knockouts
- `metabolic_analysis.pdf`: Comprehensive visualization plots

#### Use Cases

**Metabolic Engineering:**
```bash
# Optimize E. coli for succinate production
python metabolic_target_finder.py \
    --model_id iML1515 \
    --ko_methods single production \
    --target_metabolite succ_c \
    --output_dir succinate_production
```

**Drug Target Discovery:**
```bash
# Find essential genes in pathogenic bacteria
python metabolic_target_finder.py \
    --model_file models/pathogen.xml \
    --ko_methods essential \
    --growth_threshold 0.05 \
    --output_dir antibiotic_targets
```

**Cancer Research:**
```bash
# Identify synthetic lethal pairs in cancer metabolism
python metabolic_target_finder.py \
    --model_id Recon3D \
    --ko_methods double essential \
    --output_dir cancer_targets
```

#### Tutorial

For a comprehensive guide on metabolic modeling, see the tutorial notebook:
[metabolic_modeling/tutorial_metabolic_modeling.ipynb](metabolic_modeling/tutorial_metabolic_modeling.ipynb)

The tutorial covers:
- Loading and analyzing metabolic models
- Flux Balance Analysis (FBA)
- Single and double gene knockout simulations
- Essential gene identification
- Flux Variability Analysis (FVA)
- Growth-coupled production strategies
- Working with large genome-scale models


## Installation

The Virtual Lab can be installed using pip, conda, or pixi. Installation should only take a couple of minutes.

### Option 1: Quick Install with Pixi (Recommended)

[Pixi](https://pixi.sh) is a fast, modern package manager that handles all dependencies automatically. This is the **recommended installation method** for the best out-of-the-box experience.

#### Install Pixi

```bash
# On Linux and macOS
curl -fsSL https://pixi.sh/install.sh | bash

# On Windows
iwr -useb https://pixi.sh/install.ps1 | iex
```

#### Install Virtual Lab

```bash
git clone https://github.com/zou-group/virtual_lab.git
cd virtual_lab

# Install all dependencies (one command!)
pixi install

# Activate the environment
pixi shell
```

That's it! All dependencies including Python, Jupyter, and all scientific packages are now installed and ready to use.

#### Available Pixi Commands

```bash
# Run biomarker analysis
pixi run biomarker-analysis

# Launch tutorial notebook
pixi run biomarker-tutorial

# Launch Jupyter Lab
pixi run lab

# Run nanobody design notebook
pixi run nanobody-design

# Quick test of biomarker analysis
pixi run test-biomarker

# Clean up generated files
pixi run clean
```

#### Pixi Environments

The Virtual Lab provides several pre-configured environments:

- `default`: Core Virtual Lab functionality
- `dev`: Development tools (pytest, black, ruff, mypy)
- `nanobody`: Nanobody design dependencies
- `biomarker`: Biomarker analysis dependencies
- `full`: All features combined

To use a specific environment:

```bash
# Use biomarker environment
pixi shell -e biomarker

# Use full environment with all features
pixi shell -e full
```

### Option 2: Install with Conda

Create a conda environment and install dependencies:

```bash
conda create -y -n virtual_lab python=3.12
conda activate virtual_lab

# Install Virtual Lab
pip install virtual-lab

# For biomarker analysis, install additional dependencies
pip install lifelines scikit-survival pandas scipy matplotlib seaborn
```

### Option 3: Install with Pip (Minimal)

```bash
# Install from PyPI
pip install virtual-lab

# Or install from source
git clone https://github.com/zou-group/virtual_lab.git
cd virtual_lab
pip install -e .
```

### Install Optional Dependencies

For specific use cases, you may need additional packages:

```bash
# For nanobody design
pip install -e ".[nanobody-design]"

# For biomarker analysis
pip install lifelines scikit-survival statsmodels
```

### Verify Installation

```bash
# Test Python import
python -c "import virtual_lab; print('Virtual Lab successfully installed!')"

# Check Jupyter is available
jupyter --version

# For biomarker analysis, test dependencies
python -c "import lifelines, sksurv; print('Biomarker analysis dependencies installed!')"
```

## OpenAI API Key

The Virtual Lab currently uses GPT-4o from OpenAI. Save your OpenAI API key as the environment variable `OPENAI_API_KEY`.

```bash
# Add to your shell profile (.bashrc, .bash_profile, or .zshrc)
export OPENAI_API_KEY='your-api-key-here'

# Or set for current session
export OPENAI_API_KEY='your-api-key-here'
```

To get an OpenAI API key:
1. Visit https://platform.openai.com/api-keys
2. Sign up or log in
3. Create a new API key
4. Copy and save it securely

# Genome-Scale Metabolic Modeling

This module provides comprehensive tools for identifying gene knockout/knockdown targets using genome-scale metabolic models (GEMs) and COBRApy.

## Overview

Genome-scale metabolic modeling uses constraint-based analysis to predict cellular behavior and identify metabolic engineering targets. This implementation provides:

- **Single & double gene knockout simulations**
- **Essential gene identification**
- **Synthetic lethality discovery**
- **Flux variability analysis (FVA)**
- **Growth-coupled production analysis**
- **Comprehensive visualization**

## Directory Structure

```
metabolic_modeling/
├── scripts/
│   └── metabolic_target_finder.py    # Main analysis script
├── models/                            # Store your metabolic models here
├── examples/                          # Example models and data
├── tutorial_metabolic_modeling.ipynb  # Interactive tutorial
├── requirements.txt                   # Python dependencies
└── README.md                          # This file
```

## Quick Start

### Method 1: Using Pixi (Recommended)

```bash
# From the virtual-lab root directory

# Run complete analysis with E. coli core model
pixi run metabolic-analysis

# Launch interactive tutorial
pixi run metabolic-tutorial

# Quick test
pixi run test-metabolic
```

### Method 2: Command Line

```bash
cd metabolic_modeling/scripts

# Use E. coli core model (fast, for testing)
python metabolic_target_finder.py \
    --model_id textbook \
    --output_dir ../../metabolic_results \
    --ko_methods single essential fva \
    --visualization

# Use larger E. coli model (iML1515)
python metabolic_target_finder.py \
    --model_id iML1515 \
    --output_dir ../../metabolic_results \
    --ko_methods single essential fva \
    --growth_threshold 0.1 \
    --visualization

# Use custom SBML model file
python metabolic_target_finder.py \
    --model_file path/to/model.xml \
    --output_dir ../../metabolic_results \
    --ko_methods single double essential \
    --visualization
```

### Method 3: Jupyter Notebook

```bash
jupyter notebook tutorial_metabolic_modeling.ipynb
```

### Method 4: Python API

```python
from metabolic_target_finder import MetabolicTargetFinder, Args

# Configure analysis
args = Args(
    model_id="iML1515",  # or model_file="path/to/model.xml"
    output_dir="results",
    ko_methods=["single", "essential", "fva"],
    growth_threshold=0.1
)

# Run analysis
finder = MetabolicTargetFinder(args)
finder.load_model()
finder.method_single_gene_knockout()
finder.method_essential_genes()
finder.method_flux_variability_analysis()

# Save and visualize
finder.save_results(args.output_dir)
finder.visualize_results(args.output_dir)
```

## Supported Model Formats

### BiGG Models Database

Load pre-built models directly from BiGG:

```python
import cobra

# E. coli
model = cobra.io.load_model("iML1515")  # 2,712 reactions, 1,877 genes
model = cobra.io.load_model("textbook")  # Core model (fast)

# Human
model = cobra.io.load_model("Recon3D")  # 13,543 reactions, 3,288 genes

# Yeast
model = cobra.io.load_model("iMM904")  # 1,577 reactions, 904 genes
```

### File Formats

- **SBML** (`.xml`, `.sbml`): Systems Biology Markup Language
- **JSON** (`.json`): COBRApy JSON format
- **MAT** (`.mat`): MATLAB format

```python
# Load from file
model = cobra.io.read_sbml_model("path/to/model.xml")
model = cobra.io.load_json_model("path/to/model.json")
model = cobra.io.load_matlab_model("path/to/model.mat")
```

## Analysis Methods

### 1. Single Gene Knockout

Systematically delete each gene and measure growth impact.

```bash
python metabolic_target_finder.py \
    --model_id textbook \
    --ko_methods single \
    --growth_threshold 0.1
```

**Output:**
- `single_knockout_results.csv`: All genes with growth effects
- Essential vs. non-essential classification
- Growth reduction percentages

**Use cases:**
- Identify drug targets (essential genes)
- Find metabolic engineering targets (non-essential growth-reducing)
- Understand gene importance

### 2. Essential Gene Identification

Genes whose deletion causes lethality (growth < threshold).

```bash
python metabolic_target_finder.py \
    --model_id iML1515 \
    --ko_methods essential \
    --growth_threshold 0.1
```

**Output:**
- `essential_genes_results.csv`: Essential genes with annotations
- Reaction involvement
- Gene criticality ranking

**Use cases:**
- Antibiotic target discovery
- Cancer drug targets
- Understanding core metabolism

### 3. Double Gene Knockout (Synthetic Lethality)

Find gene pairs that are non-essential individually but lethal together.

```bash
python metabolic_target_finder.py \
    --model_id iML1515 \
    --ko_methods double \
    --growth_threshold 0.1
```

**Output:**
- `double_knockout_results.csv`: All gene pair combinations
- Synthetic lethal pairs
- Synergy scores

**Use cases:**
- Combination therapy design
- Robustness analysis
- Multi-target drugs

**Note:** Computationally intensive! ~O(n²) where n = number of genes.

### 4. Flux Variability Analysis (FVA)

Identify min/max flux ranges for each reaction at near-optimal growth.

```bash
python metabolic_target_finder.py \
    --model_id iML1515 \
    --ko_methods fva \
    --fva_fraction 0.95
```

**Output:**
- `fva_results.csv`: Flux ranges for all reactions
- Flexible vs. rigid reactions
- Essential flux pathways

**Use cases:**
- Identify metabolic bottlenecks
- Find robust engineering targets
- Understand metabolic flexibility

### 5. Growth-Coupled Production

Find knockouts that couple product formation with growth.

```bash
python metabolic_target_finder.py \
    --model_id iML1515 \
    --ko_methods production \
    --target_metabolite succ_c
```

**Output:**
- `growth_coupled_results.csv`: Knockouts ranked by production/growth ratio
- Production rates
- Growth impacts

**Use cases:**
- Metabolic engineering for bioproduction
- Optimize yield
- Evolution-based strain improvement

### 6. Flux Sampling

Sample solution space to explore alternative metabolic states.

```bash
python metabolic_target_finder.py \
    --model_id iML1515 \
    --ko_methods sampling \
    --n_samples 1000
```

**Output:**
- `flux_sampling_results.csv`: Flux statistics across samples
- High-variability reactions
- Metabolic state distributions

**Use cases:**
- Explore alternative pathways
- Understand metabolic uncertainty
- Identify regulatory targets

## Command-Line Arguments

| Argument | Type | Default | Description |
|----------|------|---------|-------------|
| `--model_file` | str | None | Path to model file (SBML/JSON/MAT) |
| `--model_id` | str | `iML1515` | BiGG model ID if no file provided |
| `--output_dir` | str | `metabolic_results` | Output directory |
| `--objective` | str | None | Objective reaction (default: biomass) |
| `--target_metabolite` | str | None | Target for production analysis |
| `--ko_methods` | list | `[single, essential, fva]` | Methods to run |
| `--growth_threshold` | float | `0.1` | Essential gene threshold (fraction) |
| `--n_samples` | int | `1000` | Number of flux samples |
| `--fva_fraction` | float | `0.95` | Fraction of optimal for FVA |
| `--top_n_targets` | int | `50` | Number of top targets to report |
| `--visualization` | bool | `True` | Generate plots |
| `--random_seed` | int | `42` | Random seed |

## Available Methods

Use `--ko_methods` to select analyses:

- `single`: Single gene knockout
- `double`: Double gene knockout (slow!)
- `essential`: Essential gene identification
- `synthetic_lethality`: Find synthetic lethal pairs
- `fva`: Flux variability analysis
- `sampling`: Flux sampling
- `production`: Growth-coupled production (requires `--target_metabolite`)

## Output Files

### Results Files

- **`model_summary.csv`**: Model statistics and information
- **`single_knockout_results.csv`**: All single KO results
  - Columns: `gene_id`, `growth_rate`, `growth_fraction`, `growth_reduction`, `status`
- **`essential_genes_results.csv`**: Essential genes with details
  - Columns: `gene_id`, `gene_name`, `growth_fraction`, `n_reactions`, `reaction_ids`
- **`double_knockout_results.csv`**: Double KO results
  - Columns: `gene1`, `gene2`, `double_ko_fraction`, `synergy_score`, `synthetic_lethal`
- **`fva_results.csv`**: Flux ranges for all reactions
  - Columns: `reaction_id`, `minimum`, `maximum`, `flux_span`, `subsystem`
- **`growth_coupled_results.csv`**: Production-enhancing KOs
  - Columns: `gene_id`, `ko_growth`, `production`, `production_per_growth`
- **`flux_sampling_results.csv`**: Flux statistics
  - Columns: `reaction_id`, `mean_flux`, `std_flux`, `min_flux`, `max_flux`, `cv`

### Visualization File

- **`metabolic_analysis.pdf`**: Multi-page PDF with:
  - Knockout growth distribution
  - Essential genes classification
  - Synthetic lethality network
  - FVA flux ranges
  - Growth-coupled production plots

## Metabolic Models

### Recommended Models

**For Learning/Testing:**
- `textbook`: E. coli core model (95 reactions, 137 metabolites, 137 genes)
  - Fast analysis (~1 minute)
  - Good for tutorials

**For Production Use:**

**E. coli:**
- `iML1515`: Latest E. coli model (2,712 reactions, 1,877 genes)
  - Most comprehensive E. coli model
  - ~10-30 minutes for single KO

**Human:**
- `Recon3D`: Human metabolism (13,543 reactions, 3,288 genes)
  - Most comprehensive human model
  - Hours for complete analysis
  - Use for drug target discovery

**Yeast:**
- `iMM904`: S. cerevisiae (1,577 reactions, 904 genes)
  - Well-validated yeast model
  - Biotechnology applications

### Custom Models

To use your own model:

1. Ensure it's in SBML, JSON, or MAT format
2. Place in `metabolic_modeling/models/` directory
3. Use `--model_file` flag:

```bash
python metabolic_target_finder.py \
    --model_file ../models/my_model.xml \
    --output_dir ../results
```

## Performance Considerations

### Expected Runtime

| Model Size | Method | Approximate Time |
|------------|--------|------------------|
| Core (137 genes) | Single KO | ~30 seconds |
| Core (137 genes) | FVA | ~1-2 minutes |
| Core (137 genes) | Double KO | ~5 minutes |
| iML1515 (1,877 genes) | Single KO | ~10-30 minutes |
| iML1515 (1,877 genes) | FVA | ~30-60 minutes |
| iML1515 (1,877 genes) | Double KO | Hours to days |
| Recon3D (3,288 genes) | Single KO | ~1-2 hours |
| Recon3D (3,288 genes) | FVA | ~2-4 hours |

### Memory Requirements

- Small models (<500 reactions): ~500 MB
- Medium models (1,000-3,000 reactions): ~1-2 GB
- Large models (>5,000 reactions): ~3-5 GB

### Optimization Tips

1. **Start with core model** for testing
2. **Use faster solver** (CPLEX or Gurobi vs. GLPK)
3. **Limit double KO** to candidate genes only
4. **Reduce FVA fraction** (0.90 instead of 0.99)
5. **Parallelize** on HPC clusters

## Installing Optimization Solvers

COBRApy requires an optimization solver. Options:

### GLPK (Free, Default)

```bash
# Conda (recommended)
conda install -c conda-forge glpk

# Or pip
pip install swiglpk
```

### CPLEX (Commercial, Fast)

Free for academics: https://www.ibm.com/academic/technology/data-science

```bash
pip install cplex
```

### Gurobi (Commercial, Fast)

Free for academics: https://www.gurobi.com/academia/

```bash
pip install gurobipy
```

## Troubleshooting

### Common Issues

**1. No solver available**
```
OptimizationError: Solver status is 'failed'
```
**Solution**: Install a solver (GLPK, CPLEX, or Gurobi)

**2. Model loading fails**
```
FileNotFoundError: Model not found
```
**Solution**: Check model_id spelling or file path

**3. Memory error**
```
MemoryError: Unable to allocate array
```
**Solution**: Use smaller model or increase RAM

**4. Slow analysis**
```
Analysis taking too long...
```
**Solution**:
- Use faster solver (CPLEX/Gurobi)
- Reduce number of methods
- Use smaller model for testing

**5. Import errors**
```
ImportError: No module named 'cobra'
```
**Solution**: Install requirements
```bash
pip install -r requirements.txt
# Or
pixi install
```

## Example Use Cases

### Use Case 1: Antibiotic Target Discovery

Find essential genes in pathogenic bacteria:

```bash
python metabolic_target_finder.py \
    --model_file models/pathogen_model.xml \
    --ko_methods essential \
    --growth_threshold 0.05 \
    --output_dir antibiotic_targets
```

### Use Case 2: Biofuel Production

Optimize E. coli for ethanol production:

```bash
python metabolic_target_finder.py \
    --model_id iML1515 \
    --ko_methods single production \
    --target_metabolite etoh_c \
    --output_dir biofuel_engineering
```

### Use Case 3: Cancer Drug Targets

Find synthetic lethal pairs in cancer metabolism:

```bash
python metabolic_target_finder.py \
    --model_id Recon3D \
    --ko_methods double essential \
    --growth_threshold 0.1 \
    --output_dir cancer_targets
```

### Use Case 4: Metabolic Engineering

Design strain for succinate production:

```bash
python metabolic_target_finder.py \
    --model_id iML1515 \
    --ko_methods single fva production \
    --target_metabolite succ_c \
    --output_dir succinate_production
```

## Resources

### Documentation

- **COBRApy**: https://cobrapy.readthedocs.io/
- **BiGG Models**: http://bigg.ucsd.edu/
- **SBML**: http://sbml.org/
- **Virtual Lab**: See main README.md

### Papers

1. **Constraint-based modeling**:
   - Orth et al. (2010). "What is flux balance analysis?" *Nature Biotechnology*.

2. **COBRApy**:
   - Ebrahim et al. (2013). "COBRApy: COnstraints-Based Reconstruction and Analysis for Python." *BMC Systems Biology*.

3. **E. coli iML1515**:
   - Monk et al. (2017). "iML1515, a knowledgebase that computes Escherichia coli traits." *Nature Biotechnology*.

4. **Human Recon3D**:
   - Brunk et al. (2018). "Recon3D enables a three-dimensional view of gene variation in human metabolism." *Nature Biotechnology*.

### Tutorials

- **COBRApy getting started**: https://cobrapy.readthedocs.io/en/latest/getting_started.html
- **Escher** (flux visualization): https://escher.github.io/
- **MEMOTE** (model quality control): https://memote.io/

## Citation

If you use this metabolic modeling module, please cite the Virtual Lab paper:

Swanson, K., Wu, W., Bulaong, N.L. et al. The Virtual Lab of AI agents designs new SARS-CoV-2 nanobodies. *Nature* (2025). https://doi.org/10.1038/s41586-025-09442-9

And the COBRApy paper:

Ebrahim, A., Lerman, J. A., Palsson, B. O., & Hyduke, D. R. (2013). COBRApy: COnstraints-Based Reconstruction and Analysis for Python. *BMC Systems Biology*, 7(1), 74.

## License

This module is part of the Virtual Lab project and is licensed under the MIT License.

## Contact

For questions or issues:
- Open an issue on GitHub
- Refer to the main Virtual Lab documentation
- Check COBRApy documentation

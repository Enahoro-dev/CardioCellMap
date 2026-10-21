# CardioCellMap 

CardioCellMap is an exploratory computational biology project investigating spatial gene expression, cardiac cell populations, and inflammatory gene patterns in post MI human heart tissue.

# Current work

Using a human myocardial infarction spatial transcriptomics dataset, the project currently explores:

Spatial expression of IL6, IL1B, CXCL8, NFKB1, and ISG15.

Estimated cell-type composition across tissue spots.

Preliminary spatial relationships between gene-expressing spots and their neighbours.

# Technologies

Python, AnnData, NumPy, pandas, and Matplotlib.

# Future direction

The long-term goal is to develop a machine learning-powered computational tool to identify spatial cellular niches associated with inflammatory gene programmes in cardiac tissue.

The current work is exploratory. Machine learning development and validation remain future objectives.

# Dataset

[Processed data on Zenodo](https://doi.org/10.5281/zenodo.6578047). This project uses `Visium_control_P1.h5ad`.


# Source study
Kuppe et al. (2022), [Spatial multi-omic map of human myocardial infarction](https://doi.org/10.1038/s41586-022-05060-x).


## Running the Analysis

### 1. Clone the Repository

```bash
git clone https://github.com/Enahoro-dev/CardioCellMap.git
cd CardioCellMap
```

### 2. Create and Activate a Virtual Environment

```bash
python3 -m venv .venv
source .venv/bin/activate
```

### 3. Install Dependencies

```bash
python -m pip install anndata numpy pandas matplotlib
```

### 4. Obtain the Dataset

Place `Visium_control_P1.h5ad` in the project's `data/` directory. Obtain the dataset from its original source and check its access and reuse terms.

### 5. Run the Analysis

```bash
python notebooks/spatial_map.py
```

The script reads the dataset from `data/Visium_control_P1.h5ad` and saves the analysis outputs in the `results/` directory.

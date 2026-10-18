
import anndata as ad
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os

# Create the results folder if it does not exist
os.makedirs("results", exist_ok=True)

# Load the dataset
adata = ad.read_h5ad(
    "data/Visium_control_P1.h5ad",
    backed="r"
)

# Genes we want to investigate
genes = ["IL6", "IL1B", "CXCL8", "NFKB1", "ISG15"]

# Get spatial coordinates for every spot
coordinates = adata.obsm["X_spatial"]


# 1. Plot gene expression across the tissue
def plot_gene(gene):
    scores = adata[:, gene].X[:].toarray().flatten()

    plt.figure(figsize=(7, 6))

    plt.scatter(
        coordinates[:, 0],
        coordinates[:, 1],
        c=scores,
        s=8
    )

    plt.gca().invert_yaxis()
    plt.colorbar(label=f"{gene} expression")
    plt.title(f"{gene} expression across cardiac tissue")
    plt.xlabel("Spatial X coordinate")
    plt.ylabel("Spatial Y coordinate")
    plt.tight_layout()

    plt.savefig(
        f"results/{gene.lower()}_spatial_map.png",
        dpi=200
    )

    plt.close()


# Generate a spatial map for every gene
for gene in genes:
    plot_gene(gene)


# 2. Calculate expression statistics for every gene
expression_results = []

for gene in genes:
    scores = adata[:, gene].X[:].toarray().flatten()

    detected = int((scores > 0).sum())
    percentage = detected / len(scores) * 100
    mean_expression = scores.mean()

    expression_results.append({
        "gene": gene,
        "detected_spots": detected,
        "percentage": round(percentage, 2),
        "mean_expression": round(float(mean_expression), 4)
    })

expression_df = pd.DataFrame(expression_results)

expression_df.to_csv(
    "results/gene_expression_summary.csv",
    index=False
)

print("\nGene expression summary:")
print(expression_df.to_string(index=False))


# 3. Preliminary neighbour analysis for every gene
neighbor_results = []

for gene in genes:
    scores = adata[:, gene].X[:].toarray().flatten()
    positive_indices = np.where(scores > 0)[0]

    neighbor_counts = []

    for spot_index in positive_indices:
        spot_coordinates = coordinates[spot_index]

        distances = np.linalg.norm(
            coordinates - spot_coordinates,
            axis=1
        )

        # Exclude the spot itself
        distances[spot_index] = np.inf

        # Select the four nearest other spots
        nearest_indices = np.argsort(distances)[:4]

        # Count neighbours with detectable expression of this gene
        count = int((scores[nearest_indices] > 0).sum())
        neighbor_counts.append(count)

    values, counts = np.unique(
        neighbor_counts,
        return_counts=True
    )

    for value, count in zip(values, counts):
        neighbor_results.append({
            "gene": gene,
            "positive_neighbors": int(value),
            "number_of_spots": int(count),
            "total_positive_spots": len(positive_indices)
        })

neighbor_df = pd.DataFrame(neighbor_results)

neighbor_df.to_csv(
    "results/gene_neighbor_summary.csv",
    index=False
)

print("\nPreliminary neighbour summary:")
print(neighbor_df.to_string(index=False))


# Plot the spatial layout of all spots
plt.figure(figsize=(7, 6))

plt.scatter(
    coordinates[:, 0],
    coordinates[:, 1],
    s=8
)

plt.gca().invert_yaxis()
plt.title("Spatial layout of all tissue spots")
plt.xlabel("Spatial X coordinate")
plt.ylabel("Spatial Y coordinate")
plt.axis("equal")
plt.tight_layout()
plt.savefig("results/tissue_spatial_layout.png", dpi=200)
plt.close()


# 4. Summarise estimated cell-type composition

cell_type_columns = [
    "Adipocyte",
    "Cardiomyocyte",
    "Endothelial",
    "Fibroblast",
    "Lymphoid",
    "Mast",
    "Myeloid",
    "Neuronal",
    "Pericyte",
    "Cycling.cells",
    "vSMCs"
]

# Calculate the mean estimated proportion of each cell type
mean_proportions = adata.obs[cell_type_columns].mean()

cell_type_df = pd.DataFrame({
    "cell_type": mean_proportions.index,
    "mean_proportion": mean_proportions.values,
    "mean_percentage": mean_proportions.values * 100
})

# Sort from the highest mean proportion to the lowest
cell_type_df = cell_type_df.sort_values(
    "mean_proportion",
    ascending=False
)

# Save the summary
cell_type_df.to_csv(
    "results/mean_cell_type_composition.csv",
    index=False
)

print("\nMean estimated cell-type composition:")
print(cell_type_df.round(4).to_string(index=False))


# Close the dataset
adata.file.close()

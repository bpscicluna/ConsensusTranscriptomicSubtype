# ConsensusTranscriptomicSubtype

Classify patients into **Consensus Transcriptomic Subtypes (CTS)** of **sepsis** using gene expression data. This package implements batch correction, supervised classification via random forests, and optional visualization of results.

---

## 🔧 Installation

From source:
```r
install.packages("ConsensusTranscriptomicSubtype_0.1.4.tar.gz", repos = NULL, type = "source")
```

Or using `devtools`:
```r
devtools::install_github("bpscicluna/ConsensusTranscriptomicSubtype")
```

---

## 🚀 Usage

### Basic classification
```r
library(ConsensusTranscriptomicSubtype)

# new_expr_data must be a data.frame or matrix
# First column = gene identifiers (optional), rownames = Ensembl IDs
new_expr <- read.csv("your_new_data.csv")

result <- run_subtype_classifier(new_expr_data = new_expr)
head(result$predictions)
```

### Save heatmap and silhouette plot
```r
result <- run_subtype_classifier(
  new_expr_data = new_expr,
  heatmap_file = "cts_heatmap.pdf",
  silhouette_file = "cts_silhouette.pdf"
)
```

---

## 📦 Internal data
The package comes with:
- `exp_core_g`: reference gene expression matrix (ENSEMBL x samples)
- `core_samples`: subtype assignments for training

---

## 📄 Vignette: Classifying Sepsis Subtypes

```r
# Load package
library(ConsensusTranscriptomicSubtype)

# Prepare your new expression matrix (same genes as training set)
expr <- read.csv("example_sepsis_expr.csv")

# Run classifier and save visual output
result <- run_subtype_classifier(
  new_expr_data = expr,
  heatmap_file = "cts_heatmap.pdf",
  silhouette_file = "cts_silhouette.pdf"
)

# View predicted classes
result$predictions
```

---

## 📫 Contact
For feedback or contributions, please open an issue or pull request.

---

## 🧪 License
This package is licensed under the **GNU General Public License v3.0 (GPL-3.0)**.

See the full license text in the [LICENSE](LICENSE) file.

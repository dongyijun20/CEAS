# CEAS Single-Cell Transcriptomics

## Overview

**Chronic enteropathy associated with SLCO2A1 (CEAS)** is a rare monogenic disorder characterized by intestinal ulceration and fibrostenosis. Despite its clinical significance, the underlying mechanisms of CEAS remain poorly understood.

In this project, we constructed the **first single-cell transcriptomic map of CEAS** using intestinal biopsy samples from six genetically confirmed patients. Our analysis revealed:

- Expansion of stromal populations, including a CEAS-specific endothelial subcluster.
- A differentiation blockade between capillary and venous endothelial states.
- Abnormal angiogenesis driven by SLCO2A1 knockdown under prostaglandin imbalance.
- Crosstalk between epithelial and stromal cells via the **HIF1–VEGFA–PTGS2** axis.
- Potential therapeutic targets including **VEGFA** and **JAK inhibitors**.

These findings identify endothelial dysfunction as a key driver of CEAS and suggest novel therapeutic avenues.

![Uploading CEAS graphical abstract.png…]()

---

## Patient Recruitment and Sample Collection

- **6 CEAS patients** from Peking Union Medical College Hospital.
- Diagnosis confirmed by clinical criteria and SLCO2A1 variants (WES + Sanger).
- **9 patient samples**: 8 lesional (stomach, duodenum, ileum, colon) + 1 non-lesional.
- **3 control samples** from 2 healthy individuals.
- Ethical approval (IRB: I-25PJ0672) and informed consent obtained.

---

## Single-Cell RNA Sequencing

- Prepared using **Chromium Single Cell 3’ v3.1 (10x Genomics)**.
- Sequenced on **Illumina NovaSeq 6000** with 100,000+ reads/cell.
- Libraries constructed with Chromium Single Cell 3’ Library Kit v3.1.

---

## Data Processing and Integration

- Processed using **Cell Ranger v6.1.1** and **Seurat v4.4.0**.
- Quality control:
  - 500–10,000 genes/cell
  - 1,000–60,000 UMIs
  - <10% mitochondrial content
- Batch correction with **Harmony v1.2.1**.
- Cell type annotation via **scIBD reference** and canonical marker genes.

---

## Downstream Analyses

### 🔹 Differential Gene Expression & Enrichment

- DEGs via **Wilcoxon test** (Seurat).
- Functional enrichment using **clusterProfiler** with KEGG and GO.

### 🔹 Pseudotime Analysis

- Performed using **Monocle v2 & v3** on endothelial cells.
- Revealed trajectory and stage-specific gene modules.

### 🔹 Cell-Cell Communication

- Analyzed with **NicheNet**.
- Identified ligand-receptor pairs impacting endothelial cells.
- Highlighted **HIF1–VEGFA–PTGS2** axis as key intercellular pathway.

### 🔹 Endothelial Subpopulation Discovery

- Used **MiloR v0.1.0** for differential abundance.
- Identified a CEAS-specific endothelial cluster (Cluster 3).

---

## Key Findings

- **Endothelial dysfunction** is a central feature of CEAS.
- **Immature angiogenesis** observed in SLCO2A1-deficient conditions.
- Potential drug targets include **VEGFA** and **JAK** inhibitors.
- This study provides a valuable resource for future CEAS research and therapy development.

---

## License

This project is licensed under the MIT License. See [LICENSE](LICENSE) for more information.

---

## Contact

For questions or collaborations, please contact:  
**Yijun Dong**  
Peking Union Medical College
dongyijunms@outlook.com

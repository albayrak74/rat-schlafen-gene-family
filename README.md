# Comparative genomic and cross-dataset transcriptomic analysis of the rat Schlafen gene family in aging and inflammation

Analysis code, derived results and figures for the manuscript by İrem Alkan and Mehmet Albayrak (submitted to *Gene*).

Corresponding author: Mehmet Albayrak, Isparta University of Applied Sciences (mehmetalbayrak@isparta.edu.tr).

## Contents

| Folder | Content |
|---|---|
| `scripts/` | Analysis scripts, numbered in the order they were run (see below) |
| `results/` | All derived tables (CSV/JSON) used in the manuscript and in the Supplementary Data file |
| `figures/` | Figures 1–7 (PNG, 200 dpi; TIFF, 600 dpi) |
| `phylogeny/` | Protein sequences, MAFFT alignment, trimmed alignment and IQ-TREE 2 output |
| `hmm/` | Pfam profile HMMs used (PF17057, PF04326, PF21026, PF09848, PF13538) and the alignments used to build the family-specific profiles |
| `manuscript/` | Manuscript source (Markdown with citation keys), verified reference metadata and the scripts that build the Word files |

## Data sources (not redistributed here)

All input data are public:

- Rat genome GRCr8, RefSeq annotation GCF_036323735.1 (GFF, protein FASTA, CDS FASTA) from NCBI.
- Mouse (GRCm39) and human (GRCh38) proteins and CDS from Ensembl; sequences of specific regions via the Ensembl REST API.
- Pfam HMMs from the InterPro API (`https://www.ebi.ac.uk/interpro/wwwapi/entry/pfam/<ACC>?annotation=hmm`).
- Expression data from GEO: GSE53960 (rat BodyMap, per-sample count files from the RAW archive), GSE149517, GSE228906, GSE261056, GSE298271, GSE192425 and GSE160153 (supplementary processed files as described in the Methods).

Place the GEO files under `geo/` (BodyMap count files under `geo/bodymap/`) and the genome files under `ref/` before running the scripts.

## Scripts

Run from the repository root. Outputs are written to `out/` (the committed copies are in `results/`).

| Script | Purpose |
|---|---|
| `b1_genes.py` | Rat Schlafen gene models, representative proteins, ProtParam metrics |
| `b2_domains.py` | Family-specific profile HMMs built from mouse/human Pfam-annotated segments; search of rat proteins |
| `b2b_pfam.py` | Pfam HMM search with gathering thresholds |
| `b3_ortho.py` | Mouse and human canonical Schlafen proteins and Pfam domains; sequences for the tree |
| `b4_promoter.py`, `b4b_enrich.py` | Promoter motif counts, genomic background and enrichment tests |
| `b5_kaks.py` | Ka/Ks (Nei–Gojobori) for rat–mouse pairs and paralogue pairs |
| `b6_completeness.py` | Proteome-wide Pfam search and tblastn search of the cluster region |
| `c1_bodymap.py` | BodyMap expression, τ, age and sex differential expression (PyDESeq2) |
| `c2_wistar.py` | Differential expression in the six Wistar datasets |
| `c3_ml.py` | Inflammation classifiers, leave-one-dataset-out, random-set null, SHAP |
| `c4_age_ml.py` | Age classifiers, leave-one-organ-out, random-set null, SHAP, external test |
| `c5_revision.py` | Sensitivity analyses: testis excluded, Ptprc adjustment, variance ranking on unscaled data |
| `d1_figs.py` | Figures 1–7 |

Phylogeny: `mafft --localpair --maxiterate 1000` on `phylogeny/slfn.fa`, trimming to columns with ≥50% occupancy, then `iqtree2 -s slfn_trim.fa -m MFP -B 1000`.

## Software versions

Python 3.13; pandas 3.0.5, NumPy 2.5.3, SciPy 1.18.1, statsmodels 0.15.0, scikit-learn 1.9.1, shap 0.54.1, PyDESeq2 0.5.4, Biopython 1.88, PyHMMER 0.12.3, matplotlib 3.11.2; MAFFT 7.505; IQ-TREE 2.0.7; BLAST+ 2.12.0.

## License

Code: MIT License. Derived tables and figures: CC BY 4.0.

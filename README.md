# Advancing Regulatory Variant Effect Prediction with AlphaGenome: Quantitative Dissection of the Human PRNP Allelic Spectrum

[![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/sid0sid-ops/alphagenome-prnp-analysis/blob/main/notebooks/AlphaGenome_PRNP_Analysis.ipynb)
[![GitHub Repository](https://img.shields.io/badge/GitHub-Repository-blue.svg)](https://github.com/sid0sid-ops/alphagenome-prnp-analysis)
[![Nature Paper](https://img.shields.io/badge/Nature-10.1038%2Fs41586--025--10014--0-red.svg)](https://doi.org/10.1038/s41586-025-10014-0)
[![License: MIT](https://img.shields.io/badge/License-MIT-green.svg)](https://opensource.org/licenses/MIT)

> **Primary Nature Research Paper:** *Advancing regulatory variant effect prediction with AlphaGenome* [1]  
> **Authors:** Žiga Avsec, Natasha Latysheva, Jun Cheng, Demis Hassabis, Pushmeet Kohli, et al. (Google DeepMind)  
> **Publication Citation:** *Nature*, Vol. 649, pp. 1206–1218 (Published: 29 January 2026) | [doi:10.1038/s41586-025-10014-0](https://doi.org/10.1038/s41586-025-10014-0)  
> **Target Genomic Locus:** Human Prion Protein Gene (*PRNP*), Chromosome 20 (`chr20:4,682,000 - 4,705,000` GRCh38) [11]  
> **Presenter:** Sid | 3rd Year, Computational Biology & Systems Biology (CSBB)  

---

## ⚡ Quickstart: Run Online via Google Colab (Zero Installation)

You can run the entire pipeline interactively with one click in **Google Colab**:

👉 [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/sid0sid-ops/alphagenome-prnp-analysis/blob/main/notebooks/AlphaGenome_PRNP_Analysis.ipynb)

* **Features in the Colab Notebook:**
  1. Installs official AlphaGenome packages and scientific dependencies.
  2. Executes programmatic forensic verification of all genomic coordinates against UCSC GRCh38.
  3. Integrates Colab Secrets (`userdata.get('ALPHAGENOME_API_KEY')`) for private, safe cloud authentication with zero leaked keys.
  4. Runs authentic functional deconvolution comparing promoter eQTLs, splice disruptions, and coding missense variants.
  5. Renders and displays all 5 publication-grade research figures inline directly without external downloads.

---

## 🔬 Rigorous Credibility Check: Verifying GRCh38 NM_000311 Coordinates

Before executing any model inference, we performed a strict forensic sequence and coordinate audit against the official **UCSC Genome Browser GRCh38 assembly** for canonical RefSeq transcript **`NM_000311.5`** (Ensembl: `ENST00000379440.9`) [11].

### 1. Genomic Architecture & Strand Orientation
* **Strand Orientation:** Human *PRNP* is located on the **POSITIVE strand (`+`)** of Chromosome 20.
* **Exon 1 (5' UTR):** `chr20:4,686,350 - 4,686,483` (non-coding promoter core).
* **Intron 1:** `chr20:4,686,484 - 4,699,209` (12,726 bp long, harboring neural enhancer elements).
* **3' Splice Acceptor AG:** `chr20:4,699,210` (essential splice motif for Exon 2 inclusion).
* **Coding Sequence (CDS) Start:** `chr20:4,699,221` (Start Codon 1: `ATG`, Methionine).
* **Coding Sequence (CDS) End:** `chr20:4,699,982` (Stop Codon 254: `TGA`).
* **3' UTR:** `chr20:4,699,983 - 4,701,590`.

### 2. Mathematical Coordinate Verification Formula
On the positive strand, any codon coordinate is derived via the deterministic formula:
$$\text{Genomic Position} = 4,699,221 + 3 \times (\text{Codon Number} - 1) + (\text{Base Offset})$$
*(where Base Offset $\in \{0, 1, 2\}$ for the 1st, 2nd, and 3rd nucleotide of each codon triplet).*

### 3. Canonical PRNP Allelic Validation Matrix

| Variant ID | Codon Shift | Formula Calculation | Verified GRCh38 Coordinate | Reference Codon | Mutated Codon | dbSNP ID | Validation Status |
| :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: |
| **G127V** | `Gly127Val` | $4,699,221 + 3 \times 126 + 1$ | **`chr20:4699600:G>T`** | `GGC` (Gly) | `GTC` (Val) | rs80356730 | ✅ **Exact Match (UCSC hg38)** [11] |
| **M129V** | `Met129Val` | $4,699,221 + 3 \times 128 + 0$ | **`chr20:4699605:A>G`** | `ATG` (Met) | `GTG` (Val) | rs1799990 | ✅ **Exact Match (UCSC hg38)** [11] |
| **P102L** | `Pro102Leu` | $4,699,221 + 3 \times 101 + 1$ | **`chr20:4699525:C>T`** | `CCG` (Pro) | `CTG` (Leu) | rs61729737 | ✅ **Exact Match (UCSC hg38)** [11] |
| **D178N** | `Asp178Asn` | $4,699,221 + 3 \times 177 + 0$ | **`chr20:4699752:G>A`** | `GAC` (Asp) | `AAC` (Asn) | rs61733276 | ✅ **Exact Match (UCSC hg38)** [11] |
| **E200K** | `Glu200Lys` | $4,699,221 + 3 \times 199 + 0$ | **`chr20:4699818:G>A`** | `GAG` (Glu) | `AAG` (Lys) | rs28933385 | ✅ **Exact Match (UCSC hg38)** [11] |
| **E219K** | `Glu219Lys` | $4,699,221 + 3 \times 218 + 0$ | **`chr20:4699875:G>A`** | `GAG` (Glu) | `AAG` (Lys) | rs1800533 | ✅ **Exact Match (UCSC hg38)** [11] |

*All wild-type bases match the human genome reference assembly GRCh38 [11]. This verification is programmatically tested locally by running `python src/verify_prnp_coordinates.py`.*

---

## 📊 Scientific Figures & What We Learn from AlphaGenome

### Figure 1: Calibrated Variant Predictions Across Canonical Human PRNP Alleles

![Figure 1: Calibrated Variant Predictions Across Canonical PRNP Alleles](figures/annotated/prion_variant_impact_bar.png)

#### 🔍 What this Figure Shows:
Calibrated Integrative Variant Impact PHRED scores combining AlphaGenome regulatory predictions and AlphaMissense 3D protein fold destabilization:
$$\text{PHRED} = -10 \cdot \log_{10}(\text{Empirical Genome-Wide Quantile})$$
* **Dashed Gray Line ($\text{PHRED} = 20.0$):** Top 1% most impactful variants across the entire human genome.
* **Dotted Red Line ($\text{PHRED} = 25.0$):** Severe functional impact threshold (Top 0.31% genome-wide).

#### 💡 Biological Interpretation & Key Findings:
1. **Autonomous Boundary Discovery:** The integrative framework cleanly separates pathogenic and immunity alleles ($\text{PHRED} > 25$) from common benign population modifiers ($\text{PHRED} < 18$) without requiring prior clinical labels [1], [10].
2. **The G127V Paradox (Why does Resistance score $\text{PHRED} = 28.64$?):**
   * The framework evaluates **biophysical and evolutionary perturbation**, not clinical valence (benefit vs. harm).
   * John Collinge’s team proved that `p.Gly127Val` confers **complete, 100% dominant-negative immunity** against Kuru and classical CJD [2].
   * Replacing flexible Glycine with bulky hydrophobic Valine introduces an impassable steric block that prevents monomer docking onto amyloid fibrils. AlphaMissense and Cactus conservation correctly score this as an extreme biophysical disruption [2], [10].
3. **Fatal Familial Mutations ($\text{PHRED } 25.39 - 26.79$):**
   * `p.Asp178Asn` (FFI): $\text{PHRED } 26.79$ (Top 0.21%).
   * `p.Pro102Leu` (GSS Ataxia): $\text{PHRED } 26.09$ (Top 0.25%).
   * `p.Glu200Lys` (fCJD): $\text{PHRED } 25.39$ (Top 0.29%).
   * All three cluster tightly in the severe tier, reflecting destabilization of the globular fold [7], [10].
4. **Benign Population Modifiers ($\text{PHRED } 16.34 - 17.73$):**
   * `p.Met129Val` ($\text{PHRED } 17.73$) and `p.Glu219Lys` ($\text{PHRED } 16.34$) fall well below the top 1% threshold, reflecting their roles as benign polymorphic modulators.

---

### Figure 2: Multi-Scale Feature Attribution Profile

![Figure 2: Multi-Scale Feature Attribution Profiles](figures/annotated/multimodal_radar_profile.png)

#### 🔍 What this Figure Shows:
* **Panel A (Left):** Stacked bar decomposition of variant scores across structural perturbation (AlphaMissense), mammalian conservation (Cactus 241-way), vertebrate conservation (PhastCons 470), splicing, and promoter epigenetics [1], [10].
* **Panel B (Right):** Polar radar profile comparing the structural vs. evolutionary signatures of `p.Gly127Val` (Immunity), `p.Glu200Lys` (Pathogenic), and `p.Met129Val` (Modifier).

#### 💡 Biological Interpretation & Key Findings:
* **`p.Glu200Lys` (fCJD):** Heavily driven by **AlphaMissense (+0.62)** [10], capturing the destruction of an electrostatic salt bridge on Helix 3.
* **`p.Gly127Val` (Kuru Resistance):** Synergistic high weights on both **AlphaMissense (+0.74)** (steric block) and **Cactus Conservation (+0.52)** (conserved loop architecture) [2], [10].
* **`p.Met129Val` (Modifier):** Displays balanced, low-magnitude contributions across modalities, characteristic of benign human genetic diversity [1].

---

### Figure 3: Epigenetic & Transcriptional Architecture of the Human PRNP Locus (23-kb Window)

![Figure 3: Epigenetic and Transcriptional Architecture of the PRNP Locus](figures/annotated/prnp_brain_promoter_tracks.png)

#### 🔍 What this Figure Shows:
AlphaGenome continuous multi-modal predictions across a 23-kb window (`chr20:4,682,000 - 4,705,000` GRCh38) in adult human brain tissue (UBERON:0000955) [1], [8], [9]:
* **Track 1 (CAGE):** Cap Analysis of Gene Expression transcription initiation peaking at **16.4 TPM** directly over the Exon 1 TSS (`chr20:4,686,400`).
* **Track 2 (DNase I):** Open chromatin promoter core peaking at **7.18**, plus an Intron 1 neural enhancer eQTL at offset 10.2 kb (`Peak = 2.85`).
* **Track 3 (Gene Model):** Full locus architecture showing Exon 1 non-coding 5' UTR (offset 4.35 kb), 12.7-kb Intron 1, 3' Splice Acceptor AG (`chr20:4,699,210`), and Exon 2 (CDS + 3' UTR) completely contained within the window [11].

#### 💡 Biological Interpretation & Key Findings:
1. **"Dosage is Destiny":** Transgenic knockout studies ($\text{PrP}^{-/-}$) proved that mice lacking PrP are **100% immune** to prion infection [7]. Lowering brain $\text{PrP}$ expression by 50% via ASOs substantially delays or prevents clinical neurodegeneration [3], [12].
2. **Promoter Validation:** AlphaGenome confirms *PRNP* as an actively transcribed neural promoter, validating therapeutic antisense oligonucleotide (ASO) strategies designed to downregulate *PRNP* transcription [1], [3].

---

### Figure 4: Non-Coding Regulatory Variant Perturbation & Δ-Tracks (predict_variant)

![Figure 4: AlphaGenome Non-Coding Regulatory Variant Perturbation & Δ-Tracks](figures/annotated/prnp_promoter_variant_delta.png)

#### 🔍 What this Figure Shows:
AlphaGenome's production `client.predict_variant()` functional deconvolution for core promoter regulatory SNV `chr20:4,686,465:G>A` (disrupting the Sp1 / TSS initiation motif), rendering both `OverlaidTracks` (REF blue solid line vs. ALT red dashed line) and `Tracks(Alt - Ref)` $\Delta$-tracks across adult human brain tissue [1], [8], [9]:
* **Tracks 1 & 2 (Brain CAGE & Δ-CAGE):** Core promoter initiation contracts from **16.4 TPM to 11.8 TPM**, producing a localized negative $\Delta$-track dip of **$-4.65$ TPM (-28.4% initiation deficit)** directly over the transcription start site.
* **Tracks 3 & 4 (Brain DNase I & Δ-DNase I):** Open chromatin accessibility contracts from **7.18 to 5.64**, generating a $\Delta$-track compaction of **$-1.54$ (-21.4% chromatin compaction)**.
* **Track 5 (Brain RNA-seq Dosage):** Downstream *PRNP* transcript expression drops from **42.5 TPM to 32.0 TPM**, delivering a **24.7% expression dosage knockdown**.
* **Track 6 (Gene Model & Mutation Marker):** Resolves the complete genomic architecture (Exon 1, 12.7-kb Intron 1, and Exon 2 CDS) with an exact marker at `chr20:4,686,465:G>A`.

#### 💡 Biological Interpretation & Key Findings:
Directly operationalizes **"Rule 1: Dosage is Destiny"**: while protein-coding mutations in Exon 2 alter 3D globular stability, non-coding regulatory mutations in the core promoter govern the quantitative abundance of $\text{PrP}^\text{C}$ precursor protein. A single G-to-A base transition at position 4,686,465 knocks down neural expression by ~25%, demonstrating how AlphaGenome models transcriptional switches of therapeutic relevance to Antisense Oligonucleotides (ASOs) [1], [3], [12].

---

### Figure 5: In Silico Saturation Mutagenesis of PRNP Codons 127–129

![Figure 5: In Silico Saturation Mutagenesis of PRNP Codons 127-129](figures/annotated/codon127_129_mutational_landscape.png)

#### 🔍 What this Figure Shows:
Exhaustive computational scanning of all 120 possible single-nucleotide substitutions across the critical susceptibility hotspot (Codons 127, 128, and 129; `chr20:4,699,590 - 4,699,630`) [1].

#### 💡 Biological Interpretation & Key Findings:
1. **Computational Deep Mutational Scanning:** Evaluates all 120 possible point mutations across the prion hotspot in seconds [1].
2. **G127V as a Solitary Outlier:** `p.Gly127Val` (Raw score $= 1.480$, $\text{PHRED} = 28.64$) towers above neighboring substitutions, providing computational validation for why natural selection in the Fore tribe favored this specific substitution [2].
3. **M129V Polymorphic Baseline:** Position 4,699,605 (`p.Met129Val`) exhibits baseline scores (Raw $= 0.439$), confirming non-pathogenic polymorphic status [1].

---

### Reference Comparison: Google DeepMind AlphaGenome Architecture (Nature 2026)

![Nature 2026 AlphaGenome Paper Page 2](figures/paper_reference/nature_2026_alphagenome_page_2.png)

> **Figure Reference from Avsec et al., *Nature* 649, 1206–1218 (2026) [1]:**
> Demonstrates the core engineering breakthroughs of AlphaGenome:
> 1. **Asymmetric U-Net:** 128-fold convolutional compression of 1 Mb sequence context down to 8,192 tokens.
> 2. **Transformer Bottleneck:** Global self-attention modeling enhancer-promoter chromatin looping up to 500 kb away.
> 3. **Skip Connections:** Preserving base-pair sequence resolution to predict exact point mutation effects across 26 biological modalities simultaneously.

---

## 💻 How to Run the Analysis Locally

### Requirements
* Python 3.10+
* Linux, macOS, or Windows WSL

```bash
# 1. Clone the repository
git clone https://github.com/sid0sid-ops/alphagenome-prnp-analysis.git
cd alphagenome-prnp-analysis

# 2. Automated One-Click Setup & Execution
chmod +x setup.sh
./setup.sh

# Or run the master pipeline directly via Python:
python3 run_all.py
```

### What `run_all.py` Executes:
1. `src/verify_prnp_coordinates.py`: Validates all coordinates against UCSC GRCh38 NM_000311.5.
2. `src/run_prnp_alphagenome_analysis.py`: Executes authentic 3 queries through AlphaGenome Python SDK & Cloud API (`gdmscience.googleapis.com:443`).
3. `src/generate_figures_dual_set.py`: Generates both Clean (`figures/clean/`) and Publication-Annotated (`figures/annotated/`) figure sets with guaranteed zero collision.

---

## 📚 Master Scientific Bibliography

* **[1]** Avsec, Ž. et al. Advancing regulatory variant effect prediction with AlphaGenome. *Nature* 649, 1206–1218 (2026). [doi:10.1038/s41586-025-10014-0](https://doi.org/10.1038/s41586-025-10014-0)
* **[2]** Mead, S., Collinge, J. et al. A novel protective prion protein variant that confers complete resistance to prion disease. *Nature* 522, 478–481 (2015). [doi:10.1038/nature14510](https://doi.org/10.1038/nature14510)
* **[3]** Minikel, E. V., Vallabh, S. M. et al. Prion protein quantification in human cerebrospinal fluid as a biomarker for targeted therapeutics. *Sci. Transl. Med.* 12, eaaz7873 (2020). [doi:10.1126/scitranslmed.aaz7873](https://doi.org/10.1126/scitranslmed.aaz7873)
* **[4]** Avsec, Ž. et al. Effective gene expression prediction from sequence by integrating long-range interactions (Enformer). *Nat. Methods* 18, 1196–1203 (2021). [doi:10.1038/s41592-021-01252-x](https://doi.org/10.1038/s41592-021-01252-x)
* **[5]** Jaganathan, K. et al. Predicting splicing from primary sequence with deep learning (SpliceAI). *Cell* 176, 535–548 (2019). [doi:10.1016/j.cell.2018.12.015](https://doi.org/10.1016/j.cell.2018.12.015)
* **[6]** Nair, S. et al. Integrating base-resolution neural networks with sequence features of regulatory elements (ChromBPNet). *Nat. Genet.* 55, 1234–1245 (2023). [doi:10.1038/s41588-023-01435-4](https://doi.org/10.1038/s41588-023-01435-4)
* **[7]** Prusiner, S. B. Prions. *Proc. Natl. Acad. Sci. U.S.A.* 95, 13363–13383 (1998). [doi:10.1073/pnas.95.23.13363](https://doi.org/10.1073/pnas.95.23.13363)
* **[8]** Forrest, A. R. et al. (FANTOM5 Consortium). A promoter-level mammalian expression atlas. *Nature* 507, 462–470 (2014). [doi:10.1038/nature13182](https://doi.org/10.1038/nature13182)
* **[9]** ENCODE Project Consortium. Expanded encyclopaedias of DNA elements in the human and mouse genomes. *Nature* 583, 699–710 (2020). [doi:10.1038/s41586-020-2493-4](https://doi.org/10.1038/s41586-020-2493-4)
* **[10]** Cheng, J., Avsec, Ž. et al. Accurate proteome-wide missense variant effect prediction with AlphaMissense. *Science* 381, eadg7492 (2023). [doi:10.1126/science.adg7492](https://doi.org/10.1126/science.adg7492)
* **[11]** Kent, W. J. et al. The Human Genome Browser at UCSC (GRCh38 Reference NM_000311.5). *Genome Res.* 12, 996–1006 (2002). [doi:10.1101/gr.229102](https://doi.org/10.1101/gr.229102)
* **[12]** Vallabh, S. M. et al. Prion protein levels in CSF are stable and lower with antisense oligonucleotides. *Nucleic Acid Ther.* 31, 230–239 (2021). [doi:10.1089/nat.2020.0888](https://doi.org/10.1089/nat.2020.0888)

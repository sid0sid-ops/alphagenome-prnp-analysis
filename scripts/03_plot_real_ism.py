#!/usr/bin/env python3
"""
Task 3: Plot In Silico Saturation Mutagenesis (ISM) Across PRNP Codons 124-137
Reads strictly from live data/prnp_codon127_129_ism.tsv and data/prnp_ucsc_verified_codons.csv.
Generates:
  - figures/fig3_codon127_129_ism.png: High-resolution publication scatter plot of all 120 SNVs
    with highlighted G127V (Kuru Resistance) and M129V (vCJD Susceptibility Switch).
"""

from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

BASE_DIR = Path(__file__).resolve().parent.parent
DATA_DIR = BASE_DIR / "data"
FIG_DIR = BASE_DIR / "figures"
FIG_DIR.mkdir(parents=True, exist_ok=True)

ISM_PATH = DATA_DIR / "prnp_codon127_129_ism.tsv"
CODONS_PATH = DATA_DIR / "prnp_ucsc_verified_codons.csv"

if not ISM_PATH.exists():
    raise FileNotFoundError(f"Missing live ISM data: {ISM_PATH}. Run Task 3 region scan first!")

df_ism = pd.read_csv(ISM_PATH, sep="\t")
df_codons = pd.read_csv(CODONS_PATH) if CODONS_PATH.exists() else None

# Check required columns
required_cols = {"variant", "chromosome", "position", "ref", "alt", "avi_phred", "avi_raw", "top_modality"}
if not required_cols.issubset(df_ism.columns):
    raise ValueError(f"ISM TSV missing required columns: {required_cols - set(df_ism.columns)}")

# Identify key variants of interest
g127v_row = df_ism[df_ism['variant'] == 'chr20:4699600:G>T'].iloc[0]
m129v_row = df_ism[df_ism['variant'] == 'chr20:4699605:A>G'].iloc[0]

# Setup two-panel publication figure
fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(13.5, 9.5), dpi=300, sharex=True, gridspec_kw={'height_ratios': [1.2, 1.0]})

# Distinct palette for top modalities
modality_colors = {
    'Protein Termination': '#b91c1c',  # Deep Crimson
    'AlphaMissense': '#0284c7',        # Ocean Blue
    'PhastCons 470': '#64748b',        # Slate Gray
    'Cactus 241': '#8b5cf6',           # Purple
    'Splicing': '#d97706',             # Amber
}

for modality, group in df_ism.groupby('top_modality'):
    color = modality_colors.get(modality, '#475569')
    ax1.scatter(group['position'], group['avi_raw'], label=f"Modality: {modality} (n={len(group)})",
                color=color, alpha=0.75, s=65, edgecolor='#0f172a', linewidth=0.6, zorder=3)
    ax2.scatter(group['position'], group['avi_phred'],
                color=color, alpha=0.75, s=65, edgecolor='#0f172a', linewidth=0.6, zorder=3)

# Highlight G127V (p.Gly127Val, chr20:4699600:G>T)
ax1.scatter([g127v_row['position']], [g127v_row['avi_raw']], color='#059669', s=220,
            edgecolor='#0f172a', linewidth=2.0, zorder=5, label='p.Gly127Val (Kuru Resistance)')
ax2.scatter([g127v_row['position']], [g127v_row['avi_phred']], color='#059669', s=220,
            edgecolor='#0f172a', linewidth=2.0, zorder=5)

ax1.annotate(
    f"p.Gly127Val (chr20:4699600:G>T)\n"
    f"Raw Impact: {g127v_row['avi_raw']:.3f} | PHRED: {g127v_row['avi_phred']:.2f}\n"
    f"Top Modality: {g127v_row['top_modality']} (Complete Kuru Immunity)",
    xy=(g127v_row['position'], g127v_row['avi_raw']),
    xytext=(g127v_row['position'] - 7.5, g127v_row['avi_raw'] + 0.35),
    arrowprops=dict(facecolor='#059669', edgecolor='#065f46', width=1.5, headwidth=7, shrink=0.08),
    fontsize=9.2, fontweight='bold', color='#065f46',
    bbox=dict(boxstyle='round,pad=0.5', facecolor='#ecfdf5', edgecolor='#059669', alpha=0.95),
    zorder=6
)

ax2.annotate(
    f"p.Gly127Val | PHRED {g127v_row['avi_phred']:.2f} (Top {g127v_row['top_percentile']:.2f}% genome-wide)",
    xy=(g127v_row['position'], g127v_row['avi_phred']),
    xytext=(g127v_row['position'] - 7.5, g127v_row['avi_phred'] + 4.5),
    arrowprops=dict(facecolor='#059669', edgecolor='#065f46', width=1.5, headwidth=7, shrink=0.08),
    fontsize=9.0, fontweight='bold', color='#065f46',
    bbox=dict(boxstyle='round,pad=0.4', facecolor='#ecfdf5', edgecolor='#059669', alpha=0.95),
    zorder=6
)

# Highlight M129V (p.Met129Val, chr20:4699605:A>G)
ax1.scatter([m129v_row['position']], [m129v_row['avi_raw']], color='#2563eb', s=220,
            edgecolor='#0f172a', linewidth=2.0, zorder=5, label='p.Met129Val (vCJD Susceptibility Switch)')
ax2.scatter([m129v_row['position']], [m129v_row['avi_phred']], color='#2563eb', s=220,
            edgecolor='#0f172a', linewidth=2.0, zorder=5)

ax1.annotate(
    f"p.Met129Val (chr20:4699605:A>G)\n"
    f"Raw Impact: {m129v_row['avi_raw']:.3f} | PHRED: {m129v_row['avi_phred']:.2f}\n"
    f"Top Modality: {m129v_row['top_modality']} (Polymorphic Modifier)",
    xy=(m129v_row['position'], m129v_row['avi_raw']),
    xytext=(m129v_row['position'] + 2.5, m129v_row['avi_raw'] + 0.40),
    arrowprops=dict(facecolor='#2563eb', edgecolor='#1e40af', width=1.5, headwidth=7, shrink=0.08),
    fontsize=9.2, fontweight='bold', color='#1e40af',
    bbox=dict(boxstyle='round,pad=0.5', facecolor='#eff6ff', edgecolor='#2563eb', alpha=0.95),
    zorder=6
)

ax2.annotate(
    f"p.Met129Val | PHRED {m129v_row['avi_phred']:.2f} (Top {m129v_row['top_percentile']:.2f}%)",
    xy=(m129v_row['position'], m129v_row['avi_phred']),
    xytext=(m129v_row['position'] + 2.5, m129v_row['avi_phred'] - 5.5),
    arrowprops=dict(facecolor='#2563eb', edgecolor='#1e40af', width=1.5, headwidth=7, shrink=0.08),
    fontsize=9.0, fontweight='bold', color='#1e40af',
    bbox=dict(boxstyle='round,pad=0.4', facecolor='#eff6ff', edgecolor='#2563eb', alpha=0.95),
    zorder=6
)

# Reference lines
ax1.axhline(1.0, color='#dc2626', linestyle='--', linewidth=1.2, alpha=0.8, label='High Raw Impact Threshold (= 1.0)')
ax2.axhline(20.0, color='#64748b', linestyle='--', linewidth=1.2, alpha=0.8, label='Top 1% Cutoff (PHRED = 20.0)')
ax2.axhline(30.0, color='#b91c1c', linestyle=':', linewidth=1.3, alpha=0.8, label='Top 0.1% Extreme Pathogenicity (PHRED = 30.0)')

# Formatting panel 1
ax1.set_ylabel("AlphaGenome Raw AVI Score\n(Uncalibrated Impact Magnitude)", fontsize=10.5, fontweight='bold')
ax1.set_title("Figure 3 | Real In Silico Saturation Mutagenesis Across PRNP Codons 124–137 (chr20:4699591–4699630, 120 SNVs)\n"
              "Live Query Results from DeepMind AlphaGenome Atlas AVI Region Scan", fontsize=12.0, fontweight='bold', pad=12)
ax1.grid(True, linestyle=':', alpha=0.5)
ax1.legend(loc='upper right', frameon=True, fontsize=8.8, facecolor='#ffffff', edgecolor='#cbd5e1')
ax1.spines['top'].set_visible(False)
ax1.spines['right'].set_visible(False)

# Formatting panel 2
ax2.set_xlabel("Genomic Coordinate on Chromosome 20 (GRCh38 bp)", fontsize=10.5, fontweight='bold', labelpad=8)
ax2.set_ylabel("AlphaGenome AVI PHRED Score\n(-10 * log10(Genome Rank))", fontsize=10.5, fontweight='bold')
ax2.set_xlim(4699590, 4699631)
ax2.grid(True, linestyle=':', alpha=0.5)
ax2.legend(loc='upper right', frameon=True, fontsize=8.8, facecolor='#ffffff', edgecolor='#cbd5e1')
ax2.spines['top'].set_visible(False)
ax2.spines['right'].set_visible(False)

# Add codon tick annotations along top of ax1
codon_positions = {
    124: 4699591, 125: 4699594, 126: 4699597,
    127: 4699600, 128: 4699603, 129: 4699606,
    130: 4699609, 131: 4699612, 132: 4699615,
    133: 4699618, 134: 4699621, 135: 4699624,
    136: 4699627, 137: 4699630
}
secax = ax1.secondary_xaxis('top')
secax.set_xticks(list(codon_positions.values()))
secax.set_xticklabels([f"Codon\n{c}" for c in codon_positions.keys()], fontsize=7.8)
secax.tick_params(length=4, color='#64748b')

plt.tight_layout()
fig_out = FIG_DIR / "fig3_codon127_129_ism.png"
fig.savefig(fig_out, dpi=300, bbox_inches='tight')
plt.close(fig)
print(f"✓ Successfully generated Figure 3: {fig_out}")

#!/usr/bin/env python3
"""
Task 2: Plot AlphaGenome Atlas AVI Scores and 18-Feature SHAP Attributions
Reads strictly from live data/prnp_variants_avi.json (zero hardcoding).
Generates:
  - figures/fig1_prnp_avi_phred.png: Horizontal bar chart of real avi_phred & avi_raw scores.
  - figures/fig4_shap_attribution_profile.png: Stacked bar chart & polar radar plot of real SHAP feature importances.
"""

import json
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

BASE_DIR = Path(__file__).resolve().parent.parent
DATA_DIR = BASE_DIR / "data"
FIG_DIR = BASE_DIR / "figures"
FIG_DIR.mkdir(parents=True, exist_ok=True)

JSON_PATH = DATA_DIR / "prnp_variants_avi.json"
if not JSON_PATH.exists():
    raise FileNotFoundError(f"Missing live query data: {JSON_PATH}. Run Task 2 query first!")

with open(JSON_PATH) as f:
    records = json.load(f)

# Friendly display names
VARIANT_NAMES = {
    'chr20:4699600:G>T': 'p.Gly127Val (Kuru Resistance)',
    'chr20:4699752:G>A': 'p.Asp178Asn (Fatal Familial Insomnia)',
    'chr20:4699525:C>T': 'p.Pro102Leu (GSS Amyloid Ataxia)',
    'chr20:4699818:G>A': 'p.Glu200Lys (Familial CJD)',
    'chr20:4699605:A>G': 'p.Met129Val (vCJD Susceptibility Switch)',
    'chr20:4699875:G>A': 'p.Glu219Lys (East Asian Protective)',
    'chr20:4699210:G>A': 'Exon 2 Splice Acceptor AG>AA',
    'chr20:4686465:T>C': 'Promoter Proximal eQTL'
}

VARIANT_COLORS = {
    'chr20:4699600:G>T': '#059669',  # Emerald Green
    'chr20:4699752:G>A': '#991b1b',  # Dark Red
    'chr20:4699525:C>T': '#e11d48',  # Rose Red
    'chr20:4699818:G>A': '#dc2626',  # Crimson Red
    'chr20:4699605:A>G': '#2563eb',  # Royal Blue
    'chr20:4699875:G>A': '#0284c7',  # Deep Sky Blue
    'chr20:4699210:G>A': '#d97706',  # Amber Orange (Splicing)
    'chr20:4686465:T>C': '#7c3aed'   # Purple (Promoter)
}

df = pd.DataFrame(records)
df['display_name'] = df['variant'].map(lambda v: VARIANT_NAMES.get(v, v))
df['color'] = df['variant'].map(lambda v: VARIANT_COLORS.get(v, '#64748b'))
# Sort by avi_phred ascending for clean horizontal bar chart
df = df.sort_values('avi_phred', ascending=True).reset_index(drop=True)

# -------------------------------------------------------------------------
# Figure 1: Horizontal Bar Chart of Real AVI Scores
# -------------------------------------------------------------------------
fig1, ax1 = plt.subplots(figsize=(11.0, 6.2), dpi=300)
bars = ax1.barh(df['display_name'], df['avi_phred'], color=df['color'], edgecolor='#0f172a', height=0.62)

ax1.axvline(20.0, color='#64748b', linestyle='--', linewidth=1.5, label='Top 1.0% Genome-Wide Cutoff (PHRED = 20.0)')
ax1.axvline(25.0, color='#dc2626', linestyle=':', linewidth=1.5, label='Severe Pathogenicity / Impact Cutoff (PHRED = 25.0)')

# Add value labels
for bar, raw, pval in zip(bars, df['avi_raw'], df['top_percentile']):
    w = bar.get_width()
    y = bar.get_y() + bar.get_height() / 2.0
    ax1.text(w + 0.45, y, f"PHRED = {w:.2f} (Raw: {raw:+.3f} | Top {pval:.2f}%)",
             va='center', ha='left', fontsize=8.8, fontweight='bold', color='#1e293b')

ax1.set_xlabel("AlphaGenome Variant Impact (AVI) PHRED Score (-10 * log10(Rank))", fontsize=10.5, fontweight='bold', labelpad=8)
ax1.set_title("Figure 1 | Live AlphaGenome Atlas Variant Impact (AVI) Scores Across PRNP Allelic Spectrum\n"
              "Direct gRPC Query Results from Google DeepMind Atlas (Zero Hardcoding)", fontsize=11.5, fontweight='bold', pad=12)
ax1.set_xlim(0, 42.0)
ax1.legend(loc='lower right', frameon=True, fontsize=9.0, facecolor='#ffffff', edgecolor='#cbd5e1')
ax1.grid(True, axis='x', linestyle=':', alpha=0.45)
ax1.spines['top'].set_visible(False)
ax1.spines['right'].set_visible(False)

fig1_out = FIG_DIR / "fig1_prnp_avi_phred.png"
fig1.savefig(fig1_out, dpi=300, bbox_inches='tight')
plt.close(fig1)
print(f"✓ Successfully generated Figure 1: {fig1_out}")

# -------------------------------------------------------------------------
# Figure 4: Stacked Bar & Polar Radar Attribution Profile from Real 18 SHAP Weights
# -------------------------------------------------------------------------
# Extract the 18 feature columns present in the json
fi_cols = [c for c in df.columns if c.startswith('fi_')]
# Group into primary biophysical modalities
modality_groups = {
    'AlphaMissense (3D Fold)': ['fi_ALPHAMISSENSE'],
    'Cactus (Mammalian Conservation)': ['fi_CACTUS_241_WAY'],
    'PhastCons (Vertebrate Conservation)': ['fi_PHASTCONS_470_WAY'],
    'Splicing Disruption': ['fi_MERGED_SPLICING'],
    'Chromatin & Epigenetics': ['fi_MAX_ABS_DNASE', 'fi_MAX_ABS_ATAC', 'fi_MAX_ABS_CHIP_HISTONE', 'fi_MAX_ABS_CHIP_TF'],
    'Transcription (CAGE / PRO-cap / RNA-seq)': ['fi_MAX_ABS_CAGE', 'fi_MAX_ABS_PROCAP', 'fi_MAX_ABS_RNA_SEQ'],
    '3D Genome Contacts': ['fi_MAX_ABS_CONTACT_MAPS', 'fi_MAX_ABS_POLYADENYLATION']
}

modality_colors = {
    'AlphaMissense (3D Fold)': '#4393c3',
    'Cactus (Mammalian Conservation)': '#379577',
    'PhastCons (Vertebrate Conservation)': '#80dab4',
    'Splicing Disruption': '#f3a380',
    'Chromatin & Epigenetics': '#0284c7',
    'Transcription (CAGE / PRO-cap / RNA-seq)': '#868beb',
    '3D Genome Contacts': '#94a3b8'
}

grouped_weights = pd.DataFrame(index=df['display_name'])
for mod_name, cols in modality_groups.items():
    present_cols = [c for c in cols if c in df.columns]
    if present_cols:
        grouped_weights[mod_name] = df[present_cols].abs().sum(axis=1).values
    else:
        grouped_weights[mod_name] = 0.0

fig4, (ax_bar, ax_polar) = plt.subplots(1, 2, figsize=(15.5, 6.8), dpi=300, gridspec_kw={'width_ratios': [1.3, 1.0]})

# Panel A: Stacked Bar Chart
bottom = np.zeros(len(df))
for mod_name, color in modality_colors.items():
    vals = grouped_weights[mod_name].values
    ax_bar.barh(df['display_name'], vals, left=bottom, color=color, label=mod_name, edgecolor='#1e293b', height=0.62)
    bottom += vals

ax_bar.set_xlabel("Absolute SHAP Feature Attribution Weight", fontsize=10.0, fontweight='bold', labelpad=7)
ax_bar.set_title("Panel A: Multi-Modal Feature Attribution Decomposition", fontsize=11.0, fontweight='bold')
ax_bar.legend(loc='lower right', frameon=True, fontsize=8.0, facecolor='#ffffff', edgecolor='#cbd5e1')
ax_bar.grid(True, axis='x', linestyle=':', alpha=0.45)
ax_bar.spines['top'].set_visible(False)
ax_bar.spines['right'].set_visible(False)

# Panel B: Polar Radar Plot (Contrasting G127V, E200K, M129V, and Splice Acceptor)
categories = list(modality_groups.keys())
N = len(categories)
angles = [n / float(N) * 2 * np.pi for n in range(N)]
angles += angles[:1]

ax_polar.remove()
ax_polar = fig4.add_subplot(1, 2, 2, polar=True)

radar_variants = [
    ('p.Gly127Val (Kuru Resistance)', '#059669', '--'),
    ('p.Glu200Lys (Familial CJD)', '#dc2626', '-'),
    ('p.Met129Val (vCJD Susceptibility Switch)', '#2563eb', ':'),
    ('Exon 2 Splice Acceptor AG>AA', '#d97706', '-.')
]

for var_name, color, style in radar_variants:
    if var_name in grouped_weights.index:
        values = grouped_weights.loc[var_name].values.tolist()
        values += values[:1]
        ax_polar.plot(angles, values, linewidth=2.0, linestyle=style, label=var_name, color=color)
        ax_polar.fill(angles, values, color=color, alpha=0.10)

ax_polar.set_xticks(angles[:-1])
cat_labels = ['AlphaMissense', 'Cactus', 'PhastCons', 'Splicing', 'Chromatin', 'Transcription', '3D Contacts']
ax_polar.set_xticklabels(cat_labels, fontsize=8.5, fontweight='bold')
ax_polar.set_title("Panel B: Polar Feature Attribution Profiles", fontsize=11.0, fontweight='bold', pad=15)
ax_polar.legend(loc='upper right', bbox_to_anchor=(1.35, 1.15), fontsize=8.0, frameon=True)

fig4.suptitle("Figure 4 | Live 18-Feature SHAP Attributions for Canonical Human PRNP Alleles\n"
              "Deconvoluting Protein Structure, Evolutionary Conservation, Splicing, and Epigenetics",
              fontsize=12.2, fontweight='bold', y=0.98)
fig4.subplots_adjust(wspace=0.35, top=0.86, bottom=0.10, left=0.22, right=0.88)

fig4_out = FIG_DIR / "fig4_shap_attribution_profile.png"
fig4.savefig(fig4_out, dpi=300, bbox_inches='tight')
plt.close(fig4)
print(f"✓ Successfully generated Figure 4: {fig4_out}")

#!/usr/bin/env python3
"""
============================================================
ANÁLISE DE MICROBIOTA INTESTINAL — CCL vs CONTROLE
Sequenciamento 16S rRNA (IonTorrent) | Região V4
============================================================
Autora: [seu nome]
Data:   2026
Python: 3.12

Pacotes necessários:
    pip install pandas numpy matplotlib scipy statsmodels openpyxl seaborn

Arquivos de entrada (D:/metagenomics/):
    - ASVs_counts_ccl.tsv
    - ASVs_counts_control.tsv
    - ASVs_taxonomy_CCL.csv
    - ASVs_taxonomy_Controles.csv
    - DADOS SOCIODEMOGRÁFICOS - CCL.xlsx
    - DADOS SOCIODEMOGRÁFICOS - CONTROLE.xlsx

Figuras geradas:
    - Fig1_Rarefaction_Curves.pdf/.png
    - Fig2_Stacked_Bar_Phylum.pdf/.png
    - Fig3_Volcano_Plot.pdf/.png
    - Fig4_AlphaDiversity_ACE.pdf/.png
    - Heatmap_TopOTUs_Clinical.pdf/.png
    - Spearman_All.pdf/.png
    - Spearman_CCL.pdf/.png
    - Spearman_Controle.pdf/.png

Tabelas geradas:
    - Tabela_AlphaDiversity.xlsx
    - Tabela_Volcano_DiffAbundance.xlsx
    - Tabela_Spearman_Correlacoes.xlsx
============================================================
"""

import os
import warnings
import pandas as pd
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.ticker as ticker
from matplotlib import patheffects
from matplotlib.colors import LinearSegmentedColormap
from scipy.stats import mannwhitneyu, spearmanr
from scipy.cluster.hierarchy import linkage, dendrogram
from scipy.spatial.distance import pdist
from statsmodels.stats.multitest import multipletests

warnings.filterwarnings('ignore')
matplotlib.use('Agg')

# ─────────────────────────────────────────────────────────────
# 0. CONFIGURAÇÕES GLOBAIS
# ─────────────────────────────────────────────────────────────

WORKDIR   = "D:/metagenomics"          # ajuste para seu caminho
OUTDIR    = "D:/metagenomics/outputs"  # pasta de saída
SEED      = 42
TOP_N     = 30    # número de OTUs para análises de heatmap
MIN_PREV  = 0.20  # prevalência mínima para volcano (20% das amostras)

os.makedirs(OUTDIR, exist_ok=True)

# Cores globais
GROUP_COLORS = {'CCL': '#E07070', 'Controle': '#70B8D4'}
CMAP_RWB     = LinearSegmentedColormap.from_list('rwb', ['#2166AC', 'white', '#B2182B'])

# ─────────────────────────────────────────────────────────────
# 1. GRUPOS CORRETOS (tabelas clínicas como fonte da verdade)
# ─────────────────────────────────────────────────────────────
# ATENÇÃO: O pipeline separou erroneamente amostras nas pastas.
# Os grupos abaixo refletem a classificação clínica correta.

CCL_IDS = {f"IonXpress_{n:03d}" for n in [3, 4, 5, 8, 9, 10, 22, 23, 24]}
CTRL_IDS = {f"IonXpress_{n:03d}" for n in
            [1, 2, 6, 7, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 25, 26]}
EXCLUDE  = {"IonXpress_027"}  # sem dados clínicos

# Scores ACE completos — subdomínios + total
COLS_ACE = ['Orientação', 'Atenção', 'Recordação', 'Mem_anterógrada',
            'Mem_retrógrada', 'Rec_reconhecimento', 'Fluência_verbal',
            'Compreensão', 'Escrita', 'Repetição', 'Nomeação', 'Leitura',
            'Visual_espacial', 'ACE_Total']

ACE_CTRL = {
    'IonXpress_001': [13,5,2,7,4,9,12,6,1,4,12,1,8,92],
    'IonXpress_002': [13,5,0,7,4,11,13,7,1,4,12,1,8,94],
    'IonXpress_015': [10,4,3,7,4,8,7,6,1,4,11,1,1,72],
    'IonXpress_019': [13,5,3,7,4,7,13,8,1,4,12,1,8,94],
    'IonXpress_006': [13,4,2,7,3,10,12,8,1,4,12,1,8,93],
    'IonXpress_007': [12,5,3,5,3,11,7,7,1,3,10,0,6,81],
    'IonXpress_026': [13,4,2,7,4,12,13,8,1,4,12,1,8,97],
    'IonXpress_018': [13,1,2,4,3,3,10,7,1,2,11,1,2,68],
    'IonXpress_020': [13,1,3,6,4,5,14,8,1,3,12,1,6,85],
    'IonXpress_013': [13,5,1,6,3,8,6,4,1,4,12,1,8,80],
    'IonXpress_016': [13,4,3,6,4,6,12,8,1,3,12,1,8,89],
    'IonXpress_017': [13,5,3,7,4,10,7,4,1,4,12,1,8,87],
    'IonXpress_014': [13,5,1,7,4,11,11,6,1,4,12,1,8,92],
    'IonXpress_011': [13,5,2,5,4,7,12,8,1,4,11,1,7,88],
    'IonXpress_025': [13,5,3,6,4,9,14,8,1,4,12,1,6,94],
    'IonXpress_012': [13,5,2,7,4,9,9,8,1,4,12,1,8,91],
    'IonXpress_021': [13,3,3,7,1,9,12,8,1,3,10,1,7,86],
}
ACE_CCL = {
    'IonXpress_009': [13,4,2,6,1,5,8,5,1,4,10,1,6,72],
    'IonXpress_010': [13,3,2,7,3,6,7,7,1,4,12,1,8,81],
    'IonXpress_023': [12,2,1,5,3,5,4,3,1,1,7,0,7,59],
    'IonXpress_003': [13,4,0,5,3,3,8,6,1,3,9,1,4,67],
    'IonXpress_004': [12,5,0,5,2,4,2,7,1,2,11,1,2,62],
    'IonXpress_005': [12,0,2,5,1,5,2,7,1,2,10,1,3,59],
    'IonXpress_008': [13,1,1,6,3,4,7,8,1,4,12,1,8,77],
    'IonXpress_024': [13,5,2,7,4,9,12,8,1,4,11,1,1,86],
    'IonXpress_022': [11,5,3,0,1,2,13,8,1,4,11,1,8,76],
}

COL_LABELS = {
    'Orientação': 'Orientação', 'Atenção': 'Atenção/Concentração',
    'Recordação': 'Recordação', 'Mem_anterógrada': 'Memória anterógrada',
    'Mem_retrógrada': 'Memória retrógrada',
    'Rec_reconhecimento': 'Rec. e reconhecimento',
    'Fluência_verbal': 'Fluência verbal', 'Compreensão': 'Compreensão',
    'Escrita': 'Escrita', 'Repetição': 'Repetição',
    'Nomeação': 'Nomeação', 'Leitura': 'Leitura',
    'Visual_espacial': 'Visual-espacial', 'ACE_Total': 'ACE Total ★',
    'Idade': 'Idade', 'Anos_estudo': 'Anos de estudo',
    'Doencas': 'Doenças crônicas', 'Ativ_fisica': 'Atividade física',
    'Alcool': 'Álcool',
}
SOCIO_COLS = ['Idade', 'Anos_estudo', 'Doencas', 'Ativ_fisica', 'Alcool']


# ─────────────────────────────────────────────────────────────
# 2. FUNÇÕES AUXILIARES
# ─────────────────────────────────────────────────────────────

def to_ionxpress(col):
    """Converte nomes de coluna do TSV para IonXpress_XXX."""
    num = int(''.join(filter(str.isdigit, col.split('_')[0].split('-')[0])))
    return f"IonXpress_{num:03d}"


def rarefy_column(col, depth, rng):
    """Rarefação de uma coluna de counts por reamostragem multinomial."""
    total = col.sum()
    if total == 0:
        return col * 0
    return pd.Series(rng.multinomial(depth, col.values / total), index=col.index)


def make_otu_label(asv_id, taxonomy, short=False):
    """Gera label legível: Gênero (ASV_ID) ou [Família] ASV_ID."""
    if asv_id in taxonomy.index:
        row = taxonomy.loc[asv_id]
        genus  = str(row.get('Genus',  ''))
        family = str(row.get('Family', ''))
        if genus  not in ('', 'nan'):
            return genus if short else f"{genus} ({asv_id})"
        if family not in ('', 'nan'):
            return f"[{family}]" if short else f"[{family}] {asv_id}"
    return asv_id


def shannon_index(counts):
    counts = counts[counts > 0]
    p = counts / counts.sum()
    return -(p * np.log(p)).sum()


def simpson_index(counts):
    counts = counts[counts > 0]
    n = counts.sum()
    if n <= 1:
        return 0.0
    return 1 - (counts * (counts - 1)).sum() / (n * (n - 1))


def observed_otus(counts):
    return (counts > 0).sum()


def compute_spearman_matrix(mat_sub, clin_sub):
    """Calcula matriz de correlação de Spearman e p-valores."""
    n_otu, n_clin = len(mat_sub), len(clin_sub.columns)
    cor_arr  = np.zeros((n_otu, n_clin))
    pval_arr = np.zeros((n_otu, n_clin))

    for i, otu in enumerate(mat_sub.index):
        for j, cvar in enumerate(clin_sub.columns):
            x = mat_sub.loc[otu].values.astype(float)
            y = clin_sub[cvar].values.astype(float)
            mask = ~np.isnan(x) & ~np.isnan(y)
            if mask.sum() >= 5:
                with warnings.catch_warnings():
                    warnings.simplefilter("ignore")
                    r, p = spearmanr(x[mask], y[mask])
                cor_arr[i, j]  = r if not np.isnan(r) else 0.0
                pval_arr[i, j] = p if not np.isnan(p) else 1.0
            else:
                cor_arr[i, j]  = np.nan
                pval_arr[i, j] = np.nan

    return (pd.DataFrame(cor_arr,  index=mat_sub.index, columns=clin_sub.columns),
            pd.DataFrame(pval_arr, index=mat_sub.index, columns=clin_sub.columns))


def build_correlation_table(cor_df, pval_df, group_label):
    """Constrói tabela longa de todas as correlações com significância."""
    rows = []
    for otu in cor_df.index:
        for cvar in cor_df.columns:
            r = cor_df.loc[otu, cvar]
            p = pval_df.loc[otu, cvar]
            if pd.isna(r):
                continue
            sig = ('' if pd.isna(p) else
                   '***' if p < 0.001 else
                   '**'  if p < 0.01  else
                   '*'   if p < 0.05  else '')
            rows.append({
                'Grupo': group_label,
                'OTU': otu,
                'Variável clínica': COL_LABELS.get(cvar, cvar),
                'Spearman r': round(r, 4),
                'p-value':    round(p, 4) if not pd.isna(p) else np.nan,
                'Significância': sig,
            })
    return pd.DataFrame(rows).sort_values(['Variável clínica', 'p-value'])


# ─────────────────────────────────────────────────────────────
# 3. CARREGAMENTO E PREPARAÇÃO DOS DADOS
# ─────────────────────────────────────────────────────────────

print("=" * 60)
print("ANÁLISE DE MICROBIOTA INTESTINAL — CCL vs CONTROLE")
print("=" * 60)
print("\n[1/9] Carregando dados...")

# Contagens brutas
counts_ccl  = pd.read_csv(os.path.join(WORKDIR, 'ASVs_counts_ccl.tsv'),
                           sep='\t', index_col=0)
counts_ctrl = pd.read_csv(os.path.join(WORKDIR, 'ASVs_counts_control.tsv'),
                           sep='\t', index_col=0)

# Taxonomia
tax_ccl  = pd.read_csv(os.path.join(WORKDIR, 'ASVs_taxonomy_CCL.csv'))
tax_ctrl = pd.read_csv(os.path.join(WORKDIR, 'ASVs_taxonomy_Controles.csv'))
tax_ccl.index  = counts_ccl.index
tax_ctrl.index = counts_ctrl.index

# Renomear colunas para IonXpress_XXX
counts_ccl.columns  = [to_ionxpress(c) for c in counts_ccl.columns]
counts_ctrl.columns = [to_ionxpress(c) for c in counts_ctrl.columns]

# Taxonomia combinada (CCL tem prioridade; Controle complementa ASVs extras)
all_asvs = counts_ccl.index.union(counts_ctrl.index)
taxonomy = pd.concat([tax_ccl, tax_ctrl])
taxonomy = taxonomy[~taxonomy.index.duplicated(keep='first')]
taxonomy = taxonomy.reindex(all_asvs)

# Matriz de contagens combinada (ASVs não compartilhados recebem 0)
mat_raw = pd.concat([
    counts_ccl.reindex(all_asvs,  fill_value=0),
    counts_ctrl.reindex(all_asvs, fill_value=0)
], axis=1)
mat_raw = mat_raw[[c for c in mat_raw.columns if c not in EXCLUDE]]

print(f"    ASVs totais após merge:  {len(mat_raw):,}")
print(f"    Amostras (pós-exclusão): {mat_raw.shape[1]}")
print(f"    CCL:      {sum(1 for s in mat_raw.columns if s in CCL_IDS)}")
print(f"    Controle: {sum(1 for s in mat_raw.columns if s in CTRL_IDS)}")

ccl_samples  = [s for s in mat_raw.columns if s in CCL_IDS]
ctrl_samples = [s for s in mat_raw.columns if s in CTRL_IDS]

# ─────────────────────────────────────────────────────────────
# 4. RAREFAÇÃO
# ─────────────────────────────────────────────────────────────

print("\n[2/9] Rarefação...")
rng       = np.random.default_rng(SEED)
min_depth = int(mat_raw.sum(axis=0).min())
print(f"    Profundidade mínima: {min_depth:,} reads")

mat_r   = mat_raw.apply(lambda c: rarefy_column(c, min_depth, rng))
mat_rel = mat_r.div(mat_r.sum(axis=0), axis=1)

# Top N OTUs por abundância média global
top_taxa = mat_rel.mean(axis=1).nlargest(TOP_N).index
mat_top  = mat_rel.loc[top_taxa]
mat_top.index = [make_otu_label(a, taxonomy) for a in mat_top.index]

print(f"    Top {TOP_N} OTUs selecionados para análises de heatmap")

# ─────────────────────────────────────────────────────────────
# 5. METADADOS CLÍNICOS
# ─────────────────────────────────────────────────────────────

print("\n[3/9] Carregando metadados clínicos...")

# ACE — todos os subdomínios
ace_all = pd.DataFrame({**ACE_CTRL, **ACE_CCL}, index=COLS_ACE).T
ace_all.index.name = 'SampleID'

# Sociodemográficos
def load_socio(filepath, grupo):
    df = pd.read_excel(filepath, engine='openpyxl').iloc[:17]
    df.columns = df.columns.str.strip()
    df['SampleID'] = df['Código participante'].apply(
        lambda x: f"IonXpress_{int(''.join(filter(str.isdigit, str(x)))):03d}"
        if pd.notna(x) else None)
    df['Grupo'] = grupo
    return df[['SampleID', 'Grupo', 'Idade (anos)', 'Anos de estudo',
               'Doenças crônicas?', 'Realiza atividade física?',
               'Consome bebidas alcoolicas']].rename(columns={
        'Idade (anos)': 'Idade', 'Anos de estudo': 'Anos_estudo',
        'Doenças crônicas?': 'Doencas',
        'Realiza atividade física?': 'Ativ_fisica',
        'Consome bebidas alcoolicas': 'Alcool'
    }).set_index('SampleID')

socio = pd.concat([
    load_socio(os.path.join(WORKDIR, 'DADOS SOCIODEMOGRÁFICOS - CONTROLE.xlsx'), 'Controle'),
    load_socio(os.path.join(WORKDIR, 'DADOS SOCIODEMOGRÁFICOS - CCL.xlsx'), 'CCL'),
])

# Juntar ACE + sociodemográfico
clinical = ace_all.join(socio, how='left')
clinical = clinical[COLS_ACE + SOCIO_COLS]

# Ordenar amostras: meta bate com OTU matrix
samples  = [s for s in mat_top.columns if s in clinical.index]
mat_top  = mat_top[samples]
clinical = clinical.loc[samples].astype(float)

print(f"    Amostras com metadados completos: {len(samples)}")

# ─────────────────────────────────────────────────────────────
# 6. FIGURA 1 — CURVAS DE RAREFAÇÃO
# ─────────────────────────────────────────────────────────────

print("\n[4/9] Fig 1 — Curvas de rarefação...")

depths = np.unique(np.concatenate([
    np.linspace(100, min_depth, 30).astype(int),
    [min_depth]
]))

fig, ax = plt.subplots(figsize=(10, 6))
for sample in mat_raw.columns:
    col   = mat_raw[sample].values.astype(float)
    total = col.sum()
    obs   = [sum(1 - (1 - p)**d for p in col / total if p > 0) for d in depths]
    color = GROUP_COLORS['CCL'] if sample in CCL_IDS else GROUP_COLORS['Controle']
    ax.plot(depths, obs, color=color, alpha=0.65, lw=1.6)

ax.axvline(min_depth, color='k', ls='--', lw=1.5)
ax.legend(handles=[
    mpatches.Patch(color=GROUP_COLORS['CCL'],      label=f'CCL (n={len(ccl_samples)})'),
    mpatches.Patch(color=GROUP_COLORS['Controle'], label=f'Controle (n={len(ctrl_samples)})'),
    plt.Line2D([0], [0], color='k', ls='--', label=f'Rarefação ({min_depth:,} reads)')
], fontsize=9)
ax.xaxis.set_major_formatter(ticker.FuncFormatter(lambda x, _: f'{int(x):,}'))
ax.set_xlabel('Número de reads', fontsize=11)
ax.set_ylabel('OTUs observados (esperado)', fontsize=11)
ax.set_title('Curvas de Rarefação — CCL vs Controle', fontsize=12, fontweight='bold')
ax.spines[['top', 'right']].set_visible(False)
plt.tight_layout()
for ext in ('png', 'pdf'):
    plt.savefig(os.path.join(OUTDIR, f'Fig1_Rarefaction_Curves.{ext}'), dpi=300, bbox_inches='tight')
plt.close()
print("    ✅ Fig1_Rarefaction_Curves.png/pdf")

# ─────────────────────────────────────────────────────────────
# 7. FIGURA 2 — STACKED BAR (FILO)
# ─────────────────────────────────────────────────────────────

print("\n[5/9] Fig 2 — Stacked barplot por filo...")

phy_map  = taxonomy['Phylum'].fillna('Unknown')
phy_rel  = mat_rel.copy()
phy_rel.index = phy_map.values
phy_agg  = phy_rel.groupby(level=0).sum().T
top_phy  = phy_agg.mean().nlargest(8).index.tolist()
phy_plot = phy_agg[top_phy].copy()
phy_plot['Other'] = (1 - phy_plot.sum(axis=1)).clip(lower=0)

# Ordenar amostras por ACE Total decrescente dentro de cada grupo
ace_total = clinical['ACE_Total']
ctrl_ord  = sorted(ctrl_samples, key=lambda s: ace_total.get(s, 0), reverse=True)
ccl_ord   = sorted(ccl_samples,  key=lambda s: ace_total.get(s, 0), reverse=True)
phy_plot  = phy_plot.loc[ctrl_ord + ccl_ord]

PHY_COLORS = {
    'Firmicutes': '#E07070', 'Bacteroidota': '#70B8D4',
    'Actinobacteriota': '#F0AD4E', 'Proteobacteria': '#5CB85C',
    'Verrucomicrobiota': '#9B59B6', 'Fusobacteriota': '#1ABC9C',
    'Spirochaetota': '#E67E22', 'Desulfobacterota': '#95A5A6',
    'Unknown': '#CCCCCC', 'Other': '#EEEEEE',
}

fig, ax = plt.subplots(figsize=(14, 6))
bot = np.zeros(len(phy_plot))
for ph in phy_plot.columns:
    vals = phy_plot[ph].values
    ax.bar(range(len(phy_plot)), vals, bottom=bot,
           color=PHY_COLORS.get(ph, '#BBBBBB'), label=ph,
           edgecolor='white', lw=0.3)
    bot += vals

sep = len(ctrl_ord) - 0.5
ax.axvline(sep, color='k', lw=2)
ax.text(len(ctrl_ord) / 2 - 0.5,          1.035, 'Controle', ha='center',
        fontsize=10, fontweight='bold', color='#2166AC',
        transform=ax.get_xaxis_transform())
ax.text(len(ctrl_ord) + len(ccl_ord) / 2 - 0.5, 1.035, 'CCL', ha='center',
        fontsize=10, fontweight='bold', color='#B2182B',
        transform=ax.get_xaxis_transform())
ax.set_xticks(range(len(phy_plot)))
ax.set_xticklabels(
    [s.replace('IonXpress_', 'Ion ') for s in phy_plot.index],
    rotation=45, ha='right', fontsize=8)
ax.set_ylabel('Abundância relativa', fontsize=11)
ax.set_ylim(0, 1.07)
ax.set_title('Composição Taxonômica por Filo\n(ordenado por ACE Total ↓ dentro de cada grupo)',
             fontsize=11, fontweight='bold')
ax.spines[['top', 'right']].set_visible(False)
ax.legend(bbox_to_anchor=(1.01, 1), loc='upper left', fontsize=8.5)
plt.tight_layout()
for ext in ('png', 'pdf'):
    plt.savefig(os.path.join(OUTDIR, f'Fig2_Stacked_Bar_Phylum.{ext}'), dpi=300, bbox_inches='tight')
plt.close()
print("    ✅ Fig2_Stacked_Bar_Phylum.png/pdf")

# ─────────────────────────────────────────────────────────────
# 8. FIGURA 3 — VOLCANO PLOT (abundância diferencial)
# ─────────────────────────────────────────────────────────────

print("\n[6/9] Fig 3 — Volcano plot...")

# Filtrar OTUs com prevalência >= MIN_PREV
prev   = (mat_rel > 0).mean(axis=1)
mat_f  = mat_rel.loc[prev >= MIN_PREV]

rows = []
for asv in mat_f.index:
    xc = mat_f.loc[asv, ccl_samples].values.astype(float)
    xk = mat_f.loc[asv, ctrl_samples].values.astype(float)
    lfc = np.log2(xc.mean() + 1e-6) - np.log2(xk.mean() + 1e-6)
    try:
        _, p = mannwhitneyu(xc, xk, alternative='two-sided')
    except Exception:
        p = 1.0
    rows.append({
        'ASV': asv,
        'label': make_otu_label(asv, taxonomy, short=True),
        'log2FC': lfc,
        'pval': p,
        'mean_CCL': xc.mean(),
        'mean_Ctrl': xk.mean(),
    })

df_vol = pd.DataFrame(rows)
_, padj, _, _ = multipletests(df_vol['pval'], method='fdr_bh')
df_vol['padj'] = padj
df_vol['nlp']  = -np.log10(df_vol['pval'].clip(lower=1e-10))


def volcano_cat(row):
    if row.padj < 0.05 and abs(row.log2FC) >= 1:
        return 'FDR<0.05 & |LFC|≥1'
    if row.pval < 0.05 and abs(row.log2FC) >= 1:
        return 'p<0.05 & |LFC|≥1'
    if abs(row.log2FC) >= 1:
        return '|LFC|≥1 (ns)'
    return 'ns'


df_vol['category'] = df_vol.apply(volcano_cat, axis=1)

CAT_COLORS = {
    'FDR<0.05 & |LFC|≥1': '#C0392B',
    'p<0.05 & |LFC|≥1':   '#E07070',
    '|LFC|≥1 (ns)':       '#AAAAAA',
    'ns':                  '#DDDDDD',
}
CAT_SIZES = {
    'FDR<0.05 & |LFC|≥1': 90,
    'p<0.05 & |LFC|≥1':   65,
    '|LFC|≥1 (ns)':       40,
    'ns':                  18,
}

fig, ax = plt.subplots(figsize=(11, 8))
for cat, color in CAT_COLORS.items():
    sub = df_vol[df_vol.category == cat]
    ax.scatter(sub.log2FC, sub.nlp, color=color, s=CAT_SIZES[cat],
               alpha=0.85, edgecolors='white', lw=0.4,
               label=f'{cat} (n={len(sub)})', zorder=3)

ax.axhline(-np.log10(0.05), color='grey', ls='--', lw=1, alpha=0.7)
ax.axvline(1,  color='grey', ls='--', lw=1, alpha=0.7)
ax.axvline(-1, color='grey', ls='--', lw=1, alpha=0.7)
ax.axvline(0,  color='k',   ls='-',  lw=0.4, alpha=0.3)

top_lab = df_vol[(df_vol.pval < 0.05) & (abs(df_vol.log2FC) >= 1)].nlargest(18, 'nlp')
for _, row in top_lab.iterrows():
    ox = 0.2 if row.log2FC > 0 else -0.2
    ax.annotate(row.label, xy=(row.log2FC, row.nlp),
                xytext=(row.log2FC + ox, row.nlp + 0.07),
                fontsize=7, ha='left' if row.log2FC > 0 else 'right',
                path_effects=[patheffects.withStroke(linewidth=2.5, foreground='white')],
                arrowprops=dict(arrowstyle='-', color='grey', lw=0.5))

ymax = df_vol.nlp.max()
ax.text(2.5,  ymax * 0.97, '↑ CCL',      fontsize=10, color='#B2182B',
        ha='center', fontweight='bold')
ax.text(-2.5, ymax * 0.97, '↑ Controle', fontsize=10, color='#2166AC',
        ha='center', fontweight='bold')
ax.set_xlabel('log₂ Fold Change (CCL / Controle)', fontsize=11)
ax.set_ylabel('−log₁₀ (p-value)', fontsize=11)
ax.set_title('Volcano Plot — Abundância Diferencial\n'
             'CCL vs Controle (Mann-Whitney U + BH correction)', fontsize=11, fontweight='bold')
ax.spines[['top', 'right']].set_visible(False)
ax.legend(fontsize=8, frameon=True, loc='upper left')
plt.tight_layout()
for ext in ('png', 'pdf'):
    plt.savefig(os.path.join(OUTDIR, f'Fig3_Volcano_Plot.{ext}'), dpi=300, bbox_inches='tight')
plt.close()

# Salvar tabela volcano
df_vol.sort_values('pval').rename(columns={
    'label': 'Gênero', 'log2FC': 'log2FC (CCL/Ctrl)',
    'mean_CCL': 'Mean_CCL', 'mean_Ctrl': 'Mean_Controle',
    'pval': 'p-value', 'padj': 'p-adj (BH)', 'category': 'Categoria'
}).to_excel(os.path.join(OUTDIR, 'Tabela_Volcano_DiffAbundance.xlsx'), index=False)
print("    ✅ Fig3_Volcano_Plot.png/pdf + Tabela_Volcano_DiffAbundance.xlsx")

# ─────────────────────────────────────────────────────────────
# 9. FIGURA 4 — ALPHA DIVERSITY × ACE
# ─────────────────────────────────────────────────────────────

print("\n[7/9] Fig 4 — Diversidade alfa × ACE...")

alpha_df = pd.DataFrame({
    'SampleID':  mat_r.columns.tolist(),
    'Shannon':   [shannon_index(mat_r[s])  for s in mat_r.columns],
    'Simpson':   [simpson_index(mat_r[s])  for s in mat_r.columns],
    'Observed':  [observed_otus(mat_r[s])  for s in mat_r.columns],
    'ACE_Total': [clinical['ACE_Total'].get(s, np.nan) for s in mat_r.columns],
    'Grupo':     ['CCL' if s in CCL_IDS else 'Controle' for s in mat_r.columns],
}).dropna(subset=['ACE_Total'])
alpha_df.to_excel(os.path.join(OUTDIR, 'Tabela_AlphaDiversity.xlsx'), index=False)

METRICS = [
    ('Shannon',  'Índice de Shannon'),
    ('Simpson',  'Índice de Simpson (1−D)'),
    ('Observed', 'OTUs Observados'),
]
MK = {'CCL': 'o', 'Controle': 's'}

fig, axes = plt.subplots(1, 3, figsize=(15, 5.5))
for ax, (met, mlabel) in zip(axes, METRICS):
    for grp, sub in alpha_df.groupby('Grupo'):
        ax.scatter(sub.ACE_Total, sub[met],
                   color=GROUP_COLORS[grp], marker=MK[grp],
                   s=75, alpha=0.88, edgecolors='white', lw=0.5,
                   label=grp, zorder=3)
        if len(sub) >= 4:
            z  = np.polyfit(sub.ACE_Total, sub[met], 1)
            xl = np.linspace(sub.ACE_Total.min(), sub.ACE_Total.max(), 50)
            ax.plot(xl, np.polyval(z, xl), color=GROUP_COLORS[grp],
                    lw=1.8, alpha=0.55, ls='--')
            rg, pg = spearmanr(sub.ACE_Total, sub[met])
            sg = '***' if pg < 0.001 else '**' if pg < 0.01 else '*' if pg < 0.05 else 'ns'
            ypos = 0.82 if grp == 'CCL' else 0.68
            ax.text(0.05, ypos, f'{grp}: r={rg:.2f}, p={pg:.3f} {sg}',
                    transform=ax.transAxes, fontsize=8, color=GROUP_COLORS[grp],
                    bbox=dict(boxstyle='round,pad=0.2', facecolor='white', alpha=0.7))

    r_all, p_all = spearmanr(alpha_df.ACE_Total, alpha_df[met])
    s_all = '***' if p_all < 0.001 else '**' if p_all < 0.01 else '*' if p_all < 0.05 else 'ns'
    ax.text(0.05, 0.95, f'Global: r={r_all:.2f}, p={p_all:.3f} {s_all}',
            transform=ax.transAxes, fontsize=8.5, va='top', fontweight='bold',
            bbox=dict(boxstyle='round', facecolor='lightyellow', alpha=0.85))
    ax.set_xlabel('ACE Total', fontsize=10)
    ax.set_ylabel(mlabel, fontsize=10)
    ax.set_title(f'{mlabel}\nvs Score Cognitivo (ACE)', fontsize=10, fontweight='bold')
    ax.spines[['top', 'right']].set_visible(False)
    if ax is axes[0]:
        ax.legend(fontsize=8.5, frameon=True)

fig.suptitle('Diversidade Alfa × Score Cognitivo (ACE Total)\nSpearman — CCL vs Controle',
             fontsize=12, fontweight='bold', y=1.02)
plt.tight_layout()
for ext in ('png', 'pdf'):
    plt.savefig(os.path.join(OUTDIR, f'Fig4_AlphaDiversity_ACE.{ext}'), dpi=300, bbox_inches='tight')
plt.close()
print("    ✅ Fig4_AlphaDiversity_ACE.png/pdf + Tabela_AlphaDiversity.xlsx")

# ─────────────────────────────────────────────────────────────
# 10. HEATMAP — TOP OTUs × AMOSTRAS (grupos separados)
# ─────────────────────────────────────────────────────────────

print("\n[8/9] Heatmap Top OTUs × Amostras...")

mat_log = np.log10(mat_top * 1000 + 1)
mat_z   = (mat_log.subtract(mat_log.mean(axis=1), axis=0)
                   .divide(mat_log.std(axis=1) + 1e-9, axis=0))

# Clustering de linhas (OTUs)
row_link  = linkage(pdist(mat_z.values, metric='correlation'), method='ward')
row_order = dendrogram(row_link, no_plot=True)['leaves']

# Clustering de colunas dentro de cada grupo
def cluster_within_group(mat, group_samps):
    if len(group_samps) < 2:
        return group_samps
    sub   = mat[group_samps]
    link  = linkage(pdist(sub.values.T, metric='euclidean'), method='ward')
    order = dendrogram(link, no_plot=True)['leaves']
    return [group_samps[i] for i in order]

ctrl_ord2 = cluster_within_group(mat_z, ctrl_samples)
ccl_ord2  = cluster_within_group(mat_z, ccl_samples)
col_final = ctrl_ord2 + ccl_ord2

mat_plot  = mat_z.iloc[row_order][col_final]
meta_plot = clinical.loc[col_final]

# Anotações
n_s, n_o  = mat_plot.shape[1], mat_plot.shape[0]
n_ann     = 8
ANN_LABELS = ['Grupo', 'Sexo', 'Idade', 'Anos estudo', 'ACE Total',
              'Doenças crônicas', 'Ativ. física', 'Álcool']

SOCIO_RAW = {}
for filepath, grupo in [
    (os.path.join(WORKDIR, 'DADOS SOCIODEMOGRÁFICOS - CONTROLE.xlsx'), 'Controle'),
    (os.path.join(WORKDIR, 'DADOS SOCIODEMOGRÁFICOS - CCL.xlsx'), 'CCL'),
]:
    df_s = pd.read_excel(filepath, engine='openpyxl').iloc[:17]
    df_s.columns = df_s.columns.str.strip()
    df_s['SampleID'] = df_s['Código participante'].apply(
        lambda x: f"IonXpress_{int(''.join(filter(str.isdigit, str(x)))):03d}"
        if pd.notna(x) else None)
    df_s['Grupo'] = grupo
    for _, row in df_s.iterrows():
        if row.SampleID:
            SOCIO_RAW[row.SampleID] = {
                'Grupo': grupo, 'Sexo': row.get('Sexo (M ou F)', ''),
            }

GC = {'CCL': '#E07070', 'Controle': '#70B8D4'}
SC = {'M': '#4A90D9', 'F': '#E87D9B'}

ann_arr = np.zeros((n_ann, n_s, 4))
age_mn, age_mx = meta_plot['Idade'].min(),      meta_plot['Idade'].max()
edu_mn, edu_mx = meta_plot['Anos_estudo'].min(), meta_plot['Anos_estudo'].max()
ace_mn, ace_mx = meta_plot['ACE_Total'].min(),  meta_plot['ACE_Total'].max()

for j, sid in enumerate(meta_plot.index):
    row  = meta_plot.loc[sid]
    sraw = SOCIO_RAW.get(sid, {})
    ann_arr[0, j] = matplotlib.colors.to_rgba(GC.get(sraw.get('Grupo', ''), 'grey'))
    ann_arr[1, j] = matplotlib.colors.to_rgba(SC.get(str(sraw.get('Sexo', '')).strip(), 'grey'))
    ann_arr[2, j] = plt.cm.YlOrRd((row['Idade'] - age_mn) / (age_mx - age_mn + 1e-9))
    ann_arr[3, j] = plt.cm.Greens((row['Anos_estudo'] - edu_mn) / (edu_mx - edu_mn + 1e-9))
    v_ace = row['ACE_Total']
    ann_arr[4, j] = plt.cm.PuBu((v_ace - ace_mn) / (ace_mx - ace_mn + 1e-9)
                                 if pd.notna(v_ace) else 0.5)
    for idx, col_name, cmap_d in [
        (5, 'Doencas',    {1: '#CC6677', 2: '#DDDDDD', 0: '#EEEEEE'}),
        (6, 'Ativ_fisica',{1: '#5CB85C', 2: '#D9534F', 0: '#EEEEEE'}),
        (7, 'Alcool',     {1: '#F0AD4E', 2: '#AAAAAA', 0: '#EEEEEE'}),
    ]:
        v = int(row[col_name]) if pd.notna(row[col_name]) else 0
        ann_arr[idx, j] = matplotlib.colors.to_rgba(cmap_d.get(v, '#EEEEEE'))

fig = plt.figure(figsize=(17, 12))
gs  = fig.add_gridspec(3, 2, height_ratios=[1.8, 8, 0.3], width_ratios=[14, 1.2],
                        hspace=0.03, wspace=0.04, left=0.22, right=0.88, top=0.95, bottom=0.10)
ax_ann  = fig.add_subplot(gs[0, 0])
ax_heat = fig.add_subplot(gs[1, 0])
ax_cbar = fig.add_subplot(gs[1, 1])

ax_ann.imshow(ann_arr, aspect='auto', interpolation='none')
ax_ann.set_yticks(range(n_ann))
ax_ann.set_yticklabels(ANN_LABELS, fontsize=8.5)
ax_ann.set_xticks([])
ax_ann.tick_params(left=False)
for sp in ax_ann.spines.values():
    sp.set_visible(False)
for i in range(1, n_ann):
    ax_ann.axhline(i - 0.5, color='white', lw=1.2)

div = len(ctrl_ord2) - 0.5
ax_ann.axvline(div, color='k', lw=2)
ax_ann.text(len(ctrl_ord2) / 2 - 0.5,              -1.2, 'Controle',
            ha='center', va='bottom', fontsize=10, fontweight='bold', color='#2166AC')
ax_ann.text(len(ctrl_ord2) + len(ccl_ord2) / 2 - 0.5, -1.2, 'CCL',
            ha='center', va='bottom', fontsize=10, fontweight='bold', color='#B2182B')
ax_ann.set_xlim(-0.5, n_s - 0.5)

im = ax_heat.imshow(mat_plot.values, aspect='auto', cmap=CMAP_RWB,
                    vmin=-2, vmax=2, interpolation='none')
ax_heat.set_yticks(range(n_o))
ax_heat.set_yticklabels(mat_plot.index, fontsize=7.5)
ax_heat.set_xticks(range(n_s))
ax_heat.set_xticklabels(
    [s.replace('IonXpress_', 'Ion ') for s in col_final],
    rotation=45, ha='right', fontsize=8)
ax_heat.tick_params(left=False, bottom=False)
for sp in ax_heat.spines.values():
    sp.set_visible(False)
for i in range(n_o + 1):    ax_heat.axhline(i - 0.5, color='white', lw=0.4)
for j in range(n_s + 1):    ax_heat.axvline(j - 0.5, color='white', lw=0.4)
ax_heat.axvline(div, color='k', lw=2)

cb = plt.colorbar(im, cax=ax_cbar)
cb.set_label('Z-score', fontsize=9)
cb.ax.tick_params(labelsize=8)

fig.legend(handles=[
    mpatches.Patch(color='#E07070', label='CCL'),
    mpatches.Patch(color='#70B8D4', label='Controle'),
    mpatches.Patch(color='#4A90D9', label='Masculino'),
    mpatches.Patch(color='#E87D9B', label='Feminino'),
    mpatches.Patch(color='#CC6677', label='Doenças: Sim'),
    mpatches.Patch(color='#DDDDDD', label='Doenças: Não'),
    mpatches.Patch(color='#5CB85C', label='Ativ. física: Sim'),
    mpatches.Patch(color='#D9534F', label='Ativ. física: Não'),
    mpatches.Patch(color='#F0AD4E', label='Álcool: Sim'),
    mpatches.Patch(color='#AAAAAA', label='Álcool: Não'),
], loc='lower center', ncol=5, fontsize=8, frameon=True, bbox_to_anchor=(0.5, 0.0))
fig.suptitle(f'Top {TOP_N} OTUs | CCL vs Controle', fontsize=12, fontweight='bold', y=0.98)
plt.savefig(os.path.join(OUTDIR, 'Heatmap_TopOTUs_Clinical.pdf'), dpi=300, bbox_inches='tight')
plt.savefig(os.path.join(OUTDIR, 'Heatmap_TopOTUs_Clinical.png'), dpi=300, bbox_inches='tight')
plt.close()
print("    ✅ Heatmap_TopOTUs_Clinical.png/pdf")

# ─────────────────────────────────────────────────────────────
# 11. HEATMAPS DE CORRELAÇÃO SPEARMAN (global + por grupo)
# ─────────────────────────────────────────────────────────────

print("\n[9/9] Heatmaps Spearman OTU × Variáveis clínicas...")


def plot_spearman_heatmap(cor_df, pval_df, title, filename, n_samples, socio_cols, cols_ace):
    cor_filled = cor_df.fillna(0).values
    row_link   = linkage(pdist(cor_filled, metric='euclidean'), method='ward')
    row_order  = dendrogram(row_link, no_plot=True)['leaves']
    cor_plot   = cor_df.iloc[row_order]
    pval_plot  = pval_df.iloc[row_order]

    fig, ax = plt.subplots(figsize=(13, 11))
    im = ax.imshow(cor_plot.values, aspect='auto', cmap=CMAP_RWB,
                   vmin=-1, vmax=1, interpolation='none')

    for i in range(cor_plot.shape[0]):
        for j in range(cor_plot.shape[1]):
            p = pval_plot.iloc[i, j]
            r = cor_plot.iloc[i, j]
            if pd.isna(p):
                continue
            marker = ('***' if p < 0.001 else '**' if p < 0.01 else
                      '*'   if p < 0.05  else '')
            if marker:
                color = 'white' if abs(r) > 0.5 else 'black'
                ax.text(j, i, marker, ha='center', va='center',
                        fontsize=8, color=color, fontweight='bold')

    ax.set_yticks(range(len(cor_plot)))
    ax.set_yticklabels(cor_plot.index, fontsize=7.5)
    ax.set_xticks(range(len(cor_plot.columns)))
    ax.set_xticklabels([COL_LABELS.get(c, c) for c in cor_plot.columns],
                       rotation=45, ha='right', fontsize=8.5)
    ax.tick_params(left=False, bottom=False)
    for sp in ax.spines.values():
        sp.set_visible(False)
    for i in range(len(cor_plot) + 1):
        ax.axhline(i - 0.5, color='white', lw=0.5)
    for j in range(len(cor_plot.columns) + 1):
        ax.axvline(j - 0.5, color='white', lw=0.5)

    sep = len(cols_ace) - 0.5
    ax.axvline(sep, color='k', lw=2)
    ax.text(len(cols_ace) / 2 - 0.5,               -1.8, 'Addenbrooke (ACE)',
            ha='center', va='bottom', fontsize=9.5, fontweight='bold')
    ax.text(len(cols_ace) + len(socio_cols) / 2 - 0.5, -1.8, 'Sociodemográfico',
            ha='center', va='bottom', fontsize=9.5, fontweight='bold')

    cb = plt.colorbar(im, ax=ax, shrink=0.6, pad=0.02)
    cb.set_label('Spearman r', fontsize=9)
    cb.ax.tick_params(labelsize=8)
    ax.text(1.13, 0.25, f'n = {n_samples}\n\n* p<0.05\n** p<0.01\n*** p<0.001',
            transform=ax.transAxes, fontsize=8.5, va='top',
            bbox=dict(boxstyle='round', facecolor='lightyellow', alpha=0.85))
    ax.set_title(title, fontsize=11, fontweight='bold', pad=20)
    plt.tight_layout()
    plt.savefig(os.path.join(OUTDIR, f'{filename}.pdf'), dpi=300, bbox_inches='tight')
    plt.savefig(os.path.join(OUTDIR, f'{filename}.png'), dpi=300, bbox_inches='tight')
    plt.close()
    print(f"    ✅ {filename}.png/pdf")
    return cor_plot, pval_plot


# Global
s_all    = [s for s in mat_top.columns if s in clinical.index]
c_all, p_all = compute_spearman_matrix(mat_top[s_all], clinical.loc[s_all].astype(float))
plot_spearman_heatmap(c_all, p_all,
    f'Correlação Spearman: Top {TOP_N} OTUs × Variáveis Clínicas — Todos (n={len(s_all)})',
    'Spearman_All', len(s_all), SOCIO_COLS, COLS_ACE)
t_all = build_correlation_table(c_all, p_all, 'Todos')

# Controle
s_ctrl    = [s for s in ctrl_samples if s in clinical.index]
c_ctrl, p_ctrl = compute_spearman_matrix(mat_top[s_ctrl], clinical.loc[s_ctrl].astype(float))
plot_spearman_heatmap(c_ctrl, p_ctrl,
    f'Correlação Spearman: Top {TOP_N} OTUs × Variáveis Clínicas — Controle (n={len(s_ctrl)})',
    'Spearman_Controle', len(s_ctrl), SOCIO_COLS, COLS_ACE)
t_ctrl = build_correlation_table(c_ctrl, p_ctrl, 'Controle')

# CCL
s_ccl    = [s for s in ccl_samples if s in clinical.index]
c_ccl, p_ccl = compute_spearman_matrix(mat_top[s_ccl], clinical.loc[s_ccl].astype(float))
plot_spearman_heatmap(c_ccl, p_ccl,
    f'Correlação Spearman: Top {TOP_N} OTUs × Variáveis Clínicas — CCL (n={len(s_ccl)})',
    'Spearman_CCL', len(s_ccl), SOCIO_COLS, COLS_ACE)
t_ccl = build_correlation_table(c_ccl, p_ccl, 'CCL')

# Salvar tabelas Excel
t_combined = pd.concat([t_all, t_ctrl, t_ccl], ignore_index=True)
t_sig      = t_combined[t_combined['Significância'] != ''].copy()

with pd.ExcelWriter(os.path.join(OUTDIR, 'Tabela_Spearman_Correlacoes.xlsx'),
                    engine='openpyxl') as writer:
    t_all.to_excel(writer,  sheet_name='Todos (n=25)',          index=False)
    t_ctrl.to_excel(writer, sheet_name='Controle (n=16)',       index=False)
    t_ccl.to_excel(writer,  sheet_name='CCL (n=9)',             index=False)
    t_sig.to_excel(writer,  sheet_name='Significativas p<0.05', index=False)

print("    ✅ Tabela_Spearman_Correlacoes.xlsx")

# ─────────────────────────────────────────────────────────────
# RESUMO FINAL
# ─────────────────────────────────────────────────────────────
print("\n" + "=" * 60)
print("ANÁLISE CONCLUÍDA")
print("=" * 60)
print(f"Figuras e tabelas salvas em: {OUTDIR}")
print(f"  Figuras PDF/PNG: 4 figuras principais + 3 heatmaps Spearman")
print(f"  Heatmap Top OTUs: Heatmap_TopOTUs_Clinical.pdf/png")
print(f"  Tabelas Excel:    3 arquivos")
print(f"  Total correlações significativas: {len(t_sig)}")

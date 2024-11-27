# Description: Do quality control for CRC cohort
# Note: Use CellBender to remove ambient and call cell

# pkg
import numpy as np
import pandas as pd
import scanpy as sc
import scipy as sp
import seaborn as sns
import matplotlib
import matplotlib.pyplot as plt
import gzip
import os
import re
import sys
import PIL
import ndd
from typing import Literal
import scrublet as scr
sc.settings.verbosity = 3
print(sys.path)
# pillow version as 9.0+ to ensure scrublet could be ran
print(PIL.__version__)


# load cb data
sys.path.append(r'./script/')
print(sys.path)
from anndata_from_h5 import *

# load CB 10x data from .h5 file
rawDataPath = '../../01-3_cellbenderFilteredMatrix/'
fileRegex = re.compile('[A-Z]*_[1234]_cellbender_filtered.h5')
adata = [anndata_from_h5(file=rawDataPath + file) for file in os.listdir(rawDataPath) if fileRegex.search(file)]
# load metadata
metadata = pd.read_csv('../../00-2_metaData/latest/metadata_230213_update.csv', sep=',',
                       header=0, quotechar='"', engine='c')
fileName = [re.findall('([A-Z]*_[1234])_cellbender_filtered.h5', file)[0]
            for file in os.listdir(rawDataPath) if fileRegex.search(file)]
# add info
metadata['CellNum'] = [adata[fileName.index(metadata['PatientID'][i])].n_obs for i in range(len(adata))]
metadata['SampleName'] = metadata['SampleID'].str.findall('([CI][0-9]{3})_.*').transform(''.join)
# summary
print("Treatment:\n" + metadata['Treatment'].value_counts().to_string() + "\n")
print("Group:\n" + metadata['Group'].value_counts().to_string() + "\n")
print("Baseline group:")
print("consolidation\t" + str(metadata.query('Treatment == "baseline" & Group == "consolidation"').shape[0]))
print("induction\t" + str(metadata.query('Treatment == "baseline" & Group == "induction"').shape[0]))

# add meta data and ensure obs&var unique
for i in range(len(adata)):
    adata[i].obs_names_make_unique()
    adata[i].var_names_make_unique()
    adata[i].obs['SampleID'] = metadata.query('PatientID == @fileName[@i]')['SampleID'].to_string(index=False)
    adata[i].obs['Treatment'] = metadata.query('PatientID == @fileName[@i]')['Treatment'].to_string(index=False)
    adata[i].obs['Group'] = metadata.query('PatientID == @fileName[@i]')['Group'].to_string(index=False)

# merge
adataAll = adata[0].concatenate(adata[1:])
adataAll

# filter CellBender fake "cell"
adataAll = adataAll[adataAll.X.sum(axis=1).A1 != 0]

# save
os.makedirs('./data', exist_ok=True)
save_file = './data/merged_CB_raw_data.h5ad.gz'
adataAll.write_h5ad(save_file, compression='gzip')

# Calculate and plot all basic QC metrics
adataAll = sc.read_h5ad('./data/merged_CB_raw_data.h5ad.gz')
# 1. basic metrics and gene ratio
# mitochondrial genes
adataAll.var['mt'] = adataAll.var_names.str.startswith('MT-')
# ribosomal genes
adataAll.var['ribo'] = adataAll.var_names.str.startswith(('RPS', 'RPL'))
# hemoglobin genes.
adataAll.var['hb'] = adataAll.var_names.str.contains(('^HB[^(P)]'))
# calculate
sc.pp.calculate_qc_metrics(adataAll, qc_vars=['mt', 'ribo', 'hb'], percent_top=None, log1p=False, inplace=True)
# get log2 precent hemoglobin gene
adataAll.obs['log_pct_counts_hb'] = np.log2(adataAll.obs['pct_counts_hb'] + 0.0001)
# 2. library complexity
adataAll.obs['lib_complexity'] = np.log10(adataAll.obs['n_genes_by_counts']) / np.log10(adataAll.obs['total_counts'])
# 3. Cell number
cellNum = adataAll.obs.groupby('SampleID').size().to_frame('CellNum')
cellNum['type'] = cellNum.index
# 4. Droplets entropy
allDropletEntropy = np.apply_along_axis(ndd.entropy, axis=1, arr=adataAll.X.todense(), k=adataAll.n_vars)
np.save("./data/allDropletEntropy.npy", allDropletEntropy)
# allDropletEntropy = np.load("./data/allDropletEntropy.npy")
adataAll.obs['droplet_entropy'] = allDropletEntropy

# QC metrics median
# mt median
mitoMedian = adataAll.obs.groupby('SampleID')['pct_counts_mt'].apply(lambda x: np.median(x)).to_frame('PctMtGeneMedian')
# geneNum median
geneNumMedian = adataAll.obs.groupby('SampleID')['n_genes_by_counts'].apply(lambda x: np.median(x)).to_frame('GeneNumMedian')
# umiNum median
umiNumMedian = adataAll.obs.groupby('SampleID')['total_counts'].apply(lambda x: np.median(x)).to_frame('UmiNumMedian')
# library complexity median
libComplexMedian = adataAll.obs.groupby('SampleID')['lib_complexity'].apply(lambda x: np.median(x)).to_frame('LibComplexityMedian')
# entropy
entropyMedian = adataAll.obs.groupby('SampleID')['droplet_entropy'].apply(lambda x: np.median(x)).to_frame('EntropyMedian')

# plot
# plot: `sns.violinplot`; `sns.barplot`; `sns.kdeplot`; `sns.stripplot`; `sns.boxplot`; `sns.scatterplot`; `sns.lmplot`

# save fig dir
os.makedirs('./fig', exist_ok=True)
# basic metrics
os.makedirs('./fig/before_QC/basic_metrics', exist_ok=True)
# metircs correlation scatter plot
os.makedirs('./fig/before_QC/metrics_cor_plot', exist_ok=True)
# others plots
os.makedirs('./fig/before_QC/others', exist_ok=True)

# Genes per cell
with plt.rc_context({'figure.figsize': (12, 5), 'font.family': 'Arial'}):
    ax = sc.pl.violin(adataAll, 'n_genes_by_counts',
                      jitter=False, groupby='SampleID', rotation=90, stripplot=False,
                      show=False, inner="quartile", saturation=0.6)
    ax.set_title('Genes Per Cell', size=15)
    ax.set_ylabel('')
    ax.set_yticks(ticks=[0, 500, 1000, 2000, 4000, 8000, 12000])
    ax.axhline(500, color='#FF6D28', linestyle='-.', alpha=0.8)
plt.savefig("./fig/before_QC/basic_metrics/GenesPerCell.pdf", bbox_inches='tight')

# UMIs per cell
with plt.rc_context({"figure.figsize": (12, 6), 'font.family': 'Arial'}):
    ax = sc.pl.violin(adataAll, 'total_counts',
                      jitter=False, groupby='SampleID', rotation=90, stripplot=False,
                      show=False, inner="quartile", saturation=0.6, alpha=0, log=True)
    ax.set_title('UMIs Per Cell', size=15)
    ax.set_yticks(ticks=[10, 100, 1000, 2000, 3000, 5000, 10000, 100000, 500000])
    ax.axhline(1000, color='#FF6D28', linestyle='-.', alpha=0.8)
    ax.set_ylabel('')
    ax.get_yaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
plt.savefig("./fig/before_QC/basic_metrics/UMIsPerCell.pdf", bbox_inches='tight')

myPalette = adataAll.uns['SampleID_colors']

# cell Counts
with plt.rc_context({"figure.figsize": (5, 12), 'font.family': 'Arial'}):
    ax = sns.barplot(data=cellNum, x="CellNum", y="type", palette=myPalette, saturation=0.6)
    ax.set_title('Cell Counts', size=15)
    ax.set_xlabel('')
    ax.set_ylabel('')
plt.savefig("./fig/before_QC/basic_metrics/CellCounts.pdf", bbox_inches='tight')

# Mitochondrial & Ribosomal & Hemoglobin Genes Percent
with plt.rc_context({"figure.figsize": (24, 6), 'font.family': 'Arial'}):
    fig, axes = plt.subplots(3, 1, sharex=True)
    ax1 = sns.violinplot(
        y=adataAll.obs['pct_counts_mt'], x=adataAll.obs['SampleID'],
        palette=myPalette, cut=0, sharey=False, size=1, show=False,
        inner="quartile", ax=axes[0], scale='width', saturation=1
    )
    ax2 = sns.violinplot(
        y=adataAll.obs['pct_counts_ribo'], x=adataAll.obs['SampleID'],
        palette=myPalette, cut=0, sharey=False, size=1, show=False,
        inner=None, ax=axes[1], sharex=axes[1], scale='width', saturation=1
    )
    ax3 = sns.violinplot(
        y=adataAll.obs['log_pct_counts_hb'], x=adataAll.obs['SampleID'],
        palette=myPalette, cut=0, sharey=False, size=1, show=False,
        inner=None, ax=axes[2], sharex=axes[1], scale='width', saturation=1
    )

    axes[0].set_xlabel('')
    axes[0].set_ylabel('percentMt')
    axes[1].set_xlabel('')
    axes[1].set_ylabel('percentRibo')
    axes[2].set_xlabel('')
    axes[2].set_ylabel('percentHb(log2)')
    # show mito fixed threshold
    axes[0].axhline(20, color='#FF6D28', linestyle='-.', alpha=0.8)
    fig.suptitle('Mitochondrial & Ribosomal & Hemoglobin Genes Ratio', size=15)
    fig.align_labels()
plt.xticks(rotation=90)
plt.subplots_adjust(hspace=.0)
# plt.show()
plt.savefig("./fig/before_QC/basic_metrics/SomeGenesRatio.pdf", bbox_inches='tight')

# Library complexity
with plt.rc_context({"figure.figsize": (12, 6)}):
    sns.set_theme(style="ticks", rc={"axes.facecolor": (0, 0, 0, 0), 'font.family': 'Arial'})
    ax = sns.FacetGrid(adataAll.obs, row='SampleID', aspect=16, height=0.4, hue='SampleID', palette=myPalette)
    ax.map_dataframe(sns.kdeplot, x="lib_complexity", fill=True, alpha=0.8)
    ax.map_dataframe(sns.kdeplot, x="lib_complexity", color='grey', alpha=0.8)
    ax.set_titles("")
    ax.set(yticks=[])
    ax.set_ylabels(label=None)
    ax.despine(left=False, right=False)

    def label(x, color, label):
        ax = plt.gca()
        ax.text(.05, .2, label, color='black',
                ha="left", va="center", transform=ax.transAxes)
    ax.map(label, 'SampleID')
    ax.fig.subplots_adjust(hspace=-.6)
    ax.axes[0, 0].spines['top'].set_visible(True)
    ax.fig.suptitle('Library Complexity', size=15)
    ax.set_axis_labels("log10GenesPerUMI")
plt.savefig("./fig/before_QC/basic_metrics/LibComplexRidge.pdf", bbox_inches='tight')

# all median
with plt.rc_context({"figure.figsize": (18, 4), 'font.family': 'Arial'}):
    np.random.seed(123)
    fig, axes = plt.subplots(1, 5)
    sns.stripplot(data=mitoMedian, y="PctMtGeneMedian", hue="SampleID", size=10, palette=myPalette, ax=axes[0], alpha=0.8, legend=None)
    sns.boxplot(data=mitoMedian, y='PctMtGeneMedian', showfliers=False, ax=axes[0], color="#E8E2E2")
    axes[0].set_ylabel('')
    axes[0].set_xlabel('PctMtGeneMedian')
    sns.stripplot(data=geneNumMedian, y="GeneNumMedian", hue="SampleID", size=10, palette=myPalette, ax=axes[1], alpha=0.8, legend=None)
    sns.boxplot(data=geneNumMedian, y='GeneNumMedian', showfliers=False, ax=axes[1], color="#E8E2E2")
    axes[1].set_ylabel('')
    axes[1].set_xlabel('GeneNumMedian')
    sns.stripplot(data=umiNumMedian, y="UmiNumMedian", hue="SampleID", size=10, palette=myPalette, ax=axes[2], alpha=0.8, legend=None)
    sns.boxplot(data=umiNumMedian, y='UmiNumMedian', showfliers=False, ax=axes[2], color="#E8E2E2")
    axes[2].set_ylabel('')
    axes[2].set_xlabel('UmiNumMedian')
    sns.stripplot(data=libComplexMedian, y="LibComplexityMedian", hue="SampleID", size=10, palette=myPalette, ax=axes[3], alpha=0.8, legend=None)
    sns.boxplot(data=libComplexMedian, y='LibComplexityMedian', showfliers=False, ax=axes[3], color="#E8E2E2")
    axes[3].set_ylabel('')
    axes[3].set_xlabel('LibComplexityMedian')
    sns.stripplot(data=entropyMedian, y="EntropyMedian", hue="SampleID", size=10, palette=myPalette, ax=axes[4], alpha=0.8, legend=None)
    sns.boxplot(data=entropyMedian, y='EntropyMedian', showfliers=False, ax=axes[4], color="#E8E2E2")
    axes[4].set_ylabel('')
    axes[4].set_xlabel('EntropyMedian')
plt.savefig("./fig/before_QC/basic_metrics/allMetricsMedianBoxPlot.pdf", bbox_inches='tight')

with plt.rc_context({"figure.figsize": (24, 2), 'font.family': 'Arial'}):
    ax = sns.scatterplot(data=geneNumMedian, x="SampleID", y="GeneNumMedian",
                         s=50, hue="SampleID", palette=myPalette)
    ax.axhline(20, color='#FF6D28', linestyle='-.', alpha=0.8)
    ax.set_title('mt gene median', size=15)
    ax.set_xlabel('')
    ax.set_ylabel('')
plt.xticks(rotation=90)
# change xlim
plt.xlim(-1, 60)
# plt.ylim(-1, 75)
plt.legend([], [], frameon=False)

# metrics correlation
# get corr value to check weird sample
adataAll.obs.groupby('SampleID')[['total_counts', 'pct_counts_mt']].corr().unstack().iloc[:, 1]

with plt.rc_context({"figure.figsize": (12, 12), 'font.family': 'Arial', 'axes.facecolor': 'white'}):
    # sns.set_style('ticks')
    fg = sns.lmplot(data=adataAll.obs, x='total_counts', y='pct_counts_mt', hue='SampleID',
                    palette=myPalette, legend=False,
                    scatter_kws={'s': 12, 'alpha': 0.6, 'linewidth': 0.5, 'edgecolor': 'k'},
                    line_kws={'lw': 1, 'linestyle': '--', 'alpha': 0.9})
    sns.despine(fig=None, ax=None, top=False, right=False, left=False, bottom=False, offset=None, trim=False)
    fg.axes[0, 0].set_xlabel('UMIsNum')
    fg.axes[0, 0].set_ylabel('percentMt')
    fg.axes[0, 0].legend(loc='center left', bbox_to_anchor=(1, 0.5), ncol=2, frameon=False)
plt.title("UMIsNum vs percentMt", {'fontsize': 15})
plt.savefig("./fig/before_QC/metrics_cor_plot/UMIsNumVsPercentMt.pdf", bbox_inches='tight')

# weird sample!
with plt.rc_context({"figure.figsize": (6, 6), 'font.family': 'Arial', 'axes.facecolor': 'white'}):
    fig, ax = plt.subplots()
    sns.scatterplot(
        data=adataAll.obs[(adataAll.obs['SampleID'] == 'I007_IO')],
        x='total_counts',
        y='pct_counts_mt',
        color="grey",
        alpha=0.6,
        ax=ax,
    )
    sns.kdeplot(
        data=adataAll.obs[(adataAll.obs['SampleID'] == 'I007_IO')],
        x='total_counts',
        y='pct_counts_mt',
        levels=12,
        fill=True,
        alpha=0.6,
        cmap=sns.color_palette("rocket_r", as_cmap=True),
        ax=ax,
    )
    sns.kdeplot(
        data=adataAll.obs[(adataAll.obs['SampleID'] == 'I007_IO')],
        x='total_counts',
        y='pct_counts_mt',
        levels=12,
        fill=False,
        linestyles="--",
        linewidths=1,
        color='black',
        ax=ax,
    )
plt.title("I007_IO UMIsNum vs percentMt", {'fontsize': 15})
plt.savefig("./fig/before_QC/others/I007_IO_UMIvsMt_kde.pdf", bbox_inches='tight')

# top gene
# CB rescue lots of high mt droplets, most of top gene is mt
with plt.rc_context({"figure.figsize": (4, 8), 'font.family': 'Arial'}):
    ax = sc.pl.highest_expr_genes(adataAll,
                                  n_top=20, show=False,
                                  palette=myPalette, saturation=0.5)
    ax.set_title('Top-20 Genes All Samples', size=15)
    ax.set_xlabel('percentCounts')
    ax.set_ylabel('')


# Doublet
adataAll = sc.read_h5ad('./data/merged_CB_raw_data.h5ad.gz')
os.makedirs('./fig/doublets', exist_ok=True)
# deconcatenate
adataList = [adataAll[adataAll.obs['SampleID'].isin([sample])] for sample in adataAll.obs['SampleID'].unique()]
adataAll.obs.groupby('SampleID').size()

# do scrublet
for i in range(len(adataList)):
    scrub = scr.Scrublet(adataList[i].X, expected_doublet_rate=0.08, random_state=0)
    adataList[i].obs['doublet_scores'], adataList[i].obs['predicted_doublets'] = scrub.scrub_doublets()
    # do not draw histogram
    adataList[i].obs['doublet_info'] = adataList[i].obs["predicted_doublets"].astype(str)

# merge
adataAll2 = adataList[0].concatenate(adataList[1:])

# scrublet automatic threshold looks bad
# check
scrub_test = scr.Scrublet(adataAll[adataAll.obs['SampleID'] == "C008_BL"].X, expected_doublet_rate=0.08, random_state=0)
test_score, test_predict = scrub_test.scrub_doublets()
# check auto threshold
scrub_test.plot_histogram()
# this sample the automatic can not find any doublet
np.unique(test2 == True, return_counts=True)
# 0.25 is an empirical threshold
np.unique(test1 > 0.25, return_counts=True)

# use combined threshold
adataAll2.obs['doublet_final'] = ((adataAll2.obs['doublet_scores'] > 0.25) | (adataAll2.obs['predicted_doublets'] == True))
totalDoubletRatio = sum(adataAll2.obs['doublet_final']) / adataAll2.n_obs

# all sample
cellCounts = pd.Series(adataAll2.obs.groupby('SampleID').size(), name='cellCounts')
doubleCounts = pd.Series(adataAll2.obs.groupby('SampleID')['doublet_final'].sum(), name='doubletCounts')
doubletRatio = pd.Series(doubleCounts / cellCounts, name='doubletRatio')
# merge
doubletInfo = pd.concat([cellCounts, doubletRatio], axis=1)
doubletInfo['SampleID'] = doubletInfo.index
doubletInfo
doubletInfo[['cellCounts', 'doubletRatio']].corr()

# plot
with plt.rc_context({"figure.figsize": (6, 6), 'font.family': 'Arial', 'axes.facecolor': 'white'}):
    fg = sns.lmplot(data=doubletInfo, x='cellCounts', y='doubletRatio', hue='SampleID', palette=myPalette,
                    fit_reg=True, legend=False, scatter_kws={'s': 20, 'alpha': 0.8, 'linewidth': 0.5, 'edgecolor': 'k'})
    sns.despine(fig=None, ax=None, top=False, right=False, left=False, bottom=False, offset=None, trim=False)
    fg.axes[0, 0].set_xlabel('Cell Counts')
    fg.axes[0, 0].set_ylabel('Doublet Ratio')
    fg.axes[0, 0].legend(loc='center left', bbox_to_anchor=(1, 0.5), ncol=3, frameon=False)
plt.title("CellCounts vs DoubletRatio (corr=0.4431)", {'fontsize': 15})
plt.savefig("./fig/doublets/allSampleCellNumVsDoubletRatio.pdf", bbox_inches='tight')

# wilcoxon rank test
stat, p = sp.stats.ranksums(adataAll2.obs.query('doublet_final == True')['n_genes_by_counts'],
                            adataAll2.obs.query('doublet_final != True')['n_genes_by_counts'])

# check if add nUMI and nGene to filter doublets
with plt.rc_context({'figure.figsize': (4, 4), 'font.family': 'Arial'}):
    ax = sc.pl.violin(adataAll2, 'n_genes_by_counts',
                      jitter=False, groupby='doublet_final', rotation=90, stripplot=False,
                      show=False, inner="box", saturation=0.5, palette=['lightgray', '#FF5100'])
    ax.set_title('Genes Counts (p < 1E-16)', size=15)
    ax.set_ylabel('geneNum')
    ax.set_xticklabels(['Singlets', 'Doublets'], rotation=360)
    # ax.set_yticks(ticks=[0, 500, 1000, 2000, 3000, 4000, 6000])
plt.savefig("./fig/doublets/SingletVsDoubletGeneNum.pdf", bbox_inches='tight')

# save the raw data with all metircs
save_file = './data/merged_CB_beforeFilter_data.h5ad.gz'
adataAll2.write_h5ad(save_file, compression='gzip')



# Do filter
# ====================================
# ========= droplet level =========
# ====================================

# adaptive TAD filter
def is_outlier(adata,
               metric: str,
               n_mads: int,
               side: Literal['both', 'left', 'right'] = 'both'):
    side_values = ['both', 'left', 'right']
    if side not in side_values:
        raise ValueError(f"Invalid value '{side}', allowed values are {side_values}")
    # check if a bad sample
    left_outlier = adata.obs.groupby('SampleID')[metric].apply(
        lambda x: (x < (np.median(x) - n_mads * sp.stats.median_abs_deviation(x)))
    )
    right_outlier = adata.obs.groupby('SampleID')[metric].apply(
        lambda x: (x > (np.median(x) + n_mads * sp.stats.median_abs_deviation(x)))
    )
    outlier = left_outlier | right_outlier
    if side == 'left':
        return left_outlier
    elif side == 'right':
        return right_outlier
    elif side == 'both':
        return outlier

# mt flexiable threshold
all_MAD_threshold = adataAll2.obs.groupby('SampleID')['pct_counts_mt'].apply(
    lambda x: (np.median(x) + 3 * sp.stats.median_abs_deviation(x)))

# according to the boxplot to set a flexiable threshold
with plt.rc_context({"figure.figsize": (4, 5), 'font.family': 'Arial'}):
    np.random.seed(123)
    sns.stripplot(data=all_MAD_threshold.to_frame('MAD_threshold'), y="MAD_threshold", hue="SampleID", size=10, palette=myPalette, alpha=0.8, )
    sns.boxplot(data=all_MAD_threshold.to_frame('MAD_threshold'), y='MAD_threshold', showfliers=False, color="#E8E2E2")
plt.legend([], [], frameon=False)
# 50, 10 is ok
bad_sample_threshold = 50
good_sample_threshold = 10
bad_sample_id = all_MAD_threshold[all_MAD_threshold > bad_sample_threshold].index.tolist()
good_sample_id = all_MAD_threshold[all_MAD_threshold < good_sample_threshold].index.tolist()
normal_sample_id = [x for x in all_MAD_threshold.index.tolist() if x not in bad_sample_id + good_sample_id]
# check
set(bad_sample_id + good_sample_id + normal_sample_id) == set(all_MAD_threshold.index.tolist())

# use normal sample median as bad_sample_filter
bad_sample_filter = all_MAD_threshold[normal_sample_id].median()
# good sample filter
good_sample_filter = 15

def fun_normal(x):
    return x > (np.median(x) + 3 * sp.stats.median_abs_deviation(x))
def fun_good(x):
    return x > good_sample_filter
def fun_bad(x):
    return x > bad_sample_filter

adataAll2.obs['is_mt_outlier'] = adataAll2.obs.groupby('SampleID')['pct_counts_mt'].apply(
    lambda x: fun_normal(x) if x.name in normal_sample_id else (
        fun_good(x) if x.name in good_sample_id else fun_bad(x)
    ))

# all droplet level filter
adataAll2.obs['all_outlier'] = (
    (adataAll2.obs['n_genes_by_counts'] > 8000) | \
    (adataAll2.obs['n_genes_by_counts'] < 200) | \
    (adataAll2.obs['total_counts'] > 50000) | \
    (adataAll2.obs['total_counts'] < 1000) | \
    is_outlier(adataAll2, 'droplet_entropy', 5, 'left') | \
    is_outlier(adataAll2, 'lib_complexity', 5, 'left') | \
    adataAll2.obs['doublet_final'] | \
    adataAll2.obs['is_mt_outlier']
)
# do droplet filter
adataFilter = adataAll2[~adataAll2.obs['all_outlier']]

# ====================================
# ========= Gene level =========
# ====================================

# Do gene-level filter
# gene express in low cell
sc.pp.filter_genes(adataFilter, min_cells=20)
print("Remaining genes %d" % adataFilter.n_vars)

# filter technical gene
# here we only filter 'MALAT1'
malat1 = adataFilter.var_names.str.startswith('MALAT1')
keep = np.invert(malat1)
adataFilter = adataFilter[:, keep]
print(adataFilter.n_obs, adataFilter.n_vars)

# save filtered data
save_file = './data/merged_CB_Filtered_data.h5ad.gz'
adataFilter.write_h5ad(save_file, compression='gzip')

---
title: "Lime amendment to chronically acidified forest soils results in shifts in prokaryotic communities"
output: html_document
date: "2026-5-7"
---

```{R, results='hide', fig.keep='all', message=FALSE}
options(warn = -1)
library(reticulate)
#library(kableExtra)
library(knitr)
library(phyloseq)
library(microbiome)
library(philr)
library(ape)
library(tidyr)
library(vegan)
library(randomcoloR)
library(gridExtra)
library(metacoder)
library("data.table")
library("plyr")
library(DESeq2)
library(ALDEx2) 
library(Maaslin2)
library(tidyverse)
library(readxl)
library(glue)
library(ggtext)
library(microbiomeutilities)
library(ggpubr)
conda_python(envname = 'r-reticulate', conda = "auto")

```

```{python, results='hide', fig.keep='all', message=FALSE}
import pandas as pd
import math
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse
from matplotlib.patches import Patch
import matplotlib.transforms as transforms
from scipy.stats import ttest_ind
import matplotlib as mpl
import matplotlib.cm as cm
from matplotlib.lines import Line2D
from matplotlib.offsetbox import AnchoredText
from skbio.diversity import get_alpha_diversity_metrics, get_beta_diversity_metrics, alpha_diversity, beta_diversity
from skbio import read
from skbio.tree import TreeNode
import numpy as np
from scipy.stats import ttest_ind
from deicode.preprocessing import rclr
from skbio.stats.composition import clr
from scipy.spatial import distance
from skbio.stats import ordination
from sklearn import preprocessing
#from sklearn.metrics import plot_roc_curve
from sklearn.ensemble import RandomForestClassifier, RandomForestRegressor
from sklearn.model_selection import train_test_split
from sklearn.metrics import roc_auc_score
from sklearn.metrics import confusion_matrix
from Bio import Phylo
import random
import pickle
from ete3 import Tree
#from bioinfokit.analys import stat

folder = "/Users/maggiehosmer/Library/CloudStorage/OneDrive-DalhousieUniversity/Acid_Rain_Project/March_23/"



def confidence_ellipse(x, y, ax, n_std=2.0, facecolor='none', **kwargs):
    x = np.array(x)
    y = np.array(y)
    if x.size != y.size:
        raise ValueError("x and y must be the same size")
    cov = np.cov(x, y)
    pearson = cov[0, 1]/np.sqrt(cov[0, 0] * cov[1, 1])
    # Using a special case to obtain the eigenvalues of this
    # two-dimensionl dataset.
    ell_radius_x = np.sqrt(1 + pearson)
    ell_radius_y = np.sqrt(1 - pearson)
    ellipse = Ellipse((0, 0), width=ell_radius_x * 2, height=ell_radius_y * 2,
                      facecolor=facecolor, **kwargs)
    # Calculating the stdandard deviation of x from
    # the squareroot of the variance and multiplying
    # with the given number of standard deviations.
    scale_x = np.sqrt(cov[0, 0]) * n_std
    mean_x = np.mean(x)
    # calculating the stdandard deviation of y ...
    scale_y = np.sqrt(cov[1, 1]) * n_std
    mean_y = np.mean(y)
    transf = transforms.Affine2D() \
        .rotate_deg(45) \
        .scale(scale_x, scale_y) \
        .translate(mean_x, mean_y)
    ellipse.set_transform(transf + ax.transData)
    return ax.add_patch(ellipse)

saving_figures = False
```

Draw tree:
```{python}
def draw_tree(tree, orient_tree='horizontal', vert_orient='down', axes=None, label_func=str, span=355, plot_labels=True, end_same=True, fs=10):
    # Arrays that store lines for the plot of clades
    horizontal_linecollections = []
    vertical_linecollections = []
    def get_x_positions(tree):
        """Create a mapping of each clade to its horizontal position.
        Dict of {clade: x-coord}
        """
        depths = tree.depths()
        # If there are no branch lengths, assume unit branch lengths
        if not max(depths.values()):
            depths = tree.depths(unit_branch_lengths=True)
        return depths
    def format_branch_label(clade):
                return None
    def get_y_positions(tree):
        """Create a mapping of each clade to its vertical position.
        Dict of {clade: y-coord}.
        Coordinates are negative, and integers for tips.
        """
        maxheight = tree.count_terminals()
        # Rows are defined by the tips
        heights = {tip: maxheight - i for i, tip in enumerate(reversed(tree.get_terminals()))}
        # Internal nodes: place at midpoint of children
        def calc_row(clade):
            for subclade in clade:
                if subclade not in heights:
                    calc_row(subclade)
            # Closure over heights
            heights[clade] = (
                heights[clade.clades[0]] + heights[clade.clades[-1]]
            ) / 2.0
        if tree.root.clades:
            calc_row(tree.root)
        return heights
    x_posns = get_x_positions(tree)
    y_posns = get_y_positions(tree)
    if axes is None:
        fig = plt.figure()
        if orient_tree == 'circular':
            axes = fig.add_subplot(1, 1, 1, orientation='polar')
        else:
            axes = fig.add_subplot(1, 1, 1)
    elif not isinstance(axes, plt.matplotlib.axes.Axes):
        raise ValueError("Invalid argument for axes: %s" % axes)
    leaves = [['Label', 'x loc', 'y loc', 'rotation', 'va', 'ha']]
    def draw_clade_lines(orientation="horizontal",y_here=0,x_start=0,x_here=0,y_bot=0,y_top=0,color="black",lw=".1", ls='-'):
        """Create a line.
        Graphical formatting of the lines representing clades in the plot can be
        customized by altering this function.
        """
        if orientation == "horizontal":
            axes.hlines(y_here, x_start, x_here, color=color, lw=lw, linestyle=ls)
        elif orientation == "vertical":
            axes.vlines(x_here, y_bot, y_top, color=color, linestyle=ls)
    def draw_clade(clade, x_start, color, lw, orient_tree='horizontal', vert_orient='up'):
        """Recursively draw a tree, down from the given clade."""
        x_here = x_posns[clade]
        y_here = y_posns[clade]
        xmax = max(x_posns.values())+max(x_posns.values())/30
        # phyloXML-only graphics annotations
        if hasattr(clade, "color") and clade.color is not None:
            color = clade.color.to_hex()
        if hasattr(clade, "width") and clade.width is not None:
            lw = clade.width * plt.rcParams["lines.linewidth"]
        if orient_tree == 'horizontal':
            # Draw a horizontal line from start to here
            draw_clade_lines(orientation='horizontal',y_here=y_here,x_start=x_start,x_here=x_here,color=color,lw=lw)
            if clade.name != None and end_same and '__' not in clade.name and clade.name not in ['', ' ']:
                draw_clade_lines(orientation='horizontal',y_here=y_here,x_start=xmax,x_here=x_here,color=color,lw=lw-1, ls='-.')
            # Add node/taxon labels
            if clade.name not in (None, clade.__class__.__name__, '', ' '):
                label = label_func(clade.name)
                if end_same: xplc = xmax
                else: xplc = x_here
                if plot_labels: axes.text(xplc, y_here, " %s" % label, verticalalignment="center", horizontalalignment='left', color='k', fontsize=fs)
                leaves.append([label, xplc, y_here, 0, 'center', 'left']) 
            if clade.clades:
                # Draw a vertical line connecting all children
                y_top = y_posns[clade.clades[0]]
                y_bot = y_posns[clade.clades[-1]]
                # Only apply widths to horizontal lines, like Archaeopteryx
                draw_clade_lines(orientation='vertical',x_here=x_here,y_bot=y_bot,y_top=y_top,color=color,lw=lw)
                # Draw descendents
                for child in clade:
                    draw_clade(child, x_here, color, lw)
        elif orient_tree == 'vertical':
                draw_clade_lines(orientation='vertical', x_here=y_here, y_bot=x_start, y_top=x_here,color=color,lw=lw)
                if clade.name not in (None, clade.__class__.__name__, '', ' '):
                    draw_clade_lines(orientation='vertical',x_here=y_here, y_bot=xmax, y_top=x_here,color=color,lw=lw-1, ls='-.')
                if clade.name not in (None, clade.__class__.__name__, '', ' '):
                    label = label_func(clade.name)
                    if end_same: xplc = xmax
                    else: xplc = x_here
                    if vert_orient == 'up':
                        if plot_labels: axes.text(y_here, xplc,  " %s" % label, verticalalignment='bottom', horizontalalignment='center', color='k', rotation=90, fontsize=fs)
                        leaves.append([label, y_here, xplc, 90, 'bottom', 'center']) 
                    elif vert_orient == 'down':
                        if plot_labels: axes.text(y_here, xplc,  " %s" % label, verticalalignment='top', horizontalalignment='center', color='k', rotation=90, fontsize=fs)
                        leaves.append([label, y_here, xplc, 90, 'top', 'center']) 
                if clade.clades:
                    y_top = y_posns[clade.clades[0]]
                    y_bot = y_posns[clade.clades[-1]]
                    draw_clade_lines(orientation='horizontal', y_here=x_here, x_start=y_bot, x_here=y_top, color=color,lw=lw)
                    for child in clade:
                        draw_clade(child, x_here, color, lw, orient_tree='vertical', vert_orient=vert_orient)
    def draw_clade_polar(clade, color, lw, x_start=0.1, y_start=0, span=360):
        ymax = max(y_posns.values())
        yang = span/ymax
        xmax = max(x_posns.values())+max(x_posns.values())/30
        x_here = x_posns[clade]
        y_here = y_posns[clade]
        rad = span*np.pi/180
        rad = rad/ymax
        if y_start == 0:
            y_start = rad*y_start
        y_here = rad*y_here
        if x_here != 0: 
            axes.plot([y_start, y_here], [x_start, x_here], color=color, lw=lw)
            if clade.name != None and end_same:
                axes.plot([y_start, y_here], [x_here, xmax], color=color, lw=lw-1, linestyle='-.')
        if clade.name not in (None, clade.__class__.__name__):
            label = label_func(clade.name)
            rot = y_here*(180/np.pi)
            if end_same: xplc = xmax
            else: xplc = x_here
            if rot <= 90: va, ha = 'center', 'left'
            elif rot <= 180: va, ha, rot = 'center', 'right', rot-180
            elif rot <= 270: va, ha, rot = 'center', 'right', rot-180
            else: va, ha = 'center', 'left'
            if plot_labels: axes.text(y_here, xplc, label, color='k', rotation=rot, rotation_mode='anchor', va=va, ha=ha, fontsize=fs)
            leaves.append([label, y_here, xplc, rot, va, ha])
        if clade.clades:
            y_top = y_posns[clade.clades[0]]
            y_bot = y_posns[clade.clades[-1]]
            y_top = y_top*yang*np.pi/180
            y_bot = y_bot*yang*np.pi/180
            curve = [[y_bot, y_top], [x_here, x_here]]
            x = np.linspace(curve[0][0], curve[0][1], 500)
            y = interp1d(curve[0], curve[1])(x)
            axes.plot(x, y, color=color, lw=lw)
            ymin, ymax = min(x), max(x)
            ydiff = ymax-ymin
            count = [1 for child in clade]
            count = sum(count)-2
            locs = [ymin]
            for a in range(count):
                locs.append(ydiff/(count+1)+ymin)
            locs.append(ymax)
            count = 0
            for child in clade:
                if child.name != None: 
                    y_start = y_posns[child]*rad
                else:
                    y_start = locs[count]
                draw_clade_polar(child, color, lw, x_start=x_here, y_start=y_start, span=span)
                count += 1
    plt.sca(axes)
    if orient_tree in ['horizontal', 'vertical']:
        draw_clade(tree.root, 0, "k", plt.rcParams["lines.linewidth"], orient_tree=orient_tree, vert_orient=vert_orient)
        if orient_tree == 'horizontal':
            xmax = max(x_posns.values())
            axes.set_xlim(-0.05 * xmax, 1.25 * xmax)
            # Also invert the y-axis (origin at the top)
            # Add a small vertical margin, but avoid including 0 and N+1 on the y axis
            axes.set_ylim(max(y_posns.values()) + 0.8, 0.2)
        elif orient_tree == 'vertical':
            axes.set_xlim(max(y_posns.values()) + 0.8, 0.2)
            xmax = max(x_posns.values())
            if vert_orient == 'up':
                axes.set_ylim(-0.05 * xmax, 1.25 * xmax)
            elif vert_orient == 'down':
                axes.set_ylim(1.25 * xmax, -0.05 * xmax)
        axes.set_xticks([]), axes.set_yticks([])
    elif orient_tree == 'circular':
        print('Note that if you provided an axes for this then it must be polar orientation or it will probably look very strange')
        x_start = 0
        y_start = 0
        draw_clade_polar(tree.root, "k", plt.rcParams["lines.linewidth"], x_start=x_start, y_start=y_start, span=span)
        axes.set_ylim([0, max(x_posns.values())])
        axes.yaxis.grid(False)
        axes.set_xticks([])
        axes.set_yticklabels([])
    return leaves
finished = True
```

# R
## Import data
```{R, results='hide', fig.keep='all', eval=FALSE}
folder = '/Users/maggiehosmer/Library/CloudStorage/OneDrive-DalhousieUniversity/Acid_Rain_Project/March_23/'
asv_table <- read.csv(paste(folder, "exports_9.22/feature-table_w_tax.csv", sep="")) 
sampledata <- read.csv(paste(folder, "exports_9.22/metadata_complete.csv", sep="")) #read in the metadata table
phy_tree <- read_tree(paste(folder, "exports_9.22/tree.nwk", sep='')) #read in the phylogenetic tree



taxonomy = asv_table[, c(1, 150)] #take only the OTU ID and taxonomy column to a new table
asv_table = asv_table[, 1:149] #take the OTU ID and the other columns to be the ASV table
dropping = c("DNA1", "DNA2", "DNA3", "DNA4", "DNA5", "DNA6")
asv_table_new = asv_table[ , !(names(asv_table) %in% dropping)]
asv_table = asv_table_new


asv_table_num = data.matrix(asv_table[,2:147]) #convert the ASV table to a numeric matrix
asv_table_num
rownames(asv_table_num) = asv_table[,1] #give the matrix row names
asv_table = as.matrix(asv_table_num) #convert it to a matrix

taxonomy <- separate(data = taxonomy, col = taxonomy, into = c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species"), sep = "\\;") #separate the taxonomy table so each phylogenetic level is its own column
taxmat <- taxonomy[,-1] #remove the OTU ID column from the taxonomy table
rownames(taxmat) <- taxonomy[,1] #and now give the taxonomy table the OTU IDs as row names

samples <- sampledata[, 2:10] #get the metadata columns
rownames(samples) = sampledata[,1] #and add the sample names as row names
samples = data.frame(samples, stringsAsFactors = FALSE) #convert this to a data frame

#convert these to phyloseq objects
ASV = otu_table(asv_table, taxa_are_rows = TRUE)
TAX = tax_table(taxmat)
SAMPLE = sample_data(samples)
taxa_names(TAX) <- taxonomy[,1]
physeq = phyloseq(ASV,phy_tree,TAX,SAMPLE)
taxa_names(physeq) <- paste('ASV', 1:ntaxa(physeq), sep="")
```

```{R, results='hide', fig.keep='all', fig.height=5}
tax = tax_table(physeq)
write.csv(tax, paste(folder, "R_objects_new/tax_table_new.csv", sep=""))
#then ran the add_tax.py script 
```

```{R, results='hide', fig.keep='all', fig.height=5}
taxonomy <- read.csv(paste(folder, "R_objects_new/tax_table_filled_new.csv", sep="")) 
rownames(taxonomy) <- taxonomy[,1] 
taxonomy = taxonomy[, 2:8]
TAX = tax_table(taxonomy)
taxa_names(TAX) <- rownames(taxonomy)
tax_table(physeq) = TAX
```

## Normalise

We'll perform all of the normalisations (as we did above) that we will want to use at some point
```{R, results='hide', fig.keep='all', eval=FALSE}
physeq_rare <- rarefy_even_depth(physeq, sample.size = min(sample_sums(physeq)), replace = TRUE, trimOTUs = TRUE, verbose = TRUE) #rarefy to the lowest sample depth
physeq_clr <- microbiome::transform(physeq, "clr") #convert to CLR
physeq_relabun  <- transform_sample_counts(physeq, function(x) (x / sum(x))*100) #convert to relative abundance
saveRDS(physeq, file= paste(folder, "R_objects_new/physeq_new.rds", sep=""))
saveRDS(physeq_rare, file= paste(folder, "R_objects_new/physeq_rare_new.rds", sep=""))
saveRDS(physeq_clr, file= paste(folder, "R_objects_new/physeq_clr_new.rds", sep=""))
saveRDS(physeq_relabun, file= paste(folder, "R_objects_new/physeq_relabun_new.rds", sep=""))
```

```{R, results='hide', fig.keep='all', eval=FALSE}
asv_table=otu_table(physeq)
asv_table=as.data.frame(asv_table)
write.csv(asv_table, paste(folder, "R_objects_new/asv_table_new.csv", sep=""))
```
##have run the rclr.py script
```{R, results='hide', fig.keep='all'}
asv_table_rclr = read.csv(paste(folder, "R_objects_new/asv_table_rclr_new.csv", sep="")) 
physeq_rclr=physeq
asv_table_rclr
asv_table_num = data.matrix(asv_table_rclr[,2:146]) #convert the ASV table to a numeric matrix (changed 147 to 146 6/29/22 don't know why columnds decrease, one sample lost?)
rownames(asv_table_num) = asv_table_rclr[,1] #give the matrix row names
asv_table_rclr = as.matrix(asv_table_num)
ASV = otu_table(asv_table_rclr, taxa_are_rows = TRUE)
otu_table(physeq_rclr)=ASV
saveRDS(physeq_rclr, file= paste(folder, "R_objects_new/physeq_rclr_new.rds", sep=""))
```

##R objects
```{R, results='hide', fig.keep='all'}
#16S
physeq_rare_16S = readRDS(paste(folder, "R_objects_new/", "physeq_rare_new.rds", sep=""))
physeq_relabun_16S = readRDS(paste(folder, "R_objects_new/", "physeq_relabun_new.rds", sep=""))
physeq_rclr_16S = readRDS(paste(folder, "R_objects_new/", "physeq_rclr_new.rds", sep=""))
physeq_clr_16S = readRDS(paste(folder, "R_objects_new/", "physeq_clr_new.rds", sep=""))
physeq_16S = readRDS(paste(folder, "R_objects_new/", "physeq.rds", sep=""))

#ITS
physeq_rare_ITS = readRDS(paste(folder, "R_objects_ITS/", "physeq_rare.rds", sep=""))
physeq_relabun_ITS = readRDS(paste(folder, "R_objects_ITS/", "physeq_relabun.rds", sep=""))
physeq_rclr_ITS = readRDS(paste(folder, "R_objects_ITS/", "physeq_rclr.rds", sep=""))
physeq_clr_ITS = readRDS(paste(folder, "R_objects_ITS/", "physeq_clr.rds", sep=""))
physeq_ITS = readRDS(paste(folder, "R_objects_ITS/", "physeq.rds", sep=""))

#16S
tax = as.data.frame(tax_table(physeq_16S))
write.csv(tax, paste(folder, "robyn_analysis/tables_convert_from_maggie/", "taxonomy_16S.csv", sep=""))

tree = phy_tree(physeq_16S)
write.tree(tree, file=paste(folder, "robyn_analysis/tables_convert_from_maggie/", "tree_16S.tree", sep=""))

physeq_rare_16S_df = as.data.frame(otu_table(physeq_rare_16S))
write.csv(physeq_rare_16S_df, paste(folder, "robyn_analysis/tables_convert_from_maggie/", "ft_rare_16S.csv", sep=""))
physeq_relabun_16S_df = as.data.frame(otu_table(physeq_relabun_16S))
write.csv(physeq_relabun_16S_df, paste(folder, "robyn_analysis/tables_convert_from_maggie/", "ft_relabun_16S.csv", sep=""))
physeq_rclr_16S_df = as.data.frame(otu_table(physeq_rclr_16S))
write.csv(physeq_rclr_16S_df, paste(folder, "robyn_analysis/tables_convert_from_maggie/", "ft_rclr_16S.csv", sep=""))
physeq_clr_16S_df = as.data.frame(otu_table(physeq_clr_16S))
write.csv(physeq_clr_16S_df, paste(folder, "robyn_analysis/tables_convert_from_maggie/", "ft_clr_16S.csv", sep=""))
physeq_16S_df = as.data.frame(otu_table(physeq_16S))
write.csv(physeq_16S_df, paste(folder, "robyn_analysis/tables_convert_from_maggie/", "ft_16S.csv", sep=""))
md = as.data.frame(sample_data(physeq_16S))
write.csv(md, paste(folder, "robyn_analysis/tables_convert_from_maggie/", "metadata.csv", sep=""))

#ITS
tax = as.data.frame(tax_table(physeq_ITS))
write.csv(tax, paste(folder, "robyn_analysis/tables_convert_from_maggie/", "taxonomy_ITS.csv", sep=""))

physeq_rare_ITS_df = as.data.frame(otu_table(physeq_rare_ITS))
write.csv(physeq_rare_ITS_df, paste(folder, "robyn_analysis/tables_convert_from_maggie/", "ft_rare_ITS.csv", sep=""))
physeq_relabun_ITS_df = as.data.frame(otu_table(physeq_relabun_ITS))
write.csv(physeq_relabun_ITS_df, paste(folder, "robyn_analysis/tables_convert_from_maggie/", "ft_relabun_ITS.csv", sep=""))
physeq_rclr_ITS_df = as.data.frame(otu_table(physeq_rclr_ITS))
write.csv(physeq_rclr_ITS_df, paste(folder, "robyn_analysis/tables_convert_from_maggie/", "ft_rclr_ITS.csv", sep=""))
physeq_clr_ITS_df = as.data.frame(otu_table(physeq_clr_ITS))
write.csv(physeq_clr_ITS_df, paste(folder, "robyn_analysis/tables_convert_from_maggie/", "ft_clr_ITS.csv", sep=""))
physeq_ITS_df = as.data.frame(otu_table(physeq_ITS))
write.csv(physeq_ITS_df, paste(folder, "robyn_analysis/tables_convert_from_maggie/", "ft_ITS.csv", sep=""))
md = as.data.frame(sample_data(physeq_ITS))
write.csv(md, paste(folder, "robyn_analysis/tables_convert_from_maggie/", "metadata_ITS.csv", sep=""))
```

##What we'll use:
```{R}
folder = '/Users/maggiehosmer/OneDrive - Dalhousie University/Acid_Rain_Project/March_23/'
physeq_relabun_16S = readRDS(paste(folder, "R_objects_new/", "physeq_relabun_new.rds", sep=""))
physeq_rclr_16S = readRDS(paste(folder, "R_objects_new/", "physeq_rclr_new.rds", sep=""))
physeq_16S = readRDS(paste(folder, "R_objects_new/", "physeq.rds", sep=""))

physeq_relabun_ITS = readRDS(paste(folder, "R_objects_ITS/", "physeq_relabun.rds", sep=""))
physeq_rclr_ITS = readRDS(paste(folder, "R_objects_ITS/", "physeq_rclr.rds", sep=""))
physeq_ITS = readRDS(paste(folder, "R_objects_ITS/", "physeq.rds", sep=""))
```

# PERMANOVA 16S

Bray-Curtis relative abundance, weighted unifrac relative abundance, unweighted unifrac relative abundance, robust Aitchison's distance, PHILR distance. 

Bray-Curtis on relative abundance:
```{R}
ps = physeq_relabun_16S
distance <- phyloseq::distance(ps, method="bray", weighted=F)
dist_mat = as.matrix(distance)
write.csv(dist_mat, paste(folder, 'robyn_analysis/distances/16S_relabun_braycurtis.csv', sep=''))
ads = adonis(distance ~ sample_data(ps)$Treatment*sample_data(ps)$Location*sample_data(ps)$Soil_Horizon*sample_data(ps)$Sample_within_site*sample_data(ps)$Season, parallel=12)
ads_tab = as.data.frame(ads$aov.tab)
write.csv(ads_tab, paste(folder, 'robyn_analysis/stats_tests/adonis_16S_relabun_braycurtis.csv', sep=''))

#p-adjust using the previous adonis table 
folder = "/Users/maggiehosmer/Library/CloudStorage/OneDrive-DalhousieUniversity/Acid_Rain_Project/March_23/"
ads <- read.csv(paste(folder, "robyn_analysis/stats_tests/adonis_16S_relabun_braycurtis.csv", sep = ""))
ads_table = as.data.frame(ads)
p_values = data.frame(ads_table$Pr..F.)
p_values_mat = as.matrix(p_values)
padj = p.adjust(p_values_mat, method = "fdr")
padj_tab = as.data.frame(padj)
#write.csv(padj_tab, paste(folder, 'robyn_analysis/stats_tests_mag/p_adjust/adonis_16S_relabun_braycurtis_padj_only.csv', sep=''))
ads_table$p_adjust <- c(padj)
write.csv(ads_table, paste(folder, 'robyn_analysis/stats_tests_mag/p_adjust/adonis_16S_relabun_braycurtis_padj.csv', sep=''))
```

Weighted UniFrac on relative abundance:
```{R}
ps = physeq_relabun_16S
distance <- phyloseq::distance(ps, method="unifrac", weighted=T)
dist_mat = as.matrix(distance)
write.csv(dist_mat, paste(folder, 'robyn_analysis/distances/16S_relabun_weightedunifrac.csv', sep=''))
ads = adonis(distance ~ sample_data(ps)$Treatment*sample_data(ps)$Location*sample_data(ps)$Soil_Horizon*sample_data(ps)$Sample_within_site*sample_data(ps)$Season, parallel=12)
ads_tab = as.data.frame(ads$aov.tab)
write.csv(ads_tab, paste(folder, 'robyn_analysis/stats_tests/adonis_16S_relabun_weightedunifrac.csv', sep=''))

#p-adjust using the previous adonis table
folder = "/Users/maggiehosmer/Library/CloudStorage/OneDrive-DalhousieUniversity/Acid_Rain_Project/March_23/"
ads <- read.csv(paste(folder, "robyn_analysis/stats_tests/adonis_16S_relabun_weightedunifrac.csv", sep = ""))
ads_table = as.data.frame(ads)
p_values = data.frame(ads_table$Pr..F.)
p_values_mat = as.matrix(p_values)
padj = p.adjust(p_values_mat, method = "fdr")
padj_tab = as.data.frame(padj)
write.csv(padj_tab, paste(folder, 'robyn_analysis/stats_tests_mag/p_adjust/adonis_16S_relabun_weightedunifrac_padj_only.csv', sep=''))
ads_table$p_adjust <- c(padj)
write.csv(ads_table, paste(folder, 'robyn_analysis/stats_tests_mag/p_adjust/adonis_16S_relabun_weightedunifrac_padj.csv', sep=''))
```

Unweighted UniFrac on relative abundance:
```{R}
ps = physeq_relabun_16S
distance <- phyloseq::distance(ps, method="unifrac", weighted=F)
dist_mat = as.matrix(distance)
write.csv(dist_mat, paste(folder, 'robyn_analysis/distances/16S_relabun_unweightedunifrac.csv', sep=''))
ads = adonis(distance ~ sample_data(ps)$Treatment*sample_data(ps)$Location*sample_data(ps)$Soil_Horizon*sample_data(ps)$Sample_within_site*sample_data(ps)$Season, parallel=12)
ads_tab = as.data.frame(ads$aov.tab)
write.csv(ads_tab, paste(folder, 'robyn_analysis/stats_tests/adonis_16S_relabun_unweightedunifrac.csv', sep=''))

#p-adjust using the previous adonis table
folder = "/Users/maggiehosmer/Library/CloudStorage/OneDrive-DalhousieUniversity/Acid_Rain_Project/March_23/"
ads <- read.csv(paste(folder, "robyn_analysis/stats_tests/adonis_16S_relabun_unweightedunifrac.csv", sep = ""))
ads_table = as.data.frame(ads)
p_values = data.frame(ads_table$Pr..F.)
p_values_mat = as.matrix(p_values)
padj = p.adjust(p_values_mat, method = "fdr")
padj_tab = as.data.frame(padj)
write.csv(padj_tab, paste(folder, 'robyn_analysis/stats_tests_mag/p_adjust/adonis_16S_relabun_unweightedunifrac_padj_only.csv', sep=''))
ads_table$p_adjust <- c(padj)
write.csv(ads_table, paste(folder, 'robyn_analysis/stats_tests_mag/p_adjust/adonis_16S_relabun_unweightedunifrac_padj.csv', sep=''))
```

Robust Aitchison's distance:
```{R}
ps = physeq_rclr_16S
distance <- phyloseq::distance(ps, method="euclidean", weighted=F)
dist_mat = as.matrix(distance)
write.csv(dist_mat, paste(folder, 'robyn_analysis/distances/16S_rclr_euclidean.csv', sep=''))
ads = adonis(distance ~ sample_data(ps)$Treatment*sample_data(ps)$Location*sample_data(ps)$Soil_Horizon*sample_data(ps)$Sample_within_site*sample_data(ps)$Season, parallel=12)
ads_tab = as.data.frame(ads$aov.tab)
write.csv(ads_tab, paste(folder, 'robyn_analysis/stats_tests/adonis_16S_rclr_euclidean.csv', sep=''))

#p-adjust using the previous adonis table
folder = "/Users/maggiehosmer/Library/CloudStorage/OneDrive-DalhousieUniversity/Acid_Rain_Project/March_23/"
ads <- read.csv(paste(folder, "robyn_analysis/stats_tests/adonis_16S_rclr_euclidean.csv", sep = ""))
ads_table = as.data.frame(ads)
p_values = data.frame(ads_table$Pr..F.)
p_values_mat = as.matrix(p_values)
padj = p.adjust(p_values_mat, method = "fdr")
padj_tab = as.data.frame(padj)
write.csv(padj_tab, paste(folder, 'robyn_analysis/stats_tests_mag/p_adjust/adonis_16S_rclr_euclidean_padj_only.csv', sep=''))
ads_table$p_adjust <- c(padj)
write.csv(ads_table, paste(folder, 'robyn_analysis/stats_tests_mag/p_adjust/adonis_16S_rclr_euclidean_padj.csv', sep=''))
```

PHILR distance:
```{R}
physeq_philr = physeq_16S
physeq_philr <- transform_sample_counts(physeq_philr, function(x) x+1)
phy_tree(physeq_philr) <- makeNodeLabel(phy_tree(physeq_philr), method="number", prefix='n')
otu.table <- t(otu_table(physeq_philr))
tree <- phy_tree(physeq_philr)
ps = physeq_philr

physeq.philr <- philr(otu.table, tree, part.weights='enorm.x.gm.counts', ilr.weights='blw.sqrt')
philr.dist <- dist(physeq.philr, method="euclidean")
dist_mat = as.matrix(philr.dist)
write.csv(dist_mat, paste(folder, 'robyn_analysis/distances/16S_philr.csv', sep=''))
ads = adonis(dist ~ sample_data(ps)$Treatment*sample_data(ps)$Location*sample_data(ps)$Soil_Horizon*sample_data(ps)$Sample_within_site*sample_data(ps)$Season, parallel=12)
ads_tab = as.data.frame(ads$aov.tab)
write.csv(ads_tab, paste(folder, 'robyn_analysis/stats_tests/adonis_16S_philr.csv', sep=''))

#p-adjust using the previous adonis table
folder = "/Users/maggiehosmer/Library/CloudStorage/OneDrive-DalhousieUniversity/Acid_Rain_Project/March_23/"
ads <- read.csv(paste(folder, "robyn_analysis/stats_tests/adonis_16S_philr.csv", sep = ""))
ads_table = as.data.frame(ads)
p_values = data.frame(ads_table$Pr..F.)
p_values_mat = as.matrix(p_values)
padj = p.adjust(p_values_mat, method = "fdr")
padj_tab = as.data.frame(padj)
write.csv(padj_tab, paste(folder, 'robyn_analysis/stats_tests_mag/p_adjust/adonis_16S_philr_padj_only.csv', sep=''))
ads_table$p_adjust <- c(padj)
write.csv(ads_table, paste(folder, 'robyn_analysis/stats_tests_mag/p_adjust/adonis_16S_relabun_philr_padj.csv', sep=''))
```

# PERMANOVA ITS

Bray-Curtis relative abundance, robust Aitchison's distance. 

Bray-Curtis on relative abundance:
```{R}
ps = physeq_relabun_ITS
distance <- phyloseq::distance(ps, method="bray", weighted=F)
dist_mat = as.matrix(distance)
write.csv(dist_mat, paste(folder, 'robyn_analysis/distances/ITS_relabun_braycurtis.csv', sep=''))
ads = adonis(distance ~ sample_data(ps)$Treatment*sample_data(ps)$Location*sample_data(ps)$Soil_Horizon*sample_data(ps)$Sample_within_site*sample_data(ps)$Season, parallel=12)
ads_tab = as.data.frame(ads$aov.tab)
write.csv(ads_tab, paste(folder, 'robyn_analysis/stats_tests/adonis_ITS_relabun_braycurtis.csv', sep=''))

#p-adjust using the previous adonis table
folder = "/Users/maggiehosmer/Library/CloudStorage/OneDrive-DalhousieUniversity/Acid_Rain_Project/March_23/"
ads <- read.csv(paste(folder, "robyn_analysis/stats_tests/adonis_ITS_relabun_braycurtis.csv", sep = ""))
ads_table = as.data.frame(ads)
p_values = data.frame(ads_table$Pr..F.)
p_values_mat = as.matrix(p_values)
padj = p.adjust(p_values_mat, method = "fdr")
padj_tab = as.data.frame(padj)
write.csv(padj_tab, paste(folder, 'robyn_analysis/stats_tests_mag/p_adjust/adonis_ITS_relabun_braycurtis_padj_only.csv', sep=''))
ads_table$p_adjust <- c(padj)
ads_table
write.csv(ads_table, paste(folder, 'robyn_analysis/stats_tests_mag/p_adjust/adonis_ITS_relabun_braycurtis_padj.csv', sep=''))
```

Robust Aitchison's distance:
```{R}
ps = physeq_rclr_ITS
distance <- phyloseq::distance(ps, method="euclidean", weighted=F)
dist_mat = as.matrix(distance)
write.csv(dist_mat, paste(folder, 'robyn_analysis/distances/ITS_rclr_euclidean.csv', sep=''))
ads = adonis(distance ~ sample_data(ps)$Treatment*sample_data(ps)$Location*sample_data(ps)$Soil_Horizon*sample_data(ps)$Sample_within_site*sample_data(ps)$Season, parallel=12)
ads_tab = as.data.frame(ads$aov.tab)
write.csv(ads_tab, paste(folder, 'robyn_analysis/stats_tests/adonis_ITS_rclr_euclidean.csv', sep=''))

#p-adjust using the previous adonis table 
folder = "/Users/maggiehosmer/Library/CloudStorage/OneDrive-DalhousieUniversity/Acid_Rain_Project/March_23/"
ads <- read.csv(paste(folder, "robyn_analysis/stats_tests/adonis_ITS_rclr_euclidean.csv", sep = ""))
ads_table = as.data.frame(ads)
p_values = data.frame(ads_table$Pr..F.)
p_values_mat = as.matrix(p_values)
padj = p.adjust(p_values_mat, method = "fdr")
padj_tab = as.data.frame(padj)
write.csv(padj_tab, paste(folder, 'robyn_analysis/stats_tests_mag/p_adjust/adonis_ITS_rclr_euclidean_padj_only.csv', sep=''))
ads_table$p_adjust <- c(padj)
write.csv(ads_table, paste(folder, 'robyn_analysis/stats_tests_mag/p_adjust/adonis_ITS_rclr_euclidean_padj.csv', sep=''))
```

# PERMANOVA ITS genus

Bray-Curtis relative abundance, robust Aitchison's distance. 

Bray-Curtis on relative abundance:
```{R}
ps = physeq_relabun_ITS
rnk = "ta6"
ps = tax_glom(ps, taxrank=rnk, NArm=TRUE)
distance <- phyloseq::distance(ps, method="bray", weighted=F)
dist_mat = as.matrix(distance)
write.csv(dist_mat, paste(folder, 'robyn_analysis/distances/ITS_genus_relabun_braycurtis.csv', sep=''))
ads = adonis(distance ~ sample_data(ps)$Treatment*sample_data(ps)$Location*sample_data(ps)$Soil_Horizon*sample_data(ps)$Sample_within_site*sample_data(ps)$Season, parallel=12)
ads_tab = as.data.frame(ads$aov.tab)
write.csv(ads_tab, paste(folder, 'robyn_analysis/stats_tests/adonis_ITS_genus_relabun_braycurtis.csv', sep=''))

#p-adjust using the previous adonis table 
folder = "/Users/maggiehosmer/Library/CloudStorage/OneDrive-DalhousieUniversity/Acid_Rain_Project/March_23/"
ads <- read.csv(paste(folder, "robyn_analysis/stats_tests/adonis_ITS_genus_relabun_braycurtis.csv", sep = ""))
ads_table = as.data.frame(ads)
p_values = data.frame(ads_table$Pr..F.)
p_values_mat = as.matrix(p_values)
padj = p.adjust(p_values_mat, method = "fdr")
padj_tab = as.data.frame(padj)
write.csv(padj_tab, paste(folder, 'robyn_analysis/stats_tests_mag/p_adjust/adonis_ITS_genus_relabun_braycurtis_padj_only.csv', sep=''))
ads_table$p_adjust <- c(padj)
write.csv(ads_table, paste(folder, 'robyn_analysis/stats_tests_mag/p_adjust/adonis_ITS_genus_relabun_braycurtis_padj.csv', sep=''))
```

Robust Aitchison's distance:
```{R}
ps = physeq_rclr_ITS
rnk = "ta6"
ps = tax_glom(ps, taxrank=rnk, NArm=TRUE)
distance <- phyloseq::distance(ps, method="euclidean", weighted=F)
dist_mat = as.matrix(distance)
write.csv(dist_mat, paste(folder, 'robyn_analysis/distances/ITS_genus_rclr_euclidean.csv', sep=''))
ads = adonis(distance ~ sample_data(ps)$Treatment*sample_data(ps)$Location*sample_data(ps)$Soil_Horizon*sample_data(ps)$Sample_within_site*sample_data(ps)$Season, parallel=12)
ads_tab = as.data.frame(ads$aov.tab)
write.csv(ads_tab, paste(folder, 'robyn_analysis/stats_tests/adonis_ITS_genus_rclr_euclidean.csv', sep=''))

#p-adjust using the previous adonis table 
folder = "/Users/maggiehosmer/Library/CloudStorage/OneDrive-DalhousieUniversity/Acid_Rain_Project/March_23/"
ads <- read.csv(paste(folder, "robyn_analysis/stats_tests/adonis_ITS_genus_rclr_euclidean.csv", sep = ""))
ads_table = as.data.frame(ads)
p_values = data.frame(ads_table$Pr..F.)
p_values_mat = as.matrix(p_values)
padj = p.adjust(p_values_mat, method = "fdr")
padj_tab = as.data.frame(padj)
write.csv(padj_tab, paste(folder, 'robyn_analysis/stats_tests_mag/p_adjust/adonis_ITS_genus_rclr_euclidean_padj_only.csv', sep=''))
ads_table$p_adjust <- c(padj)
write.csv(ads_table, paste(folder, 'robyn_analysis/stats_tests_mag/p_adjust/adonis_ITS_genus_rclr_euclidean_padj.csv', sep=''))
```

# Stacked bar plots

Function:
```{python}
def barplot(ft, tax, md):
  rename = {}
  for row in tax.index:
    rename[row] = tax.loc[row, 'ta3']+' '+tax.loc[row, 'ta6']
  samples = {}
  samples_treat_loc_only = {}
  for col in ft.columns:
    samples[col] = md.loc[col, 'Treatment_Location_Horizon']
    new_sample = md.loc[col, 'Soil_Horizon'].strip()+' '+md.loc[col, 'Treatment'].strip()
    samples_treat_loc_only[md.loc[col, 'Treatment_Location_Horizon']] = new_sample
  level = ft.copy(deep=True).rename(index=rename, columns=samples)
  level = level.groupby(by=level.index, axis=0).sum()
  level = level.groupby(by=level.columns, axis=1).mean().rename(columns=samples_treat_loc_only)
  level_group = level.groupby(by=level.columns, axis=1).mean()
  level_group['Mean'] = level_group.mean(axis=1)
  level_group = level_group.sort_values(by=['Mean'], ascending=False)
  level_group = level_group.head(30)
  level = level.loc[level_group.index, :]
  order = list(set(sorted(list(level.columns))))
  order = ['Upper Forest Floor Control Site1', 'Upper Forest Floor Control Site2', 'Upper Forest Floor Control Site3', 'Upper Forest Floor Control Site4', 'Upper Forest Floor Control Site5', 'Upper Forest Floor Treatment Site1', 'Upper Forest Floor Treatment Site2', 'Upper Forest Floor Treatment Site3', 'Upper Forest Floor Treatment Site4', 'Upper Forest Floor Treatment Site5', 'Lower Forest Floor Control Site1', 'Lower Forest Floor Control Site2', 'Lower Forest Floor Control Site3', 'Lower Forest Floor Control Site4', 'Lower Forest Floor Control Site5', 'Lower Forest Floor Treatment Site1', 'Lower Forest Floor Treatment Site2', 'Lower Forest Floor Treatment Site3', 'Lower Forest Floor Treatment Site4', 'Lower Forest Floor Treatment Site5', 'Upper B Horizon Control Site1', 'Upper B Horizon Control Site2', 'Upper B Horizon Control Site3', 'Upper B Horizon Control Site4', 'Upper B Horizon Control Site5','Upper B Horizon Treatment Site1', 'Upper B Horizon Treatment Site2', 'Upper B Horizon Treatment Site3', 'Upper B Horizon Treatment Site4', 'Upper B Horizon Treatment Site5']
  level = level.loc[:, order]
  ax = level.transpose().plot.bar(stacked=True, width=0.8, edgecolor='k')
  lg = ax.legend(loc='upper left', bbox_to_anchor=(1.02, 1.02), ncol=1, fontsize=8)  
  yl = ax.set_ylabel('Relative abundance (%)')
  return
```

16S:
```{python}
ft = pd.read_csv(folder+'robyn_analysis/tables_convert_from_maggie/ft_relabun_16S.csv', index_col=0, header=0)
tax = pd.read_csv(folder+'robyn_analysis/tables_convert_from_maggie/taxonomy_16S.csv', index_col=0, header=0)
md = pd.read_csv(folder+'robyn_analysis/tables_convert_from_maggie/metadata_16S.csv', index_col=0, header=0)

barplot(ft, tax, md)
plt.savefig(folder+'robyn_analysis/figures/stacked_bar_16S.png', dpi=600, bbox_inches='tight')
```

ITS:
```{python}
ft = pd.read_csv(folder+'robyn_analysis/tables_convert_from_maggie/ft_relabun_ITS.csv', index_col=0, header=0)
tax = pd.read_csv(folder+'robyn_analysis/tables_convert_from_maggie/taxonomy_ITS.csv', index_col=0, header=0)
md = pd.read_csv(folder+'robyn_analysis/tables_convert_from_maggie/metadata_ITS.csv', index_col=0, header=0)

barplot(ft, tax, md)
plt.savefig(folder+'robyn_analysis/figures/stacked_bar_ITS.png', dpi=600, bbox_inches='tight')
```

# Heatmaps relative abundance

Function:
```{python}
def make_heatmap(ft, tax, md):
  rename = {}
  for row in tax.index:
    rename[row] = tax.loc[row, 'ta3']+' '+tax.loc[row, 'ta6']
  samples = {}
  samples_treat_loc_only = {}
  for col in ft.columns:
    samples[col] = md.loc[col, 'Treatment_Location_Horizon']
    new_sample = md.loc[col, 'Soil_Horizon'].strip()+' '+md.loc[col, 'Treatment'].strip()
    samples_treat_loc_only[md.loc[col, 'Treatment_Location_Horizon']] = new_sample
  level = ft.copy(deep=True).rename(index=rename, columns=samples)
  level = level.groupby(by=level.index, axis=0).sum()
  level = level.groupby(by=level.columns, axis=1).mean().rename(columns=samples_treat_loc_only)
  level_group = level.groupby(by=level.columns, axis=1).mean()
  level_group['Mean'] = level_group.mean(axis=1)
  level_group = level_group.sort_values(by=['Mean'], ascending=False)
  level_group = level_group.head(30)
  level_group = level_group.sort_values(by=['Mean'], ascending=True)
  level = level.loc[level_group.index, :]
  order = list(set(sorted(list(level.columns))))
  order = ['Upper Forest Floor Control', 'Upper Forest Floor Treatment', 'Lower Forest Floor Control', 'Lower Forest Floor Treatment', 'Upper B Horizon Control', 'Upper B Horizon Treatment']
  level = level.loc[:, order]
  ma, mi = max(level.max(axis=1)), min(level.min(axis=1))
  fig = plt.figure(figsize=(13,10))
  ax = plt.subplot(111)
  pc = plt.pcolor(level, axes=ax, edgecolor='k', vmax=10)
  xl = plt.xticks([a+0.5 for a in range(len(level.columns))], level.columns, rotation=90)
  yl = plt.yticks([a+0.5 for a in range(len(level.index))], level.index)
  mid = np.mean([ma, mi])
  for x in range(len(level.columns)):
    for y in range(len(level.index)):
      num = level.iloc[y, x]
      if num > 5: color = 'k'
      else: color = 'w'
      ax.text(x+0.5, y+0.5, str(round(num, 2)), color=color, fontsize=8, ha='center', va='center')
  return
```

16S:
```{python}
ft = pd.read_csv(folder+'robyn_analysis/tables_convert_from_maggie/ft_relabun_16S.csv', index_col=0, header=0)
tax = pd.read_csv(folder+'robyn_analysis/tables_convert_from_maggie/taxonomy_16S.csv', index_col=0, header=0)
md = pd.read_csv(folder+'robyn_analysis/tables_convert_from_maggie/metadata_16S.csv', index_col=0, header=0)

make_heatmap(ft, tax, md)
plt.savefig(folder+'robyn_analysis/figures/heatmap_16S.png', dpi=600, bbox_inches='tight')
```

ITS:
```{python}
ft = pd.read_csv(folder+'robyn_analysis/tables_convert_from_maggie/ft_relabun_ITS.csv', index_col=0, header=0)
tax = pd.read_csv(folder+'robyn_analysis/tables_convert_from_maggie/taxonomy_ITS.csv', index_col=0, header=0)
md = pd.read_csv(folder+'robyn_analysis/tables_convert_from_maggie/metadata_ITS.csv', index_col=0, header=0)

make_heatmap(ft, tax, md)
plt.savefig(folder+'robyn_analysis/figures/heatmap_ITS.png', dpi=600, bbox_inches='tight')
```

# Tree and heatmap

16S get trees:
```{python}
tax = pd.read_csv(folder+'robyn_analysis/tables_convert_from_maggie/taxonomy_16S.csv', index_col=0, header=0)
ft = pd.read_csv(folder+'robyn_analysis/tables_convert_from_maggie/ft_relabun_16S.csv', index_col=0, header=0)
# tax = pd.read_csv(folder+'robyn_analysis/tables_convert_from_maggie/taxonomy_ITS.csv', index_col=0, header=0)
# ft = pd.read_csv(folder+'robyn_analysis/tables_convert_from_maggie/ft_relabun_ITS.csv', index_col=0, header=0)
```

Get trees for all levels:
```{R}
tax_table = py$tax
TAX = tax_table(tax_table)
rownames(TAX) = rownames(tax_table)
phy_tree <- read_tree(paste(py$folder, 'robyn_analysis/tables_convert_from_maggie/tree_ITS.tree', sep=''))
ft = py$ft
table_num = data.matrix(ft[,1:146])
rownames(table_num) = rownames(ft)
table = otu_table(table_num, taxa_are_rows = TRUE)
physeq_all = phyloseq(table, phy_tree, TAX)

ranks = c("ta6", "ta5", "ta4", "ta3", "ta2")
names = c("genus", "family", "order", "class", "phylum")
for (a in 1:length(ranks)) {
  physeq_level = tax_glom(physeq_all, taxrank=ranks[a])
  level_tree = phy_tree(physeq_level)
  write.tree(level_tree, paste(py$folder, 'robyn_analysis/processing/', names[a], '_tree.nwk', sep=''))
}
```

## 16S

```{python, eval=FALSE}
plt.figure(figsize=(50,30))
tree_names = ['phylum', 'class', 'order', 'family', 'genus']
letters = ['A', 'B', 'C', 'D', 'E']
tax_dict = {}
tax = pd.read_csv(folder+'robyn_analysis/tables_convert_from_maggie/taxonomy_16S.csv', index_col=0, header=0)
for row in tax.index:
  tax_dict[row] = tax.loc[row, :].values
ft = pd.read_csv(folder+'robyn_analysis/tables_convert_from_maggie/ft_relabun_16S.csv', index_col=0, header=0)
ft_clr = pd.read_csv(folder+'robyn_analysis/tables_convert_from_maggie/ft_rclr_16S.csv', index_col=0, header=0)
md = pd.read_csv(folder+'robyn_analysis/tables_convert_from_maggie/metadata_16S.csv', index_col=0, header=0)

col_count, row_count = 0, 0
width = [2, 2, 2, 2, 2]
rc, cc = 60, 32 # number of rows, number of columns
rs = 25
for a in range(5):
  if a == 3: col_count, row_count = 0, 32
  ax_tree = plt.subplot2grid((rc,cc),(row_count,col_count), colspan=2, rowspan=rs, frameon=False)
  plt.title('   '+letters[a]+'   '+tree_names[a].capitalize()+'\n\n', loc='left', fontweight='bold', fontsize=16)
  col_count += 2
  ax_labels = plt.subplot2grid((rc,cc),(row_count,col_count), colspan=width[a], rowspan=rs, frameon=False)
  col_count += width[a]
  plt.xticks([]), plt.yticks([])
  ax_heat_prev = plt.subplot2grid((rc,cc),(row_count,col_count), colspan=2, rowspan=rs)
  tti = plt.xticks([]), plt.yticks([]), plt.title('Prevalence', fontweight='bold')
  ax_heat_abun = plt.subplot2grid((rc,cc),(row_count,col_count+2), colspan=2, rowspan=rs)
  tti = plt.xticks([]), plt.yticks([]), plt.title('Relative\nabundance (%)', fontweight='bold')
  ax_heat_clr = plt.subplot2grid((rc,cc),(row_count,col_count+4), colspan=2, rowspan=rs)
  tti, plt.xticks([]), plt.yticks([]), plt.title('CLR\nabundance', fontweight='bold')
  col_count += 7
  ft_level = ft.copy(deep=True)
  ft_clr_level = ft_clr.copy(deep=True)
  rename_level, rename_level_opposite = {}, {}
  tree_name = folder+'robyn_analysis/processing/'+tree_names[a]+'_tree.nwk'
  for row in ft_level.index:
    rename_level[row] = tax_dict[row][a+1]
    rename_level_opposite[tax_dict[row][a+1]] = row
  ft_level = ft_level.rename(index=rename_level)
  ft_level = ft_level.groupby(by=ft_level.index, axis=0).sum()
  ft_clr_level = ft_clr_level.rename(index=rename_level)
  ft_clr_level = ft_clr_level.groupby(by=ft_clr_level.index, axis=0).sum()
  if ft_level.shape[0] > 30:
    ft_group = ft_level.copy(deep=True).transpose()
    rename_samples = {}
    for sample in ft_group.index:
      rename_samples[sample] = md.loc[sample, 'Treatment_Horizon']
    ft_group = ft_group.rename(index=rename_samples)
    ft_group = ft_group.groupby(by=ft_group.index, axis=0).mean().transpose()
    ft_group['Mean'] = ft_group.mean(axis=1)
    ft_group = ft_group.sort_values(by=['Mean'], ascending=False)
    ft_group = ft_group.iloc[:30, :]
    ft_level = ft_level.loc[ft_group.index, :]
    ft_clr_level = ft_clr_level.loc[ft_group.index, :]
    tree = Tree(tree_name, format=1)
    rename_tree_level = {}
    keeping = []
    for node in tree.traverse("postorder"):
      if node.name in rename_level:
        rename_tree_level[rename_level[node.name]] = node.name
        if rename_level[node.name] in ft_level.index:
          keeping.append(node.name)
    tree.prune(keeping)
    tree_name = folder+'robyn_analysis/processing/'+tree_names[a]+'_reduced_tree.txt'
    tree.write(outfile=tree_name, format=1)
  tree = Phylo.read(tree_name, "newick")
  leaves = draw_tree(tree, axes=ax_tree, end_same=True, plot_labels=False)
  order = []
  for leaf in leaves[1:]:
    if leaf[0] in rename_level:
      tx = ax_tree.text(leaf[1], leaf[2], '  '+rename_level[leaf[0]], va=leaf[4], ha=leaf[5])
      order.append(rename_level[leaf[0]])
  ft_level = ft_level.loc[order, :]
  ft_clr_level = ft_clr_level.loc[order, :]
  plt.ylim([0.5, leaves[-1][2]+0.5])
  rename_sample = {}
  for col in ft_level.columns:
    rn = md.loc[col, 'Soil_Horizon']+' '+md.loc[col, 'Treatment']
    rename_sample[col] = rn.replace('  ', ' ').strip()
  ft_level = ft_level.rename(columns=rename_sample)
  ft_clr_level = ft_clr_level.rename(columns=rename_sample)
  ft_level_prev = ft_level.copy(deep=True).transpose()
  ft_level_abun = ft_level.copy(deep=True).transpose()
  ft_clr_level_abun = ft_clr_level.copy(deep=True).transpose()
  ft_level_prev[ft_level_prev > 0] = 1
  treat_order = ['Upper Forest Floor Control', 'Upper Forest Floor Treatment', 'Lower Forest Floor Control', 'Lower Forest Floor Treatment', 'Upper B Horizon Control', 'Upper B Horizon Treatment']
  ft_level_prev = ft_level_prev.groupby(by=ft_level_prev.index, axis=0).mean().transpose().loc[:, treat_order]
  ft_level_abun = ft_level_abun.groupby(by=ft_level_abun.index, axis=0).mean().transpose().loc[:, treat_order]
  ft_clr_level_abun = ft_clr_level_abun.groupby(by=ft_clr_level_abun.index, axis=0).mean().transpose().loc[:, treat_order]
  min_prev, max_prev, min_abun, max_abun, min_clr, max_clr = min(ft_level_prev.min(axis=1)), max(ft_level_prev.max(axis=1)), min(ft_level_abun.min(axis=1)), max(ft_level_abun.max(axis=1)), min(ft_clr_level_abun.min(axis=1)), max(ft_clr_level_abun.max(axis=1))
  mid_prev, mid_abun, mid_clr = np.mean([min_prev, max_prev]), np.mean([min_abun, max_abun]), np.mean([min_clr, max_clr])
  print(min_prev, max_prev, min_abun, max_abun, min_clr, max_clr)
  plt.sca(ax_heat_prev)
  plt.pcolor(ft_level_prev, cmap='PuBu', edgecolor='k')
  plt.sca(ax_heat_abun)
  plt.pcolor(ft_level_abun, cmap='RdPu', edgecolor='k', vmax=10)
  plt.sca(ax_heat_clr)
  plt.pcolor(ft_clr_level_abun, cmap='bwr', edgecolor='k', vmin=-10, vmax=10)
  axes, dfs, mids = [ax_heat_prev, ax_heat_abun, ax_heat_clr], [ft_level_prev, ft_level_abun, ft_clr_level_abun], [0.5, 5, 5]
  x = [0.5, 1.5, 2.5, 3.5, 4.5, 5.5]
  for z in range(len(axes)):
    for y in range(len(dfs[z].index.values)):
      tax = dfs[z].index.values[y]
      for t in range(len(treat_order)):
        col = 'k'
        val = dfs[z].loc[tax, treat_order[t]]
        rnd = 2
        if abs(val) >= mids[z]: col = 'w'
        if abs(val) >= 1: rnd = 1
        if abs(val) >= 10:
          tx = axes[z].text(x[t], y+0.5, str(int(val)), ha='center', va='center', color=col)
        else:
          tx = axes[z].text(x[t], y+0.5, str(round(val,rnd)), ha='center', va='center', color=col)
  for ax in [ax_heat_prev, ax_heat_abun, ax_heat_clr]:
    plt.sca(ax)
    xt = plt.xticks([0.5, 1.5, 2.5, 3.5, 4.5, 5.5], [t.replace('Treat', '\nTreat').replace('Cont', '\nCont') for t in treat_order], rotation=90)

plt.savefig(folder+'robyn_analysis/figures/heatmap_all_levels_16S.png', dpi=600, bbox_inches='tight')
```

## ITS

```{python, eval=FALSE}
plt.figure(figsize=(50,30))
tree_names = ['phylum', 'class', 'order', 'family', 'genus']
letters = ['A', 'B', 'C', 'D', 'E']
tax_dict = {}
tax = pd.read_csv(folder+'robyn_analysis/tables_convert_from_maggie/taxonomy_ITS.csv', index_col=0, header=0)
for row in tax.index:
  tax_dict[row] = tax.loc[row, :].values
ft = pd.read_csv(folder+'robyn_analysis/tables_convert_from_maggie/ft_relabun_ITS.csv', index_col=0, header=0)
ft_clr = pd.read_csv(folder+'robyn_analysis/tables_convert_from_maggie/ft_rclr_ITS.csv', index_col=0, header=0)
md = pd.read_csv(folder+'robyn_analysis/tables_convert_from_maggie/metadata_ITS.csv', index_col=0, header=0)

col_count, row_count = 0, 0
width = [2, 2, 2, 2, 2]
rc, cc = 60, 32 # number of rows, number of columns
rs = 25
for a in range(5):
  if a == 3: col_count, row_count = 0, 32
  ax_labels = plt.subplot2grid((rc,cc),(row_count,col_count), colspan=width[a], rowspan=rs, frameon=False)
  plt.title('   '+letters[a]+'   '+tree_names[a].capitalize()+'\n\n', loc='left', fontweight='bold', fontsize=16)
  col_count += width[a]
  plt.xticks([]), plt.yticks([])
  ax_heat_prev = plt.subplot2grid((rc,cc),(row_count,col_count), colspan=2, rowspan=rs)
  tti = plt.xticks([]), plt.yticks([]), plt.title('Prevalence', fontweight='bold')
  ax_heat_abun = plt.subplot2grid((rc,cc),(row_count,col_count+2), colspan=2, rowspan=rs)
  tti = plt.xticks([]), plt.yticks([]), plt.title('Relative\nabundance (%)', fontweight='bold')
  ax_heat_clr = plt.subplot2grid((rc,cc),(row_count,col_count+4), colspan=2, rowspan=rs)
  tti, plt.xticks([]), plt.yticks([]), plt.title('CLR\nabundance', fontweight='bold')
  col_count += 7
  ft_level = ft.copy(deep=True)
  ft_clr_level = ft_clr.copy(deep=True)
  rename_level, rename_level_opposite = {}, {}
  tree_name = folder+'robyn_analysis/processing/'+tree_names[a]+'_tree.nwk'
  for row in ft_level.index:
    rename_level[row] = tax_dict[row][a+1]
    rename_level_opposite[tax_dict[row][a+1]] = row
  ft_level = ft_level.rename(index=rename_level)
  ft_level = ft_level.groupby(by=ft_level.index, axis=0).sum()
  ft_clr_level = ft_clr_level.rename(index=rename_level)
  ft_clr_level = ft_clr_level.groupby(by=ft_clr_level.index, axis=0).sum()
  if ft_level.shape[0] > 30:
    ft_group = ft_level.copy(deep=True).transpose()
    rename_samples = {}
    for sample in ft_group.index:
      rename_samples[sample] = md.loc[sample, 'Treatment_Horizon']
    ft_group = ft_group.rename(index=rename_samples)
    ft_group = ft_group.groupby(by=ft_group.index, axis=0).mean().transpose()
    ft_group['Mean'] = ft_group.mean(axis=1)
    ft_group = ft_group.sort_values(by=['Mean'], ascending=False)
    ft_group = ft_group.iloc[:30, :]
    ft_level = ft_level.loc[ft_group.index, :]
    ft_clr_level = ft_clr_level.loc[ft_group.index, :]
  order = list(ft_level.index.values)
  ft_level = ft_level.loc[order, :]
  ft_clr_level = ft_clr_level.loc[order, :]
  #plt.ylim([0.5, len(order)+0.5])
  rename_sample = {}
  for col in ft_level.columns:
    rn = md.loc[col, 'Soil_Horizon']+' '+md.loc[col, 'Treatment']
    rename_sample[col] = rn.replace('  ', ' ').strip()
  ft_level = ft_level.rename(columns=rename_sample)
  ft_clr_level = ft_clr_level.rename(columns=rename_sample)
  ft_level_prev = ft_level.copy(deep=True).transpose()
  ft_level_abun = ft_level.copy(deep=True).transpose()
  ft_clr_level_abun = ft_clr_level.copy(deep=True).transpose()
  ft_level_prev[ft_level_prev > 0] = 1
  treat_order = ['Upper Forest Floor Control', 'Upper Forest Floor Treatment', 'Lower Forest Floor Control', 'Lower Forest Floor Treatment', 'Upper B Horizon Control', 'Upper B Horizon Treatment']
  ft_level_prev = ft_level_prev.groupby(by=ft_level_prev.index, axis=0).mean().transpose().loc[:, treat_order]
  ft_level_abun = ft_level_abun.groupby(by=ft_level_abun.index, axis=0).mean().transpose().loc[:, treat_order]
  ft_clr_level_abun = ft_clr_level_abun.groupby(by=ft_clr_level_abun.index, axis=0).mean().transpose().loc[:, treat_order]
  min_prev, max_prev, min_abun, max_abun, min_clr, max_clr = min(ft_level_prev.min(axis=1)), max(ft_level_prev.max(axis=1)), min(ft_level_abun.min(axis=1)), max(ft_level_abun.max(axis=1)), min(ft_clr_level_abun.min(axis=1)), max(ft_clr_level_abun.max(axis=1))
  mid_prev, mid_abun, mid_clr = np.mean([min_prev, max_prev]), np.mean([min_abun, max_abun]), np.mean([min_clr, max_clr])
  plt.sca(ax_heat_prev)
  plt.pcolor(ft_level_prev, cmap='PuBu', edgecolor='k')
  yt = plt.yticks([a+0.5 for a in range(len(order))], order)
  plt.sca(ax_heat_abun)
  plt.pcolor(ft_level_abun, cmap='RdPu', edgecolor='k', vmax=10)
  plt.sca(ax_heat_clr)
  plt.pcolor(ft_clr_level_abun, cmap='bwr', edgecolor='k', vmin=-10, vmax=10)
  axes, dfs, mids = [ax_heat_prev, ax_heat_abun, ax_heat_clr], [ft_level_prev, ft_level_abun, ft_clr_level_abun], [0.5, 5, 5]
  x = [0.5, 1.5, 2.5, 3.5, 4.5, 5.5]
  for z in range(len(axes)):
    for y in range(len(dfs[z].index.values)):
      tax = dfs[z].index.values[y]
      for t in range(len(treat_order)):
        col = 'k'
        val = dfs[z].loc[tax, treat_order[t]]
        rnd = 2
        rnd = 2
        if abs(val) >= mids[z]: col = 'w'
        if abs(val) >= 1: rnd = 1
        if abs(val) >= 10:
          tx = axes[z].text(x[t], y+0.5, str(int(val)), ha='center', va='center', color=col)
        else:
          tx = axes[z].text(x[t], y+0.5, str(round(val,rnd)), ha='center', va='center', color=col)
  for ax in [ax_heat_prev, ax_heat_abun, ax_heat_clr]:
    plt.sca(ax)
    xt = plt.xticks([0.5, 1.5, 2.5, 3.5, 4.5, 5.5], [t.replace('Treat', '\nTreat').replace('Cont', '\nCont') for t in treat_order], rotation=90)

plt.savefig(folder+'robyn_analysis/figures/heatmap_all_levels_ITS.png', dpi=600, bbox_inches='tight')
```


# Alpha diversity R

Treatment 16S:
```{R}
plot_richness(physeq_rare_16S, x="Treatment", measures=c("Observed", "Chao1", "Simpson", "Shannon")) + geom_boxplot()
```

Horizon 16S:
```{R}
plot_richness(physeq_rare_16S, x="Soil_Horizon", measures=c("Observed", "Chao1", "Simpson", "Shannon")) + geom_boxplot()
```

Treatment ITS:
```{R}
plot_richness(physeq_rare_ITS, x="Treatment", measures=c("Observed", "Chao1", "Simpson", "Shannon")) + geom_boxplot()
```

Horizon ITS:
```{R}
plot_richness(physeq_rare_ITS, x="Soil_Horizon", measures=c("Observed", "Chao1", "Simpson", "Shannon")) + geom_boxplot()
```

# Alpha diversity python

16S:
```{python}
#Treatment, Horizon, Treatment and horizon #2, 3, 6
ft = pd.read_csv(folder+'robyn_analysis/tables_convert_from_maggie/ft_rare_16S.csv', index_col=0, header=0)
md = pd.read_csv(folder+'robyn_analysis/tables_convert_from_maggie/metadata.csv', index_col=0, header=0)
tree = read(folder+'robyn_analysis/tables_convert_from_maggie/tree_16S.tree', format="newick", into=TreeNode)
groups = [['Control', 'Treatment'], ['Upper Forest Floor', 'Lower Forest Floor', 'Upper B Horizon'], ['Control_Upper Forest Floor', 'Treatment_Upper Forest Floor', 'Control_Lower Forest Floor', 'Treatment_Lower Forest Floor', 'Control_Upper B Horizon', 'Treatment_Upper B Horizon']]
color_dict = {'Control':'#16A085', 'Treatment':'#F1C40F', 'Upper Forest Floor':'#E74C3C', 'Lower Forest Floor':'#8E44AD', 'Upper B Horizon':'#3498DB'}
colors = [[color_dict['Control'], color_dict['Treatment']], [color_dict['Upper Forest Floor'], color_dict['Lower Forest Floor'], color_dict['Upper B Horizon']], [color_dict['Control'], color_dict['Treatment'], color_dict['Control'], color_dict['Treatment'], color_dict['Control'], color_dict['Treatment']]]
x_labs = [['Control', 'Treatment'], ['Upper Forest Floor', 'Lower Forest Floor', 'Upper B Horizon'], ['Control', 'Treatment', 'Control', 'Treatment', 'Control', 'Treatment']]
x = [[0,1], [0,1,2], [0,1,2.5,3.5,5,6]]
compare = [[[0,1]], [[0,1], [0,2], [1,2]], [[0,1], [2,3], [4,5]]]
x_locs_hor, x_labs_hor = [0.5, 3, 5.5], groups[1]

a_div_measures = ['observed_otus', 'chao1', 'shannon', 'simpson', 'faith_pd']
a_div_labels = ['Number of ASVs', 'Chao1 richness', 'Shannon diversity', "Simpson's diversity", "Faith's phylogenetic diversity"]
fig = plt.figure(figsize=(10,12.5))
for a in range(len(a_div_measures)):
  axes = [plt.subplot2grid((len(a_div_measures),11), (a,0), colspan=2), plt.subplot2grid((len(a_div_measures),11), (a,2), colspan=3), plt.subplot2grid((len(a_div_measures),11), (a,5), colspan=6)]
  if a_div_measures[a] != 'faith_pd': alpha_div = alpha_diversity(a_div_measures[a], ft.transpose())
  else: alpha_div = alpha_diversity(a_div_measures[a], ft.transpose(), otu_ids=ft.index.values, tree=tree, validate=False)
  alpha_div = alpha_div.to_frame().rename(columns={0:a_div_measures[a]})
  alpha_div.index = ft.columns
  for b in range(len(groups)):
    plt.sca(axes[b])
    this_group = {}
    for c in range(len(groups[b])): this_group[groups[b][c]] = []
    for sample in alpha_div.index:
      if b == 0: md_group = md.loc[sample, 'Treatment'].strip()
      elif b == 1: md_group = md.loc[sample, 'Soil_Horizon'].strip()
      else: md_group = md.loc[sample, 'Treatment'].strip()+'_'+md.loc[sample, 'Soil_Horizon'].strip()
      this_group[md_group].append(alpha_div.loc[sample, a_div_measures[a]])
    for d in range(len(groups[b])):
      box = axes[b].boxplot(this_group[groups[b][d]], positions=[x[b][d]], widths=0.6, showfliers=False)
      for item in ['boxes', 'whiskers', 'fliers', 'medians', 'caps']: bi = plt.setp(box[item], color='k')
      scat = axes[b].scatter(np.random.normal(x[b][d], scale=0.075, size=len(this_group[groups[b][d]])), this_group[groups[b][d]], color=colors[b][d], alpha=0.5, edgecolors='gray')
    if a < len(a_div_measures)-1: ti = plt.xticks(x[b], [])
    else: 
      ti = plt.xticks(x[b], x_labs[b], rotation=90)
      if b == 2:
        for e in range(len(x_locs_hor)):
          xl = plt.text(x_locs_hor[e], -75, x_labs_hor[e], ha='center', va='top')
    high, low = float(axes[b].get_ylim()[1]), float(axes[b].get_ylim()[0])
    span = high-low
    diff = span*0.01
    high = high+diff*8
    sig = 0
    for o in range(len(compare[b])):
      n1, n2 = compare[b][o][0], compare[b][o][1]
      t1, t2 = groups[b][n1], groups[b][n2]
      stat, p = ttest_ind(this_group[t1], this_group[t2])
      if p > 0.05: continue
      if sig > 0: high = high+diff*16
      li = plt.plot([x[b][n1], x[b][n1], x[b][n2], x[b][n2]], [high-diff, high, high, high-diff], color='k')
      string = '*'
      if p <= 0.01: string = '**'
      if p <= 0.005: string = '***'
      tx = plt.text(np.mean([x[b][n1], x[b][n2]]), high-diff*4, string, va='bottom', ha='center', fontsize=10)
      if b == 1: sig += 1
    yl = plt.ylim([low,high+span*0.2])
  axes[0].set_ylabel(a_div_labels[a], fontweight='bold')
  if a == 0:
    axes[0].set_title('Treatment', fontweight='bold')
    axes[1].set_title('Soil Horizon', fontweight='bold')
    axes[2].set_title('Treatment at each Soil Horizon', fontweight='bold')

plt.subplots_adjust(wspace=1.4)
plt.savefig(folder+'robyn_analysis/figures/alpha_diversity_16S.png', dpi=600, bbox_inches='tight')
```

ITS:
```{python}
#Treatment, Horizon, Treatment and horizon #2, 3, 6
ft = pd.read_csv(folder+'robyn_analysis/tables_convert_from_maggie/ft_rare_ITS.csv', index_col=0, header=0)
md = pd.read_csv(folder+'robyn_analysis/tables_convert_from_maggie/metadata_ITS.csv', index_col=0, header=0)
groups = [['Control', 'Treatment'], ['Upper Forest Floor', 'Lower Forest Floor', 'Upper B Horizon'], ['Control_Upper Forest Floor', 'Treatment_Upper Forest Floor', 'Control_Lower Forest Floor', 'Treatment_Lower Forest Floor', 'Control_Upper B Horizon', 'Treatment_Upper B Horizon']]
color_dict = {'Control':'#16A085', 'Treatment':'#F1C40F', 'Upper Forest Floor':'#E74C3C', 'Lower Forest Floor':'#8E44AD', 'Upper B Horizon':'#3498DB'}
colors = [[color_dict['Control'], color_dict['Treatment']], [color_dict['Upper Forest Floor'], color_dict['Lower Forest Floor'], color_dict['Upper B Horizon']], [color_dict['Control'], color_dict['Treatment'], color_dict['Control'], color_dict['Treatment'], color_dict['Control'], color_dict['Treatment']]]
x_labs = [['Control', 'Treatment'], ['Upper Forest Floor', 'Lower Forest Floor', 'Upper B Horizon'], ['Control', 'Treatment', 'Control', 'Treatment', 'Control', 'Treatment']]
x = [[0,1], [0,1,2], [0,1,2.5,3.5,5,6]]
compare = [[[0,1]], [[0,1], [0,2], [1,2]], [[0,1], [2,3], [4,5]]]
x_locs_hor, x_labs_hor = [0.5, 3, 5.5], groups[1]

a_div_measures = ['observed_otus', 'chao1', 'shannon', 'simpson']
a_div_labels = ['Number of ASVs', 'Chao1 richness', 'Shannon diversity', "Simpson's diversity"]
fig = plt.figure(figsize=(10,10))
for a in range(len(a_div_measures)):
  axes = [plt.subplot2grid((len(a_div_measures),11), (a,0), colspan=2), plt.subplot2grid((len(a_div_measures),11), (a,2), colspan=3), plt.subplot2grid((len(a_div_measures),11), (a,5), colspan=6)]
  if a_div_measures[a] != 'faith_pd': alpha_div = alpha_diversity(a_div_measures[a], ft.transpose())
  else: alpha_div = alpha_diversity(a_div_measures[a], ft.transpose(), otu_ids=ft.index.values, tree=tree, validate=False)
  alpha_div = alpha_div.to_frame().rename(columns={0:a_div_measures[a]})
  alpha_div.index = ft.columns
  for b in range(len(groups)):
    plt.sca(axes[b])
    this_group = {}
    for c in range(len(groups[b])): this_group[groups[b][c]] = []
    for sample in alpha_div.index:
      if b == 0: md_group = md.loc[sample, 'Treatment'].strip()
      elif b == 1: md_group = md.loc[sample, 'Soil_Horizon'].strip()
      else: md_group = md.loc[sample, 'Treatment'].strip()+'_'+md.loc[sample, 'Soil_Horizon'].strip()
      this_group[md_group].append(alpha_div.loc[sample, a_div_measures[a]])
    for d in range(len(groups[b])):
      box = axes[b].boxplot(this_group[groups[b][d]], positions=[x[b][d]], widths=0.6, showfliers=False)
      for item in ['boxes', 'whiskers', 'fliers', 'medians', 'caps']: bi = plt.setp(box[item], color='k')
      scat = axes[b].scatter(np.random.normal(x[b][d], scale=0.075, size=len(this_group[groups[b][d]])), this_group[groups[b][d]], color=colors[b][d], alpha=0.5, edgecolors='gray')
    if a < len(a_div_measures)-1: ti = plt.xticks(x[b], [])
    else: 
      ti = plt.xticks(x[b], x_labs[b], rotation=90)
      if b == 2:
        for e in range(len(x_locs_hor)):
          xl = plt.text(x_locs_hor[e], -0.55, x_labs_hor[e], ha='center', va='top')
    high, low = float(axes[b].get_ylim()[1]), float(axes[b].get_ylim()[0])
    span = high-low
    diff = span*0.01
    high = high+diff*8
    sig = 0
    for o in range(len(compare[b])):
      n1, n2 = compare[b][o][0], compare[b][o][1]
      t1, t2 = groups[b][n1], groups[b][n2]
      stat, p = ttest_ind(this_group[t1], this_group[t2])
      if p > 0.05: continue
      if sig > 0: high = high+diff*16
      li = plt.plot([x[b][n1], x[b][n1], x[b][n2], x[b][n2]], [high-diff, high, high, high-diff], color='k')
      string = '*'
      if p <= 0.01: string = '**'
      if p <= 0.005: string = '***'
      tx = plt.text(np.mean([x[b][n1], x[b][n2]]), high-diff*4, string, va='bottom', ha='center', fontsize=10)
      if b == 1: sig += 1
    yl = plt.ylim([low,high+span*0.2])
  axes[0].set_ylabel(a_div_labels[a], fontweight='bold')
  if a == 0:
    axes[0].set_title('Treatment', fontweight='bold')
    axes[1].set_title('Soil Horizon', fontweight='bold')
    axes[2].set_title('Treatment at each Soil Horizon', fontweight='bold')

plt.subplots_adjust(wspace=1.4)
plt.savefig(folder+'robyn_analysis/figures/alpha_diversity_ITS.png', dpi=600, bbox_inches='tight')
```

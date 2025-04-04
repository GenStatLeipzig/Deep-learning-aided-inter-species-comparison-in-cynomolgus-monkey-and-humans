############### Module containing helper functions for CCC noteboks ###############
# Some functions were taken from other sources and modified as needed             #
# In these cases it is noted in the function where the originial was taken from.  #
# If not mentioned otherwise code by:                                             #
# Author: kmlr81, Kristina Müller                                                 #
###################################################################################

# imports
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import scanpy as sc


############### Cell-cell communication ###############

def get_expressed_genes(adata, celltype, cluster_col, expr_prop=0.1):
    """
    Modified from chapter 21. Cell-cell communication of: Heumos, L., Schaar, A.C., Lance, C. et al. 
    Best practices for single-cell analysis across modalities.Nat Rev Genet (2023). 
    https://doi.org/10.1038/s41576-023-00586-w 
    
    A function to calculate proportions of cells of a certain celltype 
    in which a specific gene is expressed, then filter with given cutoff
    to determine whether a gene can be considered expressed in the celltpye
    or not.
    
    Parameters:
    ----------
    adata: AnnData object for experiment.
    celltype: String with celltype name.
    expr_prop: Float between 0 - 1. Percentage of cells of a celltype in which gene 
               is expressed must be >= expr_prop. Default = 0.1.
    ----------
    """
    # Get adat for specific celltype
    adata_celltype = adata[adata.obs[cluster_col] == celltype].copy()
    # number of cells in which a gene is expressed / number of cells total for all genes in dataset
    percentage_expressed = adata_celltype.X.getnnz(axis=0) / adata_celltype.X.shape[0]
    
    # Get data frame connecting gene names with percentages
    gene_stats = (
        pd.DataFrame({"genes": adata_celltype.var_names, "props": percentage_expressed})
        .sort_values("genes")
    )

    # get expressed genes
    filtered_stats = gene_stats[gene_stats["props"] >= expr_prop]
    expressed_genes = filtered_stats["genes"].values

    return expressed_genes

def prep_heatmap_df(df, directed=True, weighted=True):
    """
    Prepares a given data frame for heatmap plotting. The heatmap to plot can be any combination of 
    directed/undirected and weighted/not weighted.

    Parameters:
    ----------
    df: A pandas.DataFrame object containing the columns "source", "target" and "magnitude_rank" when weighted=True
        or "number" when weighted=False.
    directed: Boolean. Default = True. When the default value is set the resulting data frame can be used to plot
              a directed heatmap. If false, the heatmap will be mirrored along the diagonal, values added up,
              to create the basis for an undirected heatmap.
    weighted: Boolean. Default = True. If the default value is given the magnitude rank is used for a weighted summation
              of potential ccc interactions between celltype pairs. If set to False the number of connections is taken
              with no weights applied.
    ---------- 
    """

    #  Rename columns, if neccessary
    if "source_cluster" not in df.columns:
        df = df.rename(columns={"source": "source_cluster", 
                                "target": "target_cluster", 
                                "magnitude_rank": "score"})
        
    if weighted:
        # Get sum of -log10 transformed magnitude values
        df = df.rename(columns={"magnitude_rank": "score"})
        df["score"] = df["score"].apply(lambda x: -np.log10(x + np.finfo(float).eps))
        df = df.groupby(["source_cluster", "target_cluster"]).sum()["score"].reset_index()
    else:
        df = df.rename(columns={"number": "score"})

    df_pivoted = df.pivot(index="source_cluster", columns="target_cluster", values="score")

    if not directed:
        upper = np.triu(df_pivoted.values) + np.tril(df_pivoted.values, k=-1).T
        lower = np.triu(df_pivoted.values, k=1).T + np.tril(df_pivoted.values, k=-1)
        df_mirrored = pd.DataFrame(upper + lower, index=df_pivoted.index, columns=df_pivoted.columns)
        df_mirrored = df_mirrored.rename_axis(None, axis=1)  
        df_mirrored = df_mirrored.rename_axis(None, axis=0)

        return df_mirrored
        
    else:
        return df_pivoted
    

def plot_heatmap(heatmap_df, 
                 title="", 
                 vmin=None, 
                 vmax=None, 
                 cmap="Reds", 
                 title_fontsize=14, 
                 y_label_rotation=0,
                 annot=False,
                 ax=None):
    """
    Wrapper function for seaborn's heatmap() function.

    Parameters:
    ----------
    heatmap_df: pandas.DataFrame object. Data frame to use for plotting the heatmap.
    title: String, default="". Default plots no title. Otherwise the given String is used as heatmap title.
    vmin, vmax: Numerical, default=None. If default is set, the min and max values found in the input data frame 
                are used for color scale boundaries. Else the given values are set.
    cmap: String matching a color map name, default="Reds". Sets the color theme of the heatmap.
    title_fontsize: Integer, default=14. Sets the font size of the heatmap title.
    y_label_rotation: Integer, default=0. Sets angle of rotation of y labels. The default sets no rotation.
    annot: Boolean, default=False. Sets whether or not to annotate the heatmap squares wit their corresponding values. 
           The default value sets no annotaiton.
    ax: matplotlib.Axes object, default=None. An axes object to plot the heatmap onto. If the default is given
        the function reaturns a fresh axes object.
    ----------  
    """
    
    heatmap = sns.heatmap(heatmap_df, cmap=cmap, vmin=vmin, vmax=vmax, annot=annot, ax=ax)
    heatmap.set_yticks(heatmap.get_yticks(), heatmap.get_yticklabels(), rotation=y_label_rotation)
    heatmap.set_title(title,fontsize=title_fontsize)

    return heatmap


############### Quality Control & Preprocessing ###############

def get_soupx_groups(adata):
    """
    Helper function that takes an AnnData object, makes a copy, normalizes, log-transforms,
    caclulates a PCA and then uses basic leiden to cluster the copy. The leiden clusters are
    extracted and returned.

    Parameters:
    ----------
    adata: AnnData object a copy of which is clustered for the SoupX analysis.
    ----------

    Returns:
    ----------
    soupx_groups: Pandas.Series object. Column with soupX groups for cells in AnnData.
    ----------
    """
    # Get a copy of the input Anndata
    adata_c = adata.copy()

    # Normalize and log-transform the copy
    sc.pp.normalize_per_cell(adata_c)
    sc.pp.log1p(adata_c)

    # Redice dimensionality and cluster with key "soupX_groups"
    sc.pp.pca(adata_c)
    sc.pp.neighbors(adata_c)
    sc.tl.leiden(adata_c, key_added="soupx_groups")

    # Get soupX groups for SoupX
    soupx_groups = adata_c.obs["soupx_groups"]

    # Delete copy of input Anndata
    del adata_c

    return soupx_groups


def is_outlier(adata, 
               metric, 
               nmads=3, 
               lower_cutoff=True):
    """
    Determines whether or not cells in a dataset are considered outliers, i.e. low quality cells, based on a 
    specific, pre-calculated QC metric and the set number of MADs (median absolute deviations).

    Copied from chapter 6.3 of the book [Single-cell best practices](https://www.sc-best-practices.org/preamble.html) 
    after Heumos, L., Schaar, A.C., Lance, C. et al. Best practices for single-cell analysis across modalities. 
    Nat Rev Genet (2023). https://doi.org/10.1038/s41576-023-00586-w 

    Parameters:
    ----------
    adata: Anndata object with cells and QC metrics
    metric: String. Name of the QC column to filter
    nmads: Integer. Number of MADS for filtering. Default = 3
    lower_cutoff: Boolean. Decides whether or not to calculate the lower MAD cutoff. 
                           Default = True.
    ----------

    Returns:
    ----------
    outlier: pandas.Series object with Boolean values denoting whether a cell is low quality or not.
    ----------
    """
    # Get QC colum for which to calculate MAD cutoffs form .obs table
    M = adata.obs[metric]

    # if lower_cutoff=True calculate upper and lower cutoff, else calculate lower cutoff only
    # Use cutoff(s) to filter QC colum and get a column True/False to show whether each cell in 
    # the data set is an outlier or not
    if lower_cutoff:
        outlier = (M < np.median(M) - nmads * median_abs_deviation(M)) | (np.median(M) + nmads * median_abs_deviation(M) < M) # Import for median_abs_deviation in notebook
    else:
        outlier = np.median(M) + nmads * median_abs_deviation(M) < M

    return outlier


def save_QC_cutoffsMAD(adata, 
                       columns, 
                       save=None, 
                       return_df=False, 
                       cols_no_lower=["pct_counts_in_top_20_genes", "pct_counts_mt"]):
    """
    Saves caclulated MAD cutoffs for a number of QC metrics as a data frame.

    Function by Vincent David Friedrich, modified and comments added by Kristina Müller

    Parameters:
    ----------
    adata: Anndata object. Contains QC information for the cells in the data set.
    columns: List of Strings. Names of QC columns for which to save MAD cutoffs.
    save: String. Path to save cutoff file to including file name and file type. Default = None. 
                  If the default is given the resulting data frame won't be saved.
    return_df: Boolean. Whether or not the resulting data frame object should be returned. Default = False.
    cols_no_lower: List of Strings. Names of QC columns with no lower cutoff. Default=["pct_counts_in_top_20_genes", "pct_counts_mt"].
                                    The value for lower_limit is set to 0 for these columns.
    ----------

    Returns:
    ----------
    If return_df = True:
        df: pandas.DataFrame object with MAD cutoffs for QC columns.
    ----------
    """
    df = pd.DataFrame(columns=['lower_limit', 'upper_limit'])
    for column in columns:
        upper_limit = np.max(adata[~adata.obs["outlier_" + column]].obs[column])

        if column in cols_no_lower:
            lower_limit = 0
        else:
            lower_limit = np.min(adata[~adata.obs["outlier_" + column]].obs[column])
        
        df = pd.concat([df, pd.DataFrame([[lower_limit, upper_limit]], columns=['lower_limit', 'upper_limit'])], ignore_index=True)

    df.index = [columns]
    if save is not None:
        df.to_csv(save)
    if return_df:
        return df
    

def calculate_clusterBased_QC(adata_in, 
                              obs_column_qc, 
                              obs_column_clustering, 
                              sort_by, 
                              id, 
                              cutoff_bad_cluster,
                              in_place=False):
    """
    Does basic cluster-based QC by taking pre-calculated QC metrics and applying them to fine grained clustering
    of the data set to determine and mark low quality clusters and high qualite clusters based on the percentage
    of low quality cells per cluster. 

    Function by Vincent David Friedrich, modified and comments added by Kristina Müller

    Paramters:
    ----------
    adata_in: Anndata Object. Contains QC columns in .obs table and columns assigning low quality cells.
    obs_columns_qc: List of Strings. List of QC column names to use.
    obs_column_clustering: String. Name of column in .obs table containing the pre-computet clustering to base QC on.
    sort_by: String. Name of column in adata.obs to sort by.
    id: String. Name of to set for new QC column in adata.obs.
    cutoff_bad_cluster: Float. Value Range 0-1. Percentage of low quality cells in a cluster must be > cutoff for the cluster
                               to be deemed low quality.
    in_place: Boolean. Decides whether to edit adata.obs in place or copy and return a new adata.
    ----------

    Returns:
    ----------
    perc_df: Pandas.DataFrame object. Data frame with percentages of good quality and bad quality cells per QC cluster
    If in_place = False:
        adata: AnnData object. Anndata with new QC column added to obs table.
    ----------
    """
    # Check whether to perform changes in adata.obs in place or return new adata with changes
    if not in_place:
        adata = adata_in.copy()
    else:
        adata = adata_in
    
    # Get data frame with percentages of high quality cells and low quality cells per QC cluster
    perc_df = adata.obs.groupby(obs_column_clustering, observed=False)[obs_column_qc].value_counts(normalize=True).unstack().sort_values(by=sort_by, ascending=False)
    
    # Get list of clusters that have a percentage of bad clles > cutoff_bad_cluster
    bad_clusters = list(perc_df[perc_df[sort_by] > cutoff_bad_cluster].index.values)

    # Add temporary QC columns to adata.obs based on bad_clusters list
    adata.obs['cell_quality_cluster' + str(id)] = adata.obs[obs_column_clustering].isin(bad_clusters)
    adata.obs['cell_quality_cluster' + str(id)] = adata.obs['cell_quality_cluster' + str(id)].map({True: 'low_quality_cluster', False: 'high_quality_cluster'})
    
    # Add final QC column to adata.obs
    adata.obs[id] = (adata.obs['cell_quality_cluster' + str(id)].astype(str) + '_' + adata.obs[obs_column_qc].astype(str)).values
    
    # Remove temporary QC columns from adata.obs
    adata.obs = adata.obs.drop(['cell_quality_cluster' + str(id)], axis=1)

    if in_place:
        return perc_df
    else:
        return adata, perc_df


def calculate_clusterBased_QC_doublets(adata_in, column_clustering, cutoff_bad_cluster=0.25, doublet_filtering='twice', in_place=False):
    """
    Does cluster-based QC based for doublet detection.

    Function by Vincent David Friedrich, modified and comments added by Kristina Müller

    Paramters:
    ----------
    adata_in: Anndata object. Anndata for which to perform cluster-based doublet QC.
    column_clustering: String. Name of column in adata.obs with clusters to use for cluster-based QC.
    cutoff_bad_cluster: Float between 0.0 and 1.0. Clusters with > cutoff percentage of cells = doublet will be marked as bad clusters.
    doublet_filtering: String. Decides how to classify doublets. 
                               Permitted values: 'once' -> a cell is a doublet if at least one tool has classified it as such
                                                 'twice' -> a cell is a doublet if at least two tools have classified it as such
    in_place: Boolean. Sets whether or not to make the changes to adata.obs in place or to make a copy and return the new anndata. Default=False. 
    ----------

    Returns:
    ----------
    If in_place=False:
        adata: Anndata object. Copy of the input Anndata with new data added.
    doub_perc_df: pandas.DataFrame object. Data frame with percentages of doublets and of singlets per cluster in the data set.
    ----------
    """
    if in_place:
        adata = adata_in 
    else:
        adata = adata_in.copy()

    # Get data frame with doublet and singlet percentages per cluster depending on filtering parameter
    if doublet_filtering == 'once':
        doub_perc_df = adata.obs.groupby(column_clustering,observed=False)['doublet_once'].value_counts(normalize=True).unstack().sort_values(by='doublet', ascending=False)
        adata.obs['doublet'] = adata.obs['doublet_once'].map({'doublet': 'bad_cell', 'singlet': 'good_cell'})
    elif doublet_filtering == 'twice':
        doub_perc_df = adata.obs.groupby(column_clustering,observed=False)['doublet_twice'].value_counts(normalize=True).unstack().sort_values(by='doublet', ascending=False)
        adata.obs['doublet'] = adata.obs['doublet_twice'].map({'doublet': 'bad_cell', 'singlet': 'good_cell'})
    elif doublet_filtering == 'thrice':
        doub_perc_df = adata.obs.groupby(column_clustering,observed=False)['doublet_thrice'].value_counts(normalize=True).unstack().sort_values(by='doublet', ascending=False)
        adata.obs['doublet'] = adata.obs['doublet_thrice'].map({'doublet': 'bad_cell', 'singlet': 'good_cell'})

    # Get list of bad clusters by cluster number using pre-set cutoff
    bad_clusters = list(doub_perc_df[doub_perc_df['doublet'] > cutoff_bad_cluster].index.values)

    # Add temporary column with cluster evaluation based on cutoff
    adata.obs['bad_cluster_' + str(cutoff_bad_cluster)] = adata.obs[column_clustering].isin(bad_clusters)
    adata.obs['bad_cluster_' + str(cutoff_bad_cluster)] = adata.obs['bad_cluster_' + str(cutoff_bad_cluster)].map({True: 'bad_cluster', False: 'good_cluster'})

    # Add doublet QC column to adata.obs
    adata.obs['doublet_QC'] = (adata.obs['bad_cluster_' + str(cutoff_bad_cluster)].astype(str) + '_' + adata.obs['doublet'].astype(str)).values

    # Remove adata.obs columns no longer needed
    adata.obs.drop(['doublet', 'bad_cluster_' + str(cutoff_bad_cluster)], axis=1) 

    if in_place:
        return doub_perc_df
    else:
        return adata, doub_perc_df


def doublet_summary(adata_in, 
                    db_col_tool1, 
                    db_col_tool2, 
                    db_col_tool3, 
                    in_place=False):
    """
    Takes information from doublet detection tools and annotates cells on whether one, two or all three 
    tools have detected a cell as doublet.

    Function by Vincent David Friedrich, modified and comments added by Kristina Müller
    
    Paramters:
    ----------
    adata_in: Anndata object. Full experiment anndata with obs table to add columns to for doublet tools.
    db_col_tool1: String. Name of column with Boolean True/False values denoting if a cell is a doublet for tool one.
    db_col_tool2: String. Name of column with Boolean True/False values denoting if a cell is a doublet for tool two.
    db_col_tool3: String. Name of column with Boolean True/False values denoting if a cell is a doublet for tool three.
    in_place: Boolean. Sets whether or not to make the changes to adata.obs in place or to make a copy and return the new anndata. Default=False.
    ----------

    Returns:
    ----------
    If in_place=False:
        adata: Anndata object with new doublet columns 
    ----------
    """
    if in_place:
        adata = adata_in
    else:
        adata = adata_in.copy()

    #m Mark cells that are detected as doublets by at least one tool out of three
    adata.obs['doublet_once'] = sum((adata.obs[db_col_tool1], adata.obs[db_col_tool2], adata.obs[db_col_tool3])) >= 1
    adata.obs['doublet_once'] = adata.obs['doublet_once'].map({True: 'doublet', False: 'singlet'})

    # Mark cells that are detected as doublets by at least two tools out of three
    adata.obs['doublet_twice'] = sum((adata.obs[db_col_tool1], adata.obs[db_col_tool2], adata.obs[db_col_tool3])) >= 2
    adata.obs['doublet_twice'] = adata.obs['doublet_twice'].map({True: 'doublet', False: 'singlet'})

    # Mark cells that are detected as doublets by all three tools
    adata.obs['doublet_thrice'] = sum((adata.obs[db_col_tool1], adata.obs[db_col_tool2], adata.obs[db_col_tool3])) == 3
    adata.obs['doublet_thrice'] = adata.obs['doublet_thrice'].map({True: 'doublet', False: 'singlet'})

    if not in_place:
        return adata


def violinplot_QC(adata, 
                  column,
                  title="" ,
                  show=True,
                  save=None,
                  has_lower=True,
                  upper=None,
                  lower=None, 
                  ax=None):
    """
    Generates a violin plot for a given QC metric with calculated upper and lower cutoffs
    for filtering out low quality cells.

    Function by Vincent David Friedrich, modified and comments added by Kristina Müller

    Parameters:
    ----------
    adata: Anndata object. Holds the data to generate the violin plot from
    column: String. Name of the QC column to visualize
    show: Boolean. Sets whether or not to display the violing plot in the notebook. Default = True
    save: String. Path to save the plot to, including plot name and file format. Default = None. If
                  the default is set, the plot is not saved.
    has_lower: Boolean. Determines whether or not a lower cutoff is drawn into theplot. Default = True.
    upper: Number value. Custom upper cutoff to be drawn in plot. Default=None and determines the cutoff from the data.
    lower: Number value. Custom lower cutoff to be drawn in plot. Default=None and determines the cutoff from the data.
    ----------
    """

    upper_limit = np.max(adata[~adata.obs["outlier_" + column]].obs[column]) if upper == None else upper
    
    ax_violin = sns.violinplot(adata.obs[column], ax=ax)
    ax_violin.axhline(y=upper_limit, color='r', linestyle='--',label = 'upper_limit')

    if has_lower:
        lower_limit = np.min(adata[~adata.obs["outlier_" + column]].obs[column]) if lower == None else lower
        ax_violin.axhline(y=lower_limit, color='k', linestyle='--',label = 'lower_limit')
        
    ax_violin.set_title(title)
    ax_violin.legend()
    if save is not None:
        ax_violin.figure.savefig(save, dpi=300, bbox_inches="tight")
    if show:
        plt.show()
    if ax is None:
        plt.close()


def pieplot_QC(QC_summary, 
               title="", 
               ax=None, 
               save=None, 
               show=True):
    """
    Generates a pie plot for QC data.

    Function by Vincent David Friedrich, modified and comments added and modified by Kristina Müller
    
    Paramters:
    ----------
    QC_summary: pandas.DataFrame object. Date frame with one colum percentages of cluster-based QC results to be plotted.
    title: String. Plot title. Default="" and will set no title.
    ax: matplotlib.Axes object. As to plot onto. Default=None. If no axes object is given a separate pyplot plot will be 
                                generated and closed after optional saving and displaying. 
    save: String. Path including file name and ending for saving the generated plot. Defalut=None and will not save the plot. 
    show: Boolean. Wheter or not to display the plot. Default=True. Default value will display the plot.
    ----------
    """
    if ax is not None:
        # Draw plot on Axis object
        QC_summary.plot(kind='pie', autopct=lambda p: '{:.0f} ({:.1f}%)'.format(p * QC_summary.sum() / 100, p), ax=ax)
        # Add y-axis title
        ax.set_ylabel('')
        # Add plot title
        ax.set_title(title)
        # Save plot to file
        if save is not None:
            ax.figure.savefig(save, dpi=300, bbox_inches="tight")
    else:
        # Generate matplotlib.figure.Figure object for plotting
        plt.figure(figsize=(10, 10))
        # Draw plot onto figue
        QC_summary.plot(kind='pie', autopct=lambda p: '{:.0f} ({:.1f}%)'.format(p * QC_summary.sum() / 100, p))
        # Set y-axis title
        plt.ylabel('')
        # Set plot title
        plt.title(title)
        # Save plot to file
        if save is not None:
            plt.savefig(save, dpi=300, bbox_inches="tight")
    if show:
        # Display plot
        plt.show()
    if ax is None:
        # Close current figure when done plotting
        plt.close()


def QC_plot_pct_clusters(df, 
                         column, 
                         cutoff=None, 
                         xlabel=None, 
                         ylabel=None, 
                         title=None, 
                         ax=None, 
                         save=None, 
                         show=True):
    """
    Generates a plot tracking percentages of marked cells per QC cluster with the option to add
    visualization of a set cut off.

    Function by Vincent David Friedrich, modified and comments added by Kristina Müller
    
    Paramters:
    ----------
    df: Pandas.DataFrame object. Data frame for plotting.
    column: String. Name of column in df to plot.
    cutoff: Float. Percentage of cells per cluster at whicht to set the visualized cutoff. Default=None.
                   When default value is given, no cutoff is visualized.
    xlabel: String. Label for x-axis of plot. Default=None. Default value sets no label.
    ylabel: String. Label for y-axis of plot. Default=None. Default value sets no label.
    title: Striing. Title for plot. Default=None. Default value sets no title.
    ax: matplotlib.axes.Axes object. Axes object to draw plot onto. Default=None. If default is given
                                     function will generate its own figure object and plot onto that.
    save: String. Path to save plot to with plot name and file ending. Default=None. Default value does
                  not save plot.
    show: Boolean. Sets whether to show plot after plotting. Default=True. Default will show plot. 
    ----------
    """
    if ax is not None:
        ax.plot(df[column])
        if cutoff is not None:
            ax.axhline(y=cutoff, color='r', linestyle='--',label = 'cutoff')
        ax.tick_params(axis='x', labelrotation=90, labelsize=6)
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
        ax.set_title(title)
        if save is not None:
                ax.figure.savefig(save, dpi=300, bbox_inches="tight")
    else:
        plt.figure()
        plt.plot(df[column])
        if cutoff is not None:
            plt.axhline(y=cutoff, color='r', linestyle='--', label='cutoff')
        plt.xticks(rotation=90,fontsize=6)
        plt.xlabel(xlabel)
        plt.ylabel(ylabel)
        plt.title(title)
        if save is not None:
                plt.savefig(save, dpi=300, bbox_inches="tight")
    if show:
        plt.show()
    if ax is None:
        plt.close()


def get_grouped_data_frame(data_frame, get_source=True):
    """
    Prepares a data frame with normalized, summed magnitude per cell type per timepoint/condition per species
    for either senders counting ligand interactions or targets counting receptor interactions.

    Function by Kristina Müller, kmlr81

    Paramters:
    ----------
    data_frame: pandas.DataFrame object. Pre-grouped and filtered liana results data frame with all species
                                         timepoints/conditions to plot for.
    get_source: Boolean. True/False deciding whether to get sources data frame or targets data frame. Default = True and
                         generates sources data frame.
    ----------

    Returns:
    ----------
    If get_source == False:
        grouped_targets: pandas.DataFrame object. Filtered and grouped data frame with normalized summed magnitude
                                                of cell-cell interactions per target cell type per condition/timepoint per species.
    else if get_source == True:
        grouped_senders: pandas.DataFrame object. Filtered and grouped data frame with normalized summed magnitude
                                                of cell-cell interactions per sender cell type per condition/timepoint per species.
    ----------
    """
    # Determine wheter to generate data frame for senders or receiver cell types
    if get_source:
        communicator = "source"
        gene_type = "ligand"
    else:
        communicator = "target"
        gene_type = "receptor"

    # Group by target, species, and group, and sum the log10magnitude_rank column
    grouped_df = data_frame.groupby([communicator, 'species', 'group']).agg({
        "n_cells_in_" + communicator + "_cluster": 'first',
        "n_" + gene_type + "s": 'first',
        'log10magnitude_rank': 'sum'
    }).reset_index()

    # Rename the summed column
    grouped_df = grouped_df.rename(columns={'log10magnitude_rank': 'sum_magnitude'})
    
    # Calculate sum_magnitude per receptor and find the maximum per species for plotting
    grouped_df.loc[:, "sum_magnitude_per_" + gene_type] = grouped_df['sum_magnitude'] / grouped_df["n_" + gene_type + "s"]
    max_sum_magnitude_per_receptor = grouped_df.groupby('species')["sum_magnitude_per_" + gene_type].transform('max')
    grouped_df.loc[:, "maxsum_magnitude_rel_" + communicator] = max_sum_magnitude_per_receptor

    return grouped_df


def prep_barplot_data_frames(data_frame, filter_magnitude=0.05, filter_specificity=0.05):
    """
    Prepares concatenated liana data frame for plotting by splitting into source and target data frames and 
    adding neccessary additional columns for normalized summed magnitude ranks per source/target.

    Function by Kristina Müller, kmlr81

    Paramters:
    ----------
    data_frame: pandas.DataFrame object. Liana results generated in 07_CCC_lf_infrence notebook.
    filter_magnitude: Float. Filter for magnitude_rank column in data_frame. The filter is applied
                             to keep only ligand-receptor pairs with a magnitude rank <= filter. Default = 0.05.
    filter_specificity: Float. Filter for specificity_rank column in data_frame. The filter is applied
                               to keep only ligand-receptor pairs with a specificity rank <= filter. Default = 0.05.
    ----------

    Returns:
    ----------
    grouped_targets: pandas.DataFrame object. Filtered and grouped data frame with normalized summed magnitude
                                              of cell-cell interactions per target cell type per condition/timepoint per species.
    grouped_senders: pandas.DataFrame object. Filtered and grouped data frame with normalized summed magnitude
                                               of cell-cell interactions per sender cell type per condition/timepoint per species.
    ----------
    """
    # Add total number of ligand-receptor interactions to each sender cell type and each receiver cell type
    data_frame.loc[:, "n_ligands"] = data_frame.groupby(["source", "species", "group"])["ligand_complex"].transform("nunique")
    data_frame.loc[:, "n_receptors"] = data_frame.groupby(["target", "species", "group"])["receptor_complex"].transform("nunique")

    # Filter all results with magnitude rank <= filter_magnitude and specificity rankm <= filter_specificity
    data_frame = data_frame[(data_frame["magnitude_rank"] <= filter_magnitude) & 
                            (data_frame["specificity_rank"] <= filter_specificity)].copy()

    data_frame.loc[:, 'log10magnitude_rank'] = -np.log10(data_frame['magnitude_rank'])

    # Group by target, species, and group, and sum the log10magnitude_rank column
    grouped_targets = get_grouped_data_frame(data_frame, get_source=False)

    # Group by source, species, and group, and sum the log10magnitude_rank column
    grouped_senders = get_grouped_data_frame(data_frame, get_source=True)

    return grouped_targets, grouped_senders
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cellTypeAssignMarkerGenes.R
\name{cellTypeAssignMarkerGenes}
\alias{cellTypeAssignMarkerGenes}
\title{\code{cellTypeAssignMarkerGenes} assigns CDSeq-identified cell types using user-provided marker gene list and plots heatmap.}
\usage{
cellTypeAssignMarkerGenes(
  cell_gep = NULL,
  marker_gene_list = NULL,
  threshold = 2,
  fig_path = getwd(),
  rowlabelsize = 1,
  collabelsize = 1,
  margins = c(3, 0),
  fig_width = 100,
  fig_height = 100,
  keysize = 1,
  srtcol = 45,
  keypar = c(3.5, 0, 3, 0),
  heatmap_name = "cellTypeAssign_heatmap.pdf",
  heatmap_name_fuzzy_assign = "cellTypeAssign_heatmap_fuzzy.pdf",
  verbose = FALSE
)
}
\arguments{
\item{cell_gep}{gene expression profile matrix with G rows (genes) and M columns (cell types).}

\item{marker_gene_list}{a G (genes) by C (cell types with known identities) matrix or dataframe that contains the marker genes for each cell type. Column names must be CellType and GeneName.}

\item{threshold}{a numeric value that provides the threshold of whether a known cell type in the marker gene list can be identified.}

\item{fig_path}{the location where the heatmap figure is saved.}

\item{rowlabelsize}{row label size}

\item{collabelsize}{column label size}

\item{margins}{a vector of length 2 indicates row and column label margins}

\item{fig_width}{figure width for pdf figure}

\item{fig_height}{figure height for pdf figure}

\item{keysize}{color key size for heatmap}

\item{srtcol}{column label angle}

\item{keypar}{color key layout}

\item{heatmap_name}{the name of heatmap figure of one-to-one assignment.}

\item{heatmap_name_fuzzy_assign}{the name of heatmap figure of fuzzy assignment.}

\item{verbose}{if TRUE, some information will be printed.}
}
\value{
cellTypeAssignMarkerGenes returns a list containing: 
GEP_markerSum (a A by B matrix where A is nrow(marker_gene_list), B is ncol(cell_gep)), 

GEP_markerSum_zscore (row-wise z score of GEP_markerSum), 

GEP_matched is cell_gep(,cell_type_idx),

cell_type_idx (column index of cell_gep that are considered matching with cell types in marker_gene_list), 

cell_type_matched stores the cell types in marker_gene_list that are considered to be matched with cell_gep, 

GEP_markerSum_zscore_matched contains only the rows of GEP_markerSum_zscore that are considered to be matched with some cell types in cell_gep. GEP_markerSum_zscore_matched and GEP_markerSum_zscore have same columns. 

cell_type_matched_fuzzy is a zero-one matrix that has the same size as GEP_markerSum_zscore_matched. If (i,j) element is one, means ith cell type in marker_gene_list is assigned to jth element in cell_gep.
}
\description{
\code{cellTypeAssignMarkerGenes} assigns CDSeq-identified cell types using user-provided marker gene list and plots heatmap.
}

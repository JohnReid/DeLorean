#' Single cell expression data and meta data from Trapnell et al. (2014).
#'
#' @docType data
#' @keywords datasets
#'
#' @name trapnell.expr
#' @aliases trapnell.cell.meta
#' @aliases trapnell.gene.meta
#'
#' @usage data(TrapnellDeLorean)
#'
#' @format There are three objects in this data:
#' \itemize{
#'   \item trapnell.expr A matrix of log expression values with
#'     no missing data. Rows are named by genes and columns are
#'     named by cells/samples.
#'   \item trapnell.gene.meta A data frame containing meta-data
#'     about the genes.
#'   \item trapnell.cell.meta A data frame containing meta-data
#'     about the cells
#' }
#'
#' @source \url{http://www.nature.com/nbt/journal/v32/n4/full/nbt.2859.html}
#'
NULL

#' Single cell expression data and meta data from McDavid et al. (2014).
#' They investigated differential expression in actively
#' cycling cells: "expression of 333 genes was interrogated in 930
#' cells, across three cell lines: H9 (HTB-176), MDA-MB-231 (HTB-26),
#' and PC3 (CRL-1435)".
#'
#' @docType data
#' @keywords datasets
#'
#' @name mcdavid.expr
#' @aliases mcdavid.cell.meta
#' @aliases mcdavid.gene.meta
#'
#' @usage data(McDavidDeLorean)
#'
#' @format There are three objects in this data:
#' \itemize{
#'   \item mcdavid.expr A matrix of log expression values with
#'     no missing data. Rows are named by genes and columns are
#'     named by cells/samples.
#'   \item mcdavid.gene.meta A data frame containing meta-data
#'     about the genes.
#'   \item mcdavid.cell.meta A data frame containing meta-data
#'     about the cells
#' }
#'
#' @source \url{http://www.ploscompbiol.org/article/info\%3Adoi\%2F10.1371\%2Fjournal.pcbi.1003696}
#'
NULL

#' Single cell expression data and meta data from Guo et al. (2012).
#' They investigated the expression of 48 genes in 500 mouse embryonic cells.
#'
#' @docType data
#' @keywords datasets
#'
#' @name guo.expr
#' @aliases guo.cell.meta
#' @aliases guo.gene.meta
#'
#' @usage data(GuoDeLorean)
#'
#' @format There are three objects in this data:
#' \itemize{
#'   \item guo.expr A matrix of log expression values with
#'     no missing data. Rows are named by genes and columns are
#'     named by cells/samples.
#'   \item guo.gene.meta A data frame containing meta-data
#'     about the genes.
#'   \item guo.cell.meta A data frame containing meta-data
#'     about the cells
#' }
#'
#' @source \url{http://www.sciencedirect.com/science/article/pii/S1534580710001103}
#'
NULL

#' Kouno et al. investigated the transcriptional network controlling how
#' THP-1 human myeloid monocytic leukemia cells differentiate into
#' macrophages. They provide expression values for 45 genes in 960 single
#' cells captured across 8 distinct time points.
#'
#' @docType data
#' @keywords datasets
#'
#' @name kouno.expr
#' @aliases kouno.cell.meta
#' @aliases kouno.gene.meta
#'
#' @usage data(KounoDeLorean)
#'
#' @format There are three objects in this data:
#' \itemize{
#'   \item kouno.expr A matrix of log expression values with
#'     no missing data. Rows are named by genes and columns are
#'     named by cells/samples.
#'   \item kouno.gene.meta A data frame containing meta-data
#'     about the genes.
#'   \item kouno.cell.meta A data frame containing meta-data
#'     about the cells
#' }
#'
#' @source \url{http://genomebiology.com/2013/14/10/R118/abstract}
#'
NULL

#' Windram et al. investigated the defense response in Arabidopsis
#' thaliana to the necrotrophic fungal pathogen Botrytis cinerea.
#' They collected data at 24 time points in two conditions for
#' 30336 genes.
#'
#' @docType data
#' @keywords datasets
#'
#' @name windram.expr
#' @aliases windram.cell.meta
#' @aliases windram.gene.meta
#'
#' @usage data(WindramDeLorean)
#'
#' @format There are three objects in this data:
#' \itemize{
#'   \item windram.expr A matrix of log expression values with
#'     no missing data. Rows are named by genes and columns are
#'     named by cells/samples.
#'   \item windram.gene.meta A data frame containing meta-data
#'     about the genes.
#'   \item windram.cell.meta A data frame containing meta-data
#'     about the cells
#' }
#'
#' @source \url{http://www.plantcell.org/content/24/9/3530.long}
#'
NULL


#' Tang et al. investigated the derivation of ESCs from the inner
#' cell mass (ICM) using single cell RNA-seq and PCR. They collected
#' expression data in several cells from ESCs, the ICM, E3.5
#' and E4.5.
#'
#' @docType data
#' @keywords datasets
#'
#' @name tang.key.genes
#' @aliases tang.pcr
#' @aliases tang.pcr.cell.meta
#' @aliases tang.pcr.gene.meta
#' @aliases tang.rna.seq
#' @aliases tang.rna.seq.cell.meta
#' @aliases tang.rna.seq.gene.meta
#'
#' @usage data(TangDeLorean)
#'
#' @format There are six objects in this data:
#' \itemize{
#'   \item tang.key.genes A vector of genes named in the publication.
#'   \item tang.rna.seq A matrix of log expression values with
#'     no missing data. Rows are named by genes and columns are
#'     named by cells/samples.
#'   \item tang.rna.seq.gene.meta A data frame containing meta-data
#'     about the genes.
#'   \item tang.rna.seq.cell.meta A data frame containing meta-data
#'     about the cells
#'   \item tang.pcr A matrix of log expression values with
#'     no missing data. Rows are named by genes and columns are
#'     named by cells/samples.
#'   \item tang.pcr.gene.meta A data frame containing meta-data
#'     about the genes.
#'   \item tang.pcr.cell.meta A data frame containing meta-data
#'     about the cells
#' }
#'
#' @source \url{http://www.plantcell.org/content/24/9/3530.long}
#'
NULL


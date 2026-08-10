# Utility functions for somatic mutation analysis.

#' Add a mutation ID to a VCF data frame
#' @param vcf_df A data frame containing VCF data with columns CHROM, POS, REF, and ALT.
#' @return The input data frame with an additional column 'mut_id' in the format "CHROM_POS:REF>ALT".
.add_mut_id <- function(vcf_df) {
    vcf_df |>
        mutate(mut_id = glue("{CHROM}_{POS}:{REF}>{ALT}"))
}

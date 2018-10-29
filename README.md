# WheatFlagLeafSenescence
Scripts for the manuscript 'Identification of transcription factors regulating senescence in wheat through gene regulatory network modelling'

## 1) Mapping
[Map samples with kallisto](scripts/01_kallisto_control.pl)

## 2) Merge samples
[Merge tpm and counts into a single table for all samples](scripts/02_tximport_summarise_counts_tpm_per_gene.R)

## 3) Differential expression
### a) With gradient tool
[Prepare data for gradient tool on cyverse](scripts/03ai_prep_data_tpms_gradienttool_cyverse.R)

[Arrange output data from gradient tool into patterns for zscore > |2|](scripts/03aii_gradient_tool_arrange_output_to_patterns.R)

### b) With ImpulseDE2.
[ImpulseDE2](scripts/03b_ImpulseDE_control.R)
This also merges together the gradient tool and ImpulseDE2 results

## 4) GO enrichment of patterns
[GO term enrichment of grouped patterns with goseq](scripts/04_impulseDE_and_gradient_tool_DE_cluster_genes_GO_enrichment_control_tpm.R)

## 5) Prepare data for CSI
[Prepare data for CSI](scripts/05_prep_data_CSI.R)

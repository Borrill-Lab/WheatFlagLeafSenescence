# WheatFlagLeafSenescence
Scripts for the manuscript 'Identification of transcription factors regulating senescence in wheat through gene regulatory network modelling'

## 1) Mapping
[Map samples with kallisto](scripts/kallisto_control.pl)

## 2) Merge samples
[Merge tpm and counts into a single table for all samples](scripts/tximport_summarise_counts_tpm_per_gene.R)

## 3) Differential expression
### a) With gradient tool
[Prepare data for gradient tool on cyverse](scripts/prep_data_tpms_gradienttool_cyverse.R)
[Arrange output data from gradient tool into patterns for zscore > |2|](scripts/gradient_tool_arrange_output_to_patterns.R)

### b) With ImpulseDE2.
[ImpulseDE2](scripts/ImpulseDE_control.R)
This also merges together the gradient tool and ImpulseDE2 results

# Pathogenic SNV filter for diagnostic
# Filtering for SNV with ethnic group frequency as threshold with gnomAD frequency

# Prediction of ethnic group of individual
Bayesian naive classification of the seven gnomAD population class.
AFR: African
AMR: Admixed american
ASJ: Askenazi Jewish
EAS: East Asian
FIN: Finnish
NFE: Non-Finnish European
SAS: South Asian
Sum of frenquency log of each individual SNV

# Filtering with appropriate gene threshold
Set ethnic group threshold with most frequent pathogenic SNV for each gene. 
The pathogenicity of SNV is given by Clinvar data base. 
Add flag for SNV with global frequency greater than threshold for the SNV gene for the population of the invidual.

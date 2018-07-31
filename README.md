
# Ethnicity filter
## About
Prediction of the individual's ethnic group is made by naive Bayesian classification of variants among the seven classes of [GnomAD](https://macarthurlab.org/2017/02/27/the-genome-aggregation-database-gnomad/).
SNVs with population frequency greater than threshold for the gene of  SNV for the population of the invidual are annotated.
## Installation
```bash
git clone URL
python3 prediction.py
```

# Results summary

## Primary conclusion

The full model shows moderate external predictive ability when all required
predictors are observed. Its strict no-imputation blind validation included 157
records and yielded Pearson r² = 0.511, RMSE = 0.268, and MAE = 0.207
genomic-LRV units.

## Internal validation

| Model | Apparent Pearson r² | Stratified 10-fold CV Pearson r² | CV RMSE | CV MAE |
|---|---:|---:|---:|---:|
| Full model, n = 243 | 0.942 | 0.454 | 0.364 | 0.242 |
| Reduced model, n = 430 | 0.812 | 0.342 | 0.386 | 0.273 |

The difference between apparent and cross-validated performance indicates that
the training fit should not be interpreted as generalization performance.

## Three blind-test strategies

| Strategy | n | Pearson r² | RMSE | MAE | Bias |
|---|---:|---:|---:|---:|---:|
| Strict no-imputation full model | 157 | 0.511 | 0.268 | 0.207 | −0.046 |
| Imputation-assisted full model | 300 | 0.339 | 0.315 | 0.232 | −0.080 |
| Reduced model without `Delta_Pi` | 300 | 0.335 | 0.307 | 0.225 | −0.011 |

The strict analysis is the primary external validation because training and test
inputs have the same completeness requirements. The other two strategies assess
coverage and applicability.

## Effect of imputation

Among 300 evaluable blind-test records, 157 were complete and 143 required at
least one imputed predictor. Complete records retained Pearson r² = 0.511,
RMSE = 0.268, and MAE = 0.207. Imputed records yielded Pearson r² = 0.157,
RMSE = 0.360, MAE = 0.259, and bias = −0.117.

Imputation therefore increased coverage but did not recover the performance
observed for complete inputs. These results are labelled
**imputation-assisted validation**, not strict external validation.

## Treatment-specific performance

In the strict blind test, UV showed the most stable group-level fit
(n = 35, Pearson r² = 0.660, RMSE = 0.223, MAE = 0.188). UV-Cl also showed a
relatively high Pearson r² but substantially larger absolute errors. The O3-Cl estimate was based
on only five complete records and should not support broad mechanism-level
inference.

## Interpretation

The model extracts shared DNA-damage signatures across pathogens, with
`Delta_Breadth` providing the dominant predictive contribution. Performance is
process-dependent and is constrained by predictor availability, treatment
mechanism, and response range. The model is best interpreted as a DNA-signal
process-screening tool rather than a surrogate for viable or infectious
pathogen removal.

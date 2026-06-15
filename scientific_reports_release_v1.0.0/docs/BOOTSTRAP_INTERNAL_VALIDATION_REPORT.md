# Bootstrap internal validation

Five hundred bootstrap resamples were requested. Penalized Cox models used the same fixed predictors in every resample.
The optimism-corrected C-index is reported as an internal-validation estimate, not external validation.

|Cohort|Model|n/events|Apparent C|Mean optimism|Corrected C|Reliability|
|---|---|---:|---:|---:|---:|---|
|TARGET|log_mrd|141/24|0.679|-0.001|0.680|exploratory|
|TARGET|log_mrd+grs|141/24|0.717|0.006|0.711|exploratory|
|single_center_historical|log_mrd|69/6|0.591|0.030|0.561|severely_underpowered|
|single_center_historical|log_mrd+grs|56/5|0.672|0.077|0.596|severely_underpowered|

Calibration curves and calibration slopes were not interpreted because both datasets contain far fewer than 100 EFS events, and the historical cohort contains only five complete-case events.
# Harmonized reanalysis audit

- One parser was applied to raw TARGET, historical single-center, and contemporary single-center karyotype strings.
- Historical raw-to-analytic row alignment passed exact checks for source, diagnosis date, sex, age, and WBC.
- The report-token burden variable is explicitly not a standard complex-karyotype classification.
- Matched complete-case C-index models used ridge penalizer 0.5 and 1000 paired bootstrap resamples.
- TARGET bootstrap results are conditional internal assessment of a fixed, post hoc candidate score; feature selection was not repeated.

- TARGET: matched n=141, events=24; apparent delta C=0.040; optimism-corrected delta C=0.034; paired bootstrap test delta median=0.043 (95% interval 0.022 to 0.048).
- Single-center historical: matched n=56, events=5; apparent delta C=0.000; optimism-corrected delta C=-0.002; paired bootstrap test delta median=0.000 (95% interval 0.000 to 0.007).
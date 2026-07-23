# Generated article fragments

`figure1_addplots.tex` is generated from the raw MacBook full-SpGEMM CSV:

```bash
python3 scripts/make_figure1.py   results/macbook/fullspgemm.csv   > generated/figure1_addplots.tex
```

The file should contain exactly three PGFPlots `\addplot` blocks and their legend entries. The generated file is derivative output; the raw CSV and script remain the primary reproducibility materials.

# Global Stratification Economics and International Extraction

This repository contains the replication materials for the manuscript **“Global Stratification Economics and International Extraction.”**

The replication package is organized around two Stata files:

- `01_clean_all_big_data.do` reconstructs the compact analysis datasets from the original WIOD, SEA, PPP, and Penn World Table source files.
- `02_estimation.do` starts from the compact posted Stata datasets and reproduces the empirical calibration, robustness checks, and figures.

The Mathematica/Wolfram Language companion file is kept separate and is not modified by the Stata workflow.

---

## Repository structure

```text
global-stratification-extraction-replication/
│
├── README.md
│
├── code/
│   ├── 01_clean_all_big_data.do
│   └── 02_estimation.do
│
├── calibration_1995_2009_aggregate/
│   ├── calib_inputs_1995_2009_table2_reg.dta
│   └── calib_inputs_1995_2009_table2_ROW.dta
│
├── calibration_2016release_2000_2014/
│   ├── calib_inputs_2000_2014_table2_reg.dta
│   └── calib_inputs_2000_2014_table2_ROW.dta
│
└── mathematica/
    ├── GSE_companion2.wl
    ├── brendan_wiod2013_1995_2009.csv
    └── brendan_wiod2016_2000_2014.csv
```

The two `calibration_*` directories contain the compact Stata datasets needed to reproduce the empirical estimates. The original large WIOD and SEA files are **not required to run `02_estimation.do`**.

> **Important:** The directory names above are the directory names currently used by `02_estimation.do`. If these folders are renamed, the corresponding paths in the do-file must also be changed.

---

## Quick replication

The quickest way to reproduce the empirical results is to use the compact Stata datasets included in this repository.

### Requirements

- Stata 17 or later
- `boottest` for wild-cluster bootstrap inference

Install `boottest` in Stata if needed:

```stata
ssc install boottest
```

### Run the estimation

1. Download or clone the repository.
2. Open Stata.
3. Change the working directory to the **root of the repository**, not the `code` folder.

For example:

```stata
cd "/path/to/global-stratification-extraction-replication"
```

4. Run:

```stata
do "code/02_estimation.do"
```

`02_estimation.do` reads only the compact `.dta` files in:

```text
calibration_1995_2009_aggregate/
calibration_2016release_2000_2014/
```

It does **not** require the original large WIOD or SEA files.

---

## Full data construction

Researchers who wish to reconstruct the compact calibration datasets from the original source data can run:

```stata
do "code/01_clean_all_big_data.do"
```

This file combines the original upstream Ricci-style data construction for both WIOD vintages with the calibration-data construction used in the paper.

Before running `01_clean_all_big_data.do`, create a folder called:

```text
raw/
```

at the root of the replication directory and place the original source files there.

The expected local structure is:

```text
global-stratification-extraction-replication/
│
├── code/
│   └── 01_clean_all_big_data.do
│
└── raw/
    ├── wiot_full.dta
    ├── Socio_Economic_Accounts_July14.xlsx
    ├── Socio_Economic_Accounts.xlsx
    ├── P_Data_Extract_From_World_Development_Indicators.xlsx
    ├── dataset_2026-08-04T19_36_35.328252551Z_DEFAULT_INTEGRATION_IMF.RES_WEO_9.0.0.csv
    ├── pwt110.dta
    ├── WIOT2000_October16_ROW.dta
    ├── WIOT2001_October16_ROW.dta
    ├── ...
    └── WIOT2014_October16_ROW.dta
```

For the WIOD 2016 annual files, the construction code also accepts filenames ending in `(1)` or `(2)`, such as:

```text
WIOT2000_October16_ROW(1).dta
```

The large original data are not included in the public repository because they are substantially larger than the compact replication inputs and remain available from their original providers.

---

## Data sources

The empirical calibration uses the following sources.

### World Input-Output Database (WIOD)

Two WIOD vintages are used:

- the earlier WIOD release underlying the 1995–2009 calibration;
- the November 2016 WIOD release underlying the 2000–2014 robustness calibration.

The first construction uses `wiot_full.dta`. The second uses annual WIOD files from 2000 through 2014.

### WIOD Socio Economic Accounts (SEA)

The data construction uses:

- `Socio_Economic_Accounts_July14.xlsx` for the earlier WIOD release;
- `Socio_Economic_Accounts.xlsx` for the November 2016 release.

These files provide the labor, compensation, employment, hours, and value-added information required by the unequal-exchange construction.

### Penn World Table 11.0

`pwt110.dta` is used to construct regional capital stocks, employment, output, labor shares, and saving/investment measures used in the calibration.

### Purchasing-power-parity data

The construction also uses:

- World Bank World Development Indicators PPP data; and
- IMF World Economic Outlook PPP information for Taiwan.

---

## Empirical construction

Ricci-style net transfers are coded so that a provider of resources to the rest of the world has a **negative** net transfer.

The model's extraction rate is therefore defined as:

```text
e = - net_transfer_pct_va / 100
```

so that extraction is positive for regions experiencing a net outflow.

The principal estimating equation is the transformed version of the model's extraction schedule:

```text
ln(m/e + λδ) = ln(χ₀) - θ ln(Ω)
```

where:

```text
m = 1 + λ(1-b)
```

and the baseline parameters are:

```text
λ = 0.5
b = 0.5
δ = 0.5
```

The paper's preferred power measure is the aggregate capital ratio:

```text
Ω = K_core / K_region
```

Capital per worker is retained as a robustness measure.

The estimation sample is restricted to non-core regions that are net transfer providers in every year of the relevant sample period.

---

## Replication targets

The main annual calibration results reported in the manuscript are:

| WIOD vintage | Years | Power measure | θ |
|---|---:|---|---:|
| WIOD 2013 | 1995–2009 | Aggregate capital | 0.495 |
| WIOD 2013 | 1995–2009 | Capital per worker | 0.349 |
| WIOD 2016 | 2000–2014 | Aggregate capital | 0.443 |
| WIOD 2016 | 2000–2014 | Capital per worker | 0.376 |

The preferred specification is the WIOD 2013, 1995–2009 estimate using aggregate capital.

`02_estimation.do` also reproduces:

- conventional and region-clustered inference;
- wild-cluster bootstrap inference using Webb weights;
- the parameter-grid robustness exercise;
- the aggregate-capital and capital-per-worker specifications;
- the empirical extraction-schedule figures; and
- the WIOD 2016 robustness exercise.

---

## Mathematica / Wolfram Language companion

The theoretical model and associated model-generated figures are reproduced by:

```text
mathematica/GSE_companion2.wl
```

The Mathematica file has been left unchanged.

For the empirical sections of the Mathematica companion, keep these files in the same directory as `GSE_companion2.wl`:

```text
brendan_wiod2013_1995_2009.csv
brendan_wiod2016_2000_2014.csv
```

Run the Mathematica companion from the `mathematica/` directory so that its relative file references resolve correctly.

---

## Which file should I run?

For reproducing the results in the paper from the posted replication data:

```stata
do "code/02_estimation.do"
```

is sufficient.

`01_clean_all_big_data.do` is provided to document and reproduce the construction of the compact analysis files from the original large source datasets. It is not required for ordinary replication of the reported estimates.

---

## Notes on reproducibility

All Stata paths are relative to the repository root. The do-files use Stata's current working directory as `ROOT`, so Stata should be opened or changed to the repository root before either do-file is executed.

The large WIOD and SEA source files should not be placed in the public GitHub repository. Researchers reconstructing the compact data should obtain the source files from the original data providers and place them in the local `raw/` directory using the filenames expected by `01_clean_all_big_data.do`.

The compact `.dta` files included in the two `calibration_*` directories are the files needed for the referee-facing replication.

---

## Citation

If using these replication materials, please cite the associated manuscript:

> *Global Stratification Economics and International Extraction.*

A permanent repository DOI can be added here after the replication package is archived.

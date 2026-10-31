# Reproduction details

## Source data

The two compressed files in `inputs/` reconstruct the European coefficient database (104,300,544 bytes; SHA-256 `f82ee0744c2eb30c6b72f5fbc993a42645ebb09264917f9b940198a02c36497d`). Individual checksums are in `inputs/checksums.sha256`. Source attribution and the original data license are recorded in `inputs/SOURCE.json`.

The database covers 23 EU countries, 87,644 modeled spatial decision units and 13 management alternatives per unit. Modeled annual land value combines crop, grazing and forestry returns less annualized transition costs, expressed in 2015 US dollars. The climate index is carbon stock multiplied by 44/12 minus annual methane multiplied by 0.02, following the source convention. It combines a stock with a 20-year flow component, not observed annual emissions. Biodiversity scores are country-specific modeled indices and are not aggregated across countries.

## Run the calculations

The README links to separate **Quick** and **Full** Colab notebooks. Both produce `outputs.zip`. Local execution from the repository root uses:

```sh
python -m pip install -r requirements.txt
python code/reproduce.py
python code/reproduce.py --full
```

The first calculation rebuilds tables and figures from released numerical certificates. The second reconstructs the source database, checks all original actions and verifies the policies against the source coefficients. It uses `verify_matrices.py` for policy responses, climate-only coverage, gains-only exclusions, participation diagnostics and resource certificates. It then evaluates the Belgian safeguard sensitivity through `safeguard_sensitivity.py` and `verify_safeguard.py`. The calculation verifies published policy constructions; it does not rerun the original policy-search optimizers.

Outputs are written to `outputs/` and `outputs.zip`. The option `--output PATH` selects another results directory. For Full, `--database PATH` accepts an uncompressed copy of the same verified SQLite database. Source files are not overwritten. `code/check.py` runs the finite mathematical examples and checksum checks; `tests/test_regression.py` compares the regenerated numerical records with the reference results.

## Implementation and numerical precision

Country target checks use a tolerance of `10^-8 × max(1, |target|)`. The response sets include exact maximizers and payoff bands of `10^-6` and `10^-4` times `max(1, max_s |E_is|)` for each unit. Complete shared-price constructions include a conservative allowance for coefficient rounding.

Climate-only coverage is calculated from exact rational upper envelopes of the stored float64 coefficients, covering open price intervals, finite breakpoints, zero and the unbounded tail. Tied actions at a breakpoint can prevent universal compliance even when adjacent intervals comply. The implementation tests therefore use the complete response set at each actual policy coefficient.

At 5% biodiversity headroom, the exact minimum number of biodiversity complements is 20, 14, 11, 9, 7, 5 and 4 for at most one through seven distinct climate coefficients. All seven constructions meet the specified choice bands. At 2.5% headroom, some policies that satisfy exact choice or the smaller band fail the `10^-4` band. Under that larger band, the national minimum remains bounded between three and four. These are calculations on a fixed source technology, not estimates of parameter uncertainty.

The original national construction uses four biodiversity complements. The original common-coefficient construction, retained for financial and resource comparisons, uses coefficient 0.1 and 21 complements. A different construction with coefficient 0.18 attains the revised common minimum of 20. The original policies reach different ecological outcomes, so their transfer difference is not a fiscal cost of harmonization. The nominal count in the shared-price frontier adds distinct climate coefficient values to positive national biodiversity coefficients; it is not an expenditure measure.

## Belgian safeguard and fixed country groups

`data/safeguard_frontier.csv` reports the frontier after reducing only Belgium's modeled land-value floor by 0%, 0.5%, 1%, 3% or 5%. All other national safeguards and ecological requirements are fixed. With a 3% or 5% reduction, the exact complement counts become 19, 13, 10, 8, 6, 5 and 4. At 1%, the first five entries retain lower and upper integer bounds. The source-level calculations also verify the Belgian climate-only coefficient ceiling of 0.269. The results are not confidence intervals.

A separate comparison assigns one climate coefficient to CEE-11 and another to the remaining twelve countries. At 5% headroom, this fixed partition requires 19 complements, against 14 for freely assigned groups. Slovenia excludes coefficients of at least 0.6937438952159174, while Belgium's safeguard restricts the other group's coefficient to less than 0.269. Together these restrictions allow at most four climate-only countries. A complete construction using coefficients 0.27 and 0.18 satisfies all 69 country-response checks. Certificates and policy assignments are in `data/fixed_groups/`.

## Figures and records

`code/figures.py` regenerates Figure 1, Figure 2, Figure 3, Figure 4 and Figure S1 from the numerical tables. SVG, PNG and TIFF share the same plot geometry; PNG and TIFF have identical RGB pixels. The default raster resolution is 600 dpi. For 1000 dpi artwork, run:

```sh
python code/figures.py --data outputs/data --output outputs/submission_figures --dpi 1000
```

This also produces vector PDFs. Figure captions are supplied with the article, not embedded in the artwork. The CSV and JSON records in `data/` and the certificates in `certificates/` retain favorable results, exclusions and unresolved bounds.

# Climate incentives and biodiversity complements

Code, source coefficients and numerical evidence for restricted climate and biodiversity incentives in the land-use model of 23 EU countries studied by Polasky et al. (2026). At 5% biodiversity headroom, the minimum number of countries requiring biodiversity complements is 20, 14, 11, 9, 7, 5 and 4 with at most one through seven distinct climate coefficients. With two coefficients, the fixed CEE-11/other-country partition requires 19 complements, compared with 14 under unrestricted grouping. The analysis also examines national economic safeguards and biodiversity payment bases.

## Reproduce

[Open in Google Colab](https://colab.research.google.com/github/calinadriancomes/climate-biodiversity-complements/blob/main/notebooks/reproduce.ipynb)

Run the notebook cells in order. The default mode checks the released scientific inputs, rebuilds the principal tables and five figures, and downloads `outputs.zip`. For calculations directly from the original European coefficients, set `RUN_FULL_SOURCE_REPLAY = True` before running the cells. The notebook retrieves the repository automatically in Colab and also supports local execution.

For local reproduction with Python 3.12, from the repository root:

```sh
python -m pip install -r requirements.txt
python code/reproduce.py --full
```

Omit `--full` for the shorter calculation from published numerical certificates. Results are saved in `outputs/` and `outputs.zip`. See [reproduction details](docs/reproduction.md) for the source technology, numerical conventions and certificate checks.

## Files

- `notebooks/reproduce.ipynb`: the single reproduction notebook.
- `code/` and `tests/`: analysis, verification and numerical tests.
- `inputs/`: the exact European coefficient database in two compressed parts, provenance and checksums.
- `certificates/`: policy assignments and mathematical certificates.
- `data/` and `results/`: numerical tables and regression references.
- `figures/`: reproducible SVG, PNG and TIFF figures.

## Sources and licenses

The European coefficient data derive from [Polasky et al. (2026)](https://doi.org/10.1126/science.aea9058) and the [Dryad dataset](https://doi.org/10.5061/dryad.b2rbnzsvt), released under CC0. Original code is licensed under MIT; original scholarly text, derived materials and figures under CC BY 4.0. See the [license map](LICENSES/README.md) and [third-party attribution](LICENSES/third-party.md).

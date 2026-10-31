# Reproduction details

## Source technology

The two compressed files in `inputs/` reconstruct the exact European coefficient database. Their concatenation expands to 104,300,544 bytes with SHA-256 `f82ee0744c2eb30c6b72f5fbc993a42645ebb09264917f9b940198a02c36497d`. The individual part identities and scientific-certificate identities appear in `inputs/checksums.sha256`. Author-editable narrative, code and license files are not frozen by this data manifest.

The database covers 23 EU countries, 87,644 modeled spatial units and 13 alternatives per unit. Each source action is retained. The source verifier reconstructs net economic value from crop, grazing and forestry values minus annualized transition cost. It reconstructs the climate index as carbon times 44/12 minus annual methane times 0.02. Monetary quantities are 2015 US dollars per year. The climate index combines a stock term and a conventional 20-year methane component; it is not an annual emissions flow. Biodiversity scores are national modeled indices and are not added across countries.

The database retains original scenario configuration and provenance. Paths within those source records are historical provenance and are not used to access files. The supplied extract is sufficient for all documented calculations. An exact uncompressed copy can optionally be supplied with `--database`.

## Calculation routes

`python code/check.py` checks the scientific input identities and runs the independent finite mathematical examples. Tests include complete response sets, persistent and isolated ties, exact interval coverage, a required low climate coefficient, regional dominance, payment recentering and a counterexample in which a wider choice band fails.

`python code/reproduce.py` recomputes the original country accounts and the numerical article tables, then draws the figures. This shorter route uses the released numerical certificates and the stored shared-price results. It is not original-coefficient replay.

`python code/reproduce.py --full` also reconstructs the source database and performs these calculations:

- `verify_coefficients.py` evaluates all original alternatives for the 382 original policy assignments and 1,146 response-band checks. It verifies dominance witnesses, support certificates, scalar hulls and matched-outcome resource bounds.
- `verify_breadth.py` independently reconstructs climate-only upper envelopes using exact rational representations of the stored float64 technology. It sums national safeguard outcomes exactly and calculates cross-country coverage at 0%, 2.5% and 5% biodiversity headroom.
- `verify_additional.py` reconstructs source-unit definitions, target-sensitive gains-only exclusions, participation diagnostics and a sufficient coefficient-neighborhood bound.
- `verify_weak_exclusions.py` evaluates optimistic coordinate maxima at zero, every open price cell, every finite breakpoint and the unbounded tail for the four primary climate-only exclusions.
- `price_groups.py` proves the Belgian exclusion for every climate coefficient at least 0.269 and every nonnegative biodiversity coefficient. It enumerates the exact shared-price coverage problem subject to at least one coefficient below that ceiling, then checks 19 complete country-policy constructions against the original actions in all three response bands.
- `belgium_economic_diagnostic.py` sweeps all biodiversity-price cells at the original climate-only coverage-maximizing climate coefficient. It calculates the maximum guaranteed Belgian land value subject to its biodiversity and climate safeguards.
- `fixed_groups.py` proves the lower bound and checks the complete 19-complement construction when CEE-11 and the other twelve countries each share one climate coefficient.

The reproduction evaluates policy witnesses and numerical certificates directly. It does not rerun the original policy-search optimizers. The new frontier is established through exact coverage enumeration and attained feasible constructions. The numerical records are compared with the original coefficient-verification reference and the released shared-price results.

All commands resolve scientific inputs relative to the repository, so they also work when invoked from another working directory. Use `--output PATH` to select a separate output directory. Source data and published reader tables are never overwritten. The notebook provides separate short and full modes and one compact downloadable ZIP. The full source-coefficient stages run only when selected.

## Exact choice and choice bands

National target checks allow 10^-8 times the larger of one and the target's absolute value. The original policy replay uses its recorded float64 arithmetic. Choice bands are 0, 10^-6 and 10^-4 times each unit's scale `max(1, max_s |E_is|)`. The new complete constructions enlarge the admitted payoff set conservatively to include rounding uncertainty.

The climate-only interval export merges adjacent compliant open cells. An interior breakpoint can fail if an extra action ties there. A breakpoint's response set contains neighboring open-cell responses, so universal compliance cannot improve there. Finite unions of open-cell coverage patterns therefore establish the coverage maxima without treating every point in a merged interval as compliant. Every actual policy coefficient is checked directly against its full admitted response set.

The primary frontier is exactly 20, 14, 11, 9, 7, 5 and 4 biodiversity complements for one through seven distinct climate coefficients. Its complete constructions meet all stated choice bands. At 2.5% headroom, some low-complement constructions pass exact choice and the 10^-6 band but fail the 10^-4 band; these adverse checks remain in `data/price_groups/price_group_policy_checks.csv`. In particular, the national minimum under the larger band remains between three and four. These results concern choice tolerance on a fixed technology, not uncertainty about the source coefficients.

## Policy accounts and figures

The original national and common constructions remain in `data/policies.csv`. The original common construction uses climate coefficient 0.1 and 21 biodiversity complements. Its financial and matched-resource accounts are retained as explicit policy comparisons. The construction attaining the revised common minimum of 20 uses coefficient 0.18 and appears in `certificates/price_groups.csv`.

The two original complete policies have different attained ecological outcomes. Their transfer difference is not a matched fiscal cost of harmonization. Resource bounds instead preserve the specified reference policy's attained national biodiversity and climate outcomes. The nominal coefficient count in the new frontier adds distinct climate coefficient values and positive national biodiversity coefficients; it excludes group assignment and has no monetary interpretation.

`code/figures.py` draws all five publication figures from reader tables. Figure 1 uses the country intervals and the certified common-price ceiling. Figure 2 uses the revised complete-policy ambition bounds. Figure 3 shows the primary shared-price frontier and its nominal coefficient count. Figure 4 and Figure S1 retain the gains-only exclusion and baseline-action diagnostics. SVG, PNG and TIFF share a canvas; PNG and TIFF have identical RGB pixels. Captions are external to the figures.

Default raster reproduction uses 600 dpi to match the stored graphics. For 1000 dpi line artwork, run `python code/figures.py --data outputs/data --output outputs/submission_figures --dpi 1000`. The same command also writes vector PDFs with embedded TrueType fonts and SVG files from the unchanged figure canvases.

## Sensitivity to the Belgian land-value safeguard

The primary frontier is evaluated again after reducing only Belgium's modeled national land-value floor by 0%, 0.5%, 1%, 3% or 5%. All source actions and national ecological targets remain fixed, as do the other countries' economic safeguards. With up to seven distinct climate coefficients across all 23 countries, the minimum biodiversity-complement counts are 20, 14, 11, 9, 7, 5 and 4 at the original floor. With a 3% or 5% reduction, the exact counts become 19, 13, 10, 8, 6, 5 and 4. At 1%, the first five entries retain distinct lower and upper bounds. These are integer optimization bounds, not confidence intervals or estimates of parameter uncertainty.

The original and relaxed calculations use the same robust exact-choice and two near-best-response conventions. The complete 3% policy assignments contain 161 country records across seven constructions, with 483 recorded response checks. The independent source verifier checks the Belgian climate-only ceiling, the exact rational resource upper bounds and all complete policy assignments with extended precision.

From the repository root, the full reproduction includes this analysis through `python code/reproduce.py --full`. A standalone calculation can be run with `python code/safeguard_sensitivity.py --repository . --output outputs/safeguard`, followed by `python code/verify_safeguard.py --repository . --results outputs/safeguard --output outputs/safeguard/independent_verification.json`. The standalone commands reconstruct the unchanged source database and exact climate-only intervals when neither is supplied.

The `data/safeguard_frontier.csv` table and `certificates/safeguard_frontier_policies.csv` provide compact published records. They distinguish complete constructions from the unresolved 1% scenario and preserve the adverse German 2.5% choice-band cases elsewhere in the repository. The minimum of q + Kq at a 3% reduction is eleven for q = 5, 6 or 7, but that nominal count is not an administrative-cost estimate.

## Two fixed country groups

At 5% biodiversity headroom, fixing one climate coefficient for CEE-11 and another for the other twelve countries raises the exact minimum from 14 to 19 biodiversity complements. All original national safeguards remain in place. This descriptive partition restricts group membership; it is not an optimized regional design.

Slovenia excludes every climate coefficient at least 0.6937438952159174 under every nonnegative biodiversity coefficient. An exact rational dominance calculation bounds retained land value at USD 89,851,432.855959 annually, below its effective floor of USD 96,917,398.173443. Below this cutoff, at most one CEE-11 country can meet its safeguards without a biodiversity payment. Belgium separately limits the other group's coefficient to less than 0.269, allowing at most three climate-only countries. Thus at least 23 - 1 - 3 = 19 complements are necessary. The complete construction uses coefficients 0.27 for CEE-11 and 0.18 for the other twelve countries and passes all 69 country-response checks.

The records in `data/fixed_groups/` contain exact rational certificates, all 23 policy coefficients, original national floors and the three-band response guarantees. Full reproduction rebuilds these results from the source database and freshly reconstructed climate-only intervals. Quick reproduction checks the released records and retains them in the downloadable archive.

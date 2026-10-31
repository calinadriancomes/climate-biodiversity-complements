"""Calculate the article's numerical tables from regenerated country records."""

from pathlib import Path
import argparse
import csv
import math

ROOT = Path(__file__).resolve().parents[1]


def rows(path):
    with Path(path).open(newline="", encoding="utf-8") as stream:
        return list(csv.DictReader(stream))


def write_table(output, name, title, header, body, note):
    with (output / (name + ".csv")).open("w", newline="", encoding="utf-8") as stream:
        writer = csv.writer(stream)
        writer.writerow(header)
        writer.writerows(body)
    text = f"### {title}\n\n"
    text += "| " + " | ".join(header) + " |\n"
    text += "| " + " | ".join(["---"] * len(header)) + " |\n"
    text += "".join("| " + " | ".join(map(str, row)) + " |\n" for row in body)
    (output / (name + ".md")).write_text(text + "\n" + note + "\n", encoding="utf-8")


def main(data, output):
    data, output = Path(data), Path(output)
    output.mkdir(parents=True, exist_ok=True)
    policies = rows(data / "policies.csv")
    financial = rows(data / "financial_accounts.csv")
    labels = ["National climate coefficients", "Common climate coefficient"]
    accounts = []
    for label in labels:
        p = [r for r in policies if r["policy"] == label]
        f = [r for r in financial if r["policy"] == label]
        if len(p) != 23 or len(f) != 23:
            raise ValueError("Each complete construction needs 23 country accounts.")
        total = lambda rr, k: math.fsum(float(r[k]) for r in rr)
        accounts.append(
            [
                sum(int(r["biodiversity_complement"]) for r in p),
                total(p, "E_USD2015_per_year"),
                total(p, "C20_index_tCO2e"),
                total(p, "resource_opportunity_cost_USD2015_per_year"),
                total(p, "net_transfer_USD2015_per_year"),
                total(f, "baseline_entitlement"),
                total(f, "recentered_gross"),
                total(f, "recentered_receipts"),
                total(f, "recentered_net"),
                total(f, "outside_option_shortfall"),
                sum(int(r["units_below_outside_option"]) for r in f),
            ]
        )
    names = [
        "Positive national biodiversity coefficients",
        "Modeled land value E",
        "C20 index, billions of tCO2e-index units",
        "Resource opportunity cost versus maximum modeled E",
        "Level-contract gross and net payments",
        "Action-independent biodiversity baseline entitlement",
        "Recentered gross component payments",
        "Recentered component receipts",
        "Recentered net transfer",
        "Recentered payoff shortfall versus free outside options",
        "Modeled units below that outside option",
    ]
    body = [
        [name] + [f"{a[i]:,}" if i in (0, 10) else f"{a[i]/1e9:.3f}" for a in accounts]
        for i, name in enumerate(names)
    ]
    write_table(
        output,
        "Table_2",
        "Financial and resource accounts of the original constructions",
        ["Quantity", "Original national construction", "Original common construction"],
        body,
        "Monetary quantities are billions of 2015 US dollars per year. "
        "The two original complete constructions have different attained ecological "
        "outcomes, so their transfer difference is not a matched fiscal cost of harmonization.",
    )
    resource = rows(data / "resource_bounds.csv")
    body = []
    for country in ["Latvia", "Poland", "Greece", "Sweden"]:
        values = []
        for policy in ["national_final", "common_final"]:
            r = next(
                r for r in resource if r["country"] == country and r["policy"] == policy
            )
            values.append(
                f"{float(r['resource_gain_lower'])/1e6:.3f} to {float(r['resource_gain_upper'])/1e6:.3f}"
            )
        body.append([country] + values)
    write_table(
        output,
        "Table_3",
        "Attained-outcome resource improvement bounds",
        ["Country", "National climate construction", "Common climate construction"],
        body,
        "Millions of 2015 US dollars per year. Each bound compares a "
        "restricted policy with resource-efficient allocations that preserve "
        "at least that policy's attained national biodiversity and climate outcomes.",
    )
    breadth = rows(data / "breadth.csv")
    frontier = rows(data / "price_groups/price_group_frontier.csv")
    body = []
    for row in breadth:
        common = next(
            r
            for r in frontier
            if float(r["ambition"]) == float(row["alpha"])
            and int(r["maximum_climate_coefficients"]) == 1
            and float(r["response_band"]) == 0
        )
        national = int(row["national_B_minimum"])
        lower, upper = int(common["complement_lower_bound"]), int(
            common["complement_upper_bound"]
        )
        bounds = str(lower) if lower == upper else f"{lower} to {upper}"
        body.append(
            [f"{100*float(row['alpha']):g}%", national, bounds, lower - national]
        )
    write_table(
        output,
        "Table_S3",
        "Biodiversity ambition and complete-policy complement reach",
        [
            "Biodiversity ambition",
            "National exact minimum",
            "Common minimum bounds",
            "Minimum increase",
        ],
        body,
        "Every exact best response must meet all national safeguards. The common minimum is identified "
        "at baseline preservation and 5% headroom. At 2.5%, the lower bound and feasible upper bound differ.",
    )
    primary = sorted(
        [
            r
            for r in frontier
            if float(r["ambition"]) == 0.05 and float(r["response_band"]) == 1e-4
        ],
        key=lambda r: int(r["maximum_climate_coefficients"]),
    )
    body = [
        [
            r["maximum_climate_coefficients"],
            r["complement_lower_bound"],
            r["complement_upper_bound"],
            str(
                int(r["maximum_climate_coefficients"])
                + int(r["complement_lower_bound"])
            ),
        ]
        for r in primary
    ]
    write_table(
        output,
        "Climate_coefficient_frontier",
        "Shared climate coefficients and biodiversity complements",
        [
            "Distinct climate coefficients",
            "Complement lower bound",
            "Complement upper bound",
            "Total nominal coefficients",
        ],
        body,
        "The primary 5% headroom frontier is attained under all specified response bands through 10^-4. "
        "The nominal count adds distinct climate coefficient values to positive national biodiversity coefficients. "
        "It does not include the country-to-group assignment or measure administrative expenditure.",
    )
    print(
        "Calculated the original financial accounts, matched resource bounds, ambition comparison and policy frontier."
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data", type=Path, default=ROOT / "data")
    parser.add_argument("--output", type=Path, default=ROOT / "outputs/tables")
    args = parser.parse_args()
    main(args.data, args.output)

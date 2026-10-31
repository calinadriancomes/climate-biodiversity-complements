"""Publication figures from full-precision reader tables. No optimization."""

from pathlib import Path
import argparse, csv, re, json, io, os, tempfile
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
from PIL import Image

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data"
BLUE = "#27627D"
ORANGE = "#A4542D"
GRAY = "#606469"
LIGHT = "#DEE1E3"
LAYOUT = []
RASTER_DPI = 600
plt.rcParams.update(
    {
        "font.family": "DejaVu Sans",
        "font.size": 9,
        "axes.labelsize": 9,
        "xtick.labelsize": 8,
        "ytick.labelsize": 8.8,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "axes.linewidth": 0.65,
        "xtick.major.width": 0.6,
        "ytick.major.width": 0.6,
        "svg.fonttype": "path",
        "pdf.fonttype": 42,
        "savefig.facecolor": "white",
    }
)


def read(n):
    with (DATA / n).open(newline="", encoding="utf-8-sig") as f:
        return list(csv.DictReader(f))


def write_complete(path, payload):
    """Make a complete raster visible only after its bytes reach the file."""
    path = Path(path)
    temporary = None
    try:
        with tempfile.NamedTemporaryFile(dir=path.parent, prefix=".figure_", delete=False) as stream:
            temporary = Path(stream.name)
            stream.write(payload)
            stream.flush()
            os.fsync(stream.fileno())
        temporary.replace(path)
        if path.read_bytes() != payload:
            raise OSError("Incomplete figure export: " + path.name)
    finally:
        if temporary is not None and temporary.exists():
            temporary.unlink()


def finish(fig, out, name):
    # One figure canvas supplies all formats. TIFF uses exactly the PNG pixels.
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    width, height = fig.get_size_inches() * fig.dpi
    omitted = set()
    for axes in fig.axes:
        for axis in [axes.xaxis, axes.yaxis]:
            lo, hi = sorted(axis.get_view_interval())
            for tick in axis.get_major_ticks() + axis.get_minor_ticks():
                if not lo <= tick.get_loc() <= hi:
                    omitted.update([tick.label1, tick.label2])
    outside = []
    for text in fig.findobj(matplotlib.text.Text):
        if not text.get_visible() or not text.get_text() or text in omitted:
            continue
        box = text.get_window_extent(renderer)
        if (
            box.x0 < -0.5
            or box.y0 < -0.5
            or box.x1 > width + 0.5
            or box.y1 > height + 0.5
        ):
            outside.append(text.get_text())
    if outside:
        raise ValueError("Figure text outside canvas: " + name + ": " + repr(outside))
    LAYOUT.append(
        {
            "figure": name,
            "width_inches": float(fig.get_size_inches()[0]),
            "height_inches": float(fig.get_size_inches()[1]),
            "raster_dpi": RASTER_DPI,
            "text_inside_canvas": True,
        }
    )
    fig.savefig(
        out / (name + ".svg"),
        metadata={"Creator": None, "Date": None},
        facecolor="white",
    )
    sv = out / (name + ".svg")
    st = sv.read_text()
    st = re.sub(r"<metadata>.*?</metadata>", "", st, flags=re.S)
    st = re.sub(r"<!--.*?-->", "", st, flags=re.S)
    sv.write_text(st)
    svgpath = out / (name + ".svg")
    svgpath.write_text(svgpath.read_text().replace("matplotlib.axis_", "axis_"))
    fig.savefig(
        out / (name + ".pdf"),
        metadata={
            "Creator": None,
            "Producer": None,
            "CreationDate": None,
            "ModDate": None,
        },
        facecolor="white",
    )
    png_buffer = io.BytesIO()
    fig.savefig(png_buffer, format="png", dpi=RASTER_DPI,
                metadata={"Software": None}, facecolor="white")
    with Image.open(io.BytesIO(png_buffer.getvalue())) as im:
        raster = im.convert("RGB")
        tiff_buffer = io.BytesIO()
        raster.save(tiff_buffer, format="TIFF", compression="tiff_lzw", dpi=(RASTER_DPI, RASTER_DPI))
    write_complete(out / (name + ".png"), png_buffer.getvalue())
    write_complete(out / (name + ".tiff"), tiff_buffer.getvalue())
    with Image.open(out / (name + ".png")) as png, Image.open(out / (name + ".tiff")) as tiff:
        if not np.array_equal(np.asarray(png.convert("RGB")), np.asarray(tiff.convert("RGB"))):
            raise AssertionError("PNG and TIFF pixels differ: " + name)
    plt.close(fig)


def main(out, data=ROOT / "data", dpi=600):
    global DATA, RASTER_DPI
    DATA = Path(data)
    if dpi <= 0:
        raise ValueError("Raster resolution must be positive")
    RASTER_DPI = int(dpi)
    out = Path(out)
    LAYOUT.clear()
    out.mkdir(exist_ok=True, parents=True)
    extension = json.loads(
        (DATA / "price_groups/price_group_verification.json").read_text()
    )
    ceiling = json.loads(
        (DATA / "price_groups/climate_ceiling_certificate.json").read_text()
    )
    common = float(extension["primary_common_climate_coefficient"])
    cutoff = float(ceiling["climate_cutoff"])
    frontier = read("price_groups/price_group_frontier.csv")
    scope = [r for r in read("scope.csv") if r["included"] == "True"]
    regs = read("price_regions.csv")
    fig = plt.figure(figsize=(7.2, 5.1))
    gs = fig.add_gridspec(
        1, 2, left=0.145, right=0.985, bottom=0.18, top=0.82, wspace=0.56
    )
    for pi, group in enumerate(["CEE-11", "Other 12"]):
        ax = fig.add_subplot(gs[pi])
        rr = sorted(
            [r for r in scope if r["group"] == group], key=lambda r: r["country"]
        )
        N = len(rr)
        ax.set_title(
            "a  CEE-11" if pi == 0 else "b  Other twelve countries",
            loc="left",
            fontweight="bold",
            fontsize=10,
            pad=13,
        )
        for k, r in enumerate(rr):
            y = 11 - k
            iv = [q for q in regs if q["country"] == r["country"]]
            ax.axhline(y, color="#F0F1F1", lw=0.55, zorder=0)
            if iv:
                for q in iv:
                    a, b = float(q["a"]), float(q["b"])
                    ax.plot(
                        [a, b],
                        [y, y],
                        color=BLUE,
                        lw=2.5,
                        solid_capstyle="butt",
                        zorder=3,
                    )
                    # Endpoint ticks remain visible even for a very short interval.
                    ax.plot(
                        [a, b],
                        [y, y],
                        linestyle="none",
                        marker="|",
                        ms=5,
                        color=BLUE,
                        mew=0.9,
                        zorder=4,
                    )
            else:
                ax.text(
                    0.33,
                    y,
                    "No region",
                    ha="left",
                    va="center",
                    fontsize=8.2,
                    color=GRAY,
                )
        ax.axvline(common, color=ORANGE, ls=(0, (3, 2)), lw=0.9, zorder=2)
        ax.axvline(cutoff, color=GRAY, ls=(0, (1, 2)), lw=0.9, zorder=2)
        ax.set_yticks([11 - i for i in range(N)], [r["country"] for r in rr])
        ax.set_ylim(-0.7, 11.7)
        ax.set_xlim(-0.06, 5.5)
        ax.set_xticks(range(6))
        ax.tick_params(axis="y", length=0, pad=5)
        ax.spines["left"].set_visible(False)
        ax.grid(axis="x", color=LIGHT, lw=0.5)
        ax.set_axisbelow(True)
    fig.supxlabel(
        "Climate coefficient (annual USD2015 per C20-index unit)",
        x=0.56,
        y=0.075,
        fontsize=9,
    )
    handles = [
        Line2D([0], [0], color=BLUE, lw=2.5, label="Climate-only outer interval"),
        Line2D(
            [0],
            [0],
            color=ORANGE,
            lw=0.9,
            ls="--",
            label=f"Common policy coefficient: {common:g}",
        ),
        Line2D(
            [0],
            [0],
            color=GRAY,
            lw=0.9,
            ls=":",
            label=f"Required common coefficient: < {cutoff:g}",
        ),
    ]
    fig.legend(
        handles=handles,
        loc="upper center",
        bbox_to_anchor=(0.53, 0.995),
        ncol=2,
        frameon=False,
        fontsize=8.4,
        columnspacing=1.8,
        handlelength=2.2,
    )
    finish(fig, out, "Figure_1")

    b = read("breadth.csv")
    fig, ax = plt.subplots(figsize=(7.2, 3.55))
    fig.subplots_adjust(left=0.16, right=0.78, bottom=0.20, top=0.72)
    for i, r in enumerate(b):
        current = next(
            z
            for z in frontier
            if float(z["ambition"]) == float(r["alpha"])
            and int(z["maximum_climate_coefficients"]) == 1
            and float(z["response_band"]) == 0
        )
        y = 2 - i
        n = int(r["national_B_minimum"])
        lo = int(current["complement_lower_bound"])
        hi = int(current["complement_upper_bound"])
        if lo == hi:
            ax.plot(lo, y, "s", color=BLUE, ms=6)
            ax.text(
                lo,
                y + 0.16,
                str(lo),
                ha="center",
                va="bottom",
                fontsize=8.5,
                color=BLUE,
            )
        else:
            ax.plot([lo, hi], [y, y], lw=3, color=BLUE, solid_capstyle="butt")
            ax.plot(lo, y, "D", color=BLUE, ms=5.5)
            ax.plot(hi, y, "o", mfc="white", mec=BLUE, ms=6, mew=1.2)
            ax.text(
                (lo + hi) / 2,
                y + 0.16,
                f"{lo} to {hi}",
                ha="center",
                va="bottom",
                fontsize=8.5,
                color=BLUE,
            )
        ax.plot(n, y, "o", color=ORANGE, ms=6)
        ax.text(n, y - 0.17, str(n), ha="center", va="top", fontsize=8.5, color=ORANGE)
        extra = (
            str(lo - n) + " more" if lo == hi else "At least " + str(lo - n) + " more"
        )
        ax.text(
            24.0,
            y,
            extra,
            ha="left",
            va="center",
            fontsize=9,
            color="#222222",
            clip_on=False,
        )
    ax.set_yticks([2, 1, 0], ["0% (baseline)", "2.5%", "5%"])
    ax.set_ylim(-0.58, 2.65)
    ax.set_xlim(0, 23)
    ax.set_xticks([0, 4, 8, 12, 16, 20, 23])
    ax.grid(axis="x", color=LIGHT, lw=0.5)
    ax.set_axisbelow(True)
    ax.spines["left"].set_visible(False)
    ax.tick_params(axis="y", length=0, pad=8)
    ax.set_xlabel("Countries requiring a biodiversity complement", labelpad=8)
    fig.legend(
        handles=[
            Line2D(
                [0],
                [0],
                marker="o",
                linestyle="none",
                color=ORANGE,
                label="National exact minimum",
            ),
            Line2D(
                [0],
                [0],
                marker="s",
                linestyle="none",
                color=BLUE,
                label="Common exact minimum",
            ),
            Line2D(
                [0],
                [0],
                marker="D",
                linestyle="none",
                color=BLUE,
                label="Common lower bound",
            ),
            Line2D(
                [0],
                [0],
                marker="o",
                linestyle="none",
                mfc="white",
                mec=BLUE,
                mew=1.2,
                label="Common feasible upper bound",
            ),
        ],
        loc="upper center",
        bbox_to_anchor=(0.51, 0.995),
        ncol=2,
        frameon=False,
        fontsize=8.3,
        columnspacing=2.0,
        handletextpad=0.7,
    )
    finish(fig, out, "Figure_2")

    primary = sorted(
        [
            r
            for r in frontier
            if float(r["ambition"]) == 0.05 and float(r["response_band"]) == 1e-4
        ],
        key=lambda r: int(r["maximum_climate_coefficients"]),
    )
    if any(r["minimum_identified"] != "True" for r in primary):
        raise ValueError("The primary frontier contains an unresolved interval.")
    xx = np.array([int(r["maximum_climate_coefficients"]) for r in primary])
    kk = np.array([int(r["complement_lower_bound"]) for r in primary])
    counts = xx + kk
    fig, axes = plt.subplots(1, 2, figsize=(7.2, 3.6))
    fig.subplots_adjust(left=0.09, right=0.98, bottom=0.22, top=0.80, wspace=0.35)
    for ax, values, title, color, limits in [
        (axes[0], kk, "a  Biodiversity complements", BLUE, (2.5, 22.5)),
        (axes[1], counts, "b  Total nominal coefficients", ORANGE, (9.5, 23.5)),
    ]:
        ax.plot(xx, values, "o-", color=color, ms=5.5, lw=1.4)
        for x, y in zip(xx, values):
            ax.annotate(
                str(y),
                (x, y),
                xytext=(0, 7),
                textcoords="offset points",
                ha="center",
                va="bottom",
                fontsize=9,
                color=color,
            )
        ax.set_title(title, loc="left", fontweight="bold", fontsize=10, pad=12)
        ax.set_xticks(xx)
        ax.set_xlim(0.65, 7.35)
        ax.set_ylim(*limits)
        ax.grid(axis="y", color=LIGHT, lw=0.5)
        ax.set_axisbelow(True)
    axes[0].set_ylabel("Minimum number of countries", labelpad=8)
    axes[0].set_yticks([4, 8, 12, 16, 20])
    axes[1].set_ylabel("Climate plus biodiversity coefficients", labelpad=8)
    axes[1].set_yticks([11, 13, 15, 17, 19, 21])
    fig.supxlabel("Distinct climate coefficients", y=0.065, fontsize=9)
    finish(fig, out, "Figure_3")

    e = read("exclusion_sensitivity.csv")
    fig, ax = plt.subplots(figsize=(7.2, 3.65))
    fig.subplots_adjust(left=0.15, right=0.965, bottom=0.29, top=0.84)
    ax.axhspan(-0.012, 0, color="#F4F4F4", zorder=0)
    ax.axhline(0, color=GRAY, lw=0.8, zorder=1)
    for country, marker, color in [("Finland", "o", ORANGE), ("Sweden", "s", BLUE)]:
        rr = [r for r in e if r["country"] == country]
        xx = np.array([100 * float(r["biodiversity_ambition_fraction"]) for r in rr])
        yy = np.array([float(r["dual_certified_shortfall_lower"]) for r in rr])
        ax.plot(xx, yy, color=color, lw=1.3, label=country)
        for x, y in zip(xx, yy):
            ax.plot(
                x,
                y,
                marker=marker,
                mfc=color if y > 0 else "white",
                mec=color,
                ms=6,
                mew=1.1,
            )
    ax.set_xlim(-0.2, 5.25)
    ax.set_ylim(-0.012, 0.027)
    ax.set_xticks([0, 2.5, 5], ["0% (baseline)", "2.5%", "5%"])
    ax.set_yticks([-0.01, 0, 0.01, 0.02], ["−0.01", "0", "0.01", "0.02"])
    ax.set_ylabel("Normalized separating lower bound", labelpad=8)
    ax.set_xlabel("Biodiversity ambition as a share of national headroom", labelpad=8)
    ax.grid(axis="y", color=LIGHT, lw=0.5)
    ax.set_axisbelow(True)
    ax.legend(
        loc="lower left",
        bbox_to_anchor=(-0.015, 1.01),
        ncol=2,
        frameon=False,
        fontsize=8.8,
    )
    fig.text(
        0.55,
        0.035,
        "Filled points: all-price exclusion     Open points: inconclusive certificate",
        ha="center",
        va="bottom",
        fontsize=8.4,
        color=GRAY,
    )
    finish(fig, out, "Figure_4")

    d = read("dominance_diagnostics.csv")
    fig, ax = plt.subplots(figsize=(7.2, 2.7))
    fig.subplots_adjust(left=0.15, right=0.92, bottom=0.25, top=0.75)
    for off, base, label, color in [
        (-0.17, "level", "Biodiversity levels", BLUE),
        (0.17, "gain", "Positive biodiversity increments", ORANGE),
    ]:
        vals = [
            100
            * float(
                next(
                    r
                    for r in d
                    if r["country"] == c
                    and r["climate_base"] == "positive"
                    and r["biodiversity_base"] == base
                )["baseline_actions_removed_share"]
            )
            for c in ["Finland", "Sweden"]
        ]
        ys = np.arange(2) + off
        ax.barh(ys, vals, height=0.27, color=color, label=label, zorder=3)
        for y, v in zip(ys, vals):
            ax.text(v + 1.2, y, f"{v:.1f}%", va="center", fontsize=8.5)
    ax.set_yticks([0, 1], ["Finland", "Sweden"])
    ax.invert_yaxis()
    ax.set_xlim(0, 105)
    ax.set_xticks([0, 25, 50, 75, 100])
    ax.set_xlabel("Reference actions excluded by dominance (%)", labelpad=7)
    ax.spines["left"].set_visible(False)
    ax.tick_params(axis="y", length=0, pad=8)
    ax.grid(axis="x", color=LIGHT, lw=0.5)
    ax.set_axisbelow(True)
    ax.legend(
        loc="lower left",
        bbox_to_anchor=(-0.02, 1.06),
        ncol=2,
        frameon=False,
        fontsize=8.5,
    )
    finish(fig, out, "Figure_S1")
    (out / "figure_layout.json").write_text(json.dumps(LAYOUT, indent=2) + "\n")
    sources = {
        "Figure_1": [
            "scope.csv",
            "price_regions.csv",
            "price_groups/climate_ceiling_certificate.json",
            "price_groups/price_group_verification.json",
        ],
        "Figure_2": ["breadth.csv", "price_groups/price_group_frontier.csv"],
        "Figure_3": ["price_groups/price_group_frontier.csv"],
        "Figure_4": ["exclusion_sensitivity.csv"],
        "Figure_S1": ["dominance_diagnostics.csv"],
    }
    (out / "figure_sources.json").write_text(json.dumps(sources, indent=2) + "\n")
    print("Five figures exported from common canvases; PNG and TIFF pixels identical.")


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--data", type=Path, default=ROOT / "data")
    p.add_argument("--output", type=Path, default=ROOT / "outputs/figures")
    p.add_argument("--dpi", type=int, default=600,
                   help="Raster export resolution; use 1000 for line-art submission")
    args = p.parse_args()
    main(args.output, args.data, args.dpi)

"""Regenerate the reviewed CSV tables and figures from the public results.json.

This is a presentation script, not a solver or a new accuracy assessment.
"""

import argparse
import csv
import json
from pathlib import Path
import re
import xml.etree.ElementTree as ET

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import FixedLocator, FormatStrFormatter, LogLocator


HERE = Path(__file__).resolve().parent
OUTPUT = HERE / "runs/presentation"
PREVIEW = OUTPUT
IDS = [
    "cosh12-qp-q6", "cosh12-qp-q7", "cosh12-qp-q8",
    "cosh12-dmrg-chi384", "cosh12-dmrg-chi512",
    "cosh20-qp-q6", "cosh24-qp-q6", "cosh26-qp-q5", "cosh26-qp-q6",
    "surrogate24-window100-qp-q6", "surrogate24-window200-qp-q6",
    "surrogate24-window400-qp-q6",
]
FINALS = [
    "cosh26-qp-q6", "surrogate24-window400-qp-q6",
    "cosh12-dmrg-chi512", "cosh12-qp-q8",
]
OBSERVABLES = {
    "signed_gap": "g=(E_D-E_S)/Delta; negative means doublet ground state",
    "singlet.energy": "E_S/Delta, lowest singlet energy in the stated reference",
    "doublet.energy": "E_D/Delta, lowest doublet energy in the stated reference",
    "singlet.P0": "Singlet impurity empty-state probability",
    "singlet.P1": "Singlet impurity single-occupation probability (both spins)",
    "singlet.P2": "Singlet impurity double-occupation probability",
    "doublet.P0": "Doublet impurity empty-state probability",
    "doublet.P1": "Doublet impurity single-occupation probability (both spins)",
    "doublet.P2": "Doublet impurity double-occupation probability",
    "doublet.moment": "Impurity <(n_up-n_down)/2> in the total S_z=+1/2 doublet",
}
# Every tuple is (kind, new endpoint, old endpoint); subtraction is always new-old.
PAIRS = [
    ("qp_cutoff_refinement", "cosh26-qp-q6", "cosh26-qp-q5"),
    ("bath_level_refinement", "cosh24-qp-q6", "cosh20-qp-q6"),
    ("bath_level_refinement", "cosh26-qp-q6", "cosh24-qp-q6"),
    ("window_refinement", "surrogate24-window200-qp-q6", "surrogate24-window100-qp-q6"),
    ("window_refinement", "surrogate24-window400-qp-q6", "surrogate24-window200-qp-q6"),
    ("bond_refinement", "cosh12-dmrg-chi512", "cosh12-dmrg-chi384"),
    ("qp_cutoff_refinement", "cosh12-qp-q7", "cosh12-qp-q6"),
    ("qp_cutoff_refinement", "cosh12-qp-q8", "cosh12-qp-q7"),
    ("same_bath_backend", "cosh12-dmrg-chi512", "cosh12-qp-q8"),
    ("cross_family_and_levels", "surrogate24-window400-qp-q6", "cosh26-qp-q6"),
] + [("nrg_minus_finite", "nrg_reference", result_id) for result_id in FINALS]
OBS_COLUMNS = [
    "result_id", "role", "backend", "bath_family", "bath_levels", "qp_cutoff",
    "bond_dimension", "frequency_window", "observable", "value",
    "empirical_uncertainty", "uncertainty_kind", "requested_observable_target",
    "finite_problem_qualified", "continuum_target_established",
]
CMP_COLUMNS = [
    "comparison_id", "kind", "new_id", "old_id", "observable", "new_value",
    "old_value", "difference", "scale", "ratio", "denominator_type",
    "normalization_meaning",
]
MEANINGS = {
    "requested_observable_target": "abs(new-old)/(absolute+relative*abs(new)); finite difference, not an error bound",
    "nrg_empirical_uncertainty": "abs(NRG-finite)/NRG empirical uncertainty; conditional, nonstatistical",
}
COLORS = ["#0072B2", "#009E73", "#D55E00", "#CC79A7"]
MARKERS = ["o", "s", "D", "^"]
SHORT_LABELS = [
    r"QP / cosh 26 / $q=6$",
    "QP / surrogate 24\n" + r"$W=400$, $q=6$",
    r"DMRG / cosh 12 / $\chi=512$",
    r"QP / cosh 12 / $q=8$",
]


def tables(data):
    with (OUTPUT / "observables.csv").open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=OBS_COLUMNS)
        writer.writeheader()
        for result_id, record in list(data["finite_results"].items()) + [("nrg_reference", data["nrg_reference"])]:
            is_nrg = result_id == "nrg_reference"
            settings = {"backend": "nrg"} if is_nrg else record["settings"]
            target = data["targets"][settings["backend"]]
            for observable, value in record["values"].items():
                writer.writerow({
                    "result_id": result_id,
                    "role": "reference" if is_nrg else "final" if result_id in FINALS else "refinement",
                    **settings, "observable": observable, "value": value,
                    "empirical_uncertainty": record["empirical_uncertainties"][observable] if is_nrg else "",
                    "uncertainty_kind": record["uncertainty_kind"] if is_nrg else "",
                    "requested_observable_target": target["absolute"] + target["relative"] * abs(value),
                    "finite_problem_qualified": "" if is_nrg else str(record["finite_problem_qualified"]).lower(),
                    "continuum_target_established": "" if is_nrg else "false",
                })
    with (OUTPUT / "comparisons.csv").open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=CMP_COLUMNS)
        writer.writeheader()
        for pair in data["comparisons"]:
            new = data["nrg_reference"] if pair["new_id"] == "nrg_reference" else data["finite_results"][pair["new_id"]]
            old = data["finite_results"][pair["old_id"]]
            for observable, values in pair["observables"].items():
                writer.writerow({
                    "comparison_id": pair["id"], "kind": pair["kind"],
                    "new_id": pair["new_id"], "old_id": pair["old_id"],
                    "observable": observable, "new_value": new["values"][observable],
                    "old_value": old["values"][observable], **values,
                    "denominator_type": pair["denominator_type"],
                    "normalization_meaning": MEANINGS[pair["denominator_type"]],
                })


def save_figure(fig, name, description):
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    for text in fig.texts:
        bounds = text.get_window_extent(renderer)
        assert fig.bbox.contains(bounds.x0, bounds.y0)
        assert fig.bbox.contains(bounds.x1, bounds.y1)
        for ax in fig.axes:
            assert not bounds.overlaps(ax.xaxis.label.get_window_extent(renderer))
    fig.savefig(OUTPUT / f"{name}.svg", metadata={"Date": None, "Creator": None, "Description": description})
    fig.savefig(PREVIEW / f"{name}.png", dpi=180, metadata={"Software": "matplotlib"})
    plt.close(fig)


def figures(data):
    plt.rcParams.update({
        "font.family": "DejaVu Sans", "font.size": 11,
        "axes.titlesize": 13, "axes.titleweight": "bold", "axes.labelsize": 11,
        "text.color": "#202B33", "axes.labelcolor": "#202B33",
        "xtick.color": "#3C4852", "ytick.color": "#3C4852",
        "axes.edgecolor": "#AAB3BA", "axes.spines.top": False,
        "axes.spines.right": False, "axes.axisbelow": True,
        "grid.color": "#E4E8EB", "grid.linewidth": 0.7,
        "svg.fonttype": "none", "svg.hashsalt": "static-results",
        "savefig.facecolor": "white", "figure.facecolor": "white",
        "mathtext.fontset": "dejavusans",
    })
    finite = data["finite_results"]
    reference = data["nrg_reference"]
    nrg_pairs = data["comparisons"][-4:]
    # A narrow, vertically stacked figure remains legible when embedded on mobile.
    fig, (ax, bx) = plt.subplots(2, 1, figsize=(8.0, 10.4), gridspec_kw={"height_ratios": [1, 1.15]})
    fig.subplots_adjust(left=0.33, right=0.96, top=0.86, bottom=0.19, hspace=0.48)
    fig.text(0.055, 0.962, "Finite-bath results against NRG", fontsize=19, weight="bold")
    fig.text(0.055, 0.930, r"Single reservoir  |  $D/\Delta=100$, $U/\Delta=2$, $\Gamma/\Delta=0.4$", fontsize=11)
    fig.text(0.055, 0.904, r"Zero detuning and field  |  $g=(E_D-E_S)/\Delta$", fontsize=11)
    gap = reference["values"]["signed_gap"]
    uncertainty = reference["empirical_uncertainties"]["signed_gap"]
    ax.axvspan(gap - uncertainty, gap + uncertainty, color="#DDE6EB", zorder=0)
    ax.axvline(gap, color="#52616B", linewidth=1.2, linestyle="--")
    ax.axhspan(1.5, 3.5, color="#F5F5F5", zorder=-1)
    for index, result_id in enumerate(FINALS):
        ax.plot(finite[result_id]["values"]["signed_gap"], index, marker=MARKERS[index], color=COLORS[index], markersize=8, linestyle="none")
    ax.set_yticks(range(4), SHORT_LABELS)
    ax.set_ylim(3.5, -0.8)
    ax.set_xlim(-0.174451, -0.174357)
    ax.xaxis.set_major_locator(FixedLocator([-0.17444, -0.17440, -0.17436]))
    ax.xaxis.set_major_formatter(FormatStrFormatter("%.5f"))
    ax.tick_params(axis="y", length=0, pad=12)
    ax.set_xlabel(r"Signed gap $g$ (negative: doublet ground state)")
    ax.set_title("a  Signed gaps", loc="left", pad=14)
    ax.text(0.98, 0.98, "NRG empirical band", transform=ax.transAxes, ha="right", va="top", fontsize=10, color="#52616B")
    ax.text(0.98, 0.04, "Shaded rows: same cosh-12 bath", transform=ax.transAxes, ha="right", va="bottom", fontsize=10)
    ax.grid(axis="x")

    ticks = [r"$g$", r"$E_S/\Delta$", r"$E_D/\Delta$", r"$P_{0,S}$", r"$P_{1,S}$", r"$P_{2,S}$", r"$P_{0,D}$", r"$P_{1,D}$", r"$P_{2,D}$", r"$\langle S_{z,d}\rangle_D$"]
    bx.axvspan(1e-5, 1, color="#EEF3F6", zorder=0)
    bx.axvline(1, color="#52616B", linewidth=1, linestyle="--")
    for index, pair in enumerate(nrg_pairs):
        ratios = [pair["observables"][observable]["ratio"] for observable in OBSERVABLES]
        bx.plot(ratios, [j + (index - 1.5) * 0.15 for j in range(10)], linestyle="none", marker=MARKERS[index], color=COLORS[index], markersize=5.8, markerfacecolor="white" if index == 3 else COLORS[index], markeredgewidth=1.2)
    bx.set_xscale("log")
    bx.set_xlim(1e-4, 25)
    bx.set_ylim(9.6, -0.7)
    bx.set_yticks(range(10), ticks)
    bx.tick_params(axis="y", length=0, pad=12)
    bx.xaxis.set_major_locator(LogLocator(base=10, numticks=7))
    bx.set_xlabel(r"$|X_{\mathrm{NRG}}-X_{\mathrm{finite}}|\,/\,u_{\mathrm{NRG}}(X)$")
    bx.set_title("b  All ten observables", loc="left", pad=14)
    bx.grid(axis="x", which="major")
    fig.text(0.055, 0.106, "Markers match panel a. Larger-bath QP: all ten ratios < 0.0731.", fontsize=10)
    fig.text(0.055, 0.081, "Cosh-12: gap 2.25; singlet probabilities 12.89 NRG empirical uncertainties.", fontsize=10)
    fig.text(0.055, 0.056, r"$u_{\mathrm{NRG}}$: conditional, nonstatistical uncertainty. Finite points have no error bars.", fontsize=10)
    fig.text(0.055, 0.031, "QP/DMRG continuum targets are not established; DMRG fails the all-root residual criterion.", fontsize=10)
    save_figure(fig, "comparison", "Final signed gaps and all ten NRG-minus-finite absolute differences divided by conditional nonstatistical NRG empirical uncertainties. The cosh-12 QP and DMRG results share a bath; the larger QP baths differ. No QP/DMRG continuum accuracy is established.")

    fig, (ax, bx) = plt.subplots(2, 1, figsize=(8.0, 10.4), gridspec_kw={"height_ratios": [1.3, 1]})
    fig.subplots_adjust(left=0.34, right=0.95, top=0.86, bottom=0.18, hspace=0.55)
    fig.text(0.055, 0.962, "Refinement and residual diagnostics", fontsize=18, weight="bold")
    fig.text(0.055, 0.931, r"Requested QP/DMRG target: $\tau(X)=10^{-9}+10^{-6}|X_{\mathrm{new}}|$", fontsize=11)
    fig.text(0.055, 0.904, "Finite differences are not error bounds; small changes do not certify convergence.", fontsize=10)
    refinement_labels = [
        "Cosh 26: cutoff 5 to 6",
        "Cosh: 20 to 24 levels (q=6)",
        "Cosh: 24 to 26 levels (q=6)",
        "Surrogate 24: W=100 to 200",
        "Surrogate 24: W=200 to 400",
        "Cosh 12: bond 384 to 512",
        "Cosh 12: cutoff 6 to 7",
        "Cosh 12: cutoff 7 to 8",
    ]
    maxima = [max(row["ratio"] for row in pair["observables"].values()) for pair in data["comparisons"][:8]]
    ladder_colors = [COLORS[0], COLORS[0], COLORS[0], COLORS[1], COLORS[1], COLORS[2], COLORS[3], COLORS[3]]
    ax.axvspan(1e-9, 1, color="#EEF3F6")
    ax.axvline(1, color="#52616B", linewidth=1.2, linestyle="--")
    for index, (maximum, color) in enumerate(zip(maxima, ladder_colors, strict=True)):
        ax.plot(maximum, index, "o", color=color, markersize=7)
        ax.annotate(f"{maximum:.3g}", (maximum, index), xytext=(7, 0), textcoords="offset points", va="center", fontsize=10)
    ax.axhline(5.5, color="#CCD4DA", linewidth=1)
    ax.set_xscale("log")
    ax.set_xlim(1e-9, 45)
    ax.set_ylim(7.6, -0.7)
    ax.set_yticks(range(8), refinement_labels)
    ax.tick_params(axis="y", length=0, pad=12, labelsize=10)
    ax.set_xticks([1e-8, 1e-6, 1e-4, 1e-2, 1])
    ax.set_xlabel(r"$\max_X\ |X_{\mathrm{new}}-X_{\mathrm{old}}|\,/\,\tau(X_{\mathrm{new}})$")
    ax.set_title("a  Observable refinement changes", loc="left", pad=14)
    ax.grid(axis="x")
    ax.text(0.02, 0.015, "Surrogate: q=6 throughout; W in units of Delta", transform=ax.transAxes, fontsize=9, va="bottom")

    residual_labels = [r"Singlet root 0", r"Singlet root 1", r"Doublet root 0", r"Doublet root 1"]
    residual_colors = [COLORS[0], COLORS[1], COLORS[2], COLORS[3]]
    bx.axhspan(1e-8, 1e-7, color="#EEF3F6", zorder=0)
    bx.axhline(1e-7, color="#52616B", linewidth=1.2, linestyle="--")
    for index, (sector, root) in enumerate([("singlet", 0), ("singlet", 1), ("doublet", 0), ("doublet", 1)]):
        residuals = [finite[f"cosh12-dmrg-chi{chi}"]["residuals"][sector][root] for chi in [384, 512]]
        bx.plot([384, 512], residuals, marker=MARKERS[index], color=residual_colors[index], markersize=6, linewidth=1.5, label=residual_labels[index])
    bx.set_yscale("log")
    bx.set_ylim(4e-8, 2e-5)
    bx.set_xlim(365, 532)
    bx.set_xticks([384, 512])
    bx.set_xlabel(r"DMRG bond dimension $\chi$ (same cosh-12 bath)")
    bx.set_ylabel(r"$\|(H-E)\psi\|/\Delta$", labelpad=9)
    bx.set_title("b  Untruncated finite-H residuals", loc="left", pad=14)
    bx.grid(axis="y", which="major")
    bx.text(0.03, 0.19, r"All-root criterion: $10^{-7}$", transform=bx.transAxes, fontsize=10)
    handles, labels = bx.get_legend_handles_labels()
    fig.legend(handles, labels, loc="center left", bbox_to_anchor=(0.05, 0.32), frameon=False, fontsize=10, handlelength=1.6, labelspacing=0.7)
    fig.text(0.055, 0.077, "At bond 512, three of four residuals exceed the criterion; only doublet root 0 passes.", fontsize=10)
    fig.text(0.055, 0.049, "QP residual checks apply to the projected cutoff Hamiltonian, not omitted QP sectors.", fontsize=10)
    fig.text(0.055, 0.023, "Neither QP nor DMRG has met the requested continuum accuracy target.", fontsize=10)
    save_figure(fig, "convergence", "Eight finite-refinement maxima normalized by the requested observable target, and all four untruncated finite-H DMRG root residuals at bond dimensions 384 and 512. Three of four roots fail at 512. Finite differences are not error bounds.")


def verify(data):
    assert len(data["observables"]) == 10 and "singlet.moment" not in data["observables"]
    assert list(data["finite_results"]) == IDS
    assert data["final_result_ids"] == FINALS
    assert len(data["comparisons"]) == 14
    for result_id, row in data["finite_results"].items():
        assert set(row["values"]) == set(OBSERVABLES)
        values = row["values"]
        assert abs(values["doublet.energy"] - values["singlet.energy"] - values["signed_gap"]) < 2e-15
        for sector in ["singlet", "doublet"]:
            assert abs(sum(values[f"{sector}.P{i}"] for i in range(3)) - 1) < 1e-11
        qualified = all(r <= row["residual_tolerance"] for rs in row["residuals"].values() for r in rs)
        assert qualified == row["finite_problem_qualified"] == (row["settings"]["backend"] == "qp")
        assert row["continuum_target_established"] is False
        if row["settings"]["bath_family"] == "surrogate":
            assert row["settings"]["frequency_window"] == float(result_id.split("window")[1].split("-")[0])
    dmrg = data["finite_results"]["cosh12-dmrg-chi512"]
    assert sum(r > dmrg["residual_tolerance"] for rs in dmrg["residuals"].values() for r in rs) == 3
    with (OUTPUT / "observables.csv").open(newline="") as stream:
        reader = csv.DictReader(stream)
        assert reader.fieldnames == OBS_COLUMNS
        rows = list(reader)
    assert len(rows) == 130
    assert len({(r["result_id"], r["observable"]) for r in rows}) == 130
    for row in rows:
        is_nrg = row["result_id"] == "nrg_reference"
        source = data["nrg_reference"] if is_nrg else data["finite_results"][row["result_id"]]
        assert float(row["value"]) == source["values"][row["observable"]]
        target = data["targets"][row["backend"]]
        assert float(row["requested_observable_target"]) == target["absolute"] + target["relative"] * abs(float(row["value"]))
        if is_nrg:
            assert float(row["empirical_uncertainty"]) == source["empirical_uncertainties"][row["observable"]]
            assert float(row["empirical_uncertainty"]) <= float(row["requested_observable_target"])
        else:
            assert row["empirical_uncertainty"] == "" and row["uncertainty_kind"] == ""
            assert row["role"] == ("final" if row["result_id"] in FINALS else "refinement")
            assert row["finite_problem_qualified"] == str(source["finite_problem_qualified"]).lower()
    with (OUTPUT / "comparisons.csv").open(newline="") as stream:
        reader = csv.DictReader(stream)
        assert reader.fieldnames == CMP_COLUMNS
        rows = list(reader)
    assert len(rows) == 140
    assert len({(r["comparison_id"], r["observable"]) for r in rows}) == 140
    pairs = {pair["id"]: pair for pair in data["comparisons"]}
    for row in rows:
        pair = pairs[row["comparison_id"]]
        observable = row["observable"]
        new_value, old_value = float(row["new_value"]), float(row["old_value"])
        assert float(row["difference"]) == new_value - old_value
        if row["denominator_type"] == "requested_observable_target":
            target = data["targets"][data["finite_results"][row["new_id"]]["settings"]["backend"]]
            assert float(row["scale"]) == target["absolute"] + target["relative"] * abs(new_value)
        else:
            assert row["new_id"] == "nrg_reference"
            assert float(row["scale"]) == data["nrg_reference"]["empirical_uncertainties"][observable]
        assert float(row["ratio"]) == abs(new_value - old_value) / float(row["scale"])
        for key in ["difference", "scale", "ratio"]:
            assert float(row[key]) == pair["observables"][observable][key]
        assert row["normalization_meaning"] == MEANINGS[row["denominator_type"]]
    for name in ["comparison", "convergence"]:
        svg = ET.parse(OUTPUT / f"{name}.svg").getroot()
        assert svg.findall(".//{http://www.w3.org/2000/svg}text")
        assert not svg.findall(".//{http://purl.org/dc/elements/1.1/}date")
        image = plt.imread(PREVIEW / f"{name}.png")
        assert image.shape[0] > 1000 and image.shape[1] > 1000
        assert image.std() > 0.05
    for name in ["results.json", "observables.csv", "comparisons.csv", "comparison.svg", "convergence.svg"]:
        text = (OUTPUT / name).read_text()
        assert not re.search(r"campaign|sha256|producer|task_id|budget|attempt-|closeout|retention2|runs/", text, re.I), name
        assert not re.search(r"\b[0-9a-f]{64}\b", text), name
        print(f"{name}: {(OUTPUT / name).stat().st_size:,} bytes")
    assert (OUTPUT / "results.json").stat().st_size < 60000
    for result_id, pair in zip(FINALS, data["comparisons"][-4:], strict=True):
        worst = max(pair["observables"], key=lambda o: pair["observables"][o]["ratio"])
        print(f"{result_id}: gap={data['finite_results'][result_id]['values']['signed_gap']:.16g}; max NRG ratio={pair['observables'][worst]['ratio']:.9g} ({worst})")
    for pair in data["comparisons"][:10]:
        print(f"{pair['id']}: max target ratio={max(row['ratio'] for row in pair['observables'].values()):.9g}")
    print("Verified: 12 finite endpoints, 10 NRG observables, 14 pairs, 130/140 CSV rows, residual flags, SVG XML/text, PNG renders.")


def main(argv=None):
    global OUTPUT, PREVIEW
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, default=HERE / "runs/presentation")
    args = parser.parse_args(argv)
    OUTPUT = PREVIEW = args.output_dir.resolve()
    if OUTPUT == (HERE / "output").resolve() or (HERE / "output").resolve() in OUTPUT.parents:
        parser.error("use a separate directory; published outputs are read-only")
    OUTPUT.mkdir(parents=True, exist_ok=True)
    source = (HERE / "output/results.json").read_bytes()
    data = json.loads(source)
    (OUTPUT / "results.json").write_bytes(source)
    tables(data)
    figures(data)
    verify(data)


if __name__ == "__main__":
    main()

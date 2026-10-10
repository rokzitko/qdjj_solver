"""Regenerate presentation tables and figures from the shipped results.json.

No solver, archive, HDF5, or installed qdjj_solver package is needed. Matplotlib
is needed only for figures; --tables-only performs standard-library validation.
Requested targets normalize differences, not uncertainties or error bounds.
"""

import argparse
import csv
import json
import math
from pathlib import Path
import re
import xml.etree.ElementTree as ET


HERE = Path(__file__).resolve().parent
FINALS = ["cosh20-q7", "surrogate24-w400-q6", "cosh16-dmrg-final", "cosh16-q8"]
NRG_ID = "NRG-L2.0-K8000-Nz8"
OBSERVABLES = {
    "signed_gap": "g=(E_D-E_S)/Delta; negative means doublet ground state",
    "singlet.energy": "E_S/Delta, lowest singlet candidate in the stated reference",
    "doublet.energy": "E_D/Delta, lowest doublet candidate in the stated reference",
    "singlet.P0": "Singlet impurity empty-state probability",
    "singlet.P1": "Singlet impurity single-occupation probability (both spins)",
    "singlet.P2": "Singlet impurity double-occupation probability",
    "doublet.P0": "Doublet impurity empty-state probability",
    "doublet.P1": "Doublet impurity single-occupation probability (both spins)",
    "doublet.P2": "Doublet impurity double-occupation probability",
    "doublet.moment": "Impurity <(n_up-n_down)/2> in the total S_z=+1/2 doublet",
}
# Tuples are (kind, new endpoint, old endpoint, plot label); subtraction is new-old.
PAIRS = [
    ("qp_cutoff_refinement", "cosh16-q7", "cosh16-q6", "Cosh 16: cutoff 6 to 7"),
    ("qp_cutoff_refinement", "cosh16-q8", "cosh16-q7", "Cosh 16: cutoff 7 to 8"),
    ("bath_level_refinement", "cosh16-q8", "cosh12-q8", "Cosh: 12 to 16 levels (q=8)"),
    ("bath_level_refinement", "cosh20-q7", "cosh16-q7", "Cosh: 16 to 20 levels (q=7)"),
    ("window_refinement", "surrogate24-w400-q6", "surrogate24-w200-q6", "Surrogate 24: W=200 to 400"),
    ("nrg_twist_sensitivity", NRG_ID, "NRG-L2.0-K8000-Nz4", "Lambda 2: Nz=4 to 8"),
    ("nrg_keep_sensitivity", NRG_ID, "NRG-L2.0-K4000-Nz8", "Lambda 2: keep 4000 to 8000"),
    ("nrg_keep_sensitivity", "NRG-L2.5-K8000-Nz4", "NRG-L2.5-K4000-Nz4", "Lambda 2.5: keep 4000 to 8000"),
    ("nrg_lambda_sensitivity", "NRG-L2.0-K8000-Nz4", "NRG-L2.5-K8000-Nz4", "Lambda 2.5 to 2 (Nz=4)"),
    ("nrg_lambda_sensitivity", "NRG-L1.8-K8000-Nz4", "NRG-L2.0-K8000-Nz4", "Lambda 2 to 1.8 (Nz=4)"),
    ("bond_refinement", "cosh16-dmrg-chi512", "cosh16-dmrg-chi256", "DMRG: bonds 256 to 512"),
    ("bond_refinement", "cosh16-dmrg-chi768", "cosh16-dmrg-chi512", "DMRG: bonds 512 to 768"),
    ("bond_refinement", "cosh16-dmrg-final", "cosh16-dmrg-chi768", "DMRG: singlet 768 to 1024"),
    ("same_bath_backend", "cosh16-dmrg-final", "cosh16-q8", "Matched cosh-16 DMRG minus QP"),
] + [("nrg_minus_finite", NRG_ID, name, name) for name in FINALS]
OBS_COLUMNS = [
    "result_id", "role", "backend", "bath_family", "bath_levels", "qp_cutoff",
    "singlet_bond_dimension", "doublet_bond_dimension", "frequency_window",
    "lambda", "keep", "nz", "keepenergy", "keepmin", "nmax", "observable", "value",
    "empirical_uncertainty", "uncertainty_kind", "requested_observable_target",
    "finite_problem_qualified", "finite_flow_completed", "continuum_target_established",
]
CMP_COLUMNS = [
    "comparison_id", "kind", "new_id", "old_id", "observable", "new_value", "old_value",
    "difference", "scale", "ratio", "denominator_type", "normalization_meaning",
]
MEANING = "abs(new-old)/(absolute+relative*abs(new)); new backend requested target, not uncertainty or error bound"
COLORS = ["#0072B2", "#009E73", "#D55E00", "#CC79A7"]
MARKERS = ["o", "s", "D", "^"]
SHORT_LABELS = [
    r"QP / cosh 20 / $q=7$",
    "QP / surrogate 24\n" + r"$W=400$, $q=6$",
    "DMRG / cosh 16\n" + r"$\chi_S=1024$, $\chi_D=768$",
    r"QP / cosh 16 / $q=8$",
]
ARTIFACTS = ["results.json", "observables.csv", "comparisons.csv", "comparison.svg", "convergence.svg"]


def verify_data(data):
    """Check the scientific contract using only the standard library."""
    assert data["format_version"] == 1
    assert data["observables"] == OBSERVABLES
    assert data["final_result_ids"] == FINALS and data["nrg_comparison_id"] == NRG_ID
    assert data["model"] == dict(geometry="single", gap=1.0, bandwidth=100.0,
                                 u=2.0, gamma=0.4, detuning=0.5, field=0.0)
    assert data["conventions"]["compress"] is False
    assert data["targets"] == {b: dict(absolute=a, relative=r) for b, a, r in
                               [("qp", 1e-9, 1e-6), ("dmrg", 1e-9, 1e-6), ("nrg", 1e-6, 1e-3)]}
    finite, nrg = data["finite_results"], data["nrg_results"]
    assert len(finite) == 14 and len(nrg) == 7
    qp_ids = {f"cosh{n}-q{q}" for n, qs in [(12, [6, 7, 8]), (16, [6, 7, 8]), (20, [6, 7])] for q in qs}
    qp_ids |= {"surrogate24-w200-q6", "surrogate24-w400-q6"}
    dmrg_ids = {"cosh16-dmrg-chi256", "cosh16-dmrg-chi512", "cosh16-dmrg-chi768", "cosh16-dmrg-final"}
    assert set(finite) == qp_ids | dmrg_ids
    records = finite | nrg
    assert set(data["symmetry_diagnostics"]) == set(records)
    for name, row in records.items():
        values = row["values"]
        assert set(values) == set(OBSERVABLES)
        assert all(type(v) in (float, int) and math.isfinite(v) for v in values.values())
        assert math.isclose(values["doublet.energy"] - values["singlet.energy"], values["signed_gap"], abs_tol=2e-15, rel_tol=0)
        assert values["signed_gap"] < 0
        for branch in ("singlet", "doublet"):
            ps = [values[f"{branch}.P{i}"] for i in range(3)]
            assert all(0 <= p <= 1 for p in ps) and abs(sum(ps) - 1) < 1e-11
            assert ps[0] > ps[2]  # Positive detuning favors the empty impurity.
        assert 0 < values["doublet.moment"] <= values["doublet.P1"] / 2
        assert abs(data["symmetry_diagnostics"][name]["singlet.moment"]) < 1e-10
        assert row["continuum_target_established"] is False
        settings = row["settings"]
        if name in nrg:
            assert settings["backend"] == "nrg" and row["finite_flow_completed"] is True
            assert row["uncertainty_kind"] == "not_established" and row["empirical_uncertainties"] is None
            assert name == f"NRG-L{settings['lambda']}-K{settings['keep']}-Nz{settings['nz']}"
            assert set(settings) == {"backend", "lambda", "keep", "nz", "keepenergy", "keepmin", "nmax"}
            continue
        qp = name in qp_ids
        assert settings["backend"] == ("qp" if qp else "dmrg")
        assert row["residual_tolerance"] == (1e-8 if qp else 1e-7)
        assert row["residual_scope"] == ("projected_cutoff_hamiltonian" if qp else "untruncated_finite_hamiltonian")
        assert set(row["residuals"]) == {"singlet", "doublet"}
        for branch, residuals in row["residuals"].items():
            assert residuals and all(math.isfinite(r) and r >= 0 for r in residuals)
            assert len(residuals) == row["checked_root_counts"][branch]
            assert row["root_separation_checked"][branch] == (len(residuals) > 1)
            assert row["seed_trial_counts"][branch] in (1, 3)
            spin = data["symmetry_diagnostics"][name]["total_spin_squared"][branch]
            assert abs(spin - (0 if branch == "singlet" else 0.75)) < 1e-6
            assert row["branch_finite_problem_qualified"][branch] == all(r <= row["residual_tolerance"] for r in residuals)
        assert row["finite_problem_qualified"] == all(row["branch_finite_problem_qualified"].values())
        assert row["finite_problem_qualified"] == (qp or name == "cosh16-dmrg-final")
        if qp:
            assert settings["qp_cutoff"] == int(name.rsplit("q", 1)[1])
            assert settings["bath_levels"] == int(re.search(r"\d+", name)[0])
        else:
            assert settings["qp_cutoff"] is None and settings["bath_levels"] == 16
            chi = int(name.rsplit("chi", 1)[1]) if "chi" in name else None
            assert settings["bond_dimensions"] == (dict(singlet=chi, doublet=chi) if chi else dict(singlet=1024, doublet=768))
        assert settings["frequency_window"] == (float(name.split("-w")[1].split("-")[0]) if "-w" in name else None)
    assert len(data["comparisons"]) == len(PAIRS) == 18
    for pair, (kind, new_id, old_id, _) in zip(data["comparisons"], PAIRS, strict=True):
        assert (pair["kind"], pair["new_id"], pair["old_id"]) == (kind, new_id, old_id)
        assert pair["id"] == f"{new_id}-minus-{old_id}"
        assert pair["denominator_type"] == "requested_observable_target"
        assert set(pair["observables"]) == set(OBSERVABLES)
        new, old = records[new_id], records[old_id]
        target = data["targets"][new["settings"]["backend"]]
        for key, metric in pair["observables"].items():
            diff = new["values"][key] - old["values"][key]
            scale = target["absolute"] + target["relative"] * abs(new["values"][key])
            assert metric == dict(difference=diff, scale=scale, ratio=abs(diff) / scale)


def table_rows(data):
    """Return both tables without rounding the stored floating-point values."""
    records = data["finite_results"] | data["nrg_results"]
    observations, comparisons = [], []
    for name, row in records.items():
        settings = dict(row["settings"])
        bonds = settings.pop("bond_dimensions", {})
        target = data["targets"][settings["backend"]]
        for observable in OBSERVABLES:
            value = row["values"][observable]
            observations.append({
                "result_id": name,
                "role": "comparison_point" if name == NRG_ID else "final" if name in FINALS else "refinement",
                **settings, "singlet_bond_dimension": bonds.get("singlet", ""),
                "doublet_bond_dimension": bonds.get("doublet", ""), "observable": observable, "value": value,
                "empirical_uncertainty": "", "uncertainty_kind": row.get("uncertainty_kind", ""),
                "requested_observable_target": target["absolute"] + target["relative"] * abs(value),
                "finite_problem_qualified": str(row["finite_problem_qualified"]).lower() if name in data["finite_results"] else "",
                "finite_flow_completed": "true" if name in data["nrg_results"] else "",
                "continuum_target_established": "false",
            })
    for pair in data["comparisons"]:
        for observable, metric in pair["observables"].items():
            comparisons.append({
                "comparison_id": pair["id"], "kind": pair["kind"], "new_id": pair["new_id"], "old_id": pair["old_id"],
                "observable": observable, "new_value": records[pair["new_id"]]["values"][observable],
                "old_value": records[pair["old_id"]]["values"][observable], **metric,
                "denominator_type": pair["denominator_type"], "normalization_meaning": MEANING,
            })
    return observations, comparisons


def tables(data, output_dir):
    for name, columns, rows in zip(["observables", "comparisons"], [OBS_COLUMNS, CMP_COLUMNS], table_rows(data), strict=True):
        with (output_dir / f"{name}.csv").open("w", newline="") as stream:
            writer = csv.DictWriter(stream, fieldnames=columns)
            writer.writeheader()
            writer.writerows(rows)


def figure_data(data):
    """Numerical payloads embedded in SVG descriptions for offline verification."""
    finite = data["finite_results"]
    comparison = {
        "caption": "Signed gaps and all ten absolute NRG-minus-finite differences divided by the requested NRG target, not uncertainty. The vertical line is one finite-Lambda twist mean. No uncertainty bands or continuum certification.",
        "nrg_gap": data["nrg_results"][NRG_ID]["values"]["signed_gap"],
        "gaps": [finite[name]["values"]["signed_gap"] for name in FINALS],
        "target_ratios": [[p["observables"][o]["ratio"] for o in OBSERVABLES] for p in data["comparisons"][-4:]],
    }
    convergence = {
        "caption": "Maximum observable changes normalized by the new backend requested target, not error bounds. NRG is finite-Lambda; zero keep sensitivity is an exact zero in the stored means. DMRG residuals show one checked root per branch; the doublet stops at chi768, the singlet continues to chi1024. One warm trial is not fresh independent-seed evidence.",
        "finite_indices": [0, 1, 2, 3, 4, 10, 11, 12], "nrg_indices": [5, 6, 7, 8, 9],
        "maxima": [max(m["ratio"] for m in p["observables"].values()) for p in data["comparisons"][:13]],
        "residual_tolerance": 1e-7,
        "residuals": {},
    }
    for branch in ("singlet", "doublet"):
        names = [f"cosh16-dmrg-chi{chi}" for chi in (256, 512, 768)]
        if branch == "singlet":
            names.append("cosh16-dmrg-final")
        assert all(len(finite[n]["residuals"][branch]) == 1 for n in names)
        convergence["residuals"][branch] = {
            "bonds": [finite[n]["settings"]["bond_dimensions"][branch] for n in names],
            "values": [finite[n]["residuals"][branch][0] for n in names],
        }
    return {"comparison": comparison, "convergence": convergence}


def figures(data, output_dir, png=False):
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.ticker import FormatStrFormatter, LogLocator

    plt.rcParams.update({
        "font.family": "DejaVu Sans", "font.size": 11, "axes.titlesize": 13,
        "axes.titleweight": "bold", "axes.labelsize": 11, "text.color": "#202B33",
        "axes.labelcolor": "#202B33", "xtick.color": "#3C4852", "ytick.color": "#3C4852",
        "axes.edgecolor": "#AAB3BA", "axes.spines.top": False, "axes.spines.right": False,
        "axes.axisbelow": True, "grid.color": "#E4E8EB", "grid.linewidth": 0.7,
        "svg.fonttype": "none", "svg.hashsalt": "test2-static-results",
        "savefig.facecolor": "white", "figure.facecolor": "white", "mathtext.fontset": "dejavusans",
    })
    payloads = figure_data(data)

    def save(fig, name):
        fig.canvas.draw()
        renderer = fig.canvas.get_renderer()
        for text in fig.texts:
            bounds = text.get_window_extent(renderer)
            assert fig.bbox.contains(bounds.x0, bounds.y0) and fig.bbox.contains(bounds.x1, bounds.y1)
            for ax in fig.axes:
                assert not bounds.overlaps(ax.xaxis.label.get_window_extent(renderer))
        fig.savefig(output_dir / f"{name}.svg", metadata={
            "Date": None, "Creator": None, "Title": f"Test2: {name}",
            "Description": json.dumps(payloads[name], sort_keys=True, allow_nan=False),
        })
        if png:
            fig.savefig(output_dir / f"{name}.png", dpi=180, metadata={"Software": "matplotlib"})
        plt.close(fig)

    p = payloads["comparison"]
    fig, (ax, bx) = plt.subplots(2, 1, figsize=(8.0, 10.4), gridspec_kw={"height_ratios": [1, 1.15]})
    fig.subplots_adjust(left=0.33, right=0.96, top=0.855, bottom=0.19, hspace=0.48)
    fig.text(0.055, 0.962, "Detuned impurity: finite-bath comparison", fontsize=17, weight="bold")
    fig.text(0.055, 0.929, r"Single reservoir  |  $D/\Delta=100$, $U/\Delta=2$, $\Gamma/\Delta=0.4$", fontsize=11)
    fig.text(0.055, 0.903, r"Detuning $0.5\Delta$; zero field  |  $g=(E_D-E_S)/\Delta$", fontsize=11)
    ax.axvline(p["nrg_gap"], color="#52616B", linewidth=1.2, linestyle="--")
    ax.axhspan(1.5, 3.5, color="#F5F5F5", zorder=-1)
    for i, gap in enumerate(p["gaps"]):
        ax.plot(gap, i, marker=MARKERS[i], color=COLORS[i], markersize=8, linestyle="none")
    ax.set_yticks(range(4), SHORT_LABELS)
    ax.set_ylim(3.5, -0.8)
    ax.set_xlim(-0.06300105, -0.0629985)
    ax.set_xticks([-0.063001, -0.063000, -0.062999])
    ax.xaxis.set_major_formatter(FormatStrFormatter("%.6f"))
    ax.tick_params(axis="y", length=0, pad=12)
    ax.set_xlabel(r"Signed gap $g$ (negative: doublet ground state)")
    ax.set_title("a  Signed gaps", loc="left", pad=14)
    ax.text(0.97, 0.98, "NRG finite-Lambda mean", transform=ax.transAxes, ha="right", va="top", fontsize=10,
            bbox={"facecolor": "white", "edgecolor": "none", "pad": 1})
    ax.text(0.97, 0.04, "Shaded rows: same cosh-16 bath", transform=ax.transAxes, ha="right", va="bottom", fontsize=10)
    ax.grid(axis="x")
    ticks = [r"$g$", r"$E_S/\Delta$", r"$E_D/\Delta$", r"$P_{0,S}$", r"$P_{1,S}$", r"$P_{2,S}$", r"$P_{0,D}$", r"$P_{1,D}$", r"$P_{2,D}$", r"$\langle S_{z,d}\rangle_D$"]
    bx.axvline(1, color="#52616B", linewidth=1, linestyle="--")
    for i, ratios in enumerate(p["target_ratios"]):
        assert all(1e-5 < v < 2 for v in ratios)
        bx.plot(ratios, [j + (i - 1.5) * 0.15 for j in range(10)], linestyle="none",
                marker=MARKERS[i], color=COLORS[i], markersize=5.8,
                markerfacecolor="white" if i == 3 else COLORS[i], markeredgewidth=1.2)
    bx.set_xscale("log")
    bx.set_xlim(1e-5, 2)
    bx.set_ylim(9.6, -0.7)
    bx.set_yticks(range(10), ticks)
    bx.tick_params(axis="y", length=0, pad=12)
    bx.xaxis.set_major_locator(LogLocator(base=10, numticks=7))
    bx.set_xlabel(r"$|X_{\mathrm{NRG}}-X_{\mathrm{finite}}|\,/\,(10^{-6}+10^{-3}|X_{\mathrm{NRG}}|)$")
    bx.set_title("b  All ten observables / requested NRG target", loc="left", pad=14)
    bx.grid(axis="x")
    fig.text(0.055, 0.107, r"NRG comparison point: $\Lambda=2$, keep=8000, $N_z=8$; not a continuum reference.", fontsize=10)
    fig.text(0.055, 0.082, "Markers match panel a. Ratios are target-normalized deviations, not uncertainties.", fontsize=10)
    fig.text(0.055, 0.057, "No uncertainty bands or error bars are established for these endpoints.", fontsize=10)
    fig.text(0.055, 0.032, "Final DMRG branch residuals pass; no backend has a certified continuum target.", fontsize=10)
    save(fig, "comparison")

    p = payloads["convergence"]
    fig, axes = plt.subplots(3, 1, figsize=(8.0, 14.4), gridspec_kw={"height_ratios": [1.3, 0.9, 1.05]})
    fig.subplots_adjust(left=0.36, right=0.95, top=0.90, bottom=0.145, hspace=0.49)
    fig.text(0.055, 0.974, "Refinement and residual diagnostics", fontsize=18, weight="bold")
    fig.text(0.055, 0.950, r"QP/DMRG: $\tau(X)=10^{-9}+10^{-6}|X|$; NRG: $10^{-6}+10^{-3}|X|$", fontsize=11)
    fig.text(0.055, 0.929, "Differences are not error bounds; denominators use the new endpoint values.", fontsize=10)
    for ax, group, title in zip(axes[:2], ["finite_indices", "nrg_indices"],
                               ["a  QP and DMRG observable changes", "b  Finite-Lambda NRG sensitivity"], strict=True):
        indices = p[group]
        ax.axvline(1, color="#52616B", linewidth=1.2, linestyle="--")
        for j, index in enumerate(indices):
            value = p["maxima"][index]
            color = COLORS[2] if index >= 10 else COLORS[1] if index == 4 else COLORS[0]
            ax.plot(value, j, "o", color=color, markersize=7, clip_on=value != 0)
            ax.annotate(f"{value:.3g}", (value, j), xytext=(7, 0), textcoords="offset points", va="center", fontsize=10,
                        bbox={"facecolor": "white", "edgecolor": "none", "pad": 0.4})
        ax.set_yticks(range(len(indices)), [PAIRS[i][3] for i in indices])
        ax.set_ylim(len(indices) - 0.4, -0.7)
        ax.tick_params(axis="y", length=0, pad=10, labelsize=9.5)
        ax.set_title(title, loc="left", pad=13)
        ax.set_xlabel(r"$\max_X\ |X_{\mathrm{new}}-X_{\mathrm{old}}|\,/\,\tau(X_{\mathrm{new}})$")
        ax.grid(axis="x")
    axes[0].set_xscale("log")
    axes[0].set_xlim(1e-9, 3000)
    axes[0].set_xticks([1e-8, 1e-5, 1e-2, 1, 1e3])
    # A linear interval at zero preserves the exactly zero keep difference.
    axes[1].set_xscale("symlog", linthresh=1e-8, linscale=0.7)
    axes[1].set_xlim(0, 2)
    axes[1].set_xticks([0, 1e-6, 1e-4, 1e-2, 1], ["0", r"$10^{-6}$", r"$10^{-4}$", r"$10^{-2}$", "1"])
    ax = axes[2]
    ax.axhline(p["residual_tolerance"], color="#52616B", linewidth=1.2, linestyle="--")
    for i, (branch, points) in enumerate(p["residuals"].items()):
        ax.plot(points["bonds"], points["values"], marker=MARKERS[i], color=COLORS[i],
                markersize=6, linewidth=1.5, label=f"{branch.capitalize()}: root 0")
    ax.set_yscale("log")
    ax.set_ylim(1e-8, 1e-3)
    ax.set_xlim(215, 1070)
    ax.set_xticks([256, 512, 768, 1024])
    ax.set_xlabel(r"Branch-specific DMRG bond dimension $\chi$ (cosh 16)")
    ax.set_ylabel(r"$\|(H-E)\psi\|/\Delta$", labelpad=9)
    ax.set_title("c  Untruncated finite-H residuals", loc="left", pad=14)
    ax.grid(axis="y")
    ax.legend(loc="upper right", frameon=False, fontsize=10)
    ax.text(0.04, 0.28, r"Criterion: $10^{-7}$", transform=ax.transAxes, fontsize=10)
    fig.text(0.055, 0.080, "One checked root per DMRG branch; final doublet stops at 768, singlet at 1024.", fontsize=10)
    fig.text(0.055, 0.061, "Continuation uses one warm trial; root separation is not established.", fontsize=10)
    fig.text(0.055, 0.042, "QP residuals qualify the projected cutoff Hamiltonian, not omitted QP sectors.", fontsize=10)
    fig.text(0.055, 0.023, "NRG zero is exact in stored means (linear scale near zero); no continuum certification.", fontsize=10)
    save(fig, "convergence")


def verify(data, output_dir, with_figures=True):
    """Validate JSON, exact CSV round trips, and self-contained SVG numerical metadata."""
    verify_data(data)
    assert json.loads((output_dir / "results.json").read_text()) == data
    for name, columns, expected, count in zip(["observables", "comparisons"], [OBS_COLUMNS, CMP_COLUMNS],
                                              table_rows(data), [210, 180], strict=True):
        with (output_dir / f"{name}.csv").open(newline="") as stream:
            reader = csv.DictReader(stream)
            assert reader.fieldnames == columns
            rows = list(reader)
        assert len(rows) == len(expected) == count
        for actual, row in zip(rows, expected, strict=True):
            assert actual == {k: "" if row.get(k) is None else str(row[k]) for k in columns}
    if with_figures:
        for name, payload in figure_data(data).items():
            root = ET.parse(output_dir / f"{name}.svg").getroot()
            assert root.findall(".//{http://www.w3.org/2000/svg}text")
            description = root.find(".//{http://purl.org/dc/elements/1.1/}description")
            assert description is not None and json.loads(description.text) == payload
            assert not root.findall(".//{http://purl.org/dc/elements/1.1/}date")
            for node in root.iter():
                assert node.tag.rsplit("}", 1)[-1] not in {"script", "image", "foreignObject"}
                for key, value in node.attrib.items():
                    if key.rsplit("}", 1)[-1] == "href":
                        assert value.startswith("#")
            assert "@import" not in (output_dir / f"{name}.svg").read_text()
    for name in ARTIFACTS if with_figures else ARTIFACTS[:3]:
        text = (output_dir / name).read_text()
        assert not re.search(r"sha256|task_id|budget|attempt-|/project/|\b[0-9a-f]{64}\b", text)
    print("Verified: 14 finite endpoints, 7 finite-Lambda NRG means, 10 observables, 18 comparisons; 210/180 CSV rows.")
    if with_figures:
        print("Verified: both SVGs contain readable text, exact numerical payloads, and no external resources.")


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, default=HERE / "runs/presentation")
    parser.add_argument("--png", action="store_true", help="also create PNG previews")
    parser.add_argument("--tables-only", action="store_true", help="regenerate and validate JSON/CSVs without matplotlib")
    args = parser.parse_args(argv)
    output = args.output_dir.resolve()
    published = (HERE / "output").resolve()
    paths = [output] + [(output / name).resolve() for name in ARTIFACTS + ["comparison.png", "convergence.png"]]
    if any(p == published or published in p.parents for p in paths):
        parser.error("use a separate directory; published outputs are read-only")
    if args.tables_only and args.png:
        parser.error("--png requires figures; omit --tables-only")
    source = (published / "results.json").read_bytes()
    data = json.loads(source)
    verify_data(data)
    output.mkdir(parents=True, exist_ok=True)
    (output / "results.json").write_bytes(source)
    tables(data, output)
    if not args.tables_only:
        figures(data, output, png=args.png)
    verify(data, output, with_figures=not args.tables_only)
    assert (output / "results.json").read_bytes() == source
    print(f"Presentation written to {output}")


if __name__ == "__main__":
    main()

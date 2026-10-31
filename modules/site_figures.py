"""Human-checkable figures for one metal cofactor site (experimental, opt-in).

The same four figures, drawn the same way for any structure and metal cofactor, so that a person can
look at what the tool reports and check it against the structure and the literature before trusting a
table. Every figure carries the question to ask of it; ``checklist.md`` repeats the questions with the
numbers. The figures use only what the analysis modules report; nothing is fitted.

1. ``1_distance_ladder``: every N, O and S atom within 6 A of the metal, by distance, against each
   cutoff the tools use. Shows which atoms a threshold includes and which sit just outside it.
2. ``2_site_geometry``: two orthogonal views of the site with metal-ligand and hydrogen-bond distances.
3. ``3_sulfur_environment``: for each ligand sulfur, the N/O/S atoms within 4 A and whether the contact
   table lists the residue pair and types it polar.
4. ``4_contact_map``: the residue-residue contact map with the sequence-neighbor band and the site residues.

Metal cofactors only (the site centre is the metal atom or atoms of the cofactor residue). The shell
membership comes from :func:`modules.coordination_api.analyze_structure`; contacts and cofactor
classes from :mod:`modules.conserved_interactions`.
"""

from pathlib import Path
from typing import Dict, List, Optional, Sequence

import numpy as np
import pandas as pd

from . import conserved_interactions as ci
from . import coordination_table as ct
from . import protonation
from .coordination_api import analyze_structure
from .io_utils import unpack_pdb_file

# Categorical slots 1 to 4 of the reference palette (light surface); identity is also carried by marker and label.
SURFACE, INK, INK2, MUTED, GRID = "#fcfcfb", "#0b0b0b", "#52514e", "#8a8984", "#e6e5e1"
ROLES = {
    "coordination": ("#2a78d6", "o", "coordination (bond tier, or within the global cutoff when the pair is not in the table)"),
    "weak": ("#eb6834", "D", "weak tier (beyond the bond cutoff, within the weak cutoff of the table)"),
    "primary": ("#1baf7a", "s", "primary shell, not a ligand"),
    "secondary": ("#6250d6", "^", "secondary shell"),
    "other": (MUTED, ".", "other N, O or S"),
}
METALS = {"LI", "NA", "MG", "K", "CA", "MN", "FE", "CO", "NI", "CU", "ZN", "MO", "W", "V", "CD", "HG", "AG", "AU", "PT", "PD", "RU", "IR", "RH", "OS", "RE", "TI", "CR"}
CUTOFFS = (  # (distance, label): the thresholds the tools apply today
    (2.6, "2.6 analyze_structure\ndirect coordination"),
    (2.8, "2.8 cofactor_contacts\ncoord"),
    (3.6, "3.6 analyze_structure\nprimary shell"),
    (4.0, "4.0 contact"),
)
BOND = 1.95                 # heavy-atom pairs closer than this are drawn as covalent bonds


def _mpl():
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 9, "axes.edgecolor": MUTED, "axes.labelcolor": INK2,
                         "xtick.color": INK2, "ytick.color": INK2, "text.color": INK, "figure.facecolor": SURFACE,
                         "axes.facecolor": SURFACE, "savefig.facecolor": SURFACE, "axes.spines.top": False, "axes.spines.right": False})
    return plt


# ----------------------------------------------------------------------------- data
def collect_site(structure_path, cofactor_resname: str, *, chain: Optional[str] = None) -> Dict[str, object]:
    """Everything the figures draw, from the analysis modules (first model, one metal cofactor)."""
    path = str(structure_path)
    analysis = analyze_structure(path, cofactor_resname, direct_coordination=True, first_model_only=True)
    residues = analysis["residues"]
    cof = residues[residues["shell"] == "Cofactor"].iloc[0]
    chain = chain or str(cof["cofactor_chain"])
    structure, _ = unpack_pdb_file(path)
    model = list(structure)[0]
    metals = []
    for residue in model[chain]:
        if residue.id[0] != " " and residue.resname.strip().upper() == cofactor_resname.upper() and residue.id[1] == int(cof["cofactor_residue_number"]):
            metals += [(a.get_id(), np.array(a.coord, float), a.element.upper()) for a in residue if a.element.upper() in METALS]
    if not metals:
        raise ValueError(f"cofactor {cofactor_resname} has no metal atom; site_figures draws metal sites only")
    metal_xyz = np.array([m[1] for m in metals])
    atoms = ci.load_atoms(path, chain, structure_id=Path(path).stem)
    xyz = atoms[["x", "y", "z"]].to_numpy(float)
    d_metal = np.linalg.norm(xyz[:, None, :] - metal_xyz[None, :, :], axis=-1).min(axis=1)
    atoms = atoms.assign(d_metal=d_metal)
    shell = {int(r.residue_number): r.shell for r in residues.itertuples() if r.shell in ("PCS", "SCS") and r.chain == chain}
    features = pd.DataFrame([(atoms["structure_id"].iloc[0], "metal", n, "M", "metal", *c) for n, c, _ in metals],
                            columns=["structure_id", "feature", "atom", "element", "kind", "x", "y", "z"])
    contacts = ci.residue_contacts(atoms)
    metal_element = metals[0][2]
    weak_rules = [r for r in ct.RULES if r.metal == metal_element and r.weak_max is not None]
    weak_window = (min(r.bond_max for r in weak_rules), max(r.weak_max for r in weak_rules)) if weak_rules else None

    def role(row):
        donor = row.element in ("N", "O", "S")
        tier = ct.coordination_tier(metal_element, row.resname, row.atom, row.d_metal, default_max=2.8).tier if donor else "none"
        if tier in ("bond", "default") and shell.get(row.position) == "PCS":
            return "coordination"
        if tier == "weak" and shell.get(row.position) == "PCS":
            return "weak"
        if shell.get(row.position) == "PCS":
            return "primary"
        if shell.get(row.position) == "SCS":
            return "secondary"
        return "other"

    atoms["role"] = [role(r) for r in atoms.itertuples()]
    return {"structure_id": atoms["structure_id"].iloc[0], "cofactor": cofactor_resname, "chain": chain, "atoms": atoms, "metal_xyz": metal_xyz,
            "metal_names": [m[0] for m in metals], "metal_element": metal_element, "shell": shell, "contacts": contacts, "features": features, "analysis": analysis,
            "weak_window": weak_window}


def _label(row) -> str:
    return f"{row.resname.capitalize()}{row.position} {row.atom}"


# ----------------------------------------------------------------------------- figure 1
def figure_distance_ladder(site, path: Path, *, max_distance: float = 5.6):
    plt = _mpl()
    atoms = site["atoms"]
    pool = atoms[atoms["element"].isin(["N", "O", "S"]) & (atoms["d_metal"] <= max_distance)].sort_values("d_metal").reset_index(drop=True)
    carbons = atoms[(atoms["element"] == "C") & (atoms["d_metal"] <= 3.6)]
    height = 0.27 * len(pool) + 3.5
    fig, ax = plt.subplots(figsize=(8.6, height))
    window = site["weak_window"]
    if window:
        ax.axvspan(*window, color=ROLES["weak"][0], alpha=0.10, lw=0)
    for d, text in CUTOFFS:
        ax.axvline(d, color=MUTED, lw=1, ls=(0, (4, 3)))
    y = np.arange(len(pool))[::-1]
    for role, (color, marker, _) in ROLES.items():
        sel = pool["role"] == role
        ax.scatter(pool.loc[sel, "d_metal"], y[sel.to_numpy()], c=color, marker=marker, s=46 if role != "other" else 30, zorder=3, linewidths=0.5, edgecolors=SURFACE)
    for yy, row in zip(y, pool.itertuples()):
        if row.role != "other":
            ax.text(row.d_metal + 0.07, yy, f"{_label(row)}  {row.d_metal:.2f}", va="center", fontsize=8.5, color=INK)
    ax.set_ylim(-1, len(pool) + 1.2)
    ax.set_xlim(1.5, max_distance + 1.4)
    ax.set_yticks([])
    ax.set_xlabel(f"distance to {'/'.join(site['metal_names'])} (A)")
    ax.spines["left"].set_visible(False)
    top = ax.get_ylim()[1]
    for d, text in CUTOFFS:
        high_row = d in (2.8, 4.0)      # alternate two rows so neighboring cutoffs do not collide
        ax.text(d + (0.04 if d in (2.8, 4.0) else -0.04), top + (1.9 if high_row else 0.2), text, ha="left" if d in (2.8, 4.0) else "right", va="bottom", fontsize=7.5, color=INK2)
    ax.set_ylim(-1.6, top)
    if window:
        rules = "; ".join(f"{r.donor_class} bond <= {r.bond_max:g}, weak <= {r.weak_max:g}" for r in ct.RULES if r.metal == site["metal_element"] and r.weak_max is not None)
        ax.text(sum(window) / 2, -1.5, f"weak range from the table: {rules}", ha="center", va="bottom", fontsize=7.5, color=ROLES["weak"][0])
    handles = [plt.Line2D([], [], marker=m, color="none", markerfacecolor=c, markersize=7, label=t) for c, m, t in ROLES.values()]
    ax.legend(handles=handles, loc="upper center", bbox_to_anchor=(0.5, -0.1), ncol=2, fontsize=7.5, frameon=False)
    fig.suptitle(f"1. Distance ladder: {site['structure_id']} {site['cofactor']} site", x=0.02, ha="left", fontsize=12, fontweight="bold", y=1 - 0.12 / height)
    note = (f"Carbon atoms omitted ({len(carbons)} within 3.6 A, nearest {carbons['d_metal'].min():.2f} A: ring and backbone carbons of the ligands).\n"
            "CHECK: are the atoms you know are ligands left of the coordination cutoff? Is any ligand just right of a cutoff?\nIs any non-ligand left of it?")
    fig.text(0.02, 0.1 / height, note, ha="left", va="bottom", fontsize=8, color=INK2)
    fig.subplots_adjust(left=0.03, right=0.98, top=1 - 1.35 / height, bottom=1.95 / height)
    fig.savefig(path, dpi=170)
    plt.close(fig)
    return pool


# ----------------------------------------------------------------------------- figure 2
def _bonds(frame: pd.DataFrame):
    xyz = frame[["x", "y", "z"]].to_numpy(float)
    d = np.linalg.norm(xyz[:, None, :] - xyz[None, :, :], axis=-1)
    i, j = np.nonzero(np.triu((d < BOND) & (d > 0.5), 1))
    return list(zip(i.tolist(), j.tolist()))


def figure_site_geometry(site, path: Path):
    plt = _mpl()
    atoms = site["atoms"]
    drawn = atoms[atoms["role"] != "other"].copy()
    drawn = drawn[drawn["position"].isin(site["shell"])].reset_index(drop=True)
    metal = site["metal_xyz"].mean(axis=0)
    ligand = drawn[drawn["role"].isin(["coordination", "weak"])]
    pts = np.vstack([metal[None, :], ligand[["x", "y", "z"]].to_numpy(float)])
    centre = pts.mean(axis=0)
    _, _, vt = np.linalg.svd(pts - centre)
    xyz = drawn[["x", "y", "z"]].to_numpy(float)
    proj = (xyz - centre) @ vt.T
    mproj = (metal - centre) @ vt.T
    bonds = _bonds(drawn)
    fig, axes = plt.subplots(1, 2, figsize=(10.4, 5.4))
    views = ((0, 1, "face-on view (looking along the least-spread axis of the ligand set)"), (0, 2, "edge-on view"))
    donors = drawn[drawn["role"].isin(["coordination", "weak"])].index
    hbonds = []
    for a in drawn.index[(drawn["role"] == "secondary") & drawn["element"].isin(["N", "O"])]:
        for b in drawn.index[drawn["role"].isin(["coordination", "weak"]) & drawn["element"].isin(["S", "N", "O"])]:
            d = float(np.linalg.norm(xyz[a] - xyz[b]))
            if d <= 3.6 and drawn.at[a, "position"] != drawn.at[b, "position"]:
                hbonds.append((a, b, d))
    for ax, (u, v, title) in zip(axes, views):
        for i, j in bonds:
            ax.plot([proj[i, u], proj[j, u]], [proj[i, v], proj[j, v]], color="#b9b8b2", lw=1.1, zorder=1)
        partners = {a for a, _, _ in hbonds}
        for k, row in drawn.iterrows():
            if k in donors:
                continue
            if k in partners:
                ax.scatter(proj[k, u], proj[k, v], c=ROLES["secondary"][0], marker=ROLES["secondary"][1], s=46, zorder=3, edgecolors=SURFACE, linewidths=0.5)
            else:
                ax.scatter(proj[k, u], proj[k, v], c="#c9c8c2", marker=".", s=12, zorder=2)

        def clear_direction(k, length=30.0):
            """Offset that points away from every line drawn at atom k (bonds, metal bond, hydrogen bonds)."""
            lines = [np.array([proj[o, u] - proj[k, u], proj[o, v] - proj[k, v]]) for i, j in bonds for o in ((j,) if i == k else (i,) if j == k else ())]
            lines += [np.array([mproj[u] - proj[k, u], mproj[v] - proj[k, v]])] if k in donors else []
            lines += [np.array([proj[o, u] - proj[k, u], proj[o, v] - proj[k, v]]) for a, b, _ in hbonds for o in ((b,) if a == k else (a,) if b == k else ())]
            angles = [np.arctan2(l[1], l[0]) for l in lines if np.linalg.norm(l) > 0]
            best, best_gap = 0.0, -1.0
            for cand in np.linspace(0, 2 * np.pi, 24, endpoint=False):
                gap = min((abs((cand - a + np.pi) % (2 * np.pi) - np.pi) for a in angles), default=np.pi)
                if gap > best_gap:
                    best, best_gap = cand, gap
            return length * np.cos(best), length * np.sin(best)

        for k in donors:
            row = drawn.loc[k]
            color, marker, _ = ROLES[row.role]
            dist = float(np.linalg.norm(xyz[k] - metal))
            ax.plot([mproj[u], proj[k, u]], [mproj[v], proj[k, v]], color=color, lw=2.0, ls="-" if row.role == "coordination" else (0, (4, 2)), zorder=2)
            ax.scatter(proj[k, u], proj[k, v], c=color, marker=marker, s=90, zorder=4, edgecolors=SURFACE, linewidths=0.8)
            dx, dy = clear_direction(k)
            ax.annotate(f"{_label(row)}\n{dist:.2f} A", (proj[k, u], proj[k, v]), xytext=(dx, dy), textcoords="offset points", fontsize=8, color=INK,
                        ha="left" if dx > 3 else "right" if dx < -3 else "center", va="bottom" if dy > 3 else "top" if dy < -3 else "center",
                        arrowprops={"arrowstyle": "-", "color": color, "lw": 0.7, "shrinkA": 1, "shrinkB": 5})
        for a, b, d in hbonds:
            ax.plot([proj[a, u], proj[b, u]], [proj[a, v], proj[b, v]], color=ROLES["secondary"][0], lw=1.2, ls=":", zorder=2)
            dx, dy = clear_direction(a, 22.0)
            ax.annotate(f"{_label(drawn.loc[a])}\n{d:.2f} A", (proj[a, u], proj[a, v]), xytext=(dx, dy), textcoords="offset points", fontsize=7.5, color=ROLES["secondary"][0],
                        ha="left" if dx > 3 else "right" if dx < -3 else "center", va="bottom" if dy > 3 else "top" if dy < -3 else "center")
        ax.scatter([mproj[u]], [mproj[v]], c=INK, s=170, zorder=5, edgecolors=SURFACE, linewidths=1.2)
        away = [np.array([proj[k, u] - mproj[u], proj[k, v] - mproj[v]]) for k in donors]
        angles = [np.arctan2(a[1], a[0]) for a in away]
        cand = max(np.linspace(0, 2 * np.pi, 24, endpoint=False), key=lambda c: min(abs((c - a + np.pi) % (2 * np.pi) - np.pi) for a in angles))
        ax.annotate(site["cofactor"], (mproj[u], mproj[v]), xytext=(15 * np.cos(cand), 15 * np.sin(cand)), textcoords="offset points", fontsize=8.5, fontweight="bold", color=INK, ha="center", va="center")
        ax.set_aspect("equal")
        ax.set_title(title, loc="left", fontsize=9, color=INK2)
        ax.set_xticks([])
        ax.set_yticks([])
        for s in ("left", "bottom"):
            ax.spines[s].set_visible(False)
    handles = [plt.Line2D([], [], marker=m, color="none", markerfacecolor=c, markersize=7, label=t) for r, (c, m, t) in ROLES.items() if r != "other"]
    handles += [plt.Line2D([], [], color=ROLES["coordination"][0], lw=2, label="metal-ligand, solid; weak, dashed"),
                plt.Line2D([], [], color=ROLES["secondary"][0], lw=1.2, ls=":", label="N/O to ligand, 3.6 A or less")]
    fig.legend(handles=handles, loc="lower center", ncol=3, fontsize=7.5, frameon=False, bbox_to_anchor=(0.5, 0.075))
    fig.suptitle(f"2. Site geometry: {site['structure_id']} {site['cofactor']} site", x=0.02, ha="left", fontsize=12, fontweight="bold")
    fig.text(0.02, 0.012, "CHECK: does the coordination sphere have the ligand count and geometry you expect? Is each dashed (weak) ligand plausibly bonded?\nAre the dotted partners real hydrogen-bond donors to the ligand? Gray dots are other atoms of the shell residues.", fontsize=8, color=INK2)
    fig.subplots_adjust(left=0.02, right=0.98, top=0.88, bottom=0.2, wspace=0.12)
    fig.savefig(path, dpi=170)
    plt.close(fig)


# ----------------------------------------------------------------------------- figure 3
def _atom_pair_flags(a, b, metals: pd.DataFrame) -> tuple:
    """(polar, polar_s) for this one atom pair, from :func:`residue_contacts` on a two-atom table (``polar_s`` is the opt-in sulfur rule)."""
    frame = pd.DataFrame([("pair", 1, 1, a.resname, a.atom, a.element, a.x, a.y, a.z), ("pair", 9, 9, b.resname, b.atom, b.element, b.x, b.y, b.z)],
                         columns=["structure_id", "position", "seq_index", "resname", "atom", "element", "x", "y", "z"])
    table = ci.residue_contacts(frame, sulfur_acceptor=True, metals=metals.assign(structure_id="pair"))
    return (bool(len(table) and table["polar"].iloc[0]), bool(len(table) and table["polar_s"].iloc[0]))


def sulfur_partners(site, cutoff: float = 4.0) -> pd.DataFrame:
    """N/O/S atoms within ``cutoff`` of each ligand sulfur and how the contact table treats the residue pair."""
    atoms = site["atoms"]
    contacts = site["contacts"].set_index(["position_a", "position_b"])
    metals = pd.DataFrame([(site["structure_id"], site["metal_element"], *xyz) for xyz in site["metal_xyz"]],
                          columns=["structure_id", "element", "x", "y", "z"])
    atoms = atoms.reset_index(drop=True)
    sulfurs = atoms[(atoms["element"] == "S") & atoms["role"].isin(["coordination", "weak"])]
    xyz = atoms[["x", "y", "z"]].to_numpy(float)
    hydrogens = ci.placed_hydrogens(atoms, xyz)
    hydrogens[ci._metal_coordinated(atoms, xyz, metals)] = np.nan                          # a metal-bound N has no H
    rows = []
    for s in sulfurs.itertuples():
        for p in atoms[atoms["element"].isin(["N", "O", "S"]) & (atoms["position"] != s.position)].itertuples():
            d = float(np.linalg.norm(np.array([s.x - p.x, s.y - p.y, s.z - p.z])))
            if d > cutoff:
                continue
            key = (min(s.position, p.position), max(s.position, p.position))
            if abs(s.seq_index - p.seq_index) < 2:
                listed, polar = "excluded: sequence neighbor", None
            elif key in contacts.index:
                listed, polar = "listed", bool(contacts.loc[key, "polar"])
            else:
                listed, polar = "not listed", None
            polar_pair, polar_s_pair = _atom_pair_flags(s, p, metals)
            geometry = [protonation.hydrogen_bond_geometry(h, xyz[p.Index], xyz[s.Index]) for h in hydrogens[p.Index] if not np.isnan(h[0])]
            h_s, angle = min(geometry) if geometry else (None, None)
            rows.append({"sulfur": _label(s), "role": s.role, "partner": _label(p), "partner_role": p.role, "distance": d, "pair_listed": listed, "pair_polar": polar,
                         "partner_backbone_nh": p.atom == "N", "atom_pair_polar": polar_pair, "atom_pair_polar_s": polar_s_pair,
                         "h_s": h_s, "nh_angle": angle,
                         "h_bond_like": bool(geometry and h_s <= protonation.SULFUR_H_MAX and angle >= protonation.SULFUR_ANGLE_MIN)})
    return pd.DataFrame(rows).sort_values(["sulfur", "distance"], ignore_index=True) if rows else pd.DataFrame(
        columns=["sulfur", "role", "partner", "partner_role", "distance", "pair_listed", "pair_polar", "partner_backbone_nh", "atom_pair_polar", "atom_pair_polar_s", "h_s", "nh_angle", "h_bond_like"])


def _typed_text(row) -> str:
    """How the rules type one sulfur-partner atom pair, with the placed-hydrogen geometry when there is one."""
    if row.atom_pair_polar:
        return "polar"
    if not row.atom_pair_polar_s:
        return "not polar" + (f"; N-H {row.h_s:.1f} A, {row.nh_angle:.0f} deg" if row.h_s is not None and not pd.isna(row.h_s) else "")
    if row.h_s is None or pd.isna(row.h_s):
        return "polar_s (capable)"
    return f"polar_s; H...S {row.h_s:.2f} A, {row.nh_angle:.0f} deg" + (" points at S" if row.h_bond_like else " does not point")


def figure_sulfur_environment(site, path: Path):
    plt = _mpl()
    table = sulfur_partners(site)
    names = list(dict.fromkeys(table["sulfur"])) if len(table) else []
    fig, axes = plt.subplots(max(1, len(names)), 1, figsize=(9.2, 1.0 + 1.55 * max(1, len(names)) + 0.3 * len(table)), squeeze=False,
                             gridspec_kw={"height_ratios": [max(1, (table["sulfur"] == n).sum()) for n in names] or [1]})
    for ax, name in zip(axes[:, 0], names or [None]):
        sub = table[table["sulfur"] == name].reset_index(drop=True) if name else table
        if name is None:
            ax.axis("off")
            ax.text(0.5, 0.5, "no ligand sulfur in this site", ha="center")
            continue
        y = np.arange(len(sub))[::-1]
        for yy, row in zip(y, sub.itertuples()):
            color, marker, _ = ROLES[row.partner_role]
            ax.plot([2.9, row.distance], [yy, yy], color="#cfcec8", lw=1.2, zorder=1)
            ax.scatter([row.distance], [yy], c=color, marker=marker, s=50, zorder=3, edgecolors=SURFACE, linewidths=0.5)
            ax.text(row.distance + 0.06, yy, f"{row.partner}  {row.distance:.2f}", va="center", fontsize=8.5)
            if row.pair_listed == "listed":
                text = "pair listed" + ("; flagged polar by another atom pair" if row.pair_polar else "; not flagged polar")
            else:
                text = {"excluded: sequence neighbor": "pair excluded (sequence neighbor)", "not listed": "pair not listed"}[row.pair_listed]
            ax.text(5.2, yy, text, va="center", fontsize=8.5, color=INK2)
            typed = _typed_text(row)
            ax.text(8.1, yy, typed, va="center", ha="left", fontsize=8.5, color=INK if not typed.startswith("not polar") else ROLES["weak"][0],
                    fontweight="bold" if row.partner_backbone_nh and typed.startswith("not polar") else "normal")
        ax.set_xlim(2.9, 11.2)
        ax.set_xticks([3.0, 3.5, 4.0])
        ax.set_ylim(-0.7, len(sub) - 0.3)
        ax.set_yticks([])
        ax.spines["left"].set_visible(False)
        ax.set_title(f"{name} ({'ligand' if sub['role'].iloc[0] == 'coordination' else 'weak or axial'})", loc="left", fontsize=9.5, color=INK, fontweight="bold")
        ax.text(5.2, len(sub) - 0.1, "residue pair in the contact table", fontsize=7.5, color=MUTED, va="bottom")
        ax.text(8.1, len(sub) - 0.1, "atom pair typed as", fontsize=7.5, color=MUTED, va="bottom")
        ax.axvline(4.0, color=MUTED, lw=0.8, ls=(0, (4, 3)))
    axes[-1, 0].set_xlabel("sulfur to partner atom (A)")
    fig.suptitle(f"3. Ligand sulfur environment: {site['structure_id']} {site['cofactor']} site", x=0.02, ha="left", fontsize=12, fontweight="bold", y=0.99)
    fig.text(0.02, 0.01, "CHECK: a backbone N-H within 3.7 A of a ligand Cys S (3.6 A for Met) is a candidate hydrogen-bond donor to it (placed-hydrogen geometry shown).\n"
                         "Is every such atom pair typed? 'polar' needs an N/O pair, so a sulfur contact is never typed by it;\n"
                         "'polar_s' is the opt-in sulfur-acceptor rule (sulfur_acceptor=True).", fontsize=8, color=INK2)
    fig.subplots_adjust(left=0.03, right=0.98, top=0.88 if len(names) > 1 else 0.82, bottom=0.2 if len(names) > 1 else 0.3, hspace=0.6)
    fig.savefig(path, dpi=170)
    plt.close(fig)
    return table


# ----------------------------------------------------------------------------- figure 4
def figure_contact_map(site, path: Path, *, window: float = 0.05):
    plt = _mpl()
    from matplotlib.colors import LinearSegmentedColormap

    contacts = site["contacts"]
    positions = sorted(site["atoms"]["position"].unique())
    index = {p: k for k, p in enumerate(positions)}
    n = len(positions)
    grid = np.full((n, n), np.nan)
    for r in contacts.itertuples():
        a, b = index[r.position_a], index[r.position_b]
        grid[a, b] = grid[b, a] = r.min_distance
    cmap = LinearSegmentedColormap.from_list("close", ["#0d366b", "#2a78d6", "#9ec5f4", "#e1ecfa"])
    fig, ax = plt.subplots(figsize=(7.4, 7.6))
    idx = np.arange(n)
    band = np.abs(idx[:, None] - idx[None, :]) < 2
    ax.imshow(np.where(band, 1.0, np.nan), cmap="Greys", vmin=0, vmax=3, origin="upper", interpolation="nearest", alpha=0.55)
    image = ax.imshow(grid, cmap=cmap, vmin=1.9, vmax=4.0, origin="upper", interpolation="nearest")
    near = contacts[contacts["min_distance"] >= 4.0 - window]
    ax.scatter([index[p] for p in near.position_a] + [index[p] for p in near.position_b], [index[p] for p in near.position_b] + [index[p] for p in near.position_a],
               s=26, facecolors="none", edgecolors=ROLES["weak"][0], linewidths=0.9, zorder=3)
    shell = site["shell"]
    last = -99.0
    for p, shell_name in sorted(shell.items()):
        if p not in index:
            continue
        color = ROLES["primary"][0] if shell_name == "PCS" else ROLES["secondary"][0]
        ax.axhline(index[p], color=color, lw=0.5, alpha=0.5)
        ax.axvline(index[p], color=color, lw=0.5, alpha=0.5)
        name = site["atoms"].loc[site["atoms"]["position"] == p, "resname"].iloc[0].capitalize()
        ty = max(index[p], last + 2.6)
        last = ty
        ax.annotate(f"{name}{p}", (n - 0.5, index[p]), xytext=(n + 1.2, ty), fontsize=7.5, color=color, va="center", annotation_clip=False,
                    arrowprops={"arrowstyle": "-", "color": color, "lw": 0.5, "alpha": 0.7})
    step = 10 if n > 60 else 5
    ax.set_xticks(range(0, n, step))
    ax.set_xticklabels([str(positions[k]) for k in range(0, n, step)])
    ax.set_yticks(range(0, n, step))
    ax.set_yticklabels([str(positions[k]) for k in range(0, n, step)])
    ax.set_xlim(-0.5, n - 0.5)
    ax.set_ylim(n - 0.5, -0.5)
    cbar = fig.colorbar(image, ax=ax, fraction=0.035, pad=0.14)
    cbar.set_label("closest heavy-atom distance (A)")
    cbar.ax.invert_yaxis()
    fig.suptitle(f"4. Contact map: {site['structure_id']} chain {site['chain']}", x=0.02, ha="left", fontsize=12, fontweight="bold", y=0.985)
    fig.text(0.02, 0.012, f"{len(contacts)} contacts (1.9-4.0 A, at least 2 apart in sequence). Gray band: pairs less than 2 apart, excluded by the rule.\n"
                          f"Orange rings: {len(near)} contacts within {window} A of the 4.0 A cutoff. Lines: primary (green) and secondary (violet) shell.\n"
                          "CHECK: do helices, strands and site residues show the pattern you expect? Does the gray band hide a contact you care about?",
             fontsize=8, color=INK2)
    fig.subplots_adjust(left=0.08, right=0.84, top=0.93, bottom=0.16)
    fig.savefig(path, dpi=170)
    plt.close(fig)


# ----------------------------------------------------------------------------- driver
def site_check_figures(structure_path, cofactor_resname: str, output_dir, *, chain: Optional[str] = None) -> Dict[str, Path]:
    """Write the four standard figures and ``checklist.md`` for one metal site; returns the paths."""
    out = Path(output_dir)
    out.mkdir(parents=True, exist_ok=True)
    site = collect_site(structure_path, cofactor_resname, chain=chain)
    paths = {name: out / f"{name}.png" for name in ("1_distance_ladder", "2_site_geometry", "3_sulfur_environment", "4_contact_map")}
    ladder = figure_distance_ladder(site, paths["1_distance_ladder"])
    figure_site_geometry(site, paths["2_site_geometry"])
    sulfur = figure_sulfur_environment(site, paths["3_sulfur_environment"])
    figure_contact_map(site, paths["4_contact_map"])
    checklist = out / "checklist.md"
    checklist.write_text(_checklist(site, ladder, sulfur))
    paths["checklist"] = checklist
    return paths


def _checklist(site, ladder: pd.DataFrame, sulfur: pd.DataFrame) -> str:
    window = site["weak_window"]
    lines = [f"# Site check: {site['structure_id']} {site['cofactor']} (chain {site['chain']}, first model)", "",
             "Look at each figure and answer the question. A person who knows the system should be able to say yes to every one.", "",
             "## 1. Distance ladder (`1_distance_ladder.png`)",
             f"- Coordination tier (table bond cutoff for listed donors, global 2.8 A otherwise) and in the primary shell: {', '.join(_label(r) for r in ladder[ladder.role == 'coordination'].itertuples()) or 'none'}.",
             f"- Weak tier ({'table: ' + ', '.join(f'{r.donor_class} {r.bond_max:g} to {r.weak_max:g} A' for r in ct.RULES if r.metal == site['metal_element'] and r.weak_max is not None) if window else 'no weak range in the table for this metal'}): {', '.join(f'{_label(r)} {r.d_metal:.2f}' for r in ladder[ladder.role == 'weak'].itertuples()) or 'none'}.",
             f"- Primary shell but not a ligand: {', '.join(f'{_label(r)} {r.d_metal:.2f}' for r in ladder[ladder.role == 'primary'].itertuples()) or 'none'}.",
             "- Question: are the ligands you know left of the cutoff, and is any known ligand just right of one?", "",
             "## 2. Site geometry (`2_site_geometry.png`)",
             "- Question: is the ligand count and arrangement what the literature gives? Is each dashed (weak) ligand plausibly bonded?", "",
             "## 3. Ligand sulfur environment (`3_sulfur_environment.png`)"]
    if len(sulfur):
        for r in sulfur.itertuples():
            lines.append(f"- {r.sulfur} ... {r.partner} {r.distance:.2f} A: pair {r.pair_listed}; this atom pair typed {_typed_text(r)}")
    else:
        lines.append("- no ligand sulfur")
    lines += ["- Question: is every backbone N-H donor to a ligand sulfur typed (polar_s with the opt-in rule)?", "",
              "## 4. Contact map (`4_contact_map.png`)",
              f"- {len(site['contacts'])} contacts; {int((site['contacts']['min_distance'] >= 3.95).sum())} within 0.05 A of the 4.0 A cutoff.",
              "- Question: does the pattern match the secondary structure you know, and does the sequence-neighbor band hide a contact you care about?", "",
              "## Rules in force (cutoffs used for every number above)", "",
              ci.format_rules(pd.concat([ci.class_rules("residue_contacts", sulfur_acceptor=True),
                                         ci.class_rules("cofactor_contacts", coordination_model="by_pair")], ignore_index=True))]
    return "\n".join(lines)


def main(argv: Optional[Sequence[str]] = None) -> int:
    import argparse

    parser = argparse.ArgumentParser(description="Write the standard human-check figures for one metal cofactor site.")
    parser.add_argument("structure")
    parser.add_argument("cofactor", help="cofactor residue name, for example CU")
    parser.add_argument("--chain")
    parser.add_argument("--output-dir", default=None, help="default site_check_output/<structure name>")
    args = parser.parse_args(argv)
    out = args.output_dir or str(Path("site_check_output") / Path(args.structure).stem)
    for name, path in site_check_figures(args.structure, args.cofactor, out, chain=args.chain).items():
        print(name, path)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

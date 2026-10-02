"""Enforce the stable-core / experimental split described in docs/stability.md.

The stable core is what a training pipeline depends on. It must not import the
plotting, comparison, wire, or substrate code, and must import without plotly or
matplotlib installed.
"""

import ast
from pathlib import Path
import subprocess
import sys
import unittest

ROOT = Path(__file__).resolve().parents[1]
MODULES = ROOT / "modules"

STABLE_CORE = {
    "chemistry", "moieties", "motif_registry", "structure_utils", "structure_processing",
    "io_utils", "cofactor_classes", "coordination_api", "assembly_policy", "ml_export",
    "sensitivity",
}
EXPERIMENTAL = {
    "plotting", "analysis", "display", "network_comparison", "network_alignment",
    "comparison_runner", "comparison_plots", "conservation_viewer", "protein_wires",
    "wire_viewer", "substrate_seeds", "substrate_viewer", "legacy_adapter",
}
HEAVY_PLOTTING = {"plotly", "matplotlib", "mpl_toolkits"}


def local_imports(name):
    tree = ast.parse((MODULES / f"{name}.py").read_text())
    found = set()
    for node in ast.walk(tree):
        if isinstance(node, ast.ImportFrom) and node.level == 1:
            if node.module:
                found.add(node.module.split(".")[0])
            else:
                found |= {alias.name for alias in node.names}
        elif isinstance(node, ast.ImportFrom) and node.module and node.module.startswith("modules."):
            found.add(node.module.split(".")[1])
    return {n for n in found if (MODULES / f"{n}.py").exists()}


def external_imports(name):
    tree = ast.parse((MODULES / f"{name}.py").read_text())
    found = set()
    for node in ast.walk(tree):
        if isinstance(node, ast.Import):
            found |= {alias.name.split(".")[0] for alias in node.names}
        elif isinstance(node, ast.ImportFrom) and node.level == 0 and node.module:
            found.add(node.module.split(".")[0])
    return found


class ModuleBoundaryTests(unittest.TestCase):
    def test_every_module_is_classified(self):
        on_disk = {p.stem for p in MODULES.glob("*.py")} - {"__init__"}
        classified = STABLE_CORE | EXPERIMENTAL
        legacy_support = {"atom_utils", "constants", "deduplicate_chains", "reporting",
                          "residue_matching", "structure_io"}
        self.assertEqual(on_disk - classified - legacy_support, set(),
                         "new module: add it to a tier here and in docs/stability.md")
        self.assertEqual(STABLE_CORE & EXPERIMENTAL, set())

    def test_stable_core_imports_only_stable_core(self):
        for name in sorted(STABLE_CORE):
            with self.subTest(module=name):
                self.assertLessEqual(local_imports(name), STABLE_CORE, name)
                self.assertEqual(external_imports(name) & HEAVY_PLOTTING, set(), name)

    def test_stable_core_imports_without_plotting_libraries(self):
        blocker = (
            "import sys\n"
            "for name in ('plotly', 'matplotlib', 'mpl_toolkits'):\n"
            "    sys.modules[name] = None\n"
            + "".join(f"import modules.{name}\n" for name in sorted(STABLE_CORE))
        )
        done = subprocess.run([sys.executable, "-c", blocker], cwd=ROOT, capture_output=True, text=True)
        self.assertEqual(done.returncode, 0, done.stderr)


if __name__ == "__main__":
    unittest.main()

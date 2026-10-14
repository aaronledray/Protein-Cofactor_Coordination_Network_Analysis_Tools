"""mmCIF files may omit the optional insertion-code column."""

from pathlib import Path
import tempfile
import unittest

from modules.coordination_api import analyze_structure

ROOT = Path(__file__).resolve().parents[1]
PLASTOCYANIN = ROOT / "reference_structures/0_Plastocyanin/1ag6.cif"


def strip_insertion_code_column(text):
    """Remove _atom_site.pdbx_PDB_ins_code and its values from the atom loop."""
    lines = text.splitlines()
    start = next(i for i, l in enumerate(lines) if l.startswith("_atom_site.group_PDB"))
    header = []
    i = start
    while lines[i].startswith("_atom_site."):
        header.append(lines[i].strip()); i += 1
    drop = header.index("_atom_site.pdbx_PDB_ins_code")
    out = lines[:start] + [h for h in header if h != "_atom_site.pdbx_PDB_ins_code"]
    while i < len(lines) and not lines[i].startswith(("#", "loop_", "_")):
        fields = lines[i].split()
        if len(fields) == len(header):
            del fields[drop]
            out.append(" ".join(fields))
        else:
            out.append(lines[i])
        i += 1
    return "\n".join(out + lines[i:]) + "\n"


class MmcifToleranceTests(unittest.TestCase):
    def test_missing_insertion_code_column_gives_same_tables(self):
        text = PLASTOCYANIN.read_text()
        self.assertIn("_atom_site.pdbx_PDB_ins_code", text)
        stripped = strip_insertion_code_column(text)
        self.assertNotIn("_atom_site.pdbx_PDB_ins_code", stripped)
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "1ag6.cif"
            path.write_text(stripped)
            tables = analyze_structure(path, "CU", shells=3)
        reference = analyze_structure(PLASTOCYANIN, "CU", shells=3)
        for name in ("residues", "atoms", "links", "contacts"):
            self.assertTrue(
                tables[name].reset_index(drop=True).equals(reference[name].reset_index(drop=True)), name
            )


if __name__ == "__main__":
    unittest.main()

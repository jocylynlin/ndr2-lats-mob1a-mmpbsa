#!/usr/bin/env python3
"""Print '1-N' where N is the last non-solvent residue of a prmtop (solvent = WAT/Na+/Cl-)."""
import re, sys
from pathlib import Path
t = Path(sys.argv[1]).read_text(errors="replace")
F = {m.group(1): m.group(3) for m in re.finditer(r"%FLAG (\S+)\s*\n%FORMAT\(([^)]+)\)\s*\n(.*?)(?=%FLAG|\Z)", t, re.S)}
s = F["RESIDUE_LABEL"].replace("\n", "")
lab = [s[i:i+4].strip() for i in range(0, len(s), 4)]
n = max(i for i, x in enumerate(lab, 1) if x not in ("WAT", "Na+", "Cl-", "HOH"))
print(f"1-{n}")

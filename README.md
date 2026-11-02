# MM-PBSA of NDR2/LATS1-MOB1A complexes | data and analysis deposit

Supporting data, inputs and analysis code for:

> *Lessons Learned from Modeling Highly Charged Protein–Protein Interactions with MM–PBSA:
> A Case Study of NDR2/LATS1–MOB1A Complexes.*
> T. Lin, Z. Huang, W. Wang and R. Luo. Submitted to *Int. J. Mol. Sci.*

Everything reported in the paper and its Supporting Information can be recomputed from the
material here plus the per-frame energies archived separately (see **Archived separately**).

---

## 1. What was simulated

Three **conditions**, each with six **species**, each with **100 independent replicas**:

| | |
|---|---|
| Conditions | `opc-salt` (reference: OPC water, 80 added Na⁺/Cl⁻ ion pairs per box), `tip3p` (TIP3P water, same ionic content), `counterion` (OPC, neutralising counterions only) |
| Species | `lats1_complex`, `lats1_kinase`, `lats1_mob1a`, `ndr2_complex`, `ndr2_kinase`, `ndr2_mob1a` |
| Replicas | 100 per species per condition (1,800 trajectories) |
| Protocol | 30,000 steps steepest descent + 30,000 conjugate gradient; 1 ns NVT heating/equilibration at 310.15 K; 50 ns NPT equilibration; **40 ns production**, coordinates every 10 ps (4,000 frames) |
| Integrator | 2 fs timestep, 10 Å real-space cutoff, default PME |
| Force field | AMBER ff19SB; OPC or TIP3P water; Li–Merz 12-6 monovalent ions; ZAFF bonded model for the MOB1A zinc site |

Two additional sets were run for the revision:

- **Window test** — 20 replicas per system in the reference condition extended 40 → 100 ns,
  with MM–PBSA recomputed over the 40–100 ns window (`extend_c1.sh`, `late_mmpbsa_parallel.sh`).
- **Solvent-retaining set** — unbound kinases under each water model, 5 replicas × 40 ns with
  solvent coordinates written every 100 ps, for hydration and ion analysis (`setup_c2.sh`,
  `c2_hydration_one.sh`).

### Construct definitions

MOB1A UniProt Q9H8S9 residues 33–212; LATS1 residues 640–698 (from PDB **5BRK**); NDR2 residues
25–83 (from PDB **5XQZ chain B**). Five NDR2 glutamates (Glu40, Glu43, Glu55, Glu61, Glu62) are
modelled neutral so that both kinase segments carry net charge +6 — a deliberate control for the
between-system comparison, discussed in Section 2.4 of the paper. All chain termini are uncapped
and charged.

---

## 2. Software

| | |
|---|---|
| MD | AMBER 2024, `pmemd.cuda` |
| MM–PBSA / normal modes | AmberTools 2024, `MMPBSA.py` |
| Trajectory analysis | `cpptraj` (AmberTools 2024) |
| Analysis / figures | Python 3.13.2, numpy 2.2.5, scipy 1.15.3, pandas 2.3.1, matplotlib 3.10.0 |

MM–PBSA settings: internal PBSA solver, `mbondi` radii, solute dielectric **20.0**, exterior
dielectric 80.0, ionic strength **0.150 M** (held fixed for all conditions, independent of the
MD salt content), `inp=2` nonpolar model. Snapshots every 100 ps from the 40 ns production run
(400 per replica). Normal modes: 10 equally spaced snapshots per replica, minimised to
DRMS < 0.001 kcal/mol/Å; entropies computed by `MMPBSA.py` at its default 298.15 K and reported
as −TΔS at 310.15 K.

---

## 3. Layout

```
scripts/simulation/   system construction, MD, MM–PBSA, normal modes, cpptraj analyses
scripts/analysis/     the scripts that produce every published figure and table
data/                 derived tables the analysis scripts read
topologies/           AMBER prmtop files for all 18 system × condition combinations
inputs/               tleap, minimisation, equilibration, production, MM–PBSA and nmode inputs
MANIFEST.sha256       checksums for every file in this deposit
```

---

## 4. Provenance: every figure and table

| Output | Script | Reads |
|---|---|---|
| **Figure 1** covalent decomposition | `analysis/plot_covalent_decomp_n100.py` | per-replica MM–PBSA energies |
| **Figures 2, 3** free-energy surfaces | `analysis/fes_physcv.py` | `physcv` CV time series, `data/cv_metrics.json` |
| **Table 1** protocol comparison | `analysis/table1.py` | per-replica MM–PBSA energies |
| **Tables 2, 3** water / ionic comparison | `analysis/table1.py`, `data/species_breakdown_n100.csv`, `data/species_breakdown_counterion.csv` | per-replica MM–PBSA energies |
| **Figures S1, S2** time-series stability | `analysis/convplot.py` | per-frame MM–PBSA energies |
| **Figure S3** replica distributions | `analysis/plot_figS3_replica_distributions.py` | per-replica MM–PBSA energies |
| **Figure S4** bootstrap sign stability | `analysis/plot_covalent_bootstrap_n100.py` | per-replica covalent components |
| **Figures S5, S6** RDFs | `analysis/plot_c2_rdf_si.py` | `cpptraj` RDFs from the solvent-retaining set |
| **Tables S1, S2** covalent decomposition | `analysis/covdecomp.py` | per-replica covalent components |
| **Table S3** normal-mode entropy | `analysis/ent_pool.py` | per-replica nmode energies |
| **Table S4** snapshot sensitivity | `analysis/nmode_snapshot_convergence.py` | `data/nmode_snapshot_convergence_n100.csv` |
| **Table S5** production-window test | `analysis/ddh_window.py` | 0–40 ns and 40–100 ns MM–PBSA on the same 20 replicas |
| **Table S6** conformational metrics | `analysis/cv_metrics.py` | `physcv` CV time series → `data/cv_metrics.json` |
| Drift test (reported in the repository, not the paper) | `analysis/drifttest.py` | `data/convergence_drift_tests.txt` |

---

## 5. Two conventions that matter when recomputing

**Pooled unbound MOB1A reference.** MOB1A is the same molecule in both complexes with
parameter-identical topologies, so within each condition the two unbound MOB1A ensembles
(100 replicas each) are **pooled into one 200-replica reference** used for both systems.
Pooling is never done across conditions. With a shared reference the term cancels exactly in the
LATS1-minus-NDR2 contrast. Using each system's own MOB1A instead changes the multi-trajectory
contrast substantially — the two sub-ensembles differ by 8.46 ± 1.28 kcal/mol under the reference
condition — so results computed with a paired reference will not match the paper.
See `scripts/simulation/mt_pooled_reference.py`.

**Replica-level uncertainties.** Every reported ± is a replica-level standard error:
SE = s_replica / √N with N = 100, where s_replica is the standard deviation of the N replica
means. Frame-level errors, which treat the 400 snapshots within a replica as independent,
are roughly three times smaller and are *not* what the paper reports. Differences between
systems use the standard error of the paired per-replica differences, in which the shared
MOB1A reference cancels before the standard error is taken.

---

## 6. Archived separately

Too large for version control, deposited at the DOI given in the paper's Data Availability
Statement:

- Per-frame, per-replica MM–PBSA energies for every system and condition (~1.6 GB)
- Per-replica normal-mode energies for the multi-trajectory reference condition (~3.6 GB;
  this is the only condition for which normal modes were computed)
- `cpptraj` analysis outputs: RDFs, hydrogen bonds, water shells, SASA, salt bridges,
  collective-variable time series (~1.3 GB)
- Equilibrated restart coordinates (`.rst`) and solvated complex topologies

## 7. Not retained

- **Solvated MD trajectories** (~324 GB) are not archived. Water-stripped solute trajectories
  for the reference condition are available from the corresponding author on request.
- A small number of **restart files and random seeds from the earliest runs** were not kept.
  Those replicas cannot be reproduced bit-for-bit, though the ensemble statistics are unaffected.
  `scripts/simulation/check_seeds.sh` reports which replicas retain seeds.

---

## 8. Running these scripts elsewhere

The scripts are deposited **exactly as executed**, so they are a faithful record of what
produced the published numbers rather than a polished package. 34 of the 49 scripts contain
absolute paths from the machine they ran on and will need editing before they run anywhere else.
The paths to change are:

| Path in the scripts | What it is |
|---|---|
| `/home3/jocelyn/projects/hippo/campaign_100x40ns` | campaign root: `runs/`, `inputs/`, `topologies/`, `mmpbsa/`, `nmode/`, `analysis/` |
| `/home3/jocelyn/projects/hippo/reviewer_comments` | analysis root for the revision scripts |
| `/home3/jocelyn/apps/amber25/bin`, `/home3/apps/amber24/bin` | AMBER and AmberTools executables |
| `/home3/jocelyn/miniconda3/bin` | Python used for the analysis scripts |

The campaign directory layout the scripts expect is:

```
campaign_100x40ns/
  topologies/<condition>/<species>.prmtop
  inputs/<condition>/<species>/            tleap, min, heat, equil, prod inputs
  runs/<condition>/<species>/replica<NNN>/ trajectories and restarts
  mmpbsa/<condition>/<protocol>/<system>_<NNN>/
  nmode/<condition>/<system>_<NNN>/
  analysis/                                 cpptraj outputs and tidy tables
```

where `<condition>` is `opc-salt`, `tip3p` or `counterion`, `<species>` is one of the six listed
in Section 1, `<protocol>` is `st` or `mt`, and `<NNN>` is the zero-padded replica index 001–100.

The SLURM submission scripts (`queue_*.sh`, `*_slurm.sh`, `rolling_submit*.sh`) encode the
partitions and resource limits of the cluster used and will need rewriting for any other
scheduler. The analysis scripts in `scripts/analysis/` have no scheduler dependency and need
only the path edits above.

---

## 9. Contact

Ray Luo — rluo@uci.edu

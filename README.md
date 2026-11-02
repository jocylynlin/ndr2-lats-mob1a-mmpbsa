# MM-PBSA of NDR2/LATS1-MOB1A complexes | data and analysis deposit

Supporting data, inputs and analysis code for:

> *Lessons Learned from Modeling Highly Charged Protein-Protein Interactions with MM-PBSA:
> A Case Study of NDR2/LATS1-MOB1A Complexes.*
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
| Force field | AMBER ff19SB; OPC or TIP3P water; Li-Merz 12-6 monovalent ions; ZAFF bonded model for the MOB1A zinc site |

Two additional sets were run for the revision:

- **Window test** - 20 replicas per system in the reference condition extended 40 → 100 ns,
  with MM-PBSA recomputed over the 40-100 ns window (`extend_c1.sh`, `late_mmpbsa_parallel.sh`).
- **Solvent-retaining set** - unbound kinases under each water model, 5 replicas × 40 ns with
  solvent coordinates written every 100 ps, for hydration and ion analysis (`setup_c2.sh`,
  `c2_hydration_one.sh`).

### Construct definitions

MOB1A UniProt Q9H8S9 residues 33-212; LATS1 residues 640-698 (from PDB **5BRK**); NDR2 residues
25-83 (from PDB **5XQZ chain B**). Five NDR2 glutamates (Glu40, Glu43, Glu55, Glu61, Glu62) are
modelled neutral so that both kinase segments carry net charge +6 - a deliberate control for the
between-system comparison, discussed in Section 2.4 of the paper. All chain termini are uncapped
and charged.

---

## 2. Software

| | |
|---|---|
| MD | AMBER 2024, `pmemd.cuda` |
| MM-PBSA / normal modes | AmberTools 2024, `MMPBSA.py` |
| Trajectory analysis | `cpptraj` (AmberTools 2024) |
| Analysis / figures | Python 3.13.2, numpy 2.2.5, scipy 1.15.3, pandas 2.3.1, matplotlib 3.10.0 |

MM-PBSA settings: internal PBSA solver, `mbondi` radii, solute dielectric **20.0**, exterior
dielectric 80.0, ionic strength **0.150 M** (held fixed for all conditions, independent of the
MD salt content), `inp=2` nonpolar model. Snapshots every 100 ps from the 40 ns production run
(400 per replica). Normal modes: 10 equally spaced snapshots per replica, minimised to
DRMS < 0.001 kcal/mol/Å; entropies computed by `MMPBSA.py` at its default 298.15 K and reported
as −TΔS at 310.15 K.

---

## 3. Layout

```
scripts/simulation/   system construction, MD, MM-PBSA, normal modes, cpptraj analyses
scripts/analysis/     the scripts that produce every published figure and table
data/                 derived summary tables the analysis scripts read
energies/             per-frame MM-PBSA and normal-mode energies, and the cpptraj
                      hydration outputs - the primary numerical record, gzipped
topologies/           dry (solute-only) AMBER topologies used for MM-PBSA
md_systems/           solvated topologies and equilibrated restarts, gzipped - rerun the MD
inputs/               tleap, minimisation, equilibration, production, MM-PBSA and nmode inputs
MANIFEST.sha256       checksums for every file in this deposit
```

---

## 4. Two topology sets - use the right one

| Directory | Contents | Use |
|---|---|---|
| `topologies/` | 3,937 atoms, **no solvent** | MM-PBSA and trajectory analysis, which operate on the stripped solute |
| `md_systems/` | 59,818 atoms, ~13,930 waters plus ions | running or restarting the MD |

`md_systems/` mirrors the `<condition>/<species>/` layout and holds, for each of the 18 system
× condition combinations, the solvated topology and the restart coordinates. **Restart from
`min2.rst`** - it is the state at the end of the 50 ns NPT equilibration, which is where
production began. The six kinase systems also carry `min1.rst`, the intermediate state between
the steepest-descent and conjugate-gradient minimisation stages; it is included for
completeness and is not the production starting point. The `build_zaff/` subdirectories under
the kinase systems hold a byte-identical copy of the parent topology, retained because the
build scripts reference that path. Both are gzipped - 213 MB becomes 37 MB - so decompress
before use:

```bash
gunzip -r md_systems/
```

To restart production for one system:

```bash
cd md_systems/opc-salt/lats1_complex
pmemd.cuda -O -i ../../../inputs/opc-salt/lats1_complex/prod.in \
           -p moblats.prmtop -c min2.rst -o prod.out -r prod.rst -x prod.nc
```

Each of the 100 replicas differs only in the random seed used for heating and thermal
equilibration; `scripts/simulation/generate_replicas.sh` shows how they were generated, and
`scripts/simulation/check_seeds.sh` reports which seeds were retained.

---

## 5. The primary numerical record: `energies/`

Every value in the paper derives from these files. All are gzipped; `gunzip -r energies/`
decompresses in place.

| Path | Contents |
|---|---|
| `energies/mmpbsa_perframe/*.csv.gz` | per-frame MM-PBSA energies, one table per condition and protocol |
| `energies/nmode/<condition>/<system>/traj_<NNN>/per_frame.csv.gz` | normal-mode energies, 10 snapshots per replica, reference condition only |
| `energies/hydration/<condition>/<species>/replica<NNN>/` | cpptraj outputs from the solvent-retaining set: RDFs, hydrogen bonds, water shells, SASA |

The MM-PBSA tables are long-format, one row per frame, with columns
`system, water_model, salt_condition, protocol, delta_method, replica, frame, traj_frame,
time_ns, time_ns_abs` followed by the energy terms (`BOND, ANGLE, DIHED, VDWAALS, EEL, EPB,
ENPOLAR, EDISPER` and the grouped contributions). File naming:

- `_st` / `_mt` - single- or multi-trajectory protocol
- `_pooledref` - uses the pooled 200-replica unbound MOB1A reference. **These are the tables
  the paper's multi-trajectory numbers come from**; the plain `_mt` tables use each system's
  own MOB1A and will not reproduce the published values (see Section 7)
- `_late` - the 40-100 ns window from the 20 extended replicas, used for Table S5
- `_rec` - receptor-only evaluation used to build the pooled reference

---

## 6. Provenance: every figure and table

| Output | Script | Reads |
|---|---|---|
| **Figure 1** covalent decomposition | `analysis/plot_covalent_decomp_n100.py` | `energies/mmpbsa_perframe/` |
| **Figures 2, 3** free-energy surfaces | `analysis/fes_physcv.py` | `physcv` CV time series, `data/cv_metrics.json` |
| **Table 1** protocol comparison | `analysis/table1.py` | `energies/mmpbsa_perframe/` |
| **Tables 2, 3** water / ionic comparison | `analysis/table1.py`, `data/species_breakdown_n100.csv`, `data/species_breakdown_counterion.csv` | `energies/mmpbsa_perframe/` |
| **Figures S1, S2** time-series stability | `analysis/convplot.py` | `energies/mmpbsa_perframe/` |
| **Figure S3** replica distributions | `analysis/plot_figS3_replica_distributions.py` | `energies/mmpbsa_perframe/` |
| **Figure S4** bootstrap sign stability | `analysis/plot_covalent_bootstrap_n100.py` | `energies/mmpbsa_perframe/` |
| **Figures S5, S6** RDFs | `analysis/plot_c2_rdf_si.py` | `energies/hydration/*/rdf_*.dat` |
| **Tables S1, S2** covalent decomposition | `analysis/covdecomp.py` | `energies/mmpbsa_perframe/` |
| **Table S3** normal-mode entropy | `analysis/ent_pool.py` | `energies/nmode/` |
| **Table S4** snapshot sensitivity | `analysis/nmode_snapshot_convergence.py` | `data/nmode_snapshot_convergence_n100.csv` |
| **Table S5** production-window test | `analysis/ddh_window.py` | 0-40 ns and 40-100 ns MM-PBSA on the same 20 replicas |
| **Table S6** conformational metrics | `analysis/cv_metrics.py` | `physcv` CV time series → `data/cv_metrics.json` |
| Drift test (reported in the repository, not the paper) | `analysis/drifttest.py` | `data/convergence_drift_tests.txt` |

---

## 7. Two conventions that matter when recomputing

**Pooled unbound MOB1A reference.** MOB1A is the same molecule in both complexes with
parameter-identical topologies, so within each condition the two unbound MOB1A ensembles
(100 replicas each) are **pooled into one 200-replica reference** used for both systems.
Pooling is never done across conditions. With a shared reference the term cancels exactly in the
LATS1-minus-NDR2 contrast. Using each system's own MOB1A instead changes the multi-trajectory
contrast substantially - the two sub-ensembles differ by 8.46 ± 1.28 kcal/mol under the reference
condition - so results computed with a paired reference will not match the paper.
See `scripts/simulation/mt_pooled_reference.py`.

**Replica-level uncertainties.** Every reported ± is a replica-level standard error:
SE = s_replica / √N with N = 100, where s_replica is the standard deviation of the N replica
means. Frame-level errors, which treat the 400 snapshots within a replica as independent,
are roughly three times smaller and are *not* what the paper reports. Differences between
systems use the standard error of the paired per-replica differences, in which the shared
MOB1A reference cancels before the standard error is taken.

---

## 8. Archived separately

Too large for version control, deposited at the DOI given in the paper's Data Availability
Statement:

- AMBER intermediate files from the normal-mode runs (`_MMPBSA_*_nm.out.0`, ~3.5 GB) and the
  force dumps written during MM-PBSA (`reference.frc`, ~1.3 GB). These are byproducts of the
  calculations, not inputs to any reported result.
- Exploratory `cpptraj` analyses not used in the paper (~700 MB)

## 9. Not retained

- **Solvated MD trajectories** (~324 GB) are not archived. Water-stripped solute trajectories
  for the reference condition are available from the corresponding author on request.
- **Production restart files are absent for 280 of the 1,800 replicas**, produced before the
  campaign was extended to 100 replicas per system. Random seeds were recovered from the
  provenance records for 180 of those 280 and are listed in `data/seeds_assigned.csv`
  (columns: condition, species, replica, and the three stage seeds for heating, equilibration
  and production). For the remaining 100 they could not be recovered, so those replicas cannot
  be reproduced bit-for-bit. The missing restart files would be needed only to extend those
  particular runs further, which no reported result required; the ensemble statistics are
  unaffected. `scripts/simulation/check_seeds.sh` regenerates this audit.

- **Stereochemical screening.** `scripts/simulation/stereo_check.py` is the script used to screen
  starting models for D-configured Calpha centres, non-proline cis peptide bonds and heavy-atom
  contacts below 2.5 A.

---

## 10. Running these scripts elsewhere

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
in Section 1, `<protocol>` is `st` or `mt`, and `<NNN>` is the zero-padded replica index 001-100.

The SLURM submission scripts (`queue_*.sh`, `*_slurm.sh`, `rolling_submit*.sh`) encode the
partitions and resource limits of the cluster used and will need rewriting for any other
scheduler. The analysis scripts in `scripts/analysis/` have no scheduler dependency and need
only the path edits above.

---

## 11. Licence

Code under `scripts/` is MIT-licensed. Data under `data/`, `topologies/`, `md_systems/` and
`inputs/` is CC-BY-4.0. See `LICENSE`. If you use this material, cite the paper - `CITATION.cff`
has the details and GitHub will format it for you.

---

## 12. Contact

Ray Luo - rluo@uci.edu

# MM-PBSA of NDR2/LATS1-MOB1A complexes | data and analysis deposit

Supporting data, inputs and analysis code for:

> *Lessons Learned from Modeling Highly Charged Protein-Protein Interactions with MM-PBSA:
> A Case Study of NDR2/LATS1-MOB1A Complexes.*
> T. Lin, Z. Huang, W. Wang and R. Luo. Submitted to *Int. J. Mol. Sci.*

Release **v2.0.0**. The DOI of this version is given in the paper's Data Availability Statement.

Every value reported in the paper and its Supporting Information can be recomputed from the
files here. `scripts/analysis/reproduce_tables.py` does this for the main tables from the
contents of `energies/` alone and compares each result with the printed value:

```bash
python3 scripts/analysis/reproduce_tables.py      # run from the deposit root
```

---

## 1. What was simulated

Three **conditions**, each with six **species**, each with **100 independent replicas**:

| | |
|---|---|
| Conditions | `opc-salt` (reference: OPC water, added NaCl), `tip3p` (TIP3P water, added NaCl), `counterion` (OPC, neutralizing counterions only) |
| Species | `lats1_complex`, `lats1_kinase`, `lats1_mob1a`, `ndr2_complex`, `ndr2_kinase`, `ndr2_mob1a` |
| Replicas | 100 per species per condition (1,800 trajectories) |
| Protocol | 30,000 steps steepest descent + 30,000 conjugate gradient; 1 ns NVT heating/equilibration at 310.15 K; 50 ns NPT equilibration; **40 ns production**, coordinates every 10 ps (4,000 frames) |
| Integrator | 2 fs time step, 10 Å real-space cutoff, default PME |
| Force field | AMBER ff19SB; OPC or TIP3P water; Li-Merz 12-6 monovalent ions; ZAFF bonded model for the MOB1A zinc site |

### Added NaCl

In the two added-NaCl conditions each system was neutralized and Na⁺/Cl⁻ ion pairs were then
added so that every box matches the NaCl content of its own complex, expressed as ion pairs per
water molecule:

| Box | Added ion pairs |
|---|---|
| LATS1-MOB1A and NDR2-MOB1A complexes | 80 |
| LATS1 kinase segment | 52 |
| NDR2 kinase segment | 51 |
| MOB1A, LATS1-side preparation | 42 |
| MOB1A, NDR2-side preparation | 40 |

The counts are the same under OPC and TIP3P. `scripts/simulation/verify_salt.py` checks them
against the topologies. The Poisson-Boltzmann ionic strength used in MM-PBSA is a separate,
fixed setting (0.150 M in every condition).

### Additional sets

- **100 ns continuation** - replicas 1-20 of each complex in the reference condition, continued
  from 40 to 100 ns, with single-trajectory MM-PBSA over the 40-100 ns window
  (`extend_c1.sh`, `late_mmpbsa_parallel.sh`).
- **Solvent-retaining analysis** - 20 replicas of each unbound kinase segment under each water
  model, analyzed with solvent every 100 ps for hydration and ion distributions
  (`run_hydration.sh`).
- **PB ionic-strength sensitivity** - reference condition rescored at 0.150 and 0.300 M
  (`pb_istrng_sensitivity.sh`, `pb_sens_analyse.py`).

### Construct definitions

MOB1A UniProt Q9H8S9 residues 33-212; LATS1 residues 640-698 (from PDB **5BRK**); NDR2 residues
25-83 (from PDB **5XQZ chain B**). Five NDR2 glutamates (Glu40, Glu43, Glu55, Glu61, Glu62) are
modeled neutral so that both kinase segments carry net charge +6 - a deliberate control for the
between-system comparison, discussed in the Limitations section of the paper. All chain termini
are uncapped and charged.

---

## 2. Software

| | |
|---|---|
| MD | AMBER 2024, `pmemd.cuda` |
| MM-PBSA / normal modes | AmberTools 2024, `MMPBSA.py` |
| Trajectory analysis | `cpptraj` (AmberTools 2024) |
| Analysis / figures | Python 3.13.2, numpy 2.2.5, scipy 1.15.3, pandas 2.3.1, matplotlib 3.10.0 |

MM-PBSA settings: internal PBSA solver, `mbondi` radii, solute dielectric **20.0**, exterior
dielectric 80.0, ionic strength **0.150 M**, `inp=2` nonpolar model. Snapshots every 100 ps from
the 40 ns production run (400 per replica). Normal modes: 10 equally spaced snapshots per
replica, minimized to DRMS < 0.001 kcal/mol/Å; entropies computed by `MMPBSA.py` at its default
298.15 K and reported as −TΔS at 310.15 K.

---

## 3. Layout

```
scripts/simulation/        system construction, MD, MM-PBSA, normal modes, cpptraj analyses
scripts/analysis/          reproduce_tables.py, build_tidy_mt.py and the energy-partition scripts
scripts/analysis/figures/  the scripts that produce the figures and SI tables
data/                      seeds_assigned.csv
energies/                  per-frame energies and per-replica analysis outputs - the primary
                           numerical record, gzipped
topologies/                dry (solute-only) AMBER topologies used for MM-PBSA
md_systems/                solvated topologies and minimized coordinates, gzipped
inputs/                    tleap, minimization, heating, equilibration and production inputs,
                           with the ZAFF parameter files each tleap input loads
MANIFEST.sha256            checksums for every file in this deposit
```

---

## 4. Two topology sets

| Directory | Contents | Use |
|---|---|---|
| `topologies/` | solute only, no solvent | MM-PBSA and trajectory analysis on the stripped solute |
| `md_systems/` | solvated systems with ions | running the MD |

`md_systems/<condition>/<species>/` holds the solvated topology and `min2.rst`, the minimized
coordinates from which the heating stage of every replica starts (`heating_and_eq.sh` reads
`-c min2.rst`). `min1.rst`, where present, is the intermediate state between the two
minimization stages. Decompress before use:

```bash
gunzip -r md_systems/
```

Each of the 100 replicas differs only in its random seeds; `scripts/simulation/generate_replicas.sh`
shows how they were generated and `data/seeds_assigned.csv` lists them (columns: condition,
species, replica, and the heating, equilibration and production seeds; no header row).

---

## 5. The primary numerical record: `energies/`

| Path | Contents |
|---|---|
| `mmpbsa_perframe/mmpbsa_perframe_<condition>_<st\|mt>.csv.gz` | per-frame MM-PBSA energies, one table per condition and protocol |
| `mmpbsa_perframe/mmpbsa_perframe_opc-salt_st_late.csv.gz` | single-trajectory energies over 40-100 ns for the 20 continued replicas |
| `nmode/opc-salt/<system>/traj_<NNN>/per_frame.csv.gz` | normal-mode entropies, 10 snapshots per replica, reference condition |
| `physcv/physcv_<condition>_<bound\|unbound>.tar.gz` | interhelical angle and radius of gyration of the kinase segment, per frame, one CSV per replica |
| `terminal_core/terminal_core_opc-salt.tar.gz` | `cpptraj` bonded energies of the whole kinase segment and of its 49-residue core, bound and unbound |
| `interface/buried_area_<condition>.tar.gz`, `saltbridge_<condition>.tar.gz` | buried interface area and interfacial salt-bridge counts of the complexes, per frame |
| `hydration/<condition>/<species>/replica<NNN>/` | `cpptraj` outputs of the solvent-retaining analysis: RDFs, hydrogen bonds, water shells, SASA |
| `pb_sensitivity/pb_istrng_opc-salt.tar.gz` | raw `MMPBSA.py` per-frame output at PB ionic strengths 0.150 and 0.300 M |

The MM-PBSA tables are long-format, one row per frame, with columns
`system, water_model, salt_condition, protocol, delta_method, replica, frame, traj_frame,
time_ns, time_ns_abs`, then the difference terms (`BOND, ANGLE, DIHED, VDWAALS, EEL, 1-4 VDW,
1-4 EEL, EPB, ENPOLAR, EDISPER, TOTAL`, ...) and the same terms for each species as `cpx_*`
(complex), `rec_*` (MOB1A) and `lig_*` (kinase segment).

In the `_mt` tables `rec_*` and `lig_*` come from the independently simulated unbound MOB1A and
kinase segment, and the unprefixed difference columns pair replicas by index. **The paper's
multi-trajectory values are not those difference columns**: they use the pooled MOB1A
reference described in Section 7, formed from the `rec_*` columns.

---

## 6. Provenance: every figure and table

All scripts below are in `scripts/analysis/`.

| Output | Script | Reads |
|---|---|---|
| **Tables 1-3**, Table S1 differences, Table S3, Table S5, SI Section S5, interface values | `reproduce_tables.py` | `energies/` |
| **Figure 1**, Tables S1, S2 covalent decomposition | `figures/fig1_covalent_decomposition.py`, `figures/covalent_components.py`, `figures/tableS2_absolute.py` | MM-PBSA per-frame energies |
| **Figures 2, 3** free-energy surfaces | `figures/fig2fig3_fes.py`, `figures/physcv_data.py` | `energies/physcv/` |
| **Table S6** and the split-half comparison | `figures/tableS6_cv_metrics.py`, `figures/fes_splithalf_null.py` | `energies/physcv/` |
| Species partition of the water-model contrast (Section 2.2) | `analyze_watermodel_origin.py`, `species_absolute.py` | MM-PBSA per-frame energies |
| Radius-of-gyration correlations (Section 2.2) | `figures/rg_correlations.py` | `energies/physcv/`, MM-PBSA per-frame energies |
| **Figures S1, S2** time-series stability | `figures/figS1S2_convergence.py`, `figures/frames_cache.py` | MM-PBSA per-frame energies |
| **Figure S3** replica distributions | `figures/figS3_replica_distributions.py`, `figures/replica_values.py` | MM-PBSA per-frame energies |
| **Figure S4** bootstrap sign stability | `figures/figS4_bootstrap_forest.py` | MM-PBSA per-frame energies |
| **Figures S5, S6** and hydration-shell values | `figures/figS5S6_rdf.py`, `figures/hydration_shells.py` | `energies/hydration/` |
| **Table S3**, **Table S4** normal-mode entropy | `figures/tableS3_entropy.py`, `figures/tableS4_nmode_snapshots.py` | `energies/nmode/` |
| SI Section S5 terminal / core partition | `figures/terminal_core.py` | `energies/terminal_core/` |
| 40 ns versus 100 ns standard errors (Section 2.1) | `figures/se_40_vs_100ns.py` | `_st` and `_st_late` tables |
| PB ionic-strength sensitivity (Methods) | `../simulation/pb_sens_analyse.py` | `energies/pb_sensitivity/` |

`reproduce_tables.py` reads the deposit directly. The scripts under `figures/` are deposited as
executed and read the same data from the original working directories (Section 10).

---

## 7. Two conventions that matter when recomputing

**Pooled unbound MOB1A reference.** LATS1 and NDR2 bind the same MOB1A construct, so the
multi-trajectory values use one unbound MOB1A reference: within each condition the 100 replicas
of each of the two MOB1A preparations are pooled into a 200-replica reference used for both
systems. Pooling is never done across conditions. The shared reference cancels exactly in the
LATS1-minus-NDR2 difference. Pairing each complex with its own MOB1A preparation instead
changes the multi-trajectory difference substantially (the two preparations differ by
5.74 ± 1.01 kcal/mol under the reference condition), so results computed that way will not
match the paper.

**Replica-level uncertainties.** Every reported ± is a replica-level standard error:
SE = s_replica / √N, where s_replica is the standard deviation of the N replica means. Per-system
multi-trajectory standard errors include the pooled reference; standard errors of
LATS1-minus-NDR2 differences exclude it. The complex, kinase-segment and MOB1A simulations are
independent, so their standard errors are combined in quadrature; replicas are never paired by
index.

---

## 8. Not retained

- **MD trajectories** are not archived because of their size. They are retained and are
  available from the corresponding author on request.
- **Restart files and random seeds for 180 of the 1,800 replicas** were not retained. These are
  the earliest replicas, run before the campaign was extended to 100 replicas per system:
  20 each of the `lats1_complex` systems in all three conditions, of `ndr2_complex` in `tip3p`
  and `counterion`, and of the two `counterion` MOB1A systems, and 40 of `ndr2_complex` in
  `opc-salt`. Their seeds are absent from `data/seeds_assigned.csv`, so those replicas cannot be
  reproduced bit-for-bit. Their trajectories and energies are complete and are included in every
  reported result.
- AMBER intermediate files from the normal-mode and MM-PBSA runs are byproducts of the
  calculations, not inputs to any reported result, and are not included.

`scripts/simulation/stereo_check.py` is the script used to screen starting models for
D-configured Cα centers, non-proline cis peptide bonds and heavy-atom contacts below 2.5 Å.

---

## 9. Running these scripts elsewhere

The scripts are deposited **exactly as executed**, so they are a faithful record of what
produced the reported numbers rather than a polished package. 37 of the 57 scripts contain
absolute paths from the machine they ran on and will need editing before they run anywhere else.
The paths to change are:

| Path in the scripts | What it is |
|---|---|
| `/home3/jocelyn/projects/hippo/campaign_100x40ns` | campaign root: `runs/`, `inputs/`, `topologies/`, `mmpbsa/`, `nmode/`, `analysis/` |
| `/home3/jocelyn/projects/hippo/campaign_saltmatched` | root of the unbound added-NaCl systems, same layout, with condition directories named `opc-salt-matched` and `tip3p-matched` (these are the `opc-salt` and `tip3p` unbound systems of this deposit) |
| `/home3/jocelyn/projects/hippo/reviewer_comments` | analysis root: `analysis_*/data` hold the `cpptraj` outputs archived here under `energies/` |
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

The SLURM submission scripts (`*_slurm.sh`, `run_*.sh`) encode the
partitions and resource limits of the cluster used and will need rewriting for any other
scheduler. The analysis scripts in `scripts/analysis/` have no scheduler dependency;
`reproduce_tables.py` needs no edits and the others need only the path edits above.

---

## 10. License

The entire deposit - code, inputs, topologies, restart coordinates and derived energy data -
is released under the MIT License. See `LICENSE`. If you use this material, cite the paper - `CITATION.cff`
has the details and GitHub will format it for you.

---

## 11. Contact

Ray Luo - rluo@uci.edu

#!/usr/bin/env python
"""Psi4 BSSE sub-step: charged two-body pairs (Na+/Cl-, Na+/H2O, Cl-/H2O)
and the Na+/Cl-/H2O trimer (N=3), through the real psi4_step.BSSE
get_input()/analyze() pipeline (not hand-written psi4 scripts) -- the
validation orca_step's M3/N3 legs already did for ORCA, now for Psi4.

Same geometries as orca_step's validate_bsse_m3.py/validate_bsse_n3.py, and
the same level of theory (B3LYP-D3BJ/def2-TZVP -- psi4_step's Energy
defaults already match: method="Kohn-Sham (KS) density functional theory
(DFT)", functional="B3LYP ...", dispersion="d3bj"), so the two codes'
results are directly comparable at the same nominal level of theory.

Run inside the seamm-dev environment:
    python validate_psi4_bsse_m3.py [workdir]
"""
import json
import math
import sys
from pathlib import Path

import molsystem
import psi4_step
import seamm
import seamm_exec
from seamm.variables import Variables

seamm.flowchart_variables = Variables()

WORKDIR = Path(sys.argv[1]) if len(sys.argv) > 1 else Path("/tmp/psi4_bsse_m3_validate")
WORKDIR.mkdir(parents=True, exist_ok=True)

executor = seamm_exec.get_executor("local")
psi4_bin = "/Users/psaxe/miniconda3/envs/seamm-psi4/bin/psi4"
config = {"code": f"{psi4_bin} -n {{NTASKS}}", "installation": "local"}
ce = seamm_exec.computational_environment()


class _Parent:
    def __init__(self, basis):
        self.basis = basis
        self.model = None


class _Flowchart:
    def __init__(self, root_directory):
        self.root_directory = root_directory


def new_configuration(db_name):
    db = molsystem.SystemDB(filename=f"file:{db_name}?mode=memory&cache=shared")
    system = db.create_system(name=db_name)
    return system.create_configuration(name=db_name)


def run_bsse(label, configuration, fragment_atoms, fragment_charges, want_gradient="no"):
    node = psi4_step.BSSE()
    node._id = ("2", label)
    node.flowchart = _Flowchart(str(WORKDIR))
    node.parent = _Parent(basis="def2-TZVP")
    node.get_system_configuration = lambda arg: (None, configuration)

    node.parameters["fragments"].value = "specified"
    node.parameters["fragment atoms"].value = fragment_atoms
    node.parameters["fragment charges"].value = fragment_charges
    node.parameters["compute gradient"].value = want_gradient

    bsse_text = node.get_input()

    run_dir = WORKDIR / f"run_{label}"
    run_dir.mkdir(parents=True, exist_ok=True)
    script = "\n".join(
        [
            "import json",
            "from pathlib import Path",
            "",
            f"set basis {node.parent.basis}",
            bsse_text,
        ]
    )
    result = executor.run(
        cmd=["{code}"],
        config=config,
        directory=str(run_dir),
        files={"input.dat": script},
        return_files=["output.dat", "*.json"],
        in_situ=True,
        shell=True,
        ce=ce,
    )
    if not result:
        raise RuntimeError(f"[{label}] Psi4 run failed.")
    output_text = (run_dir / "output.dat").read_text()
    if "Psi4 exiting successfully" not in output_text:
        print(output_text[-3000:])
        raise RuntimeError(f"[{label}] Psi4 did not exit successfully.")

    node.analyze()

    with (Path(node.directory) / "bsse.json").open() as fd:
        return json.load(fd)


# ---------------------------------------------------------------------
# Na+...Cl-, R = 2.44 A (the M3-validated near-equilibrium point).
# ---------------------------------------------------------------------
print("=" * 70)
print("Na+ ... Cl-  (R = 2.44 A)")
print("=" * 70)
configuration = new_configuration("nacl")
configuration.atoms.append(symbol=["Na", "Cl"], x=[0.0, 0.0], y=[0.0, 0.0], z=[0.0, 2.44])
configuration.charge = 0
configuration.spin_multiplicity = 1
data = run_bsse("nacl", configuration, "1; 2", "1, -1")
kcal = data["interaction energy"] * 627.509474
print(f"CP interaction energy = {kcal:.3f} kcal/mol")
print("ORCA M3 (B3LYP-D3BJ/def2-TZVP, same R): -136.289 kcal/mol")
print("Approx. literature reference: ~ -133 kcal/mol at R_e ~ 2.36 A")

# ---------------------------------------------------------------------
# Na+...H2O and Cl-...H2O -- same geometries as orca_step's M3 script.
# ---------------------------------------------------------------------
r_oh = 0.9572
half_angle = math.radians(104.5 / 2.0)
h1 = (r_oh * math.sin(half_angle), 0.0, r_oh * math.cos(half_angle))
h2 = (-r_oh * math.sin(half_angle), 0.0, r_oh * math.cos(half_angle))
o = (0.0, 0.0, 0.0)

print()
print("=" * 70)
print("Na+ ... H2O  (Na+ on the O lone-pair side, Na-O = 2.30 A)")
print("=" * 70)
na = (0.0, 0.0, -2.30)
configuration = new_configuration("na_water")
configuration.atoms.append(
    symbol=["O", "H", "H", "Na"],
    x=[o[0], h1[0], h2[0], na[0]],
    y=[o[1], h1[1], h2[1], na[1]],
    z=[o[2], h1[2], h2[2], na[2]],
)
configuration.charge = 1
configuration.spin_multiplicity = 1
data = run_bsse("na_water", configuration, "1-3; 4", "0, 1")
kcal = data["interaction energy"] * 627.509474
print(f"CP interaction energy = {kcal:.3f} kcal/mol")
print("ORCA M3 (B3LYP-D3BJ/def2-TZVP, same geometry): -26.110 kcal/mol")
print("Approx. literature reference: ~ -24 kcal/mol")

print()
print("=" * 70)
print("Cl- ... H2O  (Cl- along the O-H1 bond extension, H...Cl = 2.20 A)")
print("=" * 70)
oh1_unit = (h1[0] / r_oh, h1[1] / r_oh, h1[2] / r_oh)
cl_distance_from_o = r_oh + 2.20
cl = tuple(oh1_unit[i] * cl_distance_from_o for i in range(3))
configuration = new_configuration("cl_water")
configuration.atoms.append(
    symbol=["O", "H", "H", "Cl"],
    x=[o[0], h1[0], h2[0], cl[0]],
    y=[o[1], h1[1], h2[1], cl[1]],
    z=[o[2], h1[2], h2[2], cl[2]],
)
configuration.charge = -1
configuration.spin_multiplicity = 1
data = run_bsse("cl_water", configuration, "1-3; 4", "0, -1")
kcal = data["interaction energy"] * 627.509474
print(f"CP interaction energy = {kcal:.3f} kcal/mol")
print("ORCA M3 (B3LYP-D3BJ/def2-TZVP, same geometry): -15.638 kcal/mol")
print("Approx. literature reference: ~ -13 kcal/mol")

# ---------------------------------------------------------------------
# N=3: Na+...Cl-...H2O trimer -- same geometry as orca_step's N3 script.
# ---------------------------------------------------------------------
print()
print("=" * 70)
print("Na+ ... Cl- ... H2O  (N=3, contact ion pair + water on Na+)")
print("=" * 70)
na3 = (0.0, 0.0, 0.0)
cl3 = (0.0, 0.0, 2.44)
approach_angle = math.radians(130.0)
approach_dir = (math.sin(approach_angle), 0.0, math.cos(approach_angle))
na_o_distance = 2.35
o3 = tuple(na3[i] + na_o_distance * approach_dir[i] for i in range(3))
local_x = (0.0, 1.0, 0.0)
local_z = approach_dir


def offset(lx, lz):
    return tuple(o3[i] + lx * local_x[i] + lz * local_z[i] for i in range(3))


h1_3 = offset(r_oh * math.sin(half_angle), r_oh * math.cos(half_angle))
h2_3 = offset(-r_oh * math.sin(half_angle), r_oh * math.cos(half_angle))

configuration = new_configuration("na_cl_water")
configuration.atoms.append(
    symbol=["Na", "Cl", "O", "H", "H"],
    x=[na3[0], cl3[0], o3[0], h1_3[0], h2_3[0]],
    y=[na3[1], cl3[1], o3[1], h1_3[1], h2_3[1]],
    z=[na3[2], cl3[2], o3[2], h1_3[2], h2_3[2]],
)
configuration.charge = 0
configuration.spin_multiplicity = 1
data = run_bsse("na_cl_water", configuration, "1; 2; 3-5", "1, -1, 0")
print(f"uncorrected energy       = {data['uncorrected energy']:.8f} Eh")
print(f"corrected energy         = {data['energy']:.8f} Eh")
print(f"bsse correction          = {data['bsse correction'] * 627.509474:.3f} kcal/mol")
kcal = data["interaction energy"] * 627.509474
print(f"CP interaction energy    = {kcal:.3f} kcal/mol")
print("ORCA N3 (B3LYP-D3BJ/def2-TZVP, same geometry): -150.204 kcal/mol")
print(
    "Sanity range: should be deeper than the Na+...Cl- pair alone but less "
    "than the naive pairwise sum (cooperative-saturation non-additivity)."
)

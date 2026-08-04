#!/usr/bin/env python
"""Real-Psi4 validation of the new psi4_step BSSE sub-step, on the water
dimer -- checks that BSSE.get_input()/analyze() actually work end to end
(the shared-script JSON round trip, the molecule-block builder, the
bsse_type=['cp','nocp'] parsing), and sanity-checks the numbers against the
already-validated ORCA M2 result for the same system.

This bypasses psi4_step.Psi4's normal subflowchart orchestration (no
Initialization node, no real Flowchart graph) since only BSSE's own
get_input()/analyze() are being tested here -- it hand-builds the same kind
of shared script Psi4.run() would (a `set basis` line + BSSE's own text) and
runs it directly through the real local executor.

Run inside the seamm-dev environment (needs the seamm-psi4 conda env's own
`psi4` binary):
    python validate_psi4_bsse.py [workdir]
"""
import json
import sys
from pathlib import Path

import molsystem
import psi4_step
import seamm
import seamm_exec
from seamm.variables import Variables

seamm.flowchart_variables = Variables()

WORKDIR = Path(sys.argv[1]) if len(sys.argv) > 1 else Path("/tmp/psi4_bsse_validate")
WORKDIR.mkdir(parents=True, exist_ok=True)
SDF = Path("/Users/psaxe/structures/dimers/H2O-H2O.sdf")

# --- A real molsystem Configuration for the water dimer ---
db = molsystem.SystemDB(filename="file:psi4_bsse_validate?mode=memory&cache=shared")
system = db.create_system(name="water-dimer")
configuration = system.create_configuration(name="water-dimer")
configuration.from_sdf(SDF)
configuration.charge = 0
configuration.spin_multiplicity = 1
print(f"Loaded {configuration.n_atoms} atoms: {list(configuration.atoms.symbols)}")


class _Parent:
    """Stand-in for the main Psi4 node -- BSSE.get_input() reads
    .basis (set by the Initialization sub-step, normally) and writes .model."""

    def __init__(self, basis):
        self.basis = basis
        self.model = None


class _Flowchart:
    def __init__(self, root_directory):
        self.root_directory = root_directory


node = psi4_step.BSSE()
node._id = ("2", "1")  # mimics "sub-step 1 of main Psi4 node 2"
node.flowchart = _Flowchart(str(WORKDIR))
node.parent = _Parent(basis="def2-svp")
node.get_system_configuration = lambda arg: (None, configuration)

# HF/def2-SVP: fast, matches the M2 ORCA regression's level of theory, so
# the two are directly comparable.
node.parameters["level"].value = "recommended"
node.parameters["method"].value = "Hartree-Fock (HF) self consistent field (SCF)"
node.parameters["fragments"].value = "auto (molecules)"
node.parameters["compute gradient"].value = "yes"

bsse_text = node.get_input()
print("\n=== Generated Psi4 input (BSSE.get_input()) ===")
print(bsse_text)

run_dir = WORKDIR / "run"
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

executor = seamm_exec.get_executor("local")
psi4_bin = "/Users/psaxe/miniconda3/envs/seamm-psi4/bin/psi4"
config = {"code": f"{psi4_bin} -n {{NTASKS}}", "installation": "local"}
ce = seamm_exec.computational_environment()

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
    raise RuntimeError("Psi4 run failed (executor returned falsy result).")
output_text = (run_dir / "output.dat").read_text()
if "Psi4 exiting successfully" not in output_text:
    print(output_text[-4000:])
    raise RuntimeError("Psi4 did not exit successfully; see output above.")
print("\nPsi4 exited successfully.")

node.analyze()
print("\n=== analyze() ran without error ===")

json_file = Path(node.directory) / "bsse.json"
with json_file.open() as fd:
    raw = json.load(fd)

print("\n=== Results (raw, from bsse.json) ===")
for key, value in raw.items():
    if key != "gradients":
        print(f"{key:30s} = {value}")
print(f"gradient rows = {len(raw['gradients'])}")
print(f"gradient[0]   = {raw['gradients'][0]}")

print("\n=== Sanity check vs. the ORCA M2 water-dimer regression ===")
print(
    "(different method -- HF/def2-SVP here for both, should be in the same "
    "ballpark as the ORCA HF/def2-SVP M2 numbers: uncorrected -151.8897296, "
    "corrected -151.8875402 Eh)"
)

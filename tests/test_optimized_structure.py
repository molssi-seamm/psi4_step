# -*- coding: utf-8 -*-
"""Optimization stores Psi4's optimized coordinates (bugfix 2026.11.8.1: it looked for
the run's final_structure.json in its own directory, never found it, and silently
kept the starting coordinates under an 'optimized' name)."""

import json
from types import SimpleNamespace

import numpy as np
import pytest

import psi4_step
from molsystem import RMSD, SystemDB


@pytest.fixture
def water():
    db = SystemDB(filename="file:psi4_opt_test?mode=memory&cache=shared")
    system = db.create_system()
    configuration = system.create_configuration()
    configuration.from_smiles("O")
    yield configuration
    db.close()


def _rotation(a, b, c):
    def rx(t):
        return np.array(
            [[1, 0, 0], [0, np.cos(t), -np.sin(t)], [0, np.sin(t), np.cos(t)]]
        )

    def rz(t):
        return np.array(
            [[np.cos(t), -np.sin(t), 0], [np.sin(t), np.cos(t), 0], [0, 0, 1]]
        )

    return rz(a) @ rx(b) @ rz(c)


def _distances(xyz):
    xyz = np.asarray(xyz)
    return np.linalg.norm(xyz[:, None, :] - xyz[None, :, :], axis=-1)


def _write(path, xyz, units="Angstrom"):
    path.write_text(
        json.dumps({"geom": np.asarray(xyz).ravel().tolist(), "units": units})
    )


def test_reads_angstrom_and_bohr_and_fails_loudly(tmp_path):
    node = psi4_step.Optimization()
    xyz = np.arange(9.0).reshape(3, 3)
    _write(tmp_path / "final_structure.json", xyz)
    assert np.allclose(node.optimized_coordinates(tmp_path), xyz)
    _write(tmp_path / "final_structure.json", xyz, units="Bohr")
    assert np.allclose(node.optimized_coordinates(tmp_path), xyz * node.BOHR)
    (tmp_path / "final_structure.json").unlink()
    with pytest.raises(RuntimeError, match="optimized structure"):
        node.optimized_coordinates(tmp_path)


def test_the_input_writes_this_steps_structure(monkeypatch):
    import seamm

    monkeypatch.setattr(seamm, "flowchart_variables", seamm.Variables())
    flowchart = seamm.Flowchart()
    flowchart.root_directory = "/tmp"
    node = psi4_step.Optimization(flowchart=flowchart)
    node._id = ("2", "3")
    node.parent = SimpleNamespace(basis="cc-pVDZ", model=None, extended_model=None)
    node.get_system_configuration = lambda *a, **k: (
        None,
        SimpleNamespace(spin_multiplicity=1, n_atoms=3),
    )
    text = node.get_input()
    assert 'with open("@3+final_structure.json", "w") as fd:' in text
    assert text.index("optimize(") < text.index("@3+final_structure.json")


def test_stored_coordinates_are_the_optimized_ones(water, tmp_path):
    """Psi4's structure, in its own frame, ends up in the configuration (aligned with
    the start) with the optimized internal geometry to 1e-8 Å."""
    initial = water.to_RDKMol()
    start = np.array(water.coordinates)
    # The "optimized" structure: O-H bonds stretched, rotated and moved
    optimized = start.copy()
    optimized[1:] = start[0] + 1.03 * (start[1:] - start[0])
    optimized = optimized @ _rotation(0.7, -1.1, 2.3).T + np.array([1.5, -2.0, 0.3])
    _write(tmp_path / "final_structure.json", optimized)

    node = psi4_step.Optimization()
    final = water.to_RDKMol()
    final.GetConformer(0).SetPositions(node.optimized_coordinates(tmp_path))
    result = RMSD(final, initial, symmetry=True, include_h=True, align=True)
    water.coordinates_from_RDKMol(final)

    stored = np.array(water.coordinates)
    assert np.abs(_distances(stored) - _distances(optimized)).max() < 1e-8
    assert result["RMSD"] > 0.01  # it did move
    # ... and is aligned with the start, not left in Psi4's frame
    assert np.abs(stored.mean(axis=0) - start.mean(axis=0)).max() < 0.05

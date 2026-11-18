#!/usr/bin/env python
"""Write psi4_step/data/basis/cc-pwCVTZ-ORCA.basis, a Psi4 basis file holding the
cc-pwCVTZ basis set exactly as ORCA 6.1.1 has it internally, with explicit
density-fitting sets.

Usage::

    python make_cc_pwcvtz_orca.py <directory of ORCA PrintBasis outputs> <output file>

The directory holds one ORCA output per element (<El>.out), from the input
``! UHF cc-pwCVTZ PrintBasis NoIter``; the section "BASIS SET IN INPUT FORMAT" is
read. Elements without an ORCA export are taken from the Basis Set Exchange
(cc-pwCVTZ, or cc-pVTZ for H and He) after checking, for every element that has
both, that ORCA's functions and the BSE's are the same. The fitting sets come from
the BSE: cc-pwCVTZ-RIFIT for the correlation (cc-pVTZ-RIFIT for H and He, which it
lacks) and def2-universal-JKFIT for the SCF.
"""

import re
import sys
from pathlib import Path

import basis_set_exchange as bse
import numpy as np

AM = "SPDFGHI"


def orca_shells(path):
    """The contracted shells of an ORCA PrintBasis output: [(am, exps, coefs)]."""
    text = Path(path).read_text()
    block = text.split("BASIS SET IN INPUT FORMAT", 1)[1]
    block = block.split("NewGTO", 1)[1].split("end;", 1)[0]
    lines = block.splitlines()[1:]
    shells = []
    i = 0
    while i < len(lines):
        fields = lines[i].split()
        i += 1
        if len(fields) != 2 or fields[0] not in AM:
            continue
        n = int(fields[1])
        exps, coefs = [], []
        for line in lines[i : i + n]:
            _, e, c = line.split()
            exps.append(float(e))
            coefs.append(float(c))
        i += n
        shells.append((AM.index(fields[0]), exps, coefs))
    return shells


def bse_shells(name, Z):
    """The contracted shells of a BSE basis for element Z, one per column."""
    data = bse.get_basis(name, elements=[Z], make_general=False)
    shells = []
    for shell in data["elements"][str(Z)]["electron_shells"]:
        exps = [float(e) for e in shell["exponents"]]
        for am, column in zip(
            shell["angular_momentum"] * len(shell["coefficients"]),
            shell["coefficients"],
        ):
            c = [float(x) for x in column]
            keep = [k for k, x in enumerate(c) if x != 0.0]
            shells.append((am, [exps[k] for k in keep], [c[k] for k in keep]))
    return shells


def _span(shells, am):
    """The contractions of angular momentum am as columns over all its exponents."""
    chosen = [(exps, coefs) for a, exps, coefs in shells if a == am]
    exps = sorted({round(e, 6) for e, _ in [(x, 0) for s in chosen for x in s[0]]})
    M = np.zeros((len(exps), len(chosen)))
    for j, (es, cs) in enumerate(chosen):
        for e, c in zip(es, cs):
            M[exps.index(round(e, 6)), j] += c
    return exps, M


def same_functions(a, b):
    """Whether two shell lists span the same functions: per angular momentum, the
    same exponents, the same number of contractions, and the same space. (ORCA
    removes from a contraction the primitives that are also separate shells, which
    changes the coefficients but not the space.)"""
    ams = {s[0] for s in a}
    if ams != {s[0] for s in b}:
        return False
    for am in ams:
        ea, A = _span(a, am)
        eb, B = _span(b, am)
        if ea != eb or A.shape != B.shape:
            return False
        A = A / np.linalg.norm(A, axis=0)
        B = B / np.linalg.norm(B, axis=0)
        r = np.linalg.matrix_rank
        if not (r(A, 1e-8) == r(B, 1e-8) == r(np.hstack([A, B]), 1e-8) == A.shape[1]):
            return False
    return True


def gbs_element(symbol, shells):
    lines = [f"{symbol.upper()}     0"]
    for am, exps, coefs in shells:
        lines.append(f"{AM[am]}   {len(exps)}   1.00")
        for e, c in zip(exps, coefs):
            lines.append(f"  {e:20.10f}  {c:16.10f}")
    lines.append("****")
    return lines


def gbs_section(name, elements):
    """A Psi4 [name] section: elements = {symbol: shells}."""
    lines = [f"[ {name} ]", "spherical", "****"]
    for symbol, shells in elements.items():
        lines.extend(gbs_element(symbol, shells))
    return lines


def main(directory, output):
    directory = Path(directory)
    orca = {p.stem: orca_shells(p) for p in sorted(directory.glob("*.out"))}

    pwcv = bse.get_metadata()["cc-pwcvtz"]
    pwcv_Z = [int(z) for z in pwcv["versions"][pwcv["latest_version"]]["elements"]]

    # Check ORCA against the BSE wherever both have the element
    report = []
    for symbol, shells in orca.items():
        Z = bse.lut.element_Z_from_sym(symbol)
        name = "cc-pwcvtz" if Z in pwcv_Z else ("cc-pvtz" if Z <= 2 else None)
        if name is None:
            report.append(f"{symbol}: ORCA only ({len(shells)} shells)")
            continue
        if not same_functions(shells, bse_shells(name, Z)):
            raise SystemExit(f"{symbol}: ORCA's {name} differs from the BSE's")
        report.append(f"{symbol}: ORCA = BSE {name}")

    # Elements: H, He; B-Ne; Al-Ar; and the ORCA-only ones (Ga-Kr exported)
    Zs = sorted({1, 2} | {z for z in pwcv_Z if z <= 18} | set(
        bse.lut.element_Z_from_sym(s) for s in orca
    ))
    orbital, ri, jk = {}, {}, {}
    sources = []
    for Z in Zs:
        symbol = bse.lut.element_sym_from_Z(Z, normalize=True)
        if symbol in orca:
            orbital[symbol] = orca[symbol]
            sources.append(f"{symbol} ORCA")
        else:
            name = "cc-pvtz" if Z <= 2 else "cc-pwcvtz"
            orbital[symbol] = bse_shells(name, Z)
            sources.append(f"{symbol} BSE {name}")
        rifit = "cc-pvtz-rifit" if Z <= 2 else "cc-pwcvtz-rifit"
        ri[symbol] = (rifit, bse_shells(rifit, Z))
        jk[symbol] = bse_shells("def2-universal-jkfit", Z)

    out = [
        "# cc-pwCVTZ exactly as ORCA 6.1.1 has it internally (cc-pVTZ on H and He),",
        "# with explicit density-fitting sets from the Basis Set Exchange "
        f"{bse.version()}:",
        "# cc-pwCVTZ-RIFIT (cc-pVTZ-RIFIT on H and He) for the correlation and",
        "# def2-universal-JKFIT for the SCF. Made by",
        "# devtools/basis/make_cc_pwcvtz_orca.py; see data/basis/README.md.",
        "# Orbital basis sources: " + ", ".join(sources),
        "",
        "basis {",
        "   assign orca-cc-pwcvtz",
    ]
    out.extend(gbs_section("orca-cc-pwcvtz", orbital))
    out.append("}")
    out.append("")
    out.append("df_basis_cc {")
    out.append("   assign bse-cc-pwcvtz-rifit")
    for symbol, (name, _) in ri.items():
        if name == "cc-pvtz-rifit":
            out.append(f"   assign {symbol.upper()} bse-cc-pvtz-rifit")
    out.extend(
        gbs_section(
            "bse-cc-pwcvtz-rifit",
            {s: sh for s, (n, sh) in ri.items() if n == "cc-pwcvtz-rifit"},
        )
    )
    out.extend(
        gbs_section(
            "bse-cc-pvtz-rifit",
            {s: sh for s, (n, sh) in ri.items() if n == "cc-pvtz-rifit"},
        )
    )
    out.append("}")
    out.append("")
    out.append("df_basis_scf {")
    out.append("   assign bse-def2-universal-jkfit")
    out.extend(gbs_section("bse-def2-universal-jkfit", jk))
    out.append("}")
    Path(output).write_text("\n".join(out) + "\n")
    print("\n".join(report))
    print(f"Wrote {output}: elements {', '.join(orbital)}")


if __name__ == "__main__":
    main(*sys.argv[1:])

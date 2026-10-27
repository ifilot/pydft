"""Regenerate independent XC references; requires optional PySCF installation.

Run from the repository root with PYTHONPATH=. python tests/generate_xc_reference.py.
PySCF supplies Libxc point values and an independent SCF/integral implementation.
PyDFT supplies only the molecular geometries, in Bohr.
"""

import json
from pathlib import Path

import numpy as np
import pyscf
from pyscf import dft, gto, lib
from pyscf.dft import libxc

from pydft import MoleculeBuilder


def main():
    lib.num_threads(1)
    rho = np.logspace(-8, 3, 12)
    reference = {
        "provenance": {
            "pyscf": pyscf.__version__,
            "libxc": libxc.__version__,
            "spin": 0,
            "basis": "sto-3g",
            "grid_level": 6,
            "scf_tolerance": 1e-11,
            "energy_unit": "Hartree (including nuclear repulsion)",
            "geometry_unit": "Bohr",
            "svwn3_convention": "Slater + VWN-RPA, spin-unpolarized",
        },
        "rho": rho.tolist(),
        "functionals": {},
        "molecules": {},
    }
    for name, correlation in [("svwn3", "LDA_C_VWN_RPA"), ("svwn5", "LDA_C_VWN")]:
        ex, vx, _, _ = libxc.eval_xc("LDA_X,", rho, spin=0, deriv=1)
        ec, vc, _, _ = libxc.eval_xc("," + correlation, rho, spin=0, deriv=1)
        reference["functionals"][name] = {
            "libxc_correlation": correlation,
            "ex": ex.tolist(), "vx": vx[0].tolist(),
            "ec": ec.tolist(), "vc": vc[0].tolist(),
        }

    for name in ["He", "H2", "CO"]:
        atoms = [(symbol, position.tolist()) for symbol, position in
                 MoleculeBuilder().from_name(name)]
        mol = gto.M(atom=atoms, unit="Bohr", basis="sto-3g", verbose=0)
        energies = {}
        for functional, data in reference["functionals"].items():
            calculation = dft.RKS(mol)
            calculation.xc = "LDA_X," + data["libxc_correlation"]
            calculation.grids.level = 6
            calculation.conv_tol = 1e-11
            energies[functional] = float(calculation.kernel())
            if not calculation.converged:
                raise RuntimeError(f"Reference SCF failed: {name}, {functional}")
        reference["molecules"][name] = {"atoms": atoms, "energies": energies}

    destination = Path(__file__).parent / "data" / "xc_reference.json"
    destination.parent.mkdir(exist_ok=True)
    destination.write_text(json.dumps(reference, indent=2) + "\n")
    print(f"Wrote {destination}")


if __name__ == "__main__":
    main()

"""Check against committed Libxc/PySCF values without requiring PySCF in CI."""

import json
from pathlib import Path

import numpy as np
import pytest

from pydft import DFT, Molecule
from pydft.xcfunctionals import Functionals


REFERENCE = json.loads((Path(__file__).parent / "data" / "xc_reference.json").read_text())


@pytest.mark.parametrize("functional", ["svwn3", "svwn5"])
def test_lda_against_libxc(functional):
    rho = np.array(REFERENCE["rho"])
    expected = REFERENCE["functionals"][functional]
    implementation = Functionals(functional)
    for method, energy, potential in [(implementation.calc_x, "ex", "vx"),
                                      (implementation.calc_c, "ec", "vc")]:
        actual_energy, actual_potential, sigma = method(rho, deriv_sigma=True)
        np.testing.assert_allclose(actual_energy, expected[energy], rtol=1e-10, atol=1e-12)
        np.testing.assert_allclose(actual_potential, expected[potential], rtol=1e-10, atol=1e-12)
        np.testing.assert_array_equal(sigma, np.zeros_like(rho))


@pytest.mark.parametrize("functional", ["svwn3", "svwn5"])
def test_correlation_potential_is_energy_derivative(functional):
    rho = np.array(REFERENCE["rho"])
    step = rho * 1e-5
    implementation = Functionals(functional)
    _, potential = implementation.calc_c(rho)
    plus = (rho + step) * implementation.calc_c(rho + step)[0]
    minus = (rho - step) * implementation.calc_c(rho - step)[0]
    np.testing.assert_allclose(potential, (plus - minus) / (2 * step),
                               rtol=1e-8, atol=1e-11)


def test_svwn_alias_selects_svwn5():
    rho = np.array(REFERENCE["rho"])
    alias = Functionals("svwn")
    canonical = Functionals("svwn5")
    assert not alias.is_gga()
    np.testing.assert_array_equal(alias.calc_x(rho), canonical.calc_x(rho))
    np.testing.assert_array_equal(alias.calc_c(rho), canonical.calc_c(rho))


@pytest.mark.parametrize("functional", ["svwn3", "svwn5"])
@pytest.mark.parametrize("name", ["He", "H2", "CO"])
def test_molecular_energy_against_pyscf(name, functional):
    expected = REFERENCE["molecules"][name]
    molecule = Molecule()
    for symbol, position in expected["atoms"]:
        molecule.add_atom(symbol, *position, unit="bohr")
    calculation = DFT(molecule, basis="sto3g", functional=functional,
                      nshells=128, nangpts=302)
    result = calculation.scf(tol=1e-7)
    # Use a finer grid than the defaults to resolve the numerical Hartree term.
    # PyQInt and PySCF also use slightly different rounded STO-3G coefficients.
    tolerance = 2e-4 if name == "CO" else 5e-6
    assert result["energy"] == pytest.approx(expected["energies"][functional], abs=tolerance)

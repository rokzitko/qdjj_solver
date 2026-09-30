"""Finite-Hamiltonian DMRG; requires the optional ``qdjj-solver[dmrg]`` extra."""
from ..common import (FermionOperator, Hamiltonian, Impurity, Sector, annihilate,
                      cosh_grid, create, make_model, number, reference_model)


def __getattr__(name):
    # Importing common models or the legacy chain transform must not require
    # TeNPy. Only requesting the numerical backend loads the optional extra.
    if name in ("PreparedModel", "prepare", "Eigenstates", "SolverOptions", "solve"):
        from importlib import import_module
        module = ".mpo" if name in ("PreparedModel", "prepare") else ".solver"
        try:
            return getattr(import_module(module, __name__), name)
        except ModuleNotFoundError as error:
            if error.name == "tenpy":
                raise ImportError("Install qdjj-solver[dmrg] to use the DMRG backend") from error
            raise
    raise AttributeError(name)

__all__ = ["FermionOperator", "Hamiltonian", "Impurity", "Sector", "annihilate", "create",
           "number", "cosh_grid", "make_model", "reference_model", "PreparedModel", "prepare",
           "Eigenstates", "SolverOptions", "solve"]

# -*- coding: utf-8 -*-
"""Global control parameters for Psi4"""

import logging

from psi4_step import methods, dft_functionals
import seamm

logger = logging.getLogger(__name__)

# The methods that are DFT, which take a functional.
DFT_METHODS = [name for name, record in methods.items() if record["method"] == "dft"]


#: How a step that needs gradients or Hessians gets them (Optimization and
#: Thermochemistry): analytic gradients where Psi4 has them, else finite
#: differences of energies; or finite differences of energies regardless.
ANALYTIC = "analytic where available"
FINITE_DIFFERENCES = "finite differences of energies"
DERIVATIVES = {
    "default": ANALYTIC,
    "kind": "enumeration",
    "default_units": "",
    "enumeration": (ANALYTIC, FINITE_DIFFERENCES),
    "format_string": "s",
    "description": "Derivatives:",
    "help_text": (
        "How the gradients (and so the Hessian) are found: analytically where "
        "Psi4 can for this method, or by finite differences of energies, which "
        "costs many more calculations but needs only energies -- a check on, or a "
        "way around, an analytic gradient."
    ),
}


class EnergyParameters(seamm.Parameters):
    """The control parameters for the energy."""

    parameters = {
        "level": {
            "default": "recommended",
            "kind": "string",
            "format_string": "s",
            "enumeration": ("recommended", "advanced"),
            "description": "The level of disclosure in the interface",
            "help_text": (
                "How much detail to show in the GUI. Currently 'recommended' "
                "or 'advanced', which shows everything."
            ),
        },
        "method": {
            "applies_when": {"level": "recommended"},
            "default": "Kohn-Sham (KS) density functional theory (DFT)",
            "kind": "enumeration",
            "default_units": "",
            "enumeration": [x for x in methods if methods[x]["level"] == "normal"],
            "format_string": "s",
            "description": "Method:",
            "help_text": ("The computational method to use."),
        },
        "advanced_method": {
            "applies_when": {"level": "advanced"},
            "default": "Kohn-Sham (KS) density functional theory (DFT)",
            "kind": "enumeration",
            "default_units": "",
            "enumeration": [x for x in methods],
            "format_string": "s",
            "description": "Method:",
            "help_text": ("The computational method to use."),
        },
        "functional": {
            "applies_when": {"method": DFT_METHODS},
            "default": "B3LYP Hyb-GGA Exchange-Correlation Functional",
            "kind": "enumeration",
            "default_units": "",
            "enumeration": [
                x for x in dft_functionals if dft_functionals[x]["level"] == "normal"
            ],
            "format_string": "s",
            "description": "DFT Functional:",
            "help_text": ("The exchange-correlation functional to use."),
        },
        "advanced_functional": {
            "applies_when": {"advanced_method": DFT_METHODS},
            "default": "B3LYP Hyb-GGA Exchange-Correlation Functional",
            "kind": "enumeration",
            "default_units": "",
            "enumeration": [x for x in dft_functionals],
            "format_string": "s",
            "description": "DFT Functional:",
            "help_text": ("The exchange-correlation functional to use."),
        },
        "dispersion": {
            "default": "d3bj",
            "kind": "enumeration",
            "default_units": "",
            "enumeration": ["none", "d3bj", "d3mbj", "nl"],
            "format_string": "s",
            "description": "Dispersion correction:",
            "help_text": ("The dispersion correction to use."),
        },
        "spin-restricted": {
            "default": "default",
            "kind": "enumeration",
            "default_units": "",
            "enumeration": ("default", "yes", "no"),
            "format_string": "s",
            "description": "Spin-restricted:",
            "help_text": (
                "Whether to restrict the spin (RHF, ROHF, RKS) or not "
                "(UHF, UKS)."
                " Default is restricted for singlets, unrestricted otherwise."
            ),
        },
        "use damping": {
            "default": "no",
            "kind": "boolean",
            "default_units": "",
            "enumeration": ("yes", "no"),
            "format_string": "s",
            "description": "Damp the iterations:",
            "help_text": (
                "Whether to damp the iterations using a fraction of the previous "
                "density to damp oscillations."
            ),
        },
        "damping percentage": {
            "applies_when": {"use damping": "yes"},
            "default": 20.0,
            "kind": "float",
            "default_units": "",
            "enumeration": tuple(),
            "format_string": "s",
            "description": "Percent damping:",
            "help_text": "Percent of previous density to use to damp oscillations.",
        },
        "damping convergence": {
            "applies_when": {"use damping": "yes"},
            "default": 0.0,
            "kind": "float",
            "default_units": "",
            "enumeration": tuple(),
            "format_string": "s",
            "description": "Damp to convergence:",
            "help_text": "Convergence level to stop damping. 0 = always damp.",
        },
        "use level shift": {
            "default": "no",
            "kind": "boolean",
            "default_units": "",
            "enumeration": ("yes", "no"),
            "format_string": "s",
            "description": "Use a level shift:",
            "help_text": ("Whether to use a level shift to help convergence."),
        },
        "level shift": {
            "applies_when": {"use level shift": "yes"},
            "default": 5.0,
            "kind": "float",
            "default_units": "E_h",
            "enumeration": tuple(),
            "format_string": "s",
            "description": "Level shift:",
            "help_text": "The amount to shift the occupied orbitals down.",
        },
        "level shift convergence": {
            "applies_when": {"use level shift": "yes"},
            "default": 0.01,
            "kind": "float",
            "default_units": "",
            "enumeration": tuple(),
            "format_string": "s",
            "description": "Stop when converged to:",
            "help_text": "Convergence level to stop level shifting.",
        },
        "use soscf": {
            "default": "no",
            "kind": "boolean",
            "default_units": "",
            "enumeration": ("yes", "no"),
            "format_string": "s",
            "description": "Use second order SCF:",
            "help_text": ("Whether to use the second order SCF to help convergence."),
        },
        "soscf starting convergence": {
            "applies_when": {"use soscf": "yes"},
            "default": 0.01,
            "kind": "float",
            "default_units": "",
            "enumeration": tuple(),
            "format_string": "s",
            "description": "Start when converged to:",
            "help_text": "Convergence level to start second order SCF.",
        },
        "soscf convergence": {
            "applies_when": {"use soscf": "yes"},
            "default": 0.001,
            "kind": "float",
            "default_units": "",
            "enumeration": tuple(),
            "format_string": "s",
            "description": "Microiteration convergence:",
            "help_text": "Convergence level for SOSCF microiterations.",
        },
        "soscf max iterations": {
            "applies_when": {"use soscf": "yes"},
            "default": 5,
            "kind": "integer",
            "default_units": "",
            "enumeration": tuple(),
            "format_string": "s",
            "description": "Maximum microiterations:",
            "help_text": "Maximum number of SOSCF microiterations.",
        },
        "soscf print iterations": {
            "applies_when": {"use soscf": "yes"},
            "default": "no",
            "kind": "boolean",
            "default_units": "",
            "enumeration": ("yes", "no"),
            "format_string": "s",
            "description": "Print microiterations:",
            "help_text": "Print information about the SOSCF microiterations.",
        },
        "density convergence": {
            "default": "default",
            "kind": "float",
            "default_units": "",
            "enumeration": ("default",),
            "format_string": "s",
            "description": "Density convergence criterion:",
            "help_text": (
                "Criterion for convergence of the density, default 10^-6 for "
                "SCF, 10^-8 optimization."
            ),
        },
        "energy convergence": {
            "default": "default",
            "kind": "float",
            "default_units": "",
            "enumeration": ("default",),
            "format_string": "s",
            "description": "Energy convergence criterion:",
            "help_text": (
                "Criterion for convergence of the energy, default 10^-6 for "
                "SCF, 10^-8 optimization."
            ),
        },
        "CC amplitude convergence": {
            "default": "default",
            "kind": "float",
            "default_units": "",
            "enumeration": ("default",),
            "format_string": "s",
            "description": "CC amplitude convergence:",
            "help_text": (
                "Convergence of the residuals of the amplitude (and, for gradients, "
                "lambda) equations of the correlated method: Psi4's r_convergence. "
                "'default' is Psi4's own (1e-7)."
            ),
        },
        "convergence error": {
            "default": "yes",
            "kind": "boolean",
            "default_units": "",
            "enumeration": ("yes", "no"),
            "format_string": "s",
            "description": "Error if not converged:",
            "help_text": "Whether to throw an error if not converged.",
        },
        "maximum iterations": {
            "default": 100,
            "kind": "integer",
            "default_units": "",
            "enumeration": tuple(),
            "format_string": "s",
            "description": "Maximum iterations:",
            "help_text": "Maximum number of SCF iterations.",
        },
        "freeze-cores": {
            "default": "yes",
            "kind": "enumeration",
            "default_units": "",
            "enumeration": ("yes", "no"),
            "format_string": "s",
            "description": "Freeze core orbitals:",
            "help_text": (
                "Whether to freeze the core orbitals in correlated " "methods"
            ),
        },
        "correlation integrals": {
            "default": "default",
            "kind": "enumeration",
            "default_units": "",
            "enumeration": (
                "default",
                "conventional",
                "density fitting",
                "Cholesky decomposition",
            ),
            "format_string": "s",
            "description": "Correlation integrals:",
            "help_text": (
                "How the integrals of the correlated method are done. 'default' is "
                "Psi4's default for the method: density fitting for MP2, MP2.5 and "
                "MP3, conventional for the coupled-cluster methods. Density fitting "
                "is much faster and gives analytic gradients with frozen core (via "
                "DFOCC), but changes the energies slightly, so the model's name "
                "starts with 'DF-'. Cholesky decomposition has no analytic "
                "gradients."
            ),
        },
        "correlation fitting basis": {
            "applies_when": {"correlation integrals": "density fitting"},
            "default": "default",
            "kind": "string",
            "default_units": "",
            "enumeration": ("default",),
            "format_string": "s",
            "description": "Fitting basis:",
            "help_text": (
                "The auxiliary basis for the density-fitted correlation integrals "
                "(Psi4's df_basis_cc). 'default' uses the RI basis that matches the "
                "orbital basis, e.g. cc-pVTZ-RI for cc-pVTZ."
            ),
        },
        "stability analysis": {
            "default": "no",
            "kind": "boolean",
            "default_units": "",
            "enumeration": ("yes", "no"),
            "format_string": "s",
            "description": "Do stability analysis:",
            "help_text": ("Analyze the stability of the SCF/DFT wavefunction."),
        },
        "results": {
            "default": {},
            "kind": "dictionary",
            "default_units": "",
            "enumeration": tuple(),
            "format_string": "",
            "description": "results",
            "help_text": ("The results to save to variables or in " "tables. "),
        },
        "create tables": {
            "default": "yes",
            "kind": "boolean",
            "default_units": "",
            "enumeration": ("yes", "no"),
            "format_string": "",
            "description": "Create tables as needed:",
            "help_text": (
                "Whether to create tables as needed for "
                "results being saved into tables."
            ),
        },
    }

    output = {
        "density": {
            "default": "no",
            "kind": "boolean",
            "default_units": "",
            "enumeration": ("yes", "no"),
            "format_string": "",
            "description": "Plot total density:",
            "help_text": "Whether to plot the total charge density.",
        },
        "orbitals": {
            "default": "no",
            "kind": "boolean",
            "default_units": "",
            "enumeration": ("yes", "no"),
            "format_string": "",
            "description": "Plot orbitals:",
            "help_text": "Whether to plot orbitals.",
        },
        "selected orbitals": {
            "applies_when": {"orbitals": "yes"},
            "default": "-1, HOMO, LUMO, +1",
            "kind": "string",
            "default_units": "",
            "enumeration": ("all", "-1, HOMO, LUMO, +1"),
            "format_string": "",
            "description": "Selected orbitals:",
            "help_text": "Which orbitals to plot.",
        },
    }

    # Rules shared by the dialog and the flowchart builder (see seamm.Parameters).
    # The simple conditions are the "applies_when" entries above: the method (and
    # functional) of the chosen level of disclosure, and the sub-controls of the
    # damping, level shift, second-order SCF and orbital plots. The rules below
    # depend on the method's and functional's metadata.

    unused = ()
    """Parameters that this kind of step never uses."""

    def _method_key(self, values):
        """The method parameter in use for the level of disclosure."""
        return "advanced_method" if values.get("level") == "advanced" else "method"

    def _functional_key(self, values):
        """The functional parameter in use for the level of disclosure."""
        if values.get("level") == "advanced":
            return "advanced_functional"
        return "functional"

    @staticmethod
    def _method_record(name):
        """The metadata of a method, given by its full or short name, or None."""
        if name in methods:
            return methods[name]
        for record in methods.values():
            if record["method"] == str(name).lower():
                return record
        return None

    def _dispersions(self, values):
        """The dispersion corrections of the chosen functional, or None if it is
        not known (e.g. a variable)."""
        record = dft_functionals.get(values.get(self._functional_key(values)))
        return None if record is None else tuple(record["dispersion"])

    def applies(self, key, values=None, _seen=None):
        """As seamm.Parameters.applies, plus: the parameters that this kind of step
        does not use never apply; the frozen-core choice applies only to
        methods with a frozen core (the correlated ones), and the dispersion
        correction only to functionals that have dispersion corrections."""
        if values is None:
            values = self.current_values()
        if key in self.unused:
            return False
        if not super().applies(key, values, _seen):
            return False
        if key == "freeze-cores":
            name = values.get(self._method_key(values))
            if self._is_expr(name):
                return True
            record = self._method_record(name)
            return record is None or bool(record.get("freeze core?", False))
        if key in ("correlation integrals", "CC amplitude convergence"):
            from .derivatives import TYPE_KEYWORD

            name = values.get(self._method_key(values))
            if self._is_expr(name):
                return True
            record = self._method_record(name)
            return record is not None and record["method"] in TYPE_KEYWORD
        if key == "dispersion":
            functional_key = self._functional_key(values)
            if not self.applies(functional_key, values, _seen):
                return False
            if self._is_expr(values.get(functional_key)):
                return True
            dispersions = self._dispersions(values)
            return dispersions is not None and len(dispersions) > 1
        return True

    def not_applicable_reason(self, key, values=None):
        """Why a parameter does not apply, for the builder's messages."""
        if values is None:
            values = self.current_values()
        if key in self.unused:
            return "this kind of step does not use it"
        reason = super().not_applicable_reason(key, values)
        if reason or self.applies(key, values):
            return reason
        if key == "freeze-cores":
            name = values.get(self._method_key(values))
            return f"the method '{name}' does not freeze core orbitals"
        if key == "correlation integrals":
            name = values.get(self._method_key(values))
            return f"the method '{name}' has no choice of correlation integrals"
        if key == "CC amplitude convergence":
            name = values.get(self._method_key(values))
            return f"the method '{name}' has no amplitude equations"
        if key == "dispersion":
            functional_key = self._functional_key(values)
            if not self.applies(functional_key, values):
                why = self.not_applicable_reason(functional_key, values)
                return f"it needs '{functional_key}', which does not apply" + (
                    f" ({why})" if why else ""
                )
            functional = values.get(functional_key)
            if self._dispersions(values) is None:
                return f"the functional '{functional}' is not known"
            return f"the functional '{functional}' has no dispersion corrections"
        return ""

    def choices(self, key, values=None):
        """The dispersion corrections are those of the functional; CCD needs
        density-fitted or Cholesky-decomposed integrals."""
        if values is None:
            values = self.current_values()
        if key == "dispersion":
            return self._dispersions(values)
        if key == "correlation integrals" and self._is_ccd(values):
            # Psi4 has CCD only with density fitting or Cholesky decomposition
            return ("density fitting", "Cholesky decomposition")
        return super().choices(key, values)

    def _is_ccd(self, values):
        """Whether the method is CCD."""
        record = self._method_record(values.get(self._method_key(values)))
        return record is not None and record["method"] == "ccd"

    def implied(self, values=None):
        """A functional implies a dispersion correction it has: if the current one
        is not available, the first real correction (as the dialog chooses). CCD
        implies density fitting unless Cholesky decomposition is chosen."""
        if values is None:
            values = self.current_values()
        result = {}
        if self.applies("dispersion", values):
            dispersions = self._dispersions(values)
            dispersion = values.get("dispersion")
            if (
                dispersions is not None
                and dispersion not in dispersions
                and not self._is_expr(dispersion)
            ):
                result["dispersion"] = dispersions[1]
        if self.applies("correlation integrals", values) and self._is_ccd(values):
            if values.get("correlation integrals") in ("default", "conventional"):
                result["correlation integrals"] = "density fitting"
        return result

    def __init__(self, defaults={}, data=None):
        """Initialize the instance, by default from the default
        parameters given in the class"""

        super().__init__(
            defaults={
                **EnergyParameters.parameters,
                **EnergyParameters.output,
                **defaults,
            },
            data=data,
        )

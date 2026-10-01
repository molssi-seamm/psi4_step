# -*- coding: utf-8 -*-

"""Smoke test of the Tk dialogs: create them and re-lay them out for every choice
that drives the layout, checking both ways that the controls shown are exactly
those that the parameters' rules say apply. Skipped when no display is available."""

import pytest

SUBSTEPS = ("Energy", "Optimization", "Thermochemistry", "BSSE")


@pytest.fixture(scope="module")
def root():
    import tkinter as tk

    try:
        root = tk.Tk()
    except tk.TclError:
        pytest.skip("no display available for Tk")
    root.withdraw()
    import Pmw

    Pmw.initialise(root)
    yield root
    root.destroy()


def make(root, substep):
    import seamm

    flowchart = seamm.Flowchart(namespace="org.molssi.seamm.psi4", directory=".")
    tk_flowchart = seamm.TkFlowchart(
        master=root, flowchart=flowchart, namespace="org.molssi.seamm.psi4.tk"
    )
    node = flowchart.create_node(substep)
    flowchart.add_node(node)
    plugin = tk_flowchart.plugin_manager.get(substep)
    tk_node = plugin.create_tk_node(
        tk_flowchart=tk_flowchart, node=node, canvas=tk_flowchart.canvas, x=100, y=100
    )
    tk_node.create_dialog()
    return tk_node


def shown(tk_node, key):
    """Whether a control is laid out, i.e. it and all the frames holding it are
    managed (gridded, packed or a notebook page)."""
    widget = tk_node[key]
    while widget is not None:
        manager = widget.winfo_manager()
        if manager == "wm":
            return True
        if manager == "":
            return False
        widget = widget.master
    return True


def check(tk_node):
    """The shown controls are exactly those that apply."""
    tk_node.reset_dialog()
    tk_node.reset_plotting()  # its own tab, laid out when 'orbitals' changes
    P = tk_node.node.parameters
    values = tk_node._widget_values()
    for key in P:
        if key in ("results", "create tables") or key not in tk_node:
            continue
        if shown(tk_node, key):
            assert P.applies(key, values), f"{key} is shown but does not apply"
        else:
            assert not P.applies(key, values), f"{key} applies but is not shown"
    return values


def set_and_check(tk_node, key, value):
    tk_node[key].set(value)
    return check(tk_node)


@pytest.mark.parametrize("substep", SUBSTEPS)
def test_layouts_follow_the_rules(root, substep):
    import psi4_step

    tk_node = make(root, substep)
    P = tk_node.node.parameters

    if substep == "Thermochemistry":
        set_and_check(tk_node, "use existing parameters", "yes")
        assert not shown(tk_node, "method")
        set_and_check(tk_node, "use existing parameters", "no")
        assert shown(tk_node, "method")

    for level, method_key, functional_key in (
        ("recommended", "method", "functional"),
        ("advanced", "advanced_method", "advanced_functional"),
    ):
        set_and_check(tk_node, "level", level)
        # Every method for Energy; one of each kind for the other sub-steps
        methods = P[method_key].enumeration
        if substep != "Energy":
            methods = {
                (
                    psi4_step.methods[m]["method"] == "dft",
                    psi4_step.methods[m].get("freeze core?", False),
                ): m
                for m in methods
            }.values()
        for method in methods:
            set_and_check(tk_node, method_key, method)
            record = psi4_step.methods[method]
            assert shown(tk_node, functional_key) == (record["method"] == "dft")
            assert shown(tk_node, "freeze-cores") == bool(
                record.get("freeze core?", False)
            )
        # DFT: the dispersion corrections offered are those of the functional.
        # One functional of each set of dispersion corrections.
        tk_node[method_key].set(psi4_step.energy_parameters.DFT_METHODS[0])
        functionals = {
            tuple(psi4_step.dft_functionals[f]["dispersion"]): f
            for f in P[functional_key].enumeration
        }.values()
        for functional in functionals:
            values = set_and_check(tk_node, functional_key, functional)
            dispersions = psi4_step.dft_functionals[functional]["dispersion"]
            assert shown(tk_node, "dispersion") == (len(dispersions) > 1)
            if len(dispersions) > 1:
                offered = tk_node["dispersion"].combobox.cget("values")
                assert list(offered) == list(dispersions)
                assert values["dispersion"] in dispersions
        # A variable for the method could be DFT and correlated
        tk_node[functional_key].set("B3LYP Hyb-GGA Exchange-Correlation Functional")
        values = set_and_check(tk_node, method_key, "$method")
        for key in (functional_key, "freeze-cores", "dispersion"):
            assert shown(tk_node, key)

    for switch in ("use damping", "use level shift", "use soscf", "orbitals"):
        for value in ("yes", "no"):
            set_and_check(tk_node, switch, value)

    if substep == "BSSE":
        for fragments in P["fragments"].enumeration:
            set_and_check(tk_node, "fragments", fragments)
            assert shown(tk_node, "fragment atoms") == (fragments == "specified")


def test_dispersion_kept_valid(root):
    """A functional without the chosen dispersion correction gets one it has."""
    import psi4_step

    tk_node = make(root, "Energy")
    functional = next(
        name
        for name, record in psi4_step.dft_functionals.items()
        if record["dispersion"] == ["none", "d3bj"]
    )
    tk_node["level"].set("advanced")
    tk_node["advanced_functional"].set(functional)
    tk_node["dispersion"].set("d3mbj")
    tk_node.reset_dialog()
    assert tk_node["dispersion"].get() == "d3bj"

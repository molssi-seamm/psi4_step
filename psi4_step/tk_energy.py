# -*- coding: utf-8 -*-

"""The graphical part of a Psi4 Energy node"""

import logging
import tkinter as tk
import tkinter.ttk as ttk

import psi4_step
import seamm
import seamm_widgets as sw

logger = logging.getLogger(__name__)


class TkEnergy(seamm.TkNode):
    def __init__(
        self,
        tk_flowchart=None,
        node=None,
        canvas=None,
        x=120,
        y=20,
        w=200,
        h=50,
        my_logger=logger,
    ):
        """Initialize the graphical Tk Psi4 energy step

        Keyword arguments:
        """
        self.results_widgets = []

        # Set the logging level for this module if requested
        # if 'psi4_tk_energy_log_level' in self.options:
        #     my_logger.setLevel(self.options.psi4_tk_energy_log_level)
        #     my_logger.critical(
        #         'Set log level to {}'.format(
        #             self.options.psi4_tk_energy_log_level
        #         )
        #     )

        super().__init__(
            tk_flowchart=tk_flowchart,
            node=node,
            canvas=canvas,
            x=x,
            y=y,
            w=w,
            h=h,
            my_logger=my_logger,
        )

    def right_click(self, event):
        """Probably need to add our dialog..."""

        super().right_click(event)
        self.popup_menu.add_command(label="Edit..", command=self.edit)

        self.popup_menu.tk_popup(event.x_root, event.y_root, 0)

    def create_dialog(self, title="Edit Psi4 Energy Step"):
        """Create the dialog!"""
        self.logger.debug("Creating the dialog")
        frame = super().create_dialog(title=title, widget="notebook", results_tab=True)

        # Create a frame for the calculation control
        self["calculation"] = ttk.LabelFrame(
            frame,
            borderwidth=4,
            relief="sunken",
            text="Calculation",
            labelanchor="n",
            padding=10,
        )
        # Create a frame for the convergence control
        self["convergence"] = ttk.LabelFrame(
            frame,
            borderwidth=4,
            relief="sunken",
            text="SCF Convergence Control",
            labelanchor="n",
            padding=10,
        )

        # Create all the widgets
        P = self.node.parameters
        for key in (
            "level",
            "method",
            "advanced_method",
            "functional",
            "advanced_functional",
            "dispersion",
            "spin-restricted",
            "freeze-cores",
            "stability analysis",
        ):
            self[key] = P[key].widget(self["calculation"])
        for key in (
            "use damping",
            "damping percentage",
            "damping convergence",
            "use level shift",
            "level shift",
            "level shift convergence",
            "use soscf",
            "soscf starting convergence",
            "soscf convergence",
            "soscf max iterations",
            "soscf print iterations",
            "density convergence",
            "energy convergence",
            "convergence error",
            "maximum iterations",
        ):
            self[key] = P[key].widget(self["convergence"])

        # bindings...
        self["level"].combobox.bind("<<ComboboxSelected>>", self.reset_calculation)
        self["level"].config(state="readonly")

        self["method"].combobox.bind("<<ComboboxSelected>>", self.reset_calculation)
        self["method"].combobox.bind("<Return>", self.reset_calculation)
        self["method"].combobox.bind("<FocusOut>", self.reset_calculation)

        self["advanced_method"].combobox.bind(
            "<<ComboboxSelected>>", self.reset_calculation
        )
        self["advanced_method"].combobox.bind("<Return>", self.reset_calculation)
        self["advanced_method"].combobox.bind("<FocusOut>", self.reset_calculation)

        # The dispersion corrections offered depend on the functional
        for key in ("functional", "advanced_functional"):
            self[key].combobox.bind("<<ComboboxSelected>>", self.reset_calculation)
            self[key].combobox.bind("<Return>", self.reset_calculation)
            self[key].combobox.bind("<FocusOut>", self.reset_calculation)

        for key in ("use damping", "use level shift", "use soscf"):
            self[key].bind("<<ComboboxSelected>>", self.reset_convergence)
            self[key].bind("<Return>", self.reset_convergence)
            self[key].bind("<FocusOut>", self.reset_convergence)

        # A tab for output -- orbitals, etc.
        notebook = self["notebook"]
        self["output frame"] = oframe = ttk.Frame(notebook)
        notebook.insert(self["results frame"], oframe, text="Output", sticky="new")

        # Frame to isolate widgets
        p_frame = self["plot frame"] = ttk.LabelFrame(
            self["output frame"],
            borderwidth=4,
            relief="sunken",
            text="Plots",
            labelanchor="n",
            padding=10,
        )

        for key in psi4_step.EnergyParameters.output:
            self[key] = P[key].widget(p_frame)

        # Set the callbacks for changes
        for widget in ("orbitals",):
            w = self[widget]
            w.combobox.bind("<<ComboboxSelected>>", self.reset_plotting)
            w.combobox.bind("<Return>", self.reset_plotting)
            w.combobox.bind("<FocusOut>", self.reset_plotting)
        p_frame.grid(row=0, column=0, sticky="new")
        oframe.columnconfigure(0, weight=1)

        self.reset_plotting()

        # Top level needs to call reset_dialog
        if self.node.calculation == "energy":
            self.reset_dialog()

        self.logger.debug("Finished creating the dialog")

        return frame

    def reset_dialog(self, row=0, widget=None):
        """Layout the widgets as needed for the current state"""

        frame = self["frame"]
        if row == 0:
            for slave in frame.grid_slaves():
                slave.grid_forget()

        self["calculation"].grid(row=row, column=0)
        self.reset_calculation()
        row += 1
        self["convergence"].grid(row=row, column=0)
        self.reset_convergence()
        row += 1
        self.reset_plotting()
        return row

    def _widget_values(self):
        """The dialog's current values, {name: value}, for the parameters' rules."""
        values = {}
        for key in self.node.parameters:
            if key == "results" or key not in self:
                continue
            try:
                value = self[key].get()
            except Exception:
                continue
            values[key] = value[0] if isinstance(value, tuple) else value
        return values

    def reset_calculation(self, widget=None):
        # Which controls to show, and what they offer, come from the parameters'
        # rules (psi4_step.EnergyParameters), which the flowchart builder uses too.
        P = self.node.parameters
        values = self._widget_values()

        def applies(key):
            return P.applies(key, values)

        if applies("method"):
            long_method = values["method"]
        else:
            long_method = values["advanced_method"]
        meta = None if self.is_expr(long_method) else P._method_record(long_method)
        self.node.method = None if meta is None else meta["method"]

        # Set up the results table because it depends on the method
        self.results_widgets = []
        self.setup_results()

        frame = self["calculation"]
        for slave in frame.grid_slaves():
            slave.grid_forget()

        widgets = []
        widgets2 = []
        row = 0
        self["level"].grid(row=row, column=0, columnspan=2, sticky=tk.EW)
        row += 1
        is_dft = False
        for method, functional in (
            ("method", "functional"),
            ("advanced_method", "advanced_functional"),
        ):
            if not applies(method):
                continue
            self[method].grid(row=row, column=0, columnspan=2, sticky=tk.EW)
            widgets.append(self[method])
            row += 1
            if applies(functional):
                is_dft = True
                self[functional].grid(row=row, column=1, sticky=tk.EW)
                widgets2.append(self[functional])
                row += 1
            if applies("freeze-cores"):
                self["freeze-cores"].grid(row=row, column=1, sticky=tk.EW)
                widgets2.append(self["freeze-cores"])
                row += 1
        if applies("dispersion"):
            # Offer only the functional's dispersion corrections, keeping the
            # choice valid.
            w = self["dispersion"]
            dispersions = P.choices("dispersion", values)
            if dispersions is None:
                dispersions = P["dispersion"].enumeration
            w.config(values=list(dispersions))
            implied = P.implied(values)
            if "dispersion" in implied:
                w.set(implied["dispersion"])
            w.grid(row=row, column=1, sticky=tk.W)
            widgets2.append(w)
            row += 1
        if is_dft:
            sw.align_labels(widgets2)
            frame.columnconfigure(0, minsize=30)
        for key in ("spin-restricted", "stability analysis"):
            if applies(key):
                self[key].grid(row=row, column=0, columnspan=2, sticky=tk.EW)
                widgets.append(self[key])
                row += 1
        sw.align_labels(widgets, sticky=tk.E)

        return row

    def reset_convergence(self, widget=None):
        """Layout the convergence widgets as needed for the current state"""

        frame = self["convergence"]
        for slave in frame.grid_slaves():
            slave.grid_forget()

        P = self.node.parameters
        values = self._widget_values()

        widgets = []
        widgets2 = []
        row = 0

        for key in (
            "maximum iterations",
            "density convergence",
            "energy convergence",
            "convergence error",
        ):
            if P.applies(key, values):
                self[key].grid(row=row, column=0, columnspan=2, sticky=tk.EW)
                widgets.append(self[key])
                row += 1

        # Each switch is followed by its indented sub-controls, which apply only
        # when it is on.
        for switch, keys in (
            ("use damping", ("damping percentage", "damping convergence")),
            ("use level shift", ("level shift", "level shift convergence")),
            (
                "use soscf",
                (
                    "soscf starting convergence",
                    "soscf convergence",
                    "soscf max iterations",
                    "soscf print iterations",
                ),
            ),
        ):
            if not P.applies(switch, values):
                continue
            self[switch].grid(row=row, column=0, columnspan=2, sticky=tk.EW)
            widgets.append(self[switch])
            row += 1
            for key in keys:
                if P.applies(key, values):
                    self[key].grid(row=row, column=1, sticky=tk.EW)
                    widgets2.append(self[key])
                    row += 1

        frame.columnconfigure(0, minsize=150)
        sw.align_labels(widgets, sticky=tk.E)
        sw.align_labels(widgets2, sticky=tk.E)

    def reset_plotting(self, widget=None):
        frame = self["plot frame"]
        for slave in frame.grid_slaves():
            slave.grid_forget()

        P = self.node.parameters
        values = self._widget_values()

        widgets = []

        # The Output tab is shown only if some of its plots apply
        notebook = self["notebook"]
        oframe = self["output frame"]
        output = psi4_step.EnergyParameters.output
        if any(P.applies(key, values) for key in output):
            notebook.add(oframe)
        else:
            notebook.hide(oframe)

        row = 0
        for key in (
            "density",
            "orbitals",
        ):
            if not P.applies(key, values):
                continue
            self[key].grid(row=row, column=0, columnspan=4, sticky=tk.EW)
            widgets.append(self[key])
            row += 1

        key = "selected orbitals"
        if P.applies(key, values):
            self[key].grid(row=row, column=1, columnspan=4, sticky=tk.EW)
            row += 1

        sw.align_labels(widgets, sticky=tk.E)
        frame.columnconfigure(0, minsize=50)
        frame.columnconfigure(4, weight=1)

        return row

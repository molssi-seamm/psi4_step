# -*- coding: utf-8 -*-

"""The graphical part of a Psi4 BSSE (counterpoise) sub-step."""

import logging
import tkinter as tk
import tkinter.ttk as ttk

import psi4_step
import seamm_widgets as sw

logger = logging.getLogger(__name__)


class TkBSSE(psi4_step.TkEnergy):
    """Graphical Psi4 BSSE sub-step: the energy dialog's level-of-theory
    controls plus the fragment definition."""

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
        self.results_widgets = []
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
        super().right_click(event)
        self.popup_menu.add_command(label="Edit..", command=self.edit)
        self.popup_menu.tk_popup(event.x_root, event.y_root, 0)

    def create_dialog(self, title="Edit Psi4 BSSE Step"):
        """Create the edit dialog: the Energy superclass builds the
        level-of-theory 'calculation'/'convergence' frames; this adds a
        'fragments' frame alongside them, reactive to the 'fragments' mode
        (the 'Fragment atoms' field only applies when it is 'specified')."""
        logger.debug("TkBSSE.create_dialog")

        frame = super().create_dialog(title=title)

        P = self.node.parameters

        fragments_frame = self["fragments frame"] = ttk.LabelFrame(
            frame,
            borderwidth=4,
            relief="sunken",
            text="Fragments",
            labelanchor="n",
            padding=10,
        )

        for key in (
            "fragments",
            "fragment atoms",
            "fragment charges",
            "compute gradient",
        ):
            self[key] = P[key].widget(fragments_frame)

        for sequence in ("<<ComboboxSelected>>", "<Return>", "<FocusOut>"):
            self["fragments"].combobox.bind(sequence, self.reset_fragments)

        if self.node.calculation == "bsse":
            self.reset_dialog()

        return frame

    def reset_dialog(self, widget=None):
        """Lay out the widgets, letting the Energy superclass go first."""
        rows = super().reset_dialog()

        self["fragments frame"].grid(row=0, column=1, rowspan=rows, sticky=tk.N)
        self.reset_fragments()

        return rows

    def reset_fragments(self, widget=None):
        """'Fragment atoms' only applies when defining the fragments by
        hand; 'Fragment charges'/'compute gradient' apply either way."""
        frame = self["fragments frame"]
        for slave in frame.grid_slaves():
            slave.grid_forget()

        widgets = []
        row = 0

        def add(key):
            nonlocal row
            self[key].grid(row=row, column=0, columnspan=2, sticky=tk.EW)
            widgets.append(self[key])
            row += 1

        add("fragments")
        if self["fragments"].get() == "specified":
            add("fragment atoms")
        add("fragment charges")
        add("compute gradient")

        sw.align_labels(widgets)

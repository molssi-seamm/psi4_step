# -*- coding: utf-8 -*-

"""The graphical part of a Thermochemistry step"""

import logging
import tkinter as tk
import tkinter.ttk as ttk

import psi4_step
import seamm_widgets as sw

logger = logging.getLogger(__name__)


class TkThermochemistry(psi4_step.TkEnergy):
    """
    The graphical part of a Thermochemistry step in a flowchart.

    Attributes
    ----------
    tk_flowchart : TkFlowchart = None
        The flowchart that we belong to.
    node : Node = None
        The corresponding node of the non-graphical flowchart
    canvas: tkCanvas = None
        The Tk Canvas to draw on
    dialog : Dialog
        The Pmw dialog object
    x : int = None
        The x-coordinate of the center of the picture of the node
    y : int = None
        The y-coordinate of the center of the picture of the node
    w : int = 200
        The width in pixels of the picture of the node
    h : int = 50
        The height in pixels of the picture of the node
    self[widget] : dict
        A dictionary of tk widgets built using the information
        contained in Thermochemistry_parameters.py

    See Also
    --------
    Thermochemistry, TkThermochemistry,
    ThermochemistryParameters,
    """

    def __init__(
        self,
        tk_flowchart=None,
        node=None,
        canvas=None,
        x=None,
        y=None,
        w=200,
        h=50,
        my_logger=logger,
    ):
        """
        Initialize a graphical node.

        Parameters
        ----------
        tk_flowchart: Tk_Flowchart
            The graphical flowchart that we are in.
        node: Node
            The non-graphical node for this step.
        namespace: str
            The stevedore namespace for finding sub-nodes.
        canvas: Canvas
           The Tk canvas to draw on.
        x: float
            The x position of the nodes center on the canvas.
        y: float
            The y position of the nodes cetner on the canvas.
        w: float
            The nodes graphical width, in pixels.
        h: float
            The nodes graphical height, in pixels.

        Returns
        -------
        None
        """
        self.dialog = None

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

    def create_dialog(self):
        """
        Create the dialog. A set of widgets will be chosen by default
        based on what is specified in the Thermochemistry_parameters
        module.

        Parameters
        ----------
        None

        Returns
        -------
        None

        See Also
        --------
        TkThermochemistry.reset_dialog
        """

        frame = super().create_dialog(title="Thermochemistry")

        # Shortcut for parameters
        P = self.node.parameters

        # Frame to isolate widgets
        thermo_frame = self["thermochemistry"] = ttk.LabelFrame(
            frame,
            borderwidth=4,
            relief="sunken",
            text="Thermochemistry",
            labelanchor="n",
            padding=10,
        )
        # Then create the widgets
        for key in psi4_step.ThermochemistryParameters.parameters:
            self[key] = P[key].widget(thermo_frame)

        # and binding to change as needed
        for key in ("use existing parameters",):
            self[key].combobox.bind("<<ComboboxSelected>>", self.reset_thermochemistry)
            self[key].combobox.bind("<Return>", self.reset_thermochemistry)
            self[key].combobox.bind("<FocusOut>", self.reset_thermochemistry)

        # Top level needs to call reset_dialog
        if self.node.calculation == "thermochemistry":
            self.reset_dialog()

        return frame

    def reset_dialog(self, widget=None, row=0):
        """Layout the widgets in the dialog.

        The widgets are chosen by default from the information in
        Thermochemistry_parameter.

        This function simply lays them out row by row with
        aligned labels. You may wish a more complicated layout that
        is controlled by values of some of the control parameters.
        If so, edit or override this method

        Parameters
        ----------
        widget : Tk Widget = None

        Returns
        -------
        None

        See Also
        --------
        TkThermochemistry.create_dialog
        """
        # Remove any widgets previously packed
        frame = self["frame"]
        if row == 0:
            for slave in frame.grid_slaves():
                slave.grid_forget()

        # Shortcut for parameters
        P = self.node.parameters

        self["thermochemistry"].grid(row=row, column=0)
        row += 1

        self.reset_thermochemistry()

        if not P["use existing parameters"]:
            row = super().reset_dialog(row=row)
        else:
            self.reset_calculation()
        return row

    def reset_thermochemistry(self, widget=None):
        frame = self["thermochemistry"]
        for slave in frame.grid_slaves():
            slave.grid_forget()

        widgets = []
        row = 0

        for key in ("use existing parameters", "T", "P"):
            self[key].grid(row=row, column=0, sticky=tk.EW)
            widgets.append(self[key])
            row += 1

        sw.align_labels(widgets)

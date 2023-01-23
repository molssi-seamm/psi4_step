# -*- coding: utf-8 -*-

"""Main module."""

import psi4_step


class AcceleratedOptimizationStep(object):
    my_description = {
        "description": "Optimization of the structure using Psi4 accelerated by MOPAC",
        "group": "Calculation",
        "name": "Accelerated Optimization",
    }

    def __init__(self, flowchart=None, gui=None):
        """Initialize this helper class, which is used by
        the application via stevedore to get information about
        and create node objects for the flowchart
        """
        pass

    def description(self):
        """Return a description of what this extension does"""
        return AcceleratedOptimizationStep.my_description

    def create_node(self, flowchart=None, **kwargs):
        """Return the new node object"""
        return psi4_step.AcceleratedOptimization(flowchart=flowchart, **kwargs)

    def create_tk_node(self, canvas=None, **kwargs):
        """Return the graphical Tk node object"""
        return psi4_step.TkAcceleratedOptimization(canvas=canvas, **kwargs)

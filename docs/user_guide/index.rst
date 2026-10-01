.. _user-guide:

**********
User Guide
**********
The Psi4 plug-in provides access to many different quantum chemical calculations usign
the Psi4 code.

..
   The following sections cover accessing and controlling this functionality.

   .. toctree::
      :maxdepth: 2
      :titlesonly:

Settings that depend on each other
==================================

Which settings apply depends on others: the functional only for DFT, the frozen core
only for correlated methods, a dispersion correction only from the functional's own
list, and the damping, level-shift and second-order SCF controls only when switched on.
The step's dialogs show only the settings that apply with the current choices, and the
same rules are used when a flowchart is built or edited without the editor
(``seamm-flowchart`` or SEAMM's MCP server): a setting that would have no effect is
refused, with the reason, and a value that contradicts another is refused too. See
"Flowcharts without the editor" in SEAMM's user guide.

BSSE uses Psi4's own SCF settings and makes no plots, so it does not offer those Energy
settings. Thermochemistry's plots apply only when it uses its own settings rather than
the previous step's; placed right after Initialization, where there is no previous
calculation, it always uses its own.


Index
=====

* :ref:`genindex`

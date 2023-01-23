# -*- coding: utf-8 -*-

"""Installer for the Psi4 plug-in.

This handles any further installation needed after installing the Python
package `psi4-step`.
"""

import logging
from pathlib import Path
import pkg_resources
import subprocess

import seamm_installer

logger = logging.getLogger(__name__)


class Installer(seamm_installer.InstallerBase):
    """Handle further installation needed after installing psi4-step.

    The Python package `psi4-step` should already be installed, using `pip`,
    `conda`, or similar. This plug-in-specific installer then checks for the
    Psi4 executable, installing it if needed, and registers its
    location in seamm.ini.

    There are a number of ways to determine which are the correct Psi4
    executables to use. The aim of this installer is to help the user locate
    the executables. There are a number of possibilities

    1. The correct executables are already available.

        #. If they are already registered in `seamm.ini` there is nothing else
           to do.

        #. They may be in the current path, in which case they need to be added
           to `seamm.ini`.

        #. If a module system is in use, a module may need to be loaded to give
           access to Psi4.

        #. They cannot be found automatically, so the user needs to locate the
           executables for the installer.

    #. Psi4 is not installed on the machine. In this case they can be
       installed in a Conda environment. There is one choice

        #. They can be installed in a separate environment, `seamm-psi4` by
           default.
    """

    def __init__(self, logger=logger):
        # Call the base class initialization, which sets up the commandline
        # parser, amongst other things.
        super().__init__(logger=logger)

        logger.debug("Initializing the Psi4 installer object.")

        self.section = "psi4-step"
        self.path_name = "psi4-path"
        self.executables = ["psi4"]
        self.resource_path = Path(pkg_resources.resource_filename(__name__, "data/"))
        # What Conda environment is the default?
        data = self.configuration.get_values(self.section)
        if "conda-environment" in data and data["conda-environment"] != "":
            self.environment = data["conda-environment"]
        else:
            self.environment = "seamm-psi4"

        # The environment.yaml file for Conda installations.
        path = Path(pkg_resources.resource_filename(__name__, "data/"))
        logger.debug(f"data directory: {path}")
        self.environment_file = path / "seamm-psi4.yml"

    def exe_version(self, path):
        """Get the version of the Psi4 executable.

        Parameters
        ----------
        path : pathlib.Path
            Path to the executable.

        Returns
        -------
        str
            The version reported by the executable, or 'unknown'.
        """
        try:
            result = subprocess.run(
                [str(path), "--version"],
                stdin=subprocess.DEVNULL,
                capture_output=True,
                text=True,
            )
        except Exception:
            version = "unknown"
        else:
            version = "unknown"
            lines = result.stdout.splitlines()
            for line in lines:
                line = line.strip()
                tmp = line.split()
                if len(tmp) == 1:
                    version = tmp[0]
                    break

        # Psi4 tends to leave timing.dat floating around. :-(
        timing = Path("timer.dat")
        if timing.exists():
            try:
                timing.unlink()
            except Exception:
                pass

        return version

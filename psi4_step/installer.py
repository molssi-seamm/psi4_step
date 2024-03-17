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
        self.executables = ["psi4"]
        self.resource_path = Path(pkg_resources.resource_filename(__name__, "data/"))

        # The environment.yaml file for Conda installations.
        logger.debug(f"data directory: {self.resource_path}")
        self.environment_file = self.resource_path / "seamm-psi4.yml"

    def check(self):
        """Check the status of the Psi4 installation."""
        print("Checking the Psi4 installation.")

        # What Conda environment is the default?
        path = self.configuration.path.parent / "psi4.ini"
        if not path.exists():
            text = (self.resource_path / "psi4.ini").read_text()
            path.write_text(text)
            print(f"    The psi4.ini file did not exist. Created {path}")

        self.exe_config.path = path

        # Get the current values
        data = self.exe_config.get_values("local")

        if "conda-environment" in data and data["conda-environment"] != "":
            self.environment = data["conda-environment"]
        else:
            self.environment = "seamm-psi4"

        super().check()

    def install(self):
        """Install Psi4 in a conda environment."""
        print("Installing Psi4.")

        # What Conda environment is the default?
        path = self.configuration.path.parent / "psi4.ini"
        if not path.exists():
            text = (self.resource_path / "psi4.ini").read_text()
            path.write_text(text)
            print(f"    The psi4.ini file did not exist. Created {path}")

        self.exe_config.path = path

        # Get the current values
        data = self.exe_config.get_values("local")

        if "conda-environment" in data and data["conda-environment"] != "":
            self.environment = data["conda-environment"]
        else:
            self.environment = "seamm-psi4"

        super().install()

    def show(self):
        """Show the status of the Psi4 installation."""
        print("Showing the Psi4 installation.")

        # What Conda environment is the default?
        path = self.configuration.path.parent / "psi4.ini"
        if not path.exists():
            text = (self.resource_path / "psi4.ini").read_text()
            path.write_text(text)
            print(f"    The psi4.ini file does not exist at {path}")
            print("    The 'check' command will create it if Psi4 is installed.")
            print("    Otherwise 'install' will install Psi4.")
            return

        self.exe_config.path = path

        if not self.exe_config.section_exists("local"):
            print(
                "    Psi4 is not configured: there is no 'local' section in "
                f"     {path}."
            )
            return

        # Get the current values
        data = self.exe_config.get_values("local")

        if "conda-environment" in data and data["conda-environment"] != "":
            self.environment = data["conda-environment"]
        else:
            self.environment = "seamm-psi4"

        super().show()

    def uninstall(self):
        """Uninstall the Psi4 installation."""
        print("Uninstall the Psi4 installation.")

        # What Conda environment is the default?
        path = self.configuration.path.parent / "psi4.ini"
        if not path.exists():
            text = (self.resource_path / "psi4.ini").read_text()
            path.write_text(text)
            print(
                f""""    The psi4.ini file does not exist at {path}
    Perhaps Psi4 is not installed, but if it is the 'check' command may locate it
    and create the ini file, after which 'uninstall' will remove it."""
            )
            return

        self.exe_config.path = path

        if not self.exe_config.section_exists("local"):
            print(
                f""""    The psi4.ini file at {path} does not have local section.
    Perhaps Psi4 not installed, but if it is the 'check' command may locate it
    and update the ini file, after which 'uninstall' will remove it."""
            )
            return

        # Get the current values
        data = self.exe_config.get_values("local")

        if "conda-environment" in data and data["conda-environment"] != "":
            self.environment = data["conda-environment"]
        else:
            self.environment = "seamm-psi4"

        super().uninstall()

    def update(self):
        """Updates the Psi4 installation."""
        print("Updating the Psi4 installation.")

        # What Conda environment is the default?
        path = self.configuration.path.parent / "psi4.ini"
        if not path.exists():
            text = (self.resource_path / "psi4.ini").read_text()
            path.write_text(text)
            print(f"    The psi4.ini file did not exist. Created {path}")

        self.exe_config.path = path

        # Get the current values
        data = self.exe_config.get_values("local")

        if "conda-environment" in data and data["conda-environment"] != "":
            self.environment = data["conda-environment"]
        else:
            self.environment = "seamm-psi4"

        super().update()

    def exe_version(self, config):
        """Get the version of the Psi4 executable.

        Parameters
        ----------
        config : dict
            Dictionary of configuration parameters from psi4.ini

        Returns
        -------
        str, str
            "Psi4" and the version reported by the executable, or 'unknown'.
        """
        environment = config["conda-environment"]
        conda = config["conda"]
        if environment[0] == "~":
            environment = str(Path(environment).expanduser())
            command = f"'{conda}' run --live-stream -p '{environment}' psi4 --version"
        elif Path(environment).is_absolute():
            command = f"'{conda}' run --live-stream -p '{environment}' psi4 --version"
        else:
            command = f"'{conda}' run --live-stream -n '{environment}' psi4 --version"

        logger.debug(f"    Running {command}")
        try:
            result = subprocess.run(
                command,
                stdin=subprocess.DEVNULL,
                capture_output=True,
                text=True,
                shell=True,
            )
        except Exception as e:
            logger.debug(f"    Failed to run {command}: {e}")
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

        return "Psi4", version

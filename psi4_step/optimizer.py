# -*- coding: utf-8 -*-

"""Acclerated optimization using MOPAC for Hessians."""

import logging

import numpy as np


def norm(V):
    return np.linalg.norm(V)


def abs_max(V):
    return max(abs(elem) for elem in V)


def abs_min(V):
    return min(abs(elem) for elem in V)


def rms(V):
    return np.sqrt(np.mean(V**2))


logger = logging.getLogger(__name__)


class Optimizer:
    def __init__(self, is_linear=False, energy_callback=None, surrogate_callback=None):
        self._E = None
        self._dE = None
        self._d2E = None
        self._iteration = 0
        self.is_linear = is_linear
        self._energy_callback = energy_callback
        self._surrogate_callback = surrogate_callback

    @property
    def E(self):
        """The energy."""
        return self._E

    @property
    def dE(self):
        """The derivatives of the energy."""
        return self._dE

    @property
    def d2E(self):
        """The second derivatives of the energy."""
        return self._d2E

    def optimize(self, _x0, n_steps=None):
        """Do the optimization!"""
        self.x0 = np.array(_x0).flatten()
        n_dof = len(self.x0)
        if self.is_linear:
            n_dof -= 5
        else:
            n_dof -= 6

        if n_steps is None:
            n_steps = 3 * n_dof

        x = np.array(self.x0)
        xlast = None
        print(
            f"{'Step':4} {'Energy':12} {'Max dE':9} {'rms dE':9} "
            f"{'Max step':9} {'rms step':9}"
        )
        print("---- ------------ --------- --------- --------- ---------")
        for step in range(0, n_steps + 1):
            if self._surrogate_callback is None:
                E, dE, d2E = self._energy_callback(x, ("E", "dE", "d2E"))
            else:
                Es, dEs, d2Es = self._surrogate_callback(x, ("E", "dE", "d2E"))
                E0, dE0 = self._energy_callback(x, ("E", "dE"))

                # Corrections to the surrogate energy and derivatives
                self.Ecorr = E0 - Es
                self.dEcorr = dE0 - dEs
                self.x0 = np.array(x)

                # Invert the Hessian carefully
                s = np.linalg.svd(d2Es, compute_uv=False, hermitian=True)
                if self.is_linear:
                    rcond = ((s[-6] + s[-5]) / 2) / s[0]
                else:
                    rcond = ((s[-7] + s[-6]) / 2) / s[0]

                A = np.linalg.pinv(d2Es, rcond=rcond)

                E = E0
                dE = np.array(dE0)

                Etmp = Es + self.Ecorr + np.dot(x - self.x0, self.dEcorr)
                dEtmp = dEs + self.dEcorr

                tmp = " ".join([f"{v:9.6f}" for v in dEtmp.tolist()])
                print("")
                print(f"MOPAC {Etmp:12.7f}  x: {tmp}")
                tmp = " ".join([f"{v:9.6f}" for v in dE.tolist()])
                print(f" Psi4 {E:12.7f}  x: {tmp}")
                print("")

            max_dE = abs_max(dE)
            rms_dE = rms(dE)
            if step == 0:
                print(f"{step:4d} {E:12.7f} {max_dE:9.6f} {rms_dE:9.6f}")
            else:
                delta = xlast - x
                max_dx = abs_max(delta)
                rms_dx = rms(delta)
                print(
                    f"{step:4d} {E:12.7f} {max_dE:9.6f} {rms_dE:9.6f} "
                    f"{max_dx:9.6f} {rms_dx:9.6f}"
                )

            xlast = np.array(x)

            if True:
                x = self._optimize(x, E, dE, A)
            else:
                direction = np.matmul(A, -dE)
                stepsize, Enew, x = self._linesearch(
                    x, E, dE, direction, self._surrogate
                )

    def _optimize(self, x, E, dE, A):
        n_dof = len(self.x0)
        if self.is_linear:
            n_dof -= 5
        else:
            n_dof -= 6
        n_steps = 3 * n_dof
        xlast = None
        for step in range(0, n_steps + 1):
            tmp = " ".join([f"{v:9.6f}" for v in x.tolist()])
            logger.debug(f"coordinates: {tmp}")

            tmp = " ".join([f"{v:9.6f}" for v in x.tolist()])
            print(f"\t\t\t\t\t         x: {tmp}")

            if step >= 0:
                E, dE, d2E = self._surrogate(x, ("E", "dE", "d2E"))

                # Invert the Hessian carefully
                s = np.linalg.svd(d2E, compute_uv=False, hermitian=True)
                if self.is_linear:
                    rcond = ((s[-6] + s[-5]) / 2) / s[0]
                else:
                    rcond = ((s[-7] + s[-6]) / 2) / s[0]

                A = np.linalg.pinv(d2E, rcond=rcond, hermitian=True)

            if False and step > 0:
                E0, dE0 = self._energy_callback(x, ("E", "dE"))
                E = E0
                dE = dE0

                tmp = " ".join([f"{v:9.6f}" for v in dE.tolist()])
                print(f"\t\t\t\t\tMOPAC {E:12.7f}  x: {tmp}")
                tmp = " ".join([f"{v:9.6f}" for v in dE0.tolist()])
                print(f"\t\t\t\t\t Psi4 {E0:12.7f}  x: {tmp}")

            tmp = " ".join([f"{v:9.6f}" for v in dE.tolist()])
            logger.debug(f"         dE: {tmp}")

            max_dE = abs_max(dE)
            rms_dE = rms(dE)
            if step == 0:
                print(f"\t{step:4d} {E:12.7f} {max_dE:9.6f} {rms_dE:9.6f}")
            else:
                delta = xlast - x
                max_dx = abs_max(delta)
                rms_dx = rms(delta)
                print(
                    f"\t{step:4d} {E:12.7f} {max_dE:9.6f} {rms_dE:9.6f} "
                    f"{max_dx:9.6f} {rms_dx:9.6f}"
                )

            if max_dE < 1.0e-04 and rms_dE < 1.0e-04:
                break

            xlast = np.array(x)

            direction = np.matmul(A, -dE)
            # direction = -dE

            tmp = " ".join([f"{v:9.6f}" for v in dE.tolist()])
            print(f"\t\t\t\t\t          dE: {tmp}")
            tmp = " ".join([f"{v:9.6f}" for v in direction.tolist()])
            print(f"\t\t\t\t\t   direction: {tmp}")

            direction, residuals, rank, s = np.linalg.lstsq(d2E, -dE, rcond=rcond)

            print(f"{rank=}")
            tmp = " ".join([f"{v:9.6f}" for v in direction.tolist()])
            print(f"\t\t\t\t\t   direction: {tmp}")

            stepsize, Enew, x = self._linesearch(x, E, dE, direction, self._surrogate)
            logger.debug(f"stepsize = {stepsize:.4f}")
        logger.debug("returning x")
        logger.debug(x)
        tmp = " ".join([f"{v:9.6f}" for v in x.tolist()])
        print(f"\t\t\t\t\t         x: {tmp}")
        return x

    def _linesearch(self, x0, E0, dE0, direction, callback):
        stepsize = 1.0
        dE0line = np.dot(dE0, direction)
        while True:
            x = x0 + direction * stepsize
            E, dE = self._surrogate(x)
            dEline = np.dot(dE, direction)
            logger.debug(
                f"{stepsize=} {E0:9.6f} --> {E:9.6f}, {dE0line:9.6f} --> {dEline:9.6f}"
            )
            if abs(dEline) < 1.0e-8:
                break

            a = abs(dE0line)
            b = abs(dEline)
            if dE0line < 0:
                if dEline > 0:
                    if a > b:
                        delta = 1 - b / a
                    else:
                        delta = a / b
                else:
                    if b > a:
                        # Stepped over a maximum?
                        raise RuntimeError("Y > X stepped over a maximum?")
                    else:
                        delta = 1 + b / a
            else:
                raise RuntimeError("Line search is backwards?")
            logger.debug(f"{delta=:9.6f}")
            if abs(delta - 1) < 0.1:
                break
            stepsize *= delta
        return stepsize, E, x

    def _surrogate(self, x, needed=("E", "dE")):
        if "d2E" in needed:
            Es, dEs, d2E = self._surrogate_callback(x, needed)
        else:
            Es, dEs = self._surrogate_callback(x, needed)

        E = Es + self.Ecorr + np.dot(x - self.x0, self.dEcorr)
        dE = dEs + self.dEcorr

        if "d2E" in needed:
            return E, dE, d2E
        else:
            return E, dE

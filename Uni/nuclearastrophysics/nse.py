# ---------------------
# Nuclear Statistical Equilibrium Solver
# Author: Gaétan Jalin
# Date: March 2025
# ---------------------

from typing import Tuple, List, Dict, Union, Optional
from pathlib import Path
import os

import numpy as np
from nse_solver.winvn_database import WinvnDatabase, Nuclide, _HTPF_TEMP_GRID


_WINVN = None

MP = 938.272046        # Proton mass MeV/c2
ME = 0.5110            # Electron mass MeV/c2
MN = 939.565379        # Neutron mass MeV/c2
AMU = 931.494061       # Atomic mass unit MeV/c2
PEXC = (MP + ME - AMU) # Proton mass excess
NEXC = (MN - AMU)      # Neutron mass excess
CLIGHT = 2.99792458e10 # Speed of light cm/s
AVO = 6.0221367e23     # Avogadro number
H = 6.6260755e-27      # Planck constant erg.s
H_EV = 4.135667696e-15 # Planck constant eV.s
K_ERG = 1.380658e-16   # Boltzmann constant erg/K
K_EV = 8.6173e-5       # Boltzmann constant eV/K
K_MEV = K_EV*1e-6      # Boltzmann constant MeV/K
SQRT2 = np.sqrt(2.0)

SAHA_CONST = 3.2047203586e3 # This is what you get after reducing all constants in saha eq.

NSE_MAX_NR_IT = 10000
NSE_NR_TOL = 1e-5
NSE_SMALL = 1e-10

# Screen parameters from Juodogalvis+2010
SC_COEFF_A     = -0.898004
SC_COEFF_B     = 0.96786
SC_COEFF_C     = 0.220703
SC_COEFF_D     = -0.86097
SC_COEFF_E     = -2.52692
SC_COEFF_BETA  = 0.295614
SC_COEFF_GAMMA = 1.98848


class NSESolver(object):
    _nuclides: Dict[str, Nuclide]
    _screening: bool

    _Yi: np.ndarray
    _Yp: float
    _Yn: float

    _As: np.ndarray
    _Zs: np.ndarray
    _Ns: np.ndarray
    _Bs: np.ndarray
    _PFs: np.ndarray
    _ln_consts: np.ndarray
    _ln_dists: np.ndarray
    _ln_saha_factors: np.ndarray
    
    def __init__(self, screening: bool = True):
        global _WINVN

        if _WINVN is None:
            wd = Path(__file__).parent
            winvn_path = str(os.path.join(wd, 'winvne_v2.0.dat'))
            htpf_path = str(os.path.join(wd, 'htpf.dat'))
            _WINVN = WinvnDatabase(winvn_path=winvn_path, htpf_path=htpf_path)

        # Setup network
        self._nuclides = _WINVN._nuclides.copy()
        self._screening = screening

        self._As    = np.atleast_1d(np.array([nuc.A    for nuc in self._nuclides.values()]))
        self._Zs    = np.atleast_1d(np.array([nuc.Z    for nuc in self._nuclides.values()]))
        self._Ns    = np.atleast_1d(np.array([nuc.N    for nuc in self._nuclides.values()]))
        self._spins = np.atleast_1d(np.array([nuc.spin for nuc in self._nuclides.values()]))
        self._Bs    = np.atleast_1d(np.array([nuc.B    for nuc in self._nuclides.values()]))
        self._PFs   = np.atleast_1d(np.stack([nuc.pf   for nuc in self._nuclides.values()]))

        self._gg = 2.*self._spins + 1.
        self._Gs = np.ones(len(self._nuclides))

        self._Yn = 0.5
        self._Yp = 0.5
        self._Yi = np.zeros(len(self._nuclides))

    def solve(self, dens: float, t9: float, Ye: float, reset_solver=True):
        # Set initial conditions
        if reset_solver or self._Yn == 0.5:
            Yn = 1.0 - Ye
            Yp = Ye
        else:
            Yn = self._Yn
            Yp = self._Yp

        # Initialise constants of saha equation
        self._init_solver(dens, t9)

        # Start Newton-Raphson
        for k in range(NSE_MAX_NR_IT):
            oldYn = Yn
            oldYp = Yp

            Yn, Yp, converged = self._solve_nse(oldYn, oldYp, dens, t9, Ye)

            if converged:
                break
        else:
            raise RuntimeError('Too many iteration in Newton-Raphson')

        # Calculate final abundances
        Yi = self._saha(Yn, Yp, dens, t9, Ye)

        self._Yn = Yn
        self._Yp = Yp
        self._Yi[0] = Yn
        self._Yi[1] = Yp
        self._Yi[2:] = Yi

    def getYi(self) -> np.ndarray:
        return self._Yi

    def getYn(self) -> float:
        return self._Yn

    def getYp(self) -> float:
        return self._Yp

    def getYalpha(self) -> float:
        return self._Yi[5] # he4

    def getYlight(self) -> float:
        return self._Yi[2] + self._Yi[3] + self._Yi[3] # d + t + he3

    def getYheavy(self) -> float:
        return np.sum(self._Yi[6:])

    def getMassFractions(self) -> np.ndarray:
        return self._Yi*self._As

    def getXi(self) -> np.ndarray:
        return self._Yi*self._As

    def getXn(self) -> float:
        return self._Yn

    def getXp(self) -> float:
        return self._Yp

    def getXalpha(self) -> float:
        return self._Yi[5]*4.0

    def getXlight(self) -> float:
        return self._Yi[2]*2.0 + self._Yi[3]*3.0 + self._Yi[3]*3.0

    def getXheavy(self) -> float:
        return np.sum(self._Yi[6:]*self._As[6:])

    def getAbar(self) -> float:
        return 1./np.sum(self._Yi)

    def getZbar(self) -> float:
        return np.sum(self._Yi*self._Zs) / np.sum(self._Yi)

    def getAbarHeavy(self) -> float:
        return np.sum(self._Yi[6:]*self._As[6:]) / np.sum(self._Yi[6:])

    def getZbarHeavy(self) -> float:
        return np.sum(self._Yi[6:]*self._Zs[6:]) / np.sum(self._Yi[6:])

    # Screening potential from Dirk Martin
    def _muC(self, dens: float, t9: float, Ye: float, Z: float) -> Union[float, np.ndarray]:
        aekBT = 32.080407637714124 * (Ye*dens)**(-1./3.) * t9
        Gamma = np.power(Z, 5./3.) * (1.44/197.3269718)/aekBT
        fC = np.where(Gamma > 1.0, \
            SC_COEFF_A*Gamma + 4.0*SC_COEFF_B*np.power(Gamma, 1./4.) - 4.0*SC_COEFF_C*np.power(Gamma, -1./4.) + SC_COEFF_D*np.log(Gamma) + SC_COEFF_E, \
            -(1./np.sqrt(3.))*np.power(Gamma, 1.5) + (SC_COEFF_BETA/SC_COEFF_GAMMA)*np.power(Gamma, SC_COEFF_GAMMA))
        return fC.item() if fC.size == 1 else fC

    def _init_solver(self, dens, t9):
        if (t9 <= _HTPF_TEMP_GRID[0]):
            self._Gs = self._gg * self._PFs[:,0]
        elif (t9 >= _HTPF_TEMP_GRID[-1]):
            self._Gs = self._gg * self._PFs[:,-1]
        else:
            t9i = np.searchsorted(_HTPF_TEMP_GRID, t9)

            t9_left, t9_right = _HTPF_TEMP_GRID[t9i-1], _HTPF_TEMP_GRID[t9i]
            pfs_left, pfs_right = self._PFs[:,t9i-1], self._PFs[:,t9i]

            # linearly interpolate partition functions
            #slopes = (pfs_right - pfs_left) / (t9_right - t9_left)
            #self._Gs = self._gg * (slopes*(t9 - t9_right) + pfs_right)

            # lin-log interpolate partition functions
            slope = (t9 - t9_left) / (t9_right - t9_left)
            self._Gs = self._gg * (pfs_left * np.power(pfs_right/pfs_left, slope))

        self._ln_consts = (self._As - 1.0)*np.log(SAHA_CONST * dens * (t9*1e9)**(-3./2.)) \
                         + 1.5*np.log(self._As) \
                         - self._As*np.log(2.0)
        self._ln_dists  = self._Bs/(K_MEV*t9*1e9)

        self._ln_saha_factors = np.log(self._Gs) + self._ln_consts + self._ln_dists

    def _solve_nse(self, Yn: float, Yp: float, dens: float, t9: float, Ye: float) -> Tuple[float, float, bool]:
        # Calculate abundances
        Yi = self._saha(Yn, Yp, dens, t9, Ye)

        # Charge conservation
        c = Yp + np.sum(Yi * self._Zs[2:])

        # Filter insane values
        if (np.log10(abs(c)) > 10.0):
            if (Ye > 0.4):
                # Need more neutrons in nuclei
                return Yn/SQRT2, Yp, False
            else:
                # Need more protons in nuclei
                return Yn, Yp/SQRT2, False

        c = Ye - c
        # Mass conservation
        m = Yn + Yp + np.sum(Yi * self._As[2:])
        m = 1.0 - m

        # Calculate derivatives and solve Jacobian
        dmdYn = np.sum(Yi*self._As[2:]*self._Ns[2:])/Yn + 1.0
        dmdYp = np.sum(Yi*self._As[2:]*self._Zs[2:])/Yp + 1.0
        dcdYn = np.sum(Yi*self._Zs[2:]*self._Ns[2:])/Yn
        dcdYp = np.sum(Yi*self._Zs[2:]*self._Zs[2:])/Yp + 1.0

        detJ = dmdYn*dcdYp - dcdYn*dmdYp
        delYn = 1./detJ * (dcdYp*m - dmdYp*c)
        delYp = 1./detJ * (dmdYn*c - dcdYn*m)

        # Convergence
        #if (max(abs(m), abs(c)) < NSE_NR_TOL and max(abs(delYn), abs(delYp)) < NSE_NR_TOL):
        if (np.sqrt(delYn**2.0 + delYp**2.0 + m**2.0 + c**2.0) < NSE_NR_TOL):
            return Yn, Yp, True

        # Constrain abundances to positive values
        if (delYn > 0.0 or Yn > abs(delYn)):
            Yn += delYn
        if (delYp > 0.0 or Yp > abs(delYp)):
            Yp += delYp

        return Yn, Yp, False

    def _saha(self, Yn: float, Yp: float, dens: float, t9: float, Ye: float):
        # Update abundances
        lnYn = np.log(Yn)
        lnYp = np.log(Yp)

        # Screening corrections
        if self._screening:
            muCp = self._muC(dens, t9, Ye, 1.0)
            muCi = self._muC(dens, t9, Ye, self._Zs[2:])
            return np.exp(self._ln_saha_factors[2:] + self._Ns[2:]*lnYn + self._Zs[2:]*lnYp + self._Zs[2:]*muCp - muCi)
        else:
            return np.exp(self._ln_saha_factors[2:] + self._Ns[2:]*lnYn + self._Zs[2:]*lnYp)

        # Prevent large numbers
        # This slows down A LOT
        #for i in range(2, len(self._Yi)):
        #    if ln_saha[i] < -80.0:
        #        self._Yi[i] = 0.0
        #    elif ln_saha[i] > 80.0:
        #        self._Yi[i] = 100.0
        #    else:
        #        self._Yi[i] = np.exp(ln_saha[i])


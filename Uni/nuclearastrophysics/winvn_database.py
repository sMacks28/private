from __future__ import annotations
from typing import Union, Optional, List, Dict
from dataclasses import dataclass, field

import numpy as np


# TODO Move constants someplace generic
MP = 938.272046        # Proton mass MeV/c2
ME = 0.5110            # Electron mass MeV/c2
MN = 939.565379        # Neutron mass MeV/c2
AMU = 931.494061       # Atomic mass unit MeV/c2
PEXC = (MP + ME - AMU) # Proton mass excess
NEXC = (MN - AMU)      # Neutron mass excess
CLIGHT = 2.99792458e10 # Speed of light cm/s
AVO = 6.0221367e23     # Avogadro number
H = 6.6260755e-27      # Planck constant erg.s
KERG = 1.380658e-16    # Boltzmann constant cgs

# Temperature grid in GK for the partition functions
_WINVN_TEMP_GRID = np.array( \
    [1e-1, 1.5e-1, 2e-1, 3e-1, 4e-1, 5e-1, 6e-1, 7e-1, 8e-1, 9e-1, \
     1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0] \
)
_HTPF_TEMP_GRID = np.array( \
    [1e-1, 1.5e-1, 2e-1, 3e-1, 4e-1, 5e-1, 6e-1, 7e-1, 8e-1, 9e-1, \
     1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, \
     12.0, 14.0, 16.0, 18.0, 20.0, 22.0, 24.0, 26.0, 28.0, 30.0, 35.0, 40.0, 45.0, 50.0, \
     55.0, 60.0, 65.0, 70.0, 75.0, 80.0, 85.0, 90.0, 95.0, 100.0, 105.0, 110.0, 115.0, 120.0, \
     125.0, 130.0, 135.0, 140.0, 145.0, 150.0, 155.0, 160.0, 165.0, 170.0, 175.0, 180.0, 190.0, \
     200.0, 210.0, 220.0, 230.0, 240.0, 250.0, 275.0] \
)


@dataclass(frozen=True)
class Nuclide(object):
    name: str
    A: float
    Z: float
    N: float
    spin: float
    mass_excess: float
    pf: np.ndarray               # Partition function
    B: float = field(init=False) # Binding energy

    def __post_init__(self):
        B = self.Z*PEXC + self.N*NEXC - self.mass_excess
        object.__setattr__(self, 'B', B)

    def update_pf(self, new_pf) -> None:
        object.__setattr__(self, 'pf', new_pf)

    def eval_pf(self, t9: float) -> float:
        # lin-log interpolation
        if len(self.pf) == len(_WINVN_TEMP_GRID):
            return np.exp(np.interp(t9, _WINVN_TEMP_GRID, np.log(self.pf)))
        elif len(self.pf) == len(_HTPF_TEMP_GRID):
            return np.exp(np.interp(t9, _HTPF_TEMP_GRID, np.log(self.pf)))
        else:
            raise RuntimeError('Unknown partition function')


class WinvnDatabase(object):
    _nuclides: Dict[str, Nuclide]

    def __init__(self, winvn_path: str, htpf_path: Optional[str] = None):
        self._read_winvn_file(winvn_path)
        if htpf_path:
            self._read_htpf_file(htpf_path)

    def __contains__(self, key: str) -> bool:
        return key in self._nuclides

    def __getitem__(self, key: str) -> Nuclide:
        return self._nuclides[key]

    def __len__(self) -> int:
        return len(self._nuclides)

    def __bool__(self) -> bool:
        return self._nuclides.__bool__()

    def _read_winvn_file(self, winvn_path: str) -> None:
        self._nuclides = {}
        with open(winvn_path, 'r') as winvn_file:
            # Skip first 2 lines
            line = next(winvn_file)
            line = next(winvn_file)

            # Skip list of nuclides
            for line in winvn_file:
                # Skip comments, empty lines, and initial list of nuclides
                if (line.startswith('#') or len(line.rstrip()) <= 6):
                    continue

                # Reached database entry
                # Nuclide properties
                header = line.split()
                name        = str(header[0])
                A           = float(header[1])
                Z           = float(header[2])
                N           = float(header[3])
                spin        = float(header[4])
                mass_excess = float(header[5])


                # Partition function
                pf = np.loadtxt([next(winvn_file) for i in range(3)])
                pf = np.concatenate(pf)

                # Store nuclide entry
                self._nuclides[name] = Nuclide(name, A, Z, N, spin, mass_excess, pf)

    def _read_htpf_file(self, htpf_path: str) -> None:
        with open(htpf_path, 'r') as htpf_file:
            for line in htpf_file:
                # Skip comments and empty lines
                if (line.startswith('#') or len(line.rstrip()) < 1):
                    continue

                # Read entry
                entry = line.split()
                name        = str(entry[0])
                #Z           = float(entry[1])
                #A           = float(entry[2])
                #N           = A - Z
                #spin        = float(entry[3])
                pf          = np.float64(entry[4:])

                # Update corresponding partition function
                pf = np.concatenate((self._nuclides[name].pf, pf))
                self._nuclides[name].update_pf(pf)

        # Extend other nuclides as well for vectorised interpolation
        for nuc in self._nuclides.values():
            if (len(nuc.pf) < len(_HTPF_TEMP_GRID)):
                pf = np.concatenate((nuc.pf, np.full(len(_HTPF_TEMP_GRID) - len(nuc.pf), nuc.pf[-1])))
                nuc.update_pf(pf)


""" An interface for the [MACE](https://github.com/ACEsuit/mace) calculator """

import json
import sys
# from .ase import ASEDriver
from .dummy import Dummy_driver
from ase.io import read
import numpy as np
from ipi.utils.units import unit_to_internal, unit_to_user
try:
    from apax.md import ASECalculator
except:
    ASECalculator = None
import time

__DRIVER_NAME__ = "apax"
__DRIVER_CLASS__ = "APAX_driver"

ERROR_MSG = """
APAX driver requires specification of a model checkpoint path,
and a template file that describes the chemical makeup of the structure.

Example: python driver.py -m apax -u -o directory/experiment,template.xyz
"""


class APAX_driver(Dummy_driver):
    def __init__(self, args=None, verbose=False):
        if ASECalculator is None:
            raise ImportError("Could not load apax bindings")
        super().__init__(args)

    def check_arguments(self):
        """Check the arguments required to run the driver

        This loads the potential.
        """

        if len(self.args) != 2:
            sys.exit(self.error_msg)

        model_dir = self.args[0]
        template = self.args[1]
        self.atoms = read(template, index="0")

        self.ase_calculator = ASECalculator(model_dir=model_dir)
        self.atoms.calc = self.ase_calculator

    def __call__(self, cell, pos):
        pos = unit_to_user("length", "angstrom", pos)
        cell = unit_to_user("length", "angstrom", cell.T)

        self.atoms.set_positions(pos)
        self.atoms.set_cell(cell)

        self.atoms.get_potential_energy()

        results = self.ase_calculator.results
        
        is_ens =  "energy_ensemble" in results.keys()
        has_stress = "stress" in results.keys()

        energy = results["energy"]
        force = results["forces"]

        extras = ""

        if is_ens:
            committee_e = results["energy_ensemble"]
            committee_e = unit_to_internal("energy", "electronvolt", committee_e)
            committee_f = results["forces_ensemble"]
            committee_f = unit_to_internal("force", "ev/ang", committee_f)
            extras = {"committee_pot":committee_e.tolist(),
                    "committee_force":committee_f.tolist(),
            }
        
        if has_stress:
            virials = results["stress"]
        else:
            virials = np.zeros((3,3))        
        
        pot_ipi = unit_to_internal("energy", "electronvolt", energy)
        force_ipi = unit_to_internal("force", "ev/ang", force.reshape(-1, 3))
        vir_ipi = unit_to_internal("energy", "electronvolt", virials.T)

        if extras != "":
            extras = json.dumps(extras)

        return pot_ipi, force_ipi, vir_ipi, extras
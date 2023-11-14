# This code is part of kartograf and is licensed under the MIT license.
# For details, see https://github.com/OpenFreeEnergy/kartograf
import mpmath as mp
from scipy import constants as const

import logging
import numpy as np

from rdkit.Chem import AllChem
from gufe.mapping import AtomMapping

from ._abstract_scorer import _AbstractAtomMappingScorer

logger = logging.getLogger(__name__)

class FFScorer(_AbstractAtomMappingScorer):

    def __init__(self, emin: bool = True, force_constant: float = 1.0e4, max_emin_iter: int = 1000, n_measurements: int = 37):
        """
        This class is a very simple input generation Variant.
            Don't use
        Parameters
        ----------
        n_measurements : int, optional
            points along the rotation (36 translates to 36 points between -180,180), by default 36
        emin : bool, optional
            remove severe clashes, by default True
        force_constant : float, optional
            constant strength for the position restraints during emin, by default 1.e4
        maxEminIt : int, optional
            maximal emin steps., by default 1000

        """

        self.emin = emin
        self.n_measurements = n_measurements
        self.force_constant = force_constant
        self.maxEminIt = max_emin_iter
        self.beta = ((298*const.Avogadro*const.k) / 4.184)

    def get_score(self, mapping: AtomMapping) -> float:
        """
            Under Development

        Parameters
        ----------
        mapping : AtomMapping
            AtomMapping to be scored

        Returns
        -------
        float
            normalized score
        """

        logger.info("FF Scorer:")
        rd_molA = mapping.componentA.to_rdkit()
        rd_molB = mapping.componentB.to_rdkit()
        m_dict = mapping.componentA_to_componentB

        #Absolute Energy:
        mpA = AllChem.MMFFGetMoleculeProperties(rd_molA)
        ffA = AllChem.MMFFGetMoleculeForceField(rd_molA, mpA)
        mpB = AllChem.MMFFGetMoleculeProperties(rd_molB)
        ffB = AllChem.MMFFGetMoleculeForceField(rd_molB, mpB)

        ffA.Minimize(maxIts=self.maxEminIt)
        VA_abs = ffA.CalcEnergy()

        ffB.Minimize(maxIts=self.maxEminIt)
        VB_abs = ffB.CalcEnergy()

        dVAB_abs = np.abs(VB_abs - VA_abs)
        logger.info(f"Absolute Diff A to B: dVAB_abs = VB_abs - VA_abs: "
                    f" {dVAB_abs}={VB_abs} - {VA_abs}")

        # Mapping Transformatoins A->B
        # Juggle Coords
        self.rd_molA_mapB = mapping.componentA.to_rdkit()
        self.rd_molB_mapA  = mapping.componentB.to_rdkit()

        rd_molA_mapB_conf = self.rd_molA_mapB.GetConformer()
        rd_molB_mapA_conf = self.rd_molB_mapA.GetConformer()

        rd_molA_mapB_pos = np.array(rd_molA_mapB_conf.GetPositions())
        rd_molB_mapA_pos = np.array(rd_molB_mapA_conf.GetPositions())

        for a_id, b_id in m_dict.items():
            rd_molB_mapA_conf.SetAtomPosition(b_id, rd_molA_mapB_pos[a_id])
            rd_molA_mapB_conf.SetAtomPosition(a_id, rd_molB_mapA_pos[b_id])

        # Mapping Transformation A->B:
        mpA = AllChem.MMFFGetMoleculeProperties(self.rd_molA_mapB)
        ffA = AllChem.MMFFGetMoleculeForceField(self.rd_molA_mapB, mpA)

        VA_map_start = ffA.CalcEnergy()
        ffA.Minimize(maxIts=self.maxEminIt)
        VA_map_end = ffA.CalcEnergy()
        dVA_map = VA_map_end - VA_map_start
        logger.info(f"Transform A to B: dVA_map = VA_map_end - VA_map_start:"
                    f" {dVA_map}={VA_map_end} - {VA_map_start}")


        # Mapping Transformation B->A:
        mpB = AllChem.MMFFGetMoleculeProperties(self.rd_molB_mapA )
        ffB = AllChem.MMFFGetMoleculeForceField(self.rd_molB_mapA , mpB)

        VB_map_start = ffB.CalcEnergy()
        ffB.Minimize(maxIts=self.maxEminIt)
        VB_map_end = ffB.CalcEnergy()

        dVB_map = VB_map_end-VB_map_start
        logger.info(f"Transform B to A: dVB_map = VB_map_end - VB_map_start:"
                    f" {dVB_map}={VB_map_end} - {VB_map_start}")

        # Final Result:
        dVAB_map = np.abs(VB_map_end - VA_map_end)

        dE = dVAB_abs + dVAB_map
        logger.info(f"dE=dVAB_abs + dVA_map + dVB_map: {dE}={dVAB_abs} + {dVAB_map}")

        p = float(mp.exp(-self.beta*dE))
        score = np.round(p, 2)

        logger.info("Result: " + str(score))

        return score

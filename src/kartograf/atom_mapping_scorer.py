# This code is part of OpenFE and is licensed under the MIT license.
# For details, see https://github.com/OpenFreeEnergy/kartograf

import logging
import numpy as np

from gufe.mapping import AtomMapping

from .mapping_metrics import MappingRMSDScorer
from .mapping_metrics import MappingShapeOverlapScorer, MappingShapeMismatchScorer
from .mapping_metrics import MappingVolumeRatioScorer

from .mapping_metrics._abstract_scorer import _AbstractAtomMappingScorer

log = logging.getLogger(__name__)


class DefaultKartografScorer(_AbstractAtomMappingScorer):
    """
    Warning this is highly experimental!
    """
    def __init__(self):
        """
        this could be a future scorer using the here derined metrics?

        don't use ;)
        """
        self.scorers = [MappingVolumeRatioScorer(),
                        MappingShapeOverlapScorer(),
                        MappingShapeMismatchScorer(),
                        MappingRMSDScorer()]

        self.weights = np.array([1,3,3,3])
        self.weights = self.weights/np.sum(self.weights)

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
        log.info("Kartograf Score:")
        scores = []
        for weight, scorer in zip(self.weights, self.scorers):
            s = scorer.get_score(mapping)
            log.info("\t"+scorer.__class__.__name__+"\t"+str(s)+"\tweight: "+str(weight))
            scores.append(s*weight)

        score = np.round(np.mean(scores),2)
        score = score if(score<1.0) else 1.0
        log.info("Result: "+str(score))
        return score

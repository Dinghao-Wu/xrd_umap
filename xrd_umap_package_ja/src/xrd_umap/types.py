from dataclasses import dataclass
import numpy as np
from typing import Dict, Any
@dataclass
class XRDRecord: two_theta: np.ndarray; intensity: np.ndarray; meta: Dict[str, Any]

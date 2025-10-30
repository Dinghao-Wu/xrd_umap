
import numpy as np
def l2_normalize(v, eps: float=1e-12):
    # L2 正規化（既存仕様のまま）
    n=float(np.linalg.norm(v))
    if not np.isfinite(n) or n<=eps: return None
    return v/n

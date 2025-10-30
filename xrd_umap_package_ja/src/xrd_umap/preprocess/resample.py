
import math, numpy as np
def build_grid(xmin: float, xmax: float, step: float)->np.ndarray:
    # 共通グリッドの生成（既存仕様のまま）
    if xmax<=xmin: raise ValueError("--xmax must be > --xmin")
    if step<=0: raise ValueError("--step must be > 0")
    n=int(math.floor((xmax-xmin)/step))+1
    return xmin+np.arange(n)*step
def resample_to_grid(x,y,grid):
    # 線形補間での再サンプリング（既存仕様のまま）
    return np.interp(grid,x,y,left=0.0,right=0.0)

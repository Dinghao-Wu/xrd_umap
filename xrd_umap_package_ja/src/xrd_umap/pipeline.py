
from __future__ import annotations
import os, sys, argparse, numpy as np, pandas as pd
from typing import Optional, List
from .io.readers import discover_files, read_any_xrd, XYLIB_OK
from .preprocess.resample import build_grid, resample_to_grid
from .preprocess.normalize import l2_normalize

def run_dir(input_dir: str, xmin: float, xmax: float, step: float)->int:
    if not os.path.isdir(input_dir):
        sys.stderr.write(f"[ERR] Not a directory: {input_dir}\n"); return 1
    try: grid=build_grid(xmin, xmax, step)
    except Exception as e: sys.stderr.write(f"[ERR] Grid error: {e}\n"); return 1

    files=discover_files(input_dir)
    if not files:
        sys.stderr.write("[ERR] No candidate XRD files.\n"); return 1

    rows=[]; labels=[]; skipped=[]
    for fp in files:
        base=os.path.basename(fp)
        try:
            x,y=read_any_xrd(fp)
            yi=resample_to_grid(x,y,grid)
            yi=l2_normalize(yi)
            if yi is None: raise ValueError("norm==0")
        except Exception as e:
            skipped.append((base, str(e))); continue
        rows.append(yi); labels.append(base)

    if not rows:
        sys.stderr.write("[ERR] No valid patterns.\n"); return 1
    X=np.vstack(rows)

    patterns_csv=os.path.join(os.getcwd(),"calcd_patterns.csv")
    np.savetxt(patterns_csv, X, delimiter=",")

    import pandas as pd
    tgt_df=pd.DataFrame({"file": labels, "label": list(range(len(labels))) })
    targets_csv=os.path.join(os.getcwd(), "targets.csv")
    tgt_df.to_csv(targets_csv, index=False)

    print(f"[OK] files={len(labels)} skipped={len(skipped)} gridpts={len(grid)}")
    if skipped:
        for n,r in skipped[:10]: print(f"   [skip] {n}: {r}")
        if len(skipped)>10: print(f"   ... and {len(skipped)-10} more")
    print(f"[OK] calcd_patterns -> {patterns_csv}  (numeric-only)")
    print(f"[OK] targets -> {targets_csv}")
    if not XYLIB_OK:
        print("[NOTE] xylib-py not installed; vendor formats (.ras/.uxd/.raw/...) may be skipped.")
    return 0

def main(argv: Optional[List[str]] = None)->int:
    ap=argparse.ArgumentParser(description="XRD -> dataset (common grid + L2 normalize)")
    ap.add_argument("input_dir")
    ap.add_argument("--xmin", type=float, default=10.0)
    ap.add_argument("--xmax", type=float, default=80.0)
    ap.add_argument("--step", type=float, default=0.02)
    args=ap.parse_args(argv)
    return run_dir(args.input_dir, args.xmin, args.xmax, args.step)

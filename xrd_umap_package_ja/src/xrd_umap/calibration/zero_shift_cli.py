
from __future__ import annotations
import sys, pandas as pd
from . import urlap_report_v3 as U

def zeroshift_entry(args)->int:
    if args.mode in ("fixed","scan"):
        argv = ["-i", args.input]
        if args.mode == "fixed":
            argv += ["--mode", "fixed"]
            if args.wavelength is not None: argv += ["--wavelength", str(args.wavelength)]
            if args.delta is not None: argv += ["--delta", str(args.delta)]
            if args.out: argv += ["--out", args.out]
        else:
            argv += ["--mode", "scan"]
            if args.wavelength is not None: argv += ["--wavelength", str(args.wavelength)]
            if args.delta_range: argv += ["--delta-range", str(args.delta_range[0]), str(args.delta_range[1])]
            argv += ["--step", str(args.step), "--select-by", args.select_by]
            if args.out: argv += ["--out", args.out]
        sys.argv = ["urlap_report_v3.py"] + argv
        try:
            U.main()
            return 0
        except SystemExit as e:
            return int(e.code)
    elif args.mode == "apply-curve":
        df = pd.read_csv(args.in_csv)
        if "2theta" not in df.columns or "intensity" not in df.columns:
            raise SystemExit("apply-curve expects a CSV with columns: 2theta,intensity")
        df["2theta"] = df["2theta"].astype(float) - float(args.delta)  # OBS_shown = OBS_raw - delta
        df.to_csv(args.out_csv, index=False)
        print(f"[OK] wrote {args.out_csv}")
        return 0
    else:
        print("Unknown zeroshift mode")
        return 2

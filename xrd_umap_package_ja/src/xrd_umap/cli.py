
from __future__ import annotations
import argparse, sys
from .pipeline import main as dataset_main
from .calibration.zero_shift_cli import zeroshift_entry

def main():
    parser = argparse.ArgumentParser(prog="xrd-umap", description="XRD ツール（データセット生成 + 任意のゼロ点補正）")
    sub = parser.add_subparsers(dest="cmd", required=True)

    p = sub.add_parser("dataset", help="ディレクトリからデータセットを生成（既存仕様と等価）")
    p.add_argument("input_dir")
    p.add_argument("--xmin", type=float, default=10.0)
    p.add_argument("--xmax", type=float, default=80.0)
    p.add_argument("--step", type=float, default=0.02)

    z = sub.add_parser("zeroshift", help="任意のゼロ点補正ツール（既知 hkl のみ）")
    z_sub = z.add_subparsers(dest="mode", required=True)

    # fixed
    zf = z_sub.add_parser("fixed", help="Δ を固定してレポート生成")
    zf.add_argument("-i","--input", required=True)
    zf.add_argument("--wavelength", type=float, default=None)
    zf.add_argument("--delta", type=float, default=None)
    zf.add_argument("--out", type=str, default=None)

    # scan
    zs = z_sub.add_parser("scan", help="Δ を走査して最適化しレポート生成")
    zs.add_argument("-i","--input", required=True)
    zs.add_argument("--wavelength", type=float, default=None)
    zs.add_argument("--delta-range", nargs=2, type=float, default=[-0.5, 0.5])
    zs.add_argument("--step", type=float, default=0.1)
    zs.add_argument("--select-by", choices=["rms","uncprod"], default="uncprod")
    zs.add_argument("--out", type=str, default=None)

    # apply-curve
    za = z_sub.add_parser("apply-curve", help="Δ を 2 列 XRD 曲線 (2theta,intensity) に適用")
    za.add_argument("--in-csv", required=True)
    za.add_argument("--out-csv", required=True)
    za.add_argument("--delta", type=float, required=True)

    args = parser.parse_args()

    if args.cmd == "dataset":
        code = dataset_main([args.input_dir, "--xmin", str(args.xmin), "--xmax", str(args.xmax), "--step", str(args.step)])
        raise SystemExit(code)
    elif args.cmd == "zeroshift":
        raise SystemExit(zeroshift_entry(args))


from __future__ import annotations
import os, re, glob
from typing import Callable, Tuple, Optional, List
import numpy as np

# --- Reader レジストリ（プラグイン） ---
ReaderFunc = Callable[[str], Tuple[np.ndarray, np.ndarray]]

class Reader:
    def __init__(self, name: str, loader: ReaderFunc, extensions: Tuple[str, ...] = (), priority: int = 100):
        self.name = name
        self.loader = loader
        self.extensions = tuple(e.lower() for e in extensions)
        self.priority = int(priority)

_REGISTRY: List[Reader] = []

def register_reader(reader: Reader):
    _REGISTRY.append(reader)
    _REGISTRY.sort(key=lambda r: r.priority)

def iter_readers()->List[Reader]:
    return list(_REGISTRY)

XYLIB_OK=False
try:
    import xylib  # type: ignore
    XYLIB_OK=True
except Exception:
    XYLIB_OK=False

COMMENT_PREFIXES=("#",";","!","*","//")

def is_number(s:str)->bool:
    try: float(s); return True
    except Exception: return False

# --- ASCII（2 列） ---
def read_ascii_two_columns(path: str):
    xs=[]; ys=[]
    with open(path,"r",encoding="utf-8",errors="ignore") as f:
        for line in f:
            s=line.strip()
            if not s: continue
            if s.startswith(COMMENT_PREFIXES): continue
            parts=re.split(r"[,\s;]+", s)
            if len(parts)<2: continue
            if not (is_number(parts[0]) and is_number(parts[1])): continue
            xs.append(float(parts[0])); ys.append(float(parts[1]))
    if not xs: raise ValueError("no numeric two-column rows")
    x=np.asarray(xs,float); y=np.asarray(ys,float)
    o=np.argsort(x); x=x[o]; y=y[o]
    uniq_x, idx=np.unique(x, return_index=True)
    return x[idx], y[idx]

# --- xylib（ベンダー形式） ---
def read_with_xylib(path: str):
    if not XYLIB_OK: raise ValueError("xylib not installed")
    lib=xylib.load_file(path)
    if lib is None or lib.get_block_count()<1: raise ValueError("xylib: no block")
    b=lib.get_block(0); n=b.get_point_count()
    xs=np.zeros(n,float); ys=np.zeros(n,float)
    for i in range(n):
        p=b.get_point(i); xs[i]=p.x; ys[i]=p.y
    return xs, ys

# --- XRDML ---
def read_xrdml(path: str):
    import xml.etree.ElementTree as ET
    tree=ET.parse(path); root=tree.getroot()
    dp=None
    for e in root.iter():
        if e.tag.endswith("dataPoints"): dp=e; break
    if dp is None: raise ValueError("XRDML: dataPoints missing")
    intens=None; xpos=None
    for c in dp.iter():
        if c.tag.endswith("intensities"):
            toks=((c.text or "").replace("\n"," ").split())
            intens=np.array([float(t) for t in toks], float)
        if c.tag.endswith("positions"):
            toks=((c.text or "").replace("\n"," ").split())
            if toks: xpos=np.array([float(t) for t in toks], float)
    if intens is None: raise ValueError("XRDML: intensities missing")
    if xpos is None:
        start=end=None
        for c in dp.iter():
            if c.tag.endswith("startPosition"): start=float(c.text)
            if c.tag.endswith("endPosition"): end=float(c.text)
        if start is not None and end is not None:
            n=len(intens); step=(end-start)/(n-1) if n>1 else 0.0
            xpos=start+np.arange(n)*step
    if xpos is None: raise ValueError("XRDML: positions missing")
    return np.asarray(xpos,float), np.asarray(intens,float)

ASCII_EXTS=(".txt",".csv",".xy",".xye",".dat")
XYLIB_EXTS=(".ras",".uxd",".raw",".rd",".cpi",".udf",".xdd")
XRDML_EXTS=(".xrdml",)

# 優先度（既存の解釈順に合わせる）
register_reader(Reader("ascii_two_columns", read_ascii_two_columns, extensions=ASCII_EXTS, priority=10))
register_reader(Reader("xrdml", read_xrdml, extensions=XRDML_EXTS, priority=20))
register_reader(Reader("xylib_vendor", read_with_xylib, extensions=XYLIB_EXTS, priority=30))

def discover_files(input_dir: str)->List[str]:
    # 登録された拡張子から候補を集める。なければディレクトリ内の全ファイルを対象にする。
    pats=[]
    for ext in ASCII_EXTS+XYLIB_EXTS+XRDML_EXTS:
        pats += [f"*{ext}", f"*{ext.upper()}"]
    files=[]
    for p in pats: files += glob.glob(os.path.join(input_dir,p))
    if not files:
        files=[os.path.join(input_dir,f) for f in os.listdir(input_dir) if os.path.isfile(os.path.join(input_dir,f))]
    return sorted(set(files))

def read_any_xrd(path: str):
    # 既存仕様に合わせて、ASCII → XRDML → xylib の順で試行し、必要に応じてフォールバック。
    ext=os.path.splitext(path)[1].lower()
    try: return read_ascii_two_columns(path)
    except Exception: pass
    if ext in XRDML_EXTS:
        try: return read_xrdml(path)
        except Exception: pass
    if XYLIB_OK:
        try: return read_with_xylib(path)
        except Exception: pass
    if XYLIB_OK and ext not in XYLIB_EXTS:
        try: return read_with_xylib(path)
        except Exception: pass
    # 最後のフォールバック：レジストリ順に試行
    for r in iter_readers():
        try:
            return r.loader(path)
        except Exception:
            continue
    raise ValueError("unrecognized or unreadable XRD format")

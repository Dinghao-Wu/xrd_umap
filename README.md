xrd-umap: XRDデータ前処理・解析ツール

1. 概要

このツールは、既存のXRDデータセット生成スクリプトをモジュール化したものです。元のスクリプトの振る舞いを維持しつつ、以下の機能を追加しています。

Readerプラグイン方式: 異なるファイル形式（ASCII, XRDML, ベンダー形式）を透過的に読み込みます。

ゼロ点補正 (オプション): 既知の hkl ピーク情報に基づき、回折パターンのゼロ点シフトを計算・補正します。

互換性

データセット生成:
生成される calcd_patterns.csv（数値のみの行列）および targets.csv（file→行番号の対応）の仕様は、元の xrd_to_dataset.py スクリプトと完全に同一です。

ゼロ点補正:
計算ロジックおよび peaks.txt の入出力フォーマットは、元の urlap_report_v3.py スクリプトの仕様をそのまま利用しています。

2. インストール

# 基本インストール
pip install -e .

# ベンダー固有形式 (.ras, .raw 等) のサポートが必要な場合:
pip install '.[vendor]'    # xylib-py を追加インストールします


3. 使用方法

3-1. データセット生成

指定されたディレクトリ内のXRDファイルを読み込み、共通の2θグリッドにリサンプリングし、L2正規化を行ったデータセットをCSVとして出力します。

xrd-umap dataset ./my_data_directory --xmin 10 --xmax 80 --step 0.02


出力ファイル:

calcd_patterns.csv: 数値のみの行列（ヘッダなし、インデックスなし）

targets.csv: file,label の形式（label は calcd_patterns.csv の行番号に対応）

（注） xylib-py がインストールされていない場合、ベンダー形式のファイルはスキップされます。ASCIIおよびXRDML形式は標準で読み込み可能です。

データフロー（データセット生成）:

flowchart TD
    A[入力ディレクトリ] --> B{Readerで読み込み}
    B -->|ASCII 2列 (.txt/.csv/.xy)| C[2θ, 強度]
    B -->|XRDML (.xrdml)| C
    B -->|xylib ベンダー (.ras/.raw)| C
    C --> D[共通グリッドを生成 (xmin, xmax, step)]
    D --> E[線形補間で再サンプリング]
    E --> F[L2 正規化]
    F --> G[calcd_patterns.csv (数値のみ)]
    F --> H[targets.csv (file→label)]


3-2. ゼロ点補正（オプション）

この機能は、既知の hkl 面指数を持つ回折ピーク（peaks.txt）が利用可能な場合に使用します。

peaks.txt のフォーマット:

urlap_report_v3.py と同一のフォーマットを使用します。

1行目: タイトル

2行目: 結晶系コード (1〜8)

3行目: 波長 [Å] と 初期Δ（度）（任意指定）

4行目以降: h k l 2theta_obs_raw（h>=1000 で読み込み終了）

A) 固定Δによる補正

指定した delta 値を用いて補正計算を行います。

xrd-umap zeroshift fixed -i peaks.txt --wavelength 1.5406 --delta 0.03 --out report.txt


B) Δを走査して最適化

指定した範囲（--delta-range）で delta 値を走査し、最適な値を探索します。

xrd-umap zeroshift scan -i peaks.txt --wavelength 1.5406 \
    --delta-range -0.5 0.5 --step 0.05 \
    --select-by uncprod --out report.txt


--select-by: 最適化の指標（rms または uncprod）を選択します。

C) 2列カーブ（CSV）にΔを適用

計算された delta 値を、実際の2列データ（2θ, 強度）に適用して補正します。

# 入出力は 2列CSV (2theta, intensity)
xrd-umap zeroshift apply-curve --in-csv raw.csv --out-csv shifted.csv --delta 0.03


データフロー（ゼロ点補正）:

flowchart TD
    A[peaks.txt (タイトル/系コード/波長/Δ/hkl+2θobs_raw)] --> B{モード}
    B -->|fixed| C[Δ を固定値で指定]
    B -->|scan| D[Δ∈[min,max] を走査]
    C --> E[Q空間での線形回帰 → 格子定数/不確かさ]
    D --> E
    E --> F[最適 Δ の決定 & レポート生成]
    F --> G[report.txt]
    F --> H[apply-curve で 2θ を補正 (任意)]


4. 対応フォーマット（Reader プラグイン）

ASCII 2列: .txt, .csv, .xy, .xye, .dat

PANalytical XRDML: .xrdml

ベンダー形式 (要 xylib-py): .ras, .uxd, .raw, .rd, .cpi, .udf, .xdd など

カスタムReaderの追加

Readerはプラグインとして簡単に追加登録できます。

from xrd_umap.io.readers import Reader, register_reader
import numpy as np

def read_my_format(path: str) -> tuple[np.ndarray, np.ndarray]:
    """ (2theta_array, intensity_array) を返すカスタムローダー """
    # ... 読み込みロジック ...
    x = np.array([...])
    y = np.array([...])
    return x, y

# 拡張子と優先度を指定して登録
register_reader(Reader(
    name="myfmt", 
    loader=read_my_format, 
    extensions=(".myfmt", ".custom"), 
    priority=25
))


5. ディレクトリ構成（主要部分）

src/xrd_umap/
  io/readers.py             # Reader レジストリ + 標準Reader (ASCII/XRDML/xylib)
  preprocess/resample.py  # build_grid / resample_to_grid (既存仕様のまま)
  preprocess/normalize.py # l2_normalize (既存仕様のまま)
  calibration/urlap_report_v3.py  # ゼロ点補正ロジック (原文のまま)
  calibration/zero_shift_cli.py   # ゼロ点補正CLIのアダプタ
  pipeline.py             # データセット生成の主処理 (仕様は元と同一)
  cli.py                  # メインCLI (`dataset` / `zeroshift` サブコマンド)
# xrd_umap

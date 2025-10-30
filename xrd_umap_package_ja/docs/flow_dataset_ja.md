
# データセット生成フロー（日本語）

```mermaid
flowchart TD
  A[入力ディレクトリ] --> B{Readerで読み込み}
  B -->|ASCII 2列| C[2θ, 強度]
  B -->|XRDML| C
  B -->|xylib ベンダー| C
  C --> D[共通グリッドを生成 (xmin,xmax,step)]
  D --> E[線形補間で再サンプリング]
  E --> F[L2 正規化]
  F --> G[calcd_patterns.csv (数値のみ)]
  F --> H[targets.csv (file→label)]
```

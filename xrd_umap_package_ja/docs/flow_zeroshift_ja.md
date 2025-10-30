
# ゼロ点補正フロー（日本語）

```mermaid
flowchart TD
  A[peaks.txt (タイトル/系コード/波長/Δ/hkl+2θobs_raw)] --> B{モード}
  B -->|fixed| C[Δ を指定]
  B -->|scan| D[Δ∈[min,max] を走査]
  C --> E[Q 空間での線形回帰 → 格子定数/不確かさ]
  D --> E
  E --> F[最適 Δ の決定 & レポート生成]
  F --> G[report.txt]
  F --> H[必要なら apply-curve で 2θ を補正]
```

# Capitalised panel text — Fig 1B, Fig 1C

교수님 요청 "패널 제목 다 대문자로" 반영분. 2026-08-03.

파일명이 배포본과 같으므로 **덱에 그대로 드롭인** 가능 (픽셀 치수 동일, 확대율 재조정 불필요).

| 파일 | 대체 대상 |
|---|---|
| `fig1A_smallmultiples_retitled.png` | `analysis/version 3/panels/` 동명 파일 (덱 slide 1) — ⚠️ 아래 참조 |
| `fig1B_trajectory_dsDNA.png` | `figures/panels/new color/optionB/` 동명 파일 (덱 slide 1) |
| `fig1B_trajectory_ssDNA.png` | 〃 |
| `fig1B_trajectory_pooled.png` | 〃 (현재 덱 미사용, 참고용) |
| `clonal_takeover_stack_pooled_sized_inset.png` | `analysis/version 3/panels/` 동명 파일 (덱 slide 1, Fig 1C) |
| `suppfig1B_heatmap_{ssDNA,dsDNA,pooled}.png` | `figures/panels/new color/optionB/` 동명 파일 (덱 slide 4) |
| `suppfig2_R{0..4}.png` | `figures/panels/` 동명 파일 (덱 slide 5) |
| `suppfig3A_morisita_horn.png` | `figures/panels/` 동명 파일 (덱 slide 5) |
| `suppfig3{B,C,D}_*.png` | `figures/panels/new color/optionB/` 동명 파일 (덱 slide 5) — 3B는 ⚠️ 아래 |

## 바뀐 문자열 — **패널 제목만.** 축 제목·범례·inset 라벨은 게재본 그대로 유지.

**Fig 1A**
- `output/input (titer)` → `Output/input (titer)`
- `clonal diversity` → `Clonal diversity`
- `read fraction of top 1% clones` → `Read fraction of top 1% clones`

**Fig 1B**
- `ssDNA` / `dsDNA` → **변경 없음**. 표기법이라 대문자화하면 `SsDNA`/`DsDNA`가 되어 틀림.
  → 이 두 장은 배포본과 **md5 동일**(`8212d7f790` / `3a5c78ae32`). 교체 불필요.
- pooled 패널만 산문 제목: `all (ssDNA + dsDNA)` → `All (ssDNA + dsDNA)`

**Fig 1C**
- 제목 `clonal composition` → `Clonal composition`

**Supp Fig 1B** (히트맵 3장)
- `ssDNA` / `dsDNA` → **변경 없음** → 배포본과 md5 동일(`ed8563940b` / `3828bd2d83`). 교체 불필요
- `all (ssDNA + dsDNA)` → `All (ssDNA + dsDNA)`
- 축 `selection round`, 컬러바 `depth-normalized (z)`, `n = 2145`, C1/C2/C3 스트립 라벨 전부 유지

**Supp Fig 2** (라운드별 산점도 5장)
- `all — R0` … `all — R4` → `All — R0` … `All — R4`

**Supp Fig 3**
- (A) `ss/ds overlap` → **변경 없음.** 제목이 표기법 `ss/ds` 로 시작해 `Ss/ds overlap` 이 되어버림.
  게재본 그대로 뒀으니 필요하면 리워딩을 따로 정할 것 (예: `Abundance overlap (ss/ds)`).
- (B) `shared clonotypes: ss vs ds cluster` → `Shared clonotypes: ss vs ds cluster`
  (둘째 줄 `concordant 788 (79%), discordant 215 (21%)` 는 그대로)
- (C) `rank correlation by cluster` → `Rank correlation by cluster`
- (D) `frequency bias by cluster` → `Frequency bias by cluster`

## ⚠️ Supp 3B — 덱에 옛 팔레트가 남아 있음

덱(`2026-08-02  Figures Junho v2 DK v1.pptx`) slide 5 의 Supp 3B 이미지는
`figures/panels/new color/_current_control/suppfig3B_cluster_confusion.png`
(= **리컬러 이전 팔레트**) 와 **픽셀 완전 일치**한다. 같은 그림의 3C·3D 는 새 팔레트로
교체돼 있어서, 한 그림 안에서 C2/C3 색이 서로 다르다.

16장 리컬러 패널 전수 대조 결과 (덱 이미지 vs 디스크 파일, 픽셀 평균차):

| 패널 | 덱 상태 |
|---|---|
| Fig 2A C2/C3 · Fig 2B C2/C3 | 새 팔레트 OK |
| Supp 1B ss/ds/pooled | 새 팔레트 OK |
| Supp 3C · Supp 3D | 새 팔레트 OK |
| **Supp 3B** | **옛 팔레트 (교체 누락)** |
| Fig 2A C1 R3 · Fig 2B C1 R4 | 덱 이미지가 optionB/control 어느 쪽과도 불일치 (평균차 7.1/7.3) — 별도 확인 필요 |

여기 있는 `suppfig3B_cluster_confusion.png` 는 **새 팔레트 + 대문자 제목** 이므로,
이걸 교체하면 색 누락도 같이 해소된다.

## 팔레트

전부 **현재 덱 팔레트(option B) C1 `#E63946` / C2 `#17375E` / C3 `#56C1B0`** 로 렌더됨.
- Fig 1B, Supp 1B → optionB 트리의 `make_panels.py` 를 import (팔레트가 그 안에 있음)
- Fig 1C → `make_fig1C_inset.py` 안에 팔레트가 하드코딩돼 있어 그대로 사용
- Fig 1A → 클러스터 색을 안 씀 (slate `#3D5A80` / orange `#E07A1F` / black). 팔레트 무관

⚠️ canonical `analysis/version 3/code/make_panels.py` 는 **아직 옛 팔레트**
(C2 `#457B9D` / C3 `#2A9D8F`) 다. 클러스터 색이 들어가는 패널을 여기서 렌더하면 안 된다.

## ⚠️ Fig 1A — 제목만 바꿀 수 없었던 이유 (본문 숫자 불일치)

덱에 들어 있는 Fig 1A 세 번째 패널이 **본문 Results 숫자와 다른 계열을 그리고 있다.**

|  | R0 | R1 | R2 | R3 | R4 |
|---|---|---|---|---|---|
| 본문 Results (ssDNA) | 19.4 | 19.6 | 23.7 | **68.3** | **62.2** |
| 본문 Results (dsDNA) | 19.1 | 19.3 | 24.4 | **44.7** | **37.5** |
| 덱 패널 ssDNA (눈금 판독) | ~7.5 | ~7.7 | ~6.8 | **~34.5** | **~36** |
| 덱 패널 dsDNA (눈금 판독) | ~7.3 | ~8.0 | ~8.5 | **~31.5** | **~23** |

- 본문 숫자는 `report_q40_ppm100_list.csv` 에서 **재계산해 정확히 재현됨**
  (라이브러리당 상위 1% = 11개 클론타입의 누적 read 비율; PPM/raw 어느 쪽으로 세도 동일).
- 같은 숫자가 `analysis/version 3/code/panelA_prototypes.py` 에 하드코딩돼 있고,
  그 스크립트가 만든 `panels/panelA_v3_smallmultiples.png` (6/17) 는 **본문과 일치**한다.
- 그런데 덱에 실린 `fig1A_smallmultiples_retitled.png` (8/2) 와 그 직전판
  `fig1A_smallmultiples_deployed_ORIG.png` 는 **둘 다 옛 계열**을 그린다.
  → 6/17에 지표를 top-1%-of-clonotypes 로 고치고 패널을 다시 만들었지만
  **덱에 갈아끼우지 않은 채로 남은 것**으로 보인다. 8/2 작업은 제목만 바꿨다.
- 그 옛 계열을 만드는 스크립트는 디스크에 없다. 따라서 **"제목만 바꾼" 재현이 불가능**하고,
  여기 있는 파일은 canonical 생성기(`panelA_prototypes.py`)의 검증된 숫자로 렌더한 것이다.
  즉 이 교체는 대문자화 + **본문/그림 숫자 불일치 해소**를 동시에 한다.
- y축 범위도 그래서 달라진다: 덱 0–40 → 새 파일 0–75 (68.3% 를 담아야 하므로).
- 치수: 덱 3241×1023 → 새 파일 3239×1024. 2 px 차이는 matplotlib tight-bbox 반올림
  (현재 3.10.8 vs 6/17 당시 버전). 종횡비 차 0.16% 로 배치에는 영향 없음.

## 재현성 검증

두 스크립트 모두 `--verify` 플래그를 지원한다. 라벨을 원래 소문자로 되돌려 렌더한 뒤
배포본과 md5를 비교하는 모드로, **전부 MATCH** 확인됨:

```
fig1B_trajectory_pooled.png                  MATCH  9427e10263
fig1B_trajectory_dsDNA.png                   MATCH  8212d7f790   (= 덱 slide1 이미지)
fig1B_trajectory_ssDNA.png                   MATCH  3a5c78ae32   (= 덱 slide1 이미지)
clonal_takeover_stack_pooled_sized_inset.png MATCH  186ee39495   (= 덱 slide1 이미지)
```

즉 글자 외에는 클러스터링·팔레트·폰트·치수·여백이 배포본과 완전히 동일하다.
클러스터 크기도 불변: pooled 491/264/1390 · dsDNA 262/191/640 · ssDNA 160/75/817.

## 실행

```bash
cd "figures/panels/capitalized/code"
/opt/anaconda3/bin/python render_fig1B_caps.py --verify
/opt/anaconda3/bin/python render_fig1C_caps.py --verify
```

기존 생성기는 **건드리지 않았다.** 두 스크립트는 각각
`figures/panels/new color/optionB/code/render_fig1B_promoted.py` 와
`analysis/version 3/code/make_fig1C_inset.py` 의 사본이며, 라벨 문자열과 출력 경로만 다르다.

## 남은 패널 (제목만 대상)

- Fig 1A — `output/input (titer)` · `clonal diversity` · `read fraction of top 1% clones`
- Fig 2A/B — `all — R3`, `C1 enriching — R3` 류 · Fig 2C
- Fig 3A — `R0`–`R4` (이미 대문자) · Fig 3B — `binders shift toward the ssDNA pool`
- Supp 1B/1C, Supp 2, Supp 3A–D — `ss/ds overlap` · `frequency bias by cluster` ·
  `shared clonotypes: ss vs ds cluster` 등

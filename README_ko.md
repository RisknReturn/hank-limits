# HANK 해석적 극한·집계 항등식 검증 세트 v0.3.0 (analytic-limit and aggregate-identity test suite)

[![DOI](https://zenodo.org/badge/1426447298.svg)](https://doi.org/10.5281/zenodo.23489176)

모든 판을 인용할 때는 concept DOI [10.5281/zenodo.23489176](https://doi.org/10.5281/zenodo.23489176) 를 씁니다. 이 릴리스(v0.3.0)의 DOI 는 Zenodo 보관 뒤 여기와 `CITATION.cff` 에 추가합니다. 이전 릴리스 v0.2.1 은 [10.5281/zenodo.23489177](https://doi.org/10.5281/zenodo.23489177) 입니다. GitHub 의 "Cite this repository" 버튼은 `CITATION.cff` 를 씁니다.

이 저장소는 이질적 경제주체 뉴케인지언 모형(heterogeneous-agent New Keynesian, HANK) 코드가 답이 알려진 극한(limit)과 집계 항등식(aggregate identity)에서 통과해야 하는 시험 모음입니다. 공개 패키지 `sequence-jacobian`(이하 SSJ) 1.0.0 위에서 구현했고, 각 시험은 한 지점의 통과/실패(pass/fail)만이 아니라 **극한 모수가 극한으로 갈 때의 오차와 관찰된 수렴 차수(observed order of convergence)** 를 함께 보고합니다.

- Zenodo 보관(concept DOI [10.5281/zenodo.23489176](https://doi.org/10.5281/zenodo.23489176)), MIT 라이선스(`LICENSE`). 비심사.
- 모형과 무관한 검사 함수(`hank_limits/checks.py`)와 SSJ 참조 구현(`hank_limits/models.py`)을 분리해, 다른 HANK 코드의 검증 장치(validation harness)로 재사용할 수 있게 했습니다.
- **이 세트가 입증하는 것과 입증하지 않는 것.** 통과는 구현 오류의 부재를 한정된 명제에 대해 확인한 수치 증거입니다. 통계적 유의성 검정과는 종류가 다르고, HANK 모형의 예측 우위나 모든 HANK 코드의 정확도를 보장하지 않습니다.

### v0.1 → v0.2 변경 요약 (GPT 6.1 Sol 의 R4 검토 반영; 17장 참조)

| 항목 | v0.2 변경 |
|---|---|
| 명칭·범위 | T1 = PIH **가계 블록** 극한(전체 NK 일반균형 RANK 극한 아님). T2 = **부과된 Bilbiie식 소득 귀착의 TANK 기준**(χ 구조 유도 아님). T4 = **개별 위험 소멸 극한(no-idiosyncratic-risk limit)**(완전시장 극한 아님). |
| T1 시계 | T=50·100 오차는 같은 T=200 야코비안의 부분행렬(관측창, observation window)입니다. 계산 시계를 따로 바꾼 야코비안은 부분행렬과 **정확히 같음**을 확인했습니다(차이 0). "절단 시계" 표현을 고쳤습니다. |
| T3 | 경기역행 귀착의 총소득 보존은 정상상태에서 **1차까지만** 성립함을 명시하고 크기를 보고합니다. 전이확률 Π 고정이므로 고용확률의 경기 변동은 복제하지 않습니다. |
| T4 | `r < 1/β−1` 은 근 탐색 구간이 강제하므로 독립 증거에서 제외. 대신 자본시장 잔차, Euler/KKT 잔차, 상단 격자 질량, 자산 상한·해상도 민감도를 추가했습니다. Aiyagari (1994) Table II 두 셀 대조를 추가했습니다. |
| T6 | 집단별 일반균형 발산은 "귀속 불가능성 정리"가 아니라 현재가치 0 조건을 어긴 forcing 의 성질임을 명시. |
| FastHetBlock | 음의 질량, 역반복 종료, shift·반복 횟수 민감도, SSJ 원래 전진 반복과의 대조를 진단으로 추가(`hank_limits/diagnostics.py`). |
| F1·F2 | SSJ 코드만 쓰는 최소 재현 스크립트, F2 배열 변경 대조와 정규화 역오차(backward error), 버전·wheel 해시. |
| 표 | 주 시험과 격자 연구의 수치(239행)를 극한 모수·격자·차분 간격·관측창·풀이기 열로 나눈 CSV 로 내보냅니다(`out/tables_v0_2/`). 진단·차분 간격·Aiyagari·F1·F2 결과는 각 JSON 에 있습니다. |

**v0.2 → v0.2.1 (R4 후속 검토 반영):** ① 구성상 항상 0 이던 진단 필드 `normalisation_err_last_sweep` 삭제, ② Aiyagari 둘째 셀을 "미재현·원인 미해결"로 고치고 원문 Table I 의 수준 CV·자기상관 대조와 격자 감사를 추가, ③ F2 의 실무 영향 서술을 이 예제로 한정, ④ T1 계산 시계 무관성을 고정 정상상태의 선형 가계 야코비안으로 한정, ⑤ SSJ 이슈 초안 문구 수정. 시험 결과는 같습니다.

### v0.2.1 → v0.3.0 변경: 정상상태·소득 정규화·명세 검사

계획과 사전 등록은 `docs/CHECKS_PLAN_v0_3.md`, 릴리스 전 검토 뒤의 수정과 그 재실행은 `docs/CHECKS_v0_3_0_AMENDMENT.md`, 오류 주입 시험(mutation testing)은 `docs/MUTATION_v0_1.md`(v0.2.1 묶음)와 `docs/MUTATION_v0_2.md`(독립 생성 변형으로 v0.2.1 과 v0.3 비교)에 있습니다.

| 항목 | v0.3.0 변경 |
|---|---|
| V1 | 공급자가 `ge_target` 으로 보고한 T3 일반균형 목표(`'asset_mkt'` 또는 `'goods_mkt'`; 없으면 `'asset_mkt'` 로 보고 "assumed" 로 기록)의 잔차를 b 로 나눠 기록하고 시험합니다. 풀이기가 보고한 값과, 보고된 야코비안과 검사가 넘긴 금리 경로로 직접 계산한 값을 모두 봅니다(자산시장은 `J_A` 로 J^{A,Y} dY + J^{A,r} dr, 재화시장은 `J_C` 로 dY − J^{C,Y} dY − J^{C,r} dr). 임계값: 자산시장은 자산시장 풀이 자신의 잔차에서 나온 2e-12(그대로), 재화시장은 참조 공급자를 재화시장 목표로 푼 실행에서 나온 9e-12(비례 소득)와 2e-11(경기역행 소득)입니다(v0.3.0). 다른 시장의 잔차는 임계값 없이 기록만 합니다. 적용 범위: 두 시장 가운데 하나로 푼 T3 에만 적용하고, 다른 목표가 보고되면 적용 조건 미충족으로 skip 합니다. v0.3.0 전에는 V1 이 자산시장 목표를 전제했으므로, 정확한 산술에서는 같은 해를 주는 재화시장 풀이(대조군 C03)에 경보를 냈습니다. 그 풀이의 자산시장 잔차는 비례 소득에서 F3 와 무관하게 최대 3.0e-6·b, 경기역행 소득에서 F3 의 O(h) 오차가 더해져 최대 2.4e-3·b 입니다(`docs/MUTATION_v0_2.md` 8절, `docs/CHECKS_v0_3_0_AMENDMENT.md`). |
| V2 | 공급자(provider)가 보고한 기본 수량으로 검사 쪽에서 직접 계산하는 정상상태 잔차: 정부 예산, 재화·자산시장(재정 HANK, TANK, Huggett), 자본시장·기업 1계 조건·재화시장(Aiyagari), 모든 모형의 집계 가계 예산. |
| V3 | 소득 정규화: Π 의 행 합이 1, 보고된 Π 의 정상분포에서 E[e] = 1, 총소득 Σ π y 가 모형별 기준과 같음. 기준은 Y(T1, T3), 유형별 소득 Y_H·Y_S 와 그 가중 합 λ Y_H + (1 − λ) Y_S = Y(T2), 문서화된 노동 L = 1 에서의 w L(T4), Y − Tax(T5)입니다. 경기역행 귀착(T3)에서는 Σ π e ζ = 0. |
| V4 | `hank_limits/spec.py` 의 문서화된 값으로 명세 고정(specification pins). `models.py` 에서 읽지 않습니다. 보고된 모수, 문서화된 β·EIS 에서의 Aiyagari 오일러 식, T5 경로에서의 재정 준칙·정부 예산·실질금리 준칙·재화시장을 확인합니다. v0.3.0: 공급자가 실제로 푼 외생 경로를 요청한 경로와 대조하고(`exog_G`, `exog_r`, `exog_rshock`), 실질금리 준칙 잔차에는 요청한 금리 충격을 씁니다. |
| V5 | `test_rank_benchmark_is_ricardian` 에는 RA 닫힌 해 assert 만 남기고, HA 조건은 `test_ha_multiplier_exceeds_one` 으로 나눴습니다. |
| V6 | F3 고정 시험(strict xfail, SSJ 코드만). F1·F2 고정과 같은 형식입니다. v0.3.0: V6a(hetinput 이 항상 한쪽 차분)와 V6b(`twosided=True` 를 요청해도 hetoutput 이 한쪽 차분)로 나눠, 상류에서 한쪽만 고쳐져도 그 시험만 XPASS 로 알립니다. |
| 공급자 | **선택** 필드 `ss_detail`, `J_A`, `J_C`, `ge_target`, `policy`, `irfs[...]` 의 `dC`·`dB`·`drshock` 추가. `hank_limits/interface.py` 에 공급자별로 검사가 읽는 모든 필드와 필수·선택 구분을 적었습니다. V 수량은 각자 자기 식에 필요한 필드로만 계산하고, 필드가 없으면 그 수량만 빠진 것으로 두고 이유를 남깁니다. `models.py`·`core.py` 의 계산은 바꾸지 않았고, Werning 공급자에 기본값이 기존 동작과 같은 선택 인자 `ge_target` 만 더했습니다. 오류 주입 보고서의 기준 실행(v0.2 의 4회, v0.3.0 재실행의 2회)에서 기존 검사의 요약값 276개가 v0.1 실행 B0 과 정확히 같았습니다. |
| 판정 | v0.3.0. V 시험은 값 하나라도 유한하지 않거나 음수이거나 임계값 이상이면, 다른 값이 빠져 있어도 실패합니다. 실패가 없는데 기대한 값이 빠져 있으면 **incomplete** 로 skip 하고(통과로 보지 않음), 공급자가 그 시험의 수량을 하나도 보고하지 않으면 이유와 함께 skip 합니다. 기대하는 값은 검사가 설정·소득 유형·요청한 충격별로 미리 정해 둡니다(`res['v0_3']` → `status`, `checks.v03_verdict`). |
| 검토 | 2026-10-22 릴리스 전 검토(P0·P1 결함 없음)를 반영했습니다. 위의 판정(R1), 문서와 코드의 필드 규약 일치(식 단위 필드, T3 β 는 `ss_detail` 에서 읽음; R2), 실제로 푼 외생 경로와 요청 경로의 대조(R3), 보고된 일반균형 목표로 V1 판정(R4), V6 나누기. 가짜 공급자(fake provider)로 검사기 자체를 시험하는 회귀 시험은 `tests/test_v0_3_0_checker.py` 에 있습니다. 설계·예측·임계값·재실행은 `docs/CHECKS_v0_3_0_AMENDMENT.md` 에 있습니다. 2026-10-23 병합 전 확인 뒤에는, 검사가 요청한 객체의 사본을 따로 보존하고 공급자에게는 다른 사본을 넘기도록 고쳤습니다. 그래서 공급자가 받은 인수를 직접 바꿔도 비교 기준은 바뀌지 않습니다(같은 문서 11절). |
| T4·T6 | 오류 주입 시험 v0.1 에서 나온 두 가지 설명: T4 의 오류 검출력은 진단에서 나오고(7·9절), T6 행렬 항등식은 선형대수 검사입니다. |

v0.3.0 의 pytest 전체: **114 passed, 4 xfailed** (F1, F2, V6a, V6b). 통과한 114개는 T1–T6, D, V 묶음의 시험 77개와 검사기 회귀 시험 37개입니다. 출력은 `out/pytest_output_v0_3_0.txt` 에 있습니다.

오류 주입 시험 v0.2(`docs/MUTATION_v0_2.md`): 새 검사를 모르는 별도 에이전트가 만든 변형 20개 가운데 v0.2.1 묶음은 14개, v0.3 묶음은 19개를 잡았습니다. 대조군 3개에서는 두 묶음 모두 경보를 내지 않았습니다. 두 묶음 모두 놓친 N20 은 어떤 시험도 엄격히 확인하지 않는 진단 값의 오류입니다. v0.1 변형 25개는 v0.3 묶음이 모두 잡지만, 그 검사들은 v0.1 결과를 본 뒤 설계했으므로 독립 결과가 아닙니다.

릴리스 전 검토 뒤의 재실행(`docs/CHECKS_v0_3_0_AMENDMENT.md` 10절): v0.3.0 으로 수정한 묶음을 같은 변형으로 설계 (ii)·(iii)에서 다시 돌렸습니다. (ii)는 결과가 바뀌지 않았습니다(20개 가운데 19개를 잡고 N20 을 놓침, 대조군 경보 없음). (iii)에서는 대조군 C03 만 바뀌어 더는 경보를 내지 않습니다. V1 이 이제 공급자가 보고한 일반균형 목표로 판정하기 때문이며, 설계상 그렇게 되는 결과이므로 독립 증거가 아닙니다. 위의 v0.2 결과는 당시 검사에 대한 기록으로 그대로 둡니다.

오류 주입 시험은 보고서(`docs/MUTATION_v0_1.md`, `docs/MUTATION_v0_2.md`, `docs/CHECKS_v0_3_0_AMENDMENT.md`)와, 설계 × 변형 × 시험 묶음마다 한 행인 CSV(`out/mutation_v0_1.csv`, `out/mutation_v0_2.csv`, `out/mutation_v0_3_0_rerun.csv`)부터 읽기를 권합니다. JSON 파일 `out/mutation_v0_*.json` 에는 모든 실행의 시험 결과와 검사 요약이 들어 있습니다.

---

## 1. 결과 요약 (summary)

오차 척도(error metric)는 **상대 최대 원소 오차(relative max-entry error)** `max|x − ref| / max|ref|` 입니다. 이것은 유도 행렬 노름(induced norm)이 아니며, 작은 꼬리 원소의 상대 정확도를 보장하지 않습니다. 기준값이 0 에 가까운 경우를 위해 절대·상대 허용치를 함께 쓰는 `metrics.check_close` 를 추가했습니다. 수렴 차수는 이웃한 두 설정 사이의 차수 `log(e_i/e_{i+1}) / log(h_i/h_{i+1})` 와 로그-로그 최소제곱 기울기(fitted slope)입니다.

| 시험 | 명제 (범위) | 극한 모수 | 가장 엄격한 설정의 오차 | 관찰 수렴 차수 (이론) | 상태 |
|---|---|---|---|---|---|
| T1 PIH 가계 블록 극한 | HA 가계 야코비안 dC/dY, dC/dr → 영구소득(permanent-income) 닫힌 해 (부분균형) | 소득위험 eps → 0 (β 재보정 경로) | dC/dY 9.5e-4, dC/dr 4.0e-4 (eps=0.0125); eps=0 에서 1.7e-9 | 1.95 / 2.03 (2) | 통과 |
| T2 부과 귀착 TANK 기준 | 총 iMPC = λχ I + (1−λχ) M_RA, 일반균형 산출 = (1−λ)/(1−λχ) × RANK | eps → 0 (두 유형) | 총 iMPC 8.0e-5, IRF (χ=1.5) 1.05e-3 | 저축가 1.95, IRF 1.79 (2) | 통과 |
| T3 Werning 사라지는 유동성 | 비례 소득이면 실질금리 경로에 대한 산출 반응 = RANK; 소득 탄력성 변화 시 일반화된 Euler 관계 | 차입한도 b → 0 | 비례 4.3e-4; 경기역행 6.0e-4 (이탈 0.319, 닫힌 해 0.320) | 0.99 / 0.99 (1) | 통과 |
| T4 Aiyagari 정상상태 | r → 1/β − 1 (개별 위험 → 0); 이 보정에서 위험 크기·지속성에 단조 (회귀 검사) | 로그소득 sd eps → 0 | 격차 6.2e-5 (eps=0.025) | 1.93 (2) | 통과 |
| T5 IKC 일관성 | SSJ 일반균형 dY = (I − M)^{-1} dC_partial (같은 블록 야코비안으로 구성한 일관성 시험) | 항등식 (T=100, 200, 300) | 정제 풀이 9.8e-12; SSJ 기본 풀이 5.9e-7 | 해당 없음 | 통과 + 예상 실패 2건 (F1, F2) |
| T6 HANK − RANK 분해 | 유한 T 리졸벤트 항등식과 같은 경로에 대한 차분 분해 | 항등식 (T=150, 300) | 항등식 2.2e-14, 분해 1.9e-10 (G) / 2.2e-10 (금리) | 해당 없음 | 통과 |

**차수 해석 주의.** T1·T2 의 O(eps²) 는 위험을 줄이면서 β 도 재보정하는 모형 경로의 차수이고, T4 의 O(eps²) 는 β = 0.96 을 고정한 채 위험을 줄일 때 격차 1/β − 1 − r 의 차수입니다. 어느 것도 공간 이산화가 2차 정확도라는 증거가 아닙니다. T3 의 O(b) 도 유동성 극한의 근사 차수입니다. 위험을 없애되 β 를 고정한 모든 모형이 T1 의 PIH 해로 간다고 읽으면 안 됩니다.

**허용치의 지위.** 시험 허용치는 이 보정의 **회귀 검사(regression check) 기준**입니다. 다른 구현에 같은 수렴 기울기나 같은 1e-10 기준을 요구하는 규격(conformance specification)이 아닙니다.

이 요약을 쓴 v0.2.1 의 pytest 전체: **41 passed, 2 xfailed, 약 65 ~ 75초** (11절). 위 T1–T6 결과는 v0.3.0 에서도 같고, v0.3.0 의 pytest 전체는 114 passed, 4 xfailed, 약 75초입니다(`out/pytest_output_v0_3_0.txt`). 보고서 스크립트 약 110초, 진단 스크립트 약 80초.

### 주요 발견 (key findings)

1. **T1: 관측창 의존(observation-window dependence).** HA 와 PIH 야코비안의 차이는 장기 원소(t=199, s=0)에 몰립니다. PIH 의 부(wealth)는 단위근을 갖지만 HA 의 부는 O(eps²) 속도로 평균회귀하기 때문입니다. eps=0.0125 에서 관측창 50, 100, 200 의 dC/dY 오차는 2.0e-4, 3.8e-4, 9.5e-4 입니다. 계산 시계를 T=50, 100 으로 따로 잡아 계산한 야코비안은 T=200 야코비안의 앞 블록과 **원소 단위로 정확히 같습니다**(최대 차이 0.0). SSJ 가계 야코비안의 원소는 계산 시계에 의존하지 않으므로, 이것은 절단(truncation) 효과가 아니라 긴 lag 를 관측창에 넣을수록 오차가 커지는 현상입니다. 고정된 관측창에서는 eps → 0 극한이 성립하지만, 고정된 eps 에서 관측창을 무한히 늘리면 성립하지 않습니다. 이 계산 시계 무관성은 **고정 정상상태의 선형 가계 야코비안**에 대한 결과이며, 일반균형 H_U 의 역행렬, 내생 가격 경로, 비선형 전이, 종단조건이 T 에 무관하다는 뜻은 아닙니다.
2. **이산화 오차의 원인.** 가장 작은 eps 에서 차수가 2 보다 약간 낮아지는 것은 Young 복권 방식(lottery method)의 수치 확산(numerical diffusion)이 eps 가 작을 때 실제 확산과 비슷해지기 때문입니다. 격자를 세밀하게 하면 차수가 2 를 향해 올라갑니다(T2 IRF: 2400, 4800, 9600 점에서 1.47, 1.72, 1.85; T1 dC/dr 9600 점에서 2.006).
3. **F1: 캐시된 야코비안과 출력 전파의 일관성 문제(SSJ 1.0.0, master eb4ac4c).** `HetBlock._impulse_linear` 는 전달된 `Js` 를 쓰지 않고 자신의 옵션으로 다시 미분합니다. SSJ 코드만 쓰는 최소 재현에서 `Js` 를 2배로 바꿔 넣어도 출력이 **변하지 않습니다**(차이 0). 일반균형 풀이(H_U)와 직접 충격·출력 전파가 서로 다른 선형 연산자를 쓰게 되어, 가계에 직접 들어가는 금리 충격의 IKC 항등식이 약 1e-4 로만 성립합니다. 옵션을 일치시키면 3e-13 입니다. SSJ 미분 이론의 문제가 아니라 재사용 규약의 문제이며, 옵션 일치는 우회책일 뿐 캐시 재사용과 같지 않습니다. SSJ 에 [#39](https://github.com/shade-econ/sequence-jacobian/issues/39) 로 보고했습니다.
4. **F2: 방정식 배열에 민감한 LU 안정성.** 금리 준칙을 더한 재정 HANK 에서, SSJ 가 변수별 블록으로 쌓은 H_U 에 대한 부분 피벗팅 LU 의 원소 성장률이 T=100, 200, 300 에서 19, 3.7e4, 1.1e7 입니다(조건수 3.6e2, 2.9e3, 2.1e4). 같은 시스템을 **날짜별로 재배열하면 성장률 1.02, 블록 순서를 뒤집으면 1.00** 이 되고 정규화 역오차가 3e-8 에서 1e-16 수준으로 내려갑니다. 조건수만으로 설명되는 문제가 아니라, 이 배열에서의 LU 원소 성장이 정확도 저하를 만듭니다. 따라서 "T 가 길면 SSJ 풀이가 부정확하다"가 아니라 **"이 블록 배열과 부분 피벗팅 LU 의 조합에서 원소 성장이 생긴다"** 입니다. 재배열·반복 정제·QR 중 어느 하나가 모든 모형을 해결한다고 일반화하지 않습니다. 이 예제에서 IKC 오차는 T=200·300 에서 2.2e-9 ~ 5.9e-7, pivoted QR 해와의 상대 차이는 약 1e-7 (이 환경 1.0e-7, 검토자 환경 1.2e-7) 입니다. 엄격한 항등식 검사(1e-10)에는 의미가 있습니다. 실무에서의 중요도는 모형·목표·허용치에 달리며 여기서는 비교하지 않았습니다. SSJ 에 [#40](https://github.com/shade-econ/sequence-jacobian/issues/40) 로 보고했습니다.
5. **T6: 집단별 일반균형 귀속의 범위.** 부 사분위별 첫 단계(first round) 격차는 잘 정의되고 T 와 무관합니다. 하지만 집단별 forcing 의 현재가치가 0 이 아니므로(소득이 (Y − Tax) × e 이고 집단마다 e 구성이 다름), 이를 무한 시계 소비 승수에 각각 통과시키면 예산제약과 맞지 않는 모드를 자극해 T 에 따라 발산합니다. 합계에서는 이 모드가 상쇄됩니다. 이것은 "집단별 귀속은 불가능하다"는 정리가 아닙니다. 집단별 일반균형 귀속에는 종단조건, 재정 귀착, 현재가치 보상 경로를 추가로 정해야 하고, 그 선택에 따라 귀속이 달라집니다. 현재가치 0 은 필요한 적합성 조건일 뿐 안정성·유일성을 보장하지 않습니다.
6. **분포 정확도는 정상성 잔차만으로 보장되지 않음(관찰, SSJ 결함 아님).** 느리게 혼합하는 사슬에서는 질량 0 부분공간에서의 역연산자가 커지므로, 작은 잔차가 작은 분포 오차를 뜻하지 않습니다. 중간 위험 PIH 가계(eps=0.1)에서 SSJ 원래 전진 반복의 기본 정지 규칙(10회마다 max|ΔDbeg| < 1e-10)은 37,910회에 멈추며, 그 분포의 평균 자산은 역반복 해와 **1.7e-3 (상대)** 다릅니다. 허용치를 1e-12 로 조이면 54,590회에 1.7e-5 로 줄어 같은 해로 수렴합니다. eps=0.0125 에서는 600,000회 안에 정지하지 않습니다. 빠르게 혼합하는 T4·T5 모형에서는 차이가 1e-6, 5e-9 입니다 (12절).
7. **F3 (v0.3): hetinput·hetoutput 의 단측 차분 (SSJ 1.0.0).** `HetBlock.jac_backward_prelim` 은 hetinput 에 양측 차분을 요청합니다("always use two-sided differentiation for hetinputs"). 그러나 `DifferentiableCombinedExtendedFunction.diff` 의 기본값이 `None` 이 아니라 `twosided=False` 라서, 호출할 때 인자를 생략하면 저장된 설정이 쓰이지 않습니다. 결과는 두 가지입니다. hetinput 은 항상 단측 차분으로, hetoutput 은 `jacobian(..., twosided=True)` 로 요청해도 단측 차분으로 계산됩니다. `stage_block.py` 도 같은 방식입니다. 오차는 O(h) 이고, 비선형 hetinput·hetoutput 에서만 의미가 있습니다. 이 세트에서는 T3 의 경기역행 소득이 해당하며, 오류 주입 시험 v0.1 에서 찾은 O(h) 예산 항등식 잔차의 원인입니다(그 보고서 8절). SSJ 코드만 쓰는 재현은 14절에 있습니다. SSJ 에 [#41](https://github.com/shade-econ/sequence-jacobian/issues/41) 로 보고했습니다.

---

## 2. 실행 방법 (how to run)

```bash
pip install --break-system-packages sequence-jacobian==1.0.0 pytest matplotlib   # numpy, scipy, numba 필요
cd hank-limits
python3 -m pytest                          # 약 75초; 요약 출력은 python3 -m pytest -s
python3 scripts/convergence_report.py      # 약 110초: out/convergence_v0_2.json, PNG 6장, out/tables_v0_2/*.csv
python3 scripts/diagnostics_v0_2.py        # 약 80초: out/diagnostics_v0_2.json (12절)
python3 scripts/aiyagari_table2_v0_2.py    # 약 20초: out/aiyagari_table2_v0_2.json (7절)
python3 scripts/repro_F1_v0_2.py           # SSJ 코드만 사용: out/repro_F1_v0_2.json
python3 scripts/repro_F2_v0_2.py           # SSJ 코드만 사용: out/repro_F2_v0_2.json, out/F2_HU_T300_v0_2.npz
python3 scripts/repro_F3_v0_1.py           # v0.3, SSJ 코드만 사용: hetinput 미분 (14절)
python3 scripts/repro_F3_hetblock_v0_1.py  # v0.3, SSJ 코드만 사용: HetBlock.jacobian 의 hetinput·hetoutput
python3 scripts/mutation_v0_3_0.py --v03-ref <커밋> --design ii iii   # 약 18분: 오류 주입 재실행 (v0.3.0)
```

- 모든 진입점에서 스레드를 1개로 제한합니다(`OMP_NUM_THREADS=1` 등).
- 사용 버전: **sequence-jacobian 1.0.0** (PyPI wheel SHA-256 `40715634f6fcf8f695d79687c471bedcecacdeff94522557f5186b630f46ef01`), Python 3.13.16, numpy 2.5.3, scipy 1.18.1, numba 0.68.0, OpenBLAS 0.3.34 (scipy-openblas), pytest 9.1.1. GitHub master (커밋 `eb4ac4c3157357a4cacbea4b5a7a281f5b48cab9`)도 F1 과 같은 코드입니다.
- 검토 모형(GPT 6.1 Sol, 17장)이 자기 환경(Python 3.12.14, numpy 2.5.4, scipy 1.17.0)에서 v0.1 을 재실행해 **32 passed, 2 xfailed** 를 재현했습니다. BLAS·환경에 따라 마지막 자릿수와 시간은 달라질 수 있습니다.

### 디렉터리 구조 (layout)

```
hank_limits/
  core.py         가계 블록, FastHetBlock(정상분포 직접 풀이, 진단 포함), 격자
  diagnostics.py  v0.2: Euler/KKT 잔차, 분포 감사(shift·반복 민감도, SSJ 전진 반복 대조)
  ra.py           RA, TANK, 영유동성 닫힌 해
  models.py       SSJ 참조 구현 (provider)
  checks.py       모형과 무관한 검사 함수
  interface.py    provider 규약
  metrics.py      오차 척도, 수렴 차수, 절대·상대 혼합 허용치
  suite.py        표준 설정 (tests 와 스크립트가 공유)
  spec.py         v0.3: V 검사에 쓰는 문서화된 명세 값
tests/            test_t1 ... test_t6, test_diagnostics_v0_2, test_v0_3_checks (v0.3),
                  test_v0_3_0_checker (v0.3.0: 가짜 공급자로 검사기 시험)
scripts/          convergence_report.py, diagnostics_v0_2.py, aiyagari_table2_v0_2.py,
                  repro_F1_v0_2.py, repro_F2_v0_2.py, repro_F3_*.py, a1_markov_checks_v0_2.py (설계 노트용, 13절),
                  mutation_v0_1.py, mutation_v0_2.py, mutation_v0_3_0.py (오류 주입 실행기), thresholds_v0_3_0.py,
                  rerun_v0_3_0_compare.py
out/              JSON, PNG, tables_v0_2/*.csv, pytest 출력
```

---

## 3. 공통 규약과 수치 기법 (conventions and numerics)

- **시점 규약.** SSJ 단일자산 가계: 예산식 c_t + a_t = (1 + r_t) a_{t−1} + y_t, 오일러 식 c_t^{−1/σ} = β (1 + r_{t+1}) E c_{t+1}^{−1/σ}. r_t 는 t−1 기에 보유한 자산에 t 기에 지급되는 사후 실질수익률입니다.
- **연간 보정.** 시험 세트는 연간(annual) 보정이며 정상상태 산출 Y = 1 입니다. 분기(quarterly) 보정인 UK 보정 노트와 섞지 않습니다.
- **야코비안.** SSJ 가짜 뉴스 알고리즘(fake-news algorithm)에 **양측 차분(two-sided differences)**, h = 1e-4. SSJ 기본값인 단측 차분은 dC/dr 에 약 5e-5 의 O(h) 오차를 남깁니다. 차분 간격 민감도는 4절에 있습니다.
- **`FastHetBlock`.** SSJ `HetBlock` 의 하위 클래스로, 정상분포를 찾는 방법만 바꿉니다. Young 복권 전이행렬의 고유값 1 문제를 희소 역반복(sparse inverse iteration, (1+τ)I − L 의 LU 한 번, τ = 1e-10, 자산 우선 순서로 띠 행렬 유지)으로 풀고, 역방향 반복을 이전 해로 시작합니다. 나머지(가짜 뉴스 알고리즘, 일반균형 풀이)는 SSJ 그대로입니다. 이것은 자체 수치 구성요소이므로 v0.2 에서 다음을 기록합니다(`FastHetBlock._ss_diag`, 12절).
  - 절댓값을 취하기 전 음의 질량(절대 질량 합 대비), 역반복 횟수와 종료 여부, 마지막 갱신 크기, 다듬기 단계의 질량 표류. (v0.2 의 `normalisation_err_last_sweep` 은 구성상 항상 0 이어서 v0.2.1 에서 삭제했습니다. 음의 질량은 `neg_mass_rel_max` 가 보고합니다.)
  - 결과 분포의 L1·최대 원소 잔차: 잘린(clipped) 사슬 기준과 SSJ 자신의 전진 연산자 기준.
  - shift τ (1e-12 ~ 1e-6)와 반복 횟수를 바꿨을 때 분포·평균 자산·소비의 변화.
  - 같은 정책에서 SSJ 원래 전진 반복(`HetBlock.forward_steady_state`)과의 대조.
- **격자 밖 정책.** 희소 전이 구성은 격자 상단 밖으로 나가는 정책을 끝점에 모으고(reflecting barrier), SSJ 전진 연산자는 그대로 평가합니다. 영향받는 질량(`mass_policy_off_grid`)과 SSJ 잔차를 `FastHetBlock._ss_diag` 에 계속 기록합니다. v0.2 진단 출력에 기록된 T4 점들에서는 기본 상한(200)에서 3.1e-30 이하, 상한 100 의 민감도 점검에서 최대 8.7e-13 입니다(`out/diagnostics_v0_2.json` B).

---

## 4. T1: PIH 가계 블록 극한 (household-block permanent-income limit)

**범위.** 무위험 PIH 가계 블록의 정확한 야코비안 대조와, 자산 목표에 맞춰 β 를 재보정하는 특정 eps 경로의 수렴 시험입니다. **부분균형(partial equilibrium) 가계 블록 시험이며, 전체 NK 일반균형의 RANK 극한 시험이 아닙니다.**

**명제.** 소득위험이 없고 β(1+r) = 1 이며 차입제약이 묶이지 않는 가계는 영구소득 소비자이고, q = 1/(1+r) 일 때
- dC_t/dY_s = (1−q) q^s
- dC_t/dr_s = (1−q) A q^s + σ C q (1{t ≥ s} − q^s)

입니다(오일러 식과 시점 간 예산식에서 유도, `hank_limits/ra.py`). 출처: 표준 영구소득 가설의 선형화. SSJ 방법론은 Auclert, Bardóczy, Rognlie, Straub (2021).

**정확한 조건.** r = 0.03, σ(EIS) = 0.5, y = Y e, Rouwenhorst 3상태, 지속성 0.9, 로그소득 무조건부 sd = eps. 차입한도 = 0.99 × 자연 차입한도 −y_min/r (제약 질량 3e-6 이하). 평균 자산 A = 5 가 되도록 β 를 보정하며, κ = 1 − β(1+r) 은 eps² 처럼 0 으로 갑니다(κ/eps² = 0.028). eps = 0 은 정확한 영구소득 경우입니다(β(1+r)=1, A=5 에 모든 질량; 이 경우 부의 정상분포가 유일하지 않으므로 점질량을 선택). T = 200, 자산 격자 4800 점(선형), 상한 250.

**수치 결과.**

| eps | 0 | 0.2 | 0.1 | 0.05 | 0.025 | 0.0125 |
|---|---|---|---|---|---|---|
| dC/dY 오차 | 1.7e-9 | 2.08e-1 | 5.52e-2 | 1.40e-2 | 3.60e-3 | **9.53e-4** |
| dC/dr 오차 | 2.9e-9 | 1.12e-1 | 2.63e-2 | 6.44e-3 | 1.60e-3 | **3.96e-4** |
| κ | 0 | 1.14e-3 | 2.74e-4 | 6.79e-5 | 1.71e-5 | 4.35e-6 |

- 이웃 차수: dC/dY 1.92, 1.98, 1.96, 1.92 (최소제곱 1.95); dC/dr 2.09, 2.03, 2.01, 2.01 (2.03).
- eps = 0 의 1e-9 오차는 유한차분의 O(h²) 오차뿐이며, SSJ 시점 규약과 닫힌 해가 일치함을 보여 줍니다.

**관측창과 계산 시계 (v0.2).** eps = 0.0125:

| 관측창 (같은 T=200 야코비안의 앞 블록) | 50 | 100 | 200 |
|---|---|---|---|
| dC/dY 오차 | 1.96e-4 | 3.79e-4 | 9.53e-4 |
| dC/dr 오차 | 9.25e-5 | 9.60e-5 | 3.96e-4 |
| 계산 시계를 그 값으로 따로 잡은 야코비안과의 최대 원소 차이 | 0.0 | 0.0 | — |

같은 결과가 eps = 0.1 에서도 성립합니다(`out/diagnostics_v0_2.json` C). 종단(terminal-horizon) 효과를 따로 감사할 필요가 있는 것은 일반균형 풀이(T5·T6)입니다(예: T6 RANK 승수의 1 + 6.4e-6 은 T=150 에서 잘린 부채 경로 꼬리 때문).

**차분 간격 h 민감도 (양측 차분, T = 200, v0.2).**

| eps | h = 1e-3 | 1e-4 (기본) | 1e-5 | 1e-6 |
|---|---|---|---|---|
| 0: dC/dY / dC/dr 오차 | 2.8e-9 / 정의 안 됨* | 1.7e-9 / 2.9e-9 | 1.2e-8 / 9.7e-10 | 2.3e-7 / 6.7e-9 |
| 0.0125: dC/dY / dC/dr 오차 | 9.527886e-4 / 정의 안 됨* | 9.527901e-4 / 3.9569e-4 | 9.527900e-4 / 3.9608e-4 | 9.527900e-4 / 3.9611e-4 |

\* r 를 1e-3 올리면 자연 차입한도 −y_min/(r+h) 가 격자 하한(0.99 × −y_min/r)보다 위로 올라가 소비가 음수가 됩니다. 근사 자연한도 격자에서는 r 방향 차분 간격이 약 3e-4 미만이어야 합니다. eps = 0.0125 의 dC/dY 오차는 h ≤ 1e-4 에서 상대 1.1e-7 이내로 같아 차분 간격과 무관하고, eps = 0 에서는 h 가 작아질수록 반올림 오차가 커집니다.

**격자 수렴.** eps = 0.0125 dC/dY 오차: 1200 점 1.41e-3, 2400 점 1.08e-3, 4800 점 9.53e-4, 9600 점 9.01e-4. 9600 점에서 eps 0.025 → 0.0125 차수는 dC/dY 1.96, dC/dr 2.006.

**알려진 한계 (v0.1 결과).** 극한 분포의 표준편차가 eps 와 무관하게 약 14(연소득 단위)로 넓게 남으므로, 근사 자연한도와 넓은 격자가 필요합니다. 한도 0.95 × 자연한도(600 점)에서는 제약 질량 1.8e-4 가 오차 바닥(약 6e-3)을 만듭니다. eps 를 0.0125 보다 작게 하면 수치 확산이 실제 확산보다 커져 더 세밀한 격자가 필요합니다.

---

## 5. T2: 부과된 Bilbiie식 소득 귀착의 TANK 기준 (TANK benchmark with imposed Bilbiie-style incidence)

**범위.** 주어진 소득 귀착 χ 에서의 TANK 집계·실질금리 전달 시험입니다. χ 를 노동시장·이윤 배분에서 유도하지 않았으므로 "Bilbiie 모형의 구조적 복제"가 아닙니다. χ 를 외생적으로 부과해도 기준(benchmark)으로서의 가치는 유지됩니다.

**명제.** 비율 λ 의 손에서 입으로 사는 가계(H)와 비율 1−λ 의 저축가(S). H 소득이 총소득에 대해 탄력성 χ 를 가지면 (i) 일반균형에 들어가는 집계 iMPC 행렬은 λχ I + (1 − λχ) M_RA 입니다(χ = 1 이면 λ I + (1−λ) M_RA). (ii) 영유동성(채권 순공급 0)에서 집계 오일러-IS 식은 c_t = c_{t+1} − σ (1−λ)/(1−λχ) r_{t+1} 이고, 고정 실질금리 경로에 대한 일반균형 산출 반응은 RANK 의 (1−λ)/(1−λχ) 배입니다.

**성립 조건.** 두 유형의 정상상태 소비와 EIS 가 같고, H 형이 극한에서 정확히 손에서 입으로 살며, 고정 실질금리 경로이고, **1 − λχ > 0** 인 영역(여기서 λ = 0.3, χ ∈ {1, 1.5} 이므로 1 − λχ = 0.55)입니다. 증폭 배수 0.7/0.55 = 1.27273 은 λ = 0.3, χ = 1.5 에서의 값입니다.

**출처.** Bilbiie (2020), "The New Keynesian cross", JME 114, 90-108 (집계 오일러-IS 식과 "χ > 1 일 때만 증폭"은 ECB 학회 2018 판본에서 확인; Bilbiie 의 r_t 는 SSJ 규약의 r_{t+1}). Bilbiie (2008), JET 140(1), 162-196.

**정확한 조건.** λ = 0.3. H 형: 차입한도 0, β_H = 0.95, 소득위험 eps, 500 점. S 형: T1 과 같은 가계(4800 점), λ A_H + (1−λ) A_S = 0 이 되도록 β_S 보정. 귀착: Y_H = Y^χ, Y_S = (Y − λ Y_H)/(1−λ). 일반균형: dr_t = −0.01 × 0.8^t, Y = C, 미지수 Y, T = 200. 두 유형 모두 같은 eps → 0.

**수치 결과.**

| eps | 0.2 | 0.1 | 0.05 | 0.025 | 0.0125 |
|---|---|---|---|---|---|
| H 형 M_H vs I | 0.596 | 0.342 | 3.9e-13 | 1.7e-13 | 3.9e-13 |
| S 형 M_S vs M_RA | 0.275 | 7.29e-2 | 1.85e-2 | 4.76e-3 | **1.26e-3** |
| 총 iMPC (χ=1) vs λ I + (1−λ) M_RA | 0.556 | 0.320 | 1.18e-3 | 3.03e-4 | **8.02e-5** |
| 일반균형 IRF vs (1−λ)/(1−λχ) × RANK (χ=1.5) | 0.131 | 5.44e-2 | 1.24e-2 | 3.46e-3 | **1.05e-3** |
| HA 증폭 배수 (χ=1.5 / χ=1) | 1.150 | 1.209 | 1.27273 | 1.27273 | 1.27273 |

χ = 1.5 에서 λχ M_H + (1−λχ) M_S 와 λχ I + (1−λχ) M_RA 의 차이는 eps = 0.0125 에서 4.3e-5 입니다(`out/tables_v0_2/` 의 `err_Magg_lamchi_vs_TANK_chi1.5`). H 형은 eps 가 유한한 문턱값보다 작아지면 정확히 손에서 입으로 사는 가계가 됩니다. 격자 수렴(eps=0.0125 IRF): 2400 점 1.50e-3, 4800 점 1.05e-3, 9600 점 8.70e-4.

**알려진 한계.** χ 는 축약형입니다. 채권 순공급 0 이므로 dC/dr 의 부의 효과 항은 여기서 검증되지 않습니다(T1 에서 검증). 총 iMPC 하나에 대해 기울기를 적합하는 것은 의미가 없어 이웃 차수만 보고합니다.

---

## 6. T3: Werning (2015) 사라지는 유동성 기준 (vanishing-liquidity benchmark)

**범위.** 사라지는 유동성 정제(vanishing-liquidity refinement), Werning 의 비례 소득 중립성, 소득 탄력성 변화에 따른 일반화된 Euler 관계를 시험합니다. 경기역행적 사례는 Werning 3.3절에 맞는 **소득 탄력성 변화**이며, 3.4절의 고용확률 변화 모형 자체를 구현한 것은 아닙니다. 전이확률 Π 가 고정이므로 고용확률의 경기 변동은 복제하지 않습니다.

**명제.** 차입 불가, 외부 자산 없음이면 균형에서 모두가 소득을 소비하고 실질금리는 저축 의향이 가장 큰 가계의 오일러 식이 정합니다(Werning 식 12). 정확히 유동성이 0 이면 금리가 구간으로 불확정할 수 있지만, 사라지는 유동성 정제가 한계 금리를 선택합니다. 거듭제곱 효용(식 14)과 가계 소득이 총소득에 비례(식 15)하면 집계 오일러 식의 금리 민감도가 RA 와 같습니다(Proposition 2). 출처: Werning (2015), NBER WP 21448, 3.1-3.4절. 구현은 Huggett (1993) 경제.

**소득위험 경기성 가정.** y_it = Y_t e_it (Y_t/Y)^{ζ_i}, e 는 Y 와 독립인 Markov 사슬(Rouwenhorst 3상태, 지속성 0.9, 로그 sd 0.3).
- 비례 소득: ζ = 0.
- 경기역행적 위험: ζ = (1, −0.185, −0.185). Σ_i π_i e_i ζ_i = 0 은 정상상태에서 총소득 보존을 **1차까지만** 보장합니다. 선형 SSJ 시험에는 충분합니다. 비선형 경로로 확장하면 분모 Σ π_i e_i (Y/Y_ss)^{ζ_i} 로 정규화하거나 총소득 편차를 따로 보고해야 합니다. 실제 편차(v0.2): Y = 1 ± 0.001 에서 9.3e-8, ± 0.01 에서 9.2e-6 ~ 9.3e-6, ± 0.05 에서 2.2e-4 ~ 2.4e-4 (O(dY²)).

영유동성 닫힌 해: 한계 가계 h = argmax_e E[(e'/e)^{−1/σ} | e] 의 오일러 식을 선형화하면 (1+ζ_h) ŷ_t = (1+ζ̄) ŷ_{t+1} − σ dr_{t+1}/(1+r), ζ̄ 는 한계효용 가중 평균입니다. 경기역행 경우 δ = 1.0174, 탄력성 σ/(1+ζ_h) = 0.614 입니다. 한계 저축 상태의 argmax 가 동률에 가까워지면 미분 가능성을 별도로 점검해야 합니다 (현재 보정에서는 상위 상태가 유일한 한계 저축가).

**정확한 조건.** σ = 0.5, β 는 r* = 1% 가 되도록(β = 0.86997, v0.1 결과). 차입한도 −b, 채권 순공급 0, 격자 −b + 20 b u² (400 점), 차분 간격 h = 1e-3 b. dr_t = −0.01 × 0.8^t, T = 150.

**수치 결과.**

| b | 1e-1 | 3e-2 | 1e-2 | 3e-3 | 1e-3 | 3e-4 | 1e-4 |
|---|---|---|---|---|---|---|---|
| 비례 소득: 오차 vs RANK | 6.0e-2 | 6.6e-3 | 1.29e-2 | 1.23e-2 | 4.23e-3 | 1.29e-3 | **4.30e-4** |
| 경기역행: 오차 vs 영유동성 해 | 0.100 | 2.3e-2 | 9.7e-3 | 1.70e-2 | 5.89e-3 | 1.79e-3 | **5.99e-4** |
| 경기역행: RANK 와의 이탈 | 0.452 | 0.350 | 0.333 | 0.297 | 0.312 | 0.317 | **0.319** |
| (r(b) − r*)/b | 0.58 | 1.15 | 1.85 | 2.07 | 2.081 | 2.087 | 2.088 |

점근 구간(b ≤ 1e-3) 차수 0.993 (두 경우). 격자 200, 400, 800 점에서 오차가 네 자리까지 같습니다. b = 0 에서 수치 야코비안을 직접 미분하지 않은 것은 꺾임점(kink) 때문이며, 정제와 같은 b → 0 극한으로만 검증합니다.

---

## 7. T4: Aiyagari 정상상태 (no-idiosyncratic-risk limit)

**범위.** (i) 개별 위험 소멸 극한: 소득위험 → 0 이면 r → 1/β − 1. 코드에 보험시장 완비화가 없으므로 "완전시장 극한"이 아닙니다. (ii) 이 보정에서 위험이 있으면 r < 1/β − 1 이고, 위험 크기와 지속성에 따라 r 이 낮아지는 **회귀 검사**. 일반 정리로 쓰지 않습니다.

**기계적 검사와 비기계적 검사.** `ssj_aiyagari_r` 는 근을 (−δ + 1e-3, 1/β − 1 − 1e-10) 에서 찾으므로 r < 1/β − 1 은 구성상 성립하며 독립 정확도 증거로 세지 않습니다(sanity check 로만 유지). v0.2 의 비기계적 검사:

| 점검 | 모든 T4 점에서의 최댓값 |
|---|---|
| 자본시장 잔차 \|A − K\|/K | 5.2e-10 |
| Euler 오차 (비제약 상태, 분포 가중 평균) | 1.1e-7 |
| Euler 오차 (비제약 상태, 질량 있는 점 최대) | 1.6e-5 |
| KKT 위반 (제약 상태: u'(c) < βR E u'(c') 인 경우만) | 0 |
| 상단 5개 격자점 질량 / 격자 밖 정책 질량 | 4.4e-29 / 3.1e-30 |

Euler 오차는 단위 없는(unit-free) c_impl/c − 1 이고, c(e', ·) 는 자산 격자에서 선형 보간합니다. 제약 상태의 양의 Euler 격차는 수치 오차로 세지 않습니다(`hank_limits/diagnostics.euler_kkt_residuals`).

**자산 상한·해상도 민감도 (sd = 0.4, 지속성 0.9, 가장 넓은 부 분포).** r (%):

| a_max \ 점 수 | 500 | 1000 | 2000 |
|---|---|---|---|
| 100 | 2.8818 | 2.8822 | 2.8823 |
| 200 (기본) | 2.8813 | 2.8821 | 2.8823 |
| 400 | 2.8801 | 2.8818 | 2.8822 |

**T4 가 오류를 잡는 곳 (v0.3 설명).** 오류 주입 시험 v0.1 에서 해석적 극한(r → 1/β − 1, O(eps²))과 단조성 회귀 검사는 T4 에 도달한 변형 7개를 하나도 잡지 못했습니다. T4 의 검출 4개는 모두 위의 진단(자본시장 잔차, Euler/KKT 잔차, 상단 격자 질량)에서 나왔습니다. 극한은 δ, 차입한도, 격자를 어떻게 고정해도 성립하므로, 그런 오류는 진단과 v0.3 의 V2·V4 검사(문서화된 α·δ 에서의 기업 1계 조건, 문서화된 차입한도, 문서화된 β·EIS 에서의 오일러 식)가 잡습니다.

**출처.** Aiyagari (1994), QJE 109(3), 659-684.

**정확한 조건.** β = 0.96, α = 0.36, δ = 0.08, σ(EIS) = 0.5 (위험회피 2), 차입한도 0, 노동 1, Rouwenhorst 7상태, [0, 200] 에 1000 점(u² 간격).

**수치 결과.**

| eps | 0.4 | 0.2 | 0.1 | 0.05 | 0.025 |
|---|---|---|---|---|---|
| r (%) | 2.8821 | 3.8119 | 4.0742 | 4.1429 | 4.1605 |
| 1/β − 1 − r | 1.28e-2 | 3.55e-3 | 9.24e-4 | 2.38e-4 | **6.17e-5** |
| 격차/eps² | 0.080 | 0.089 | 0.092 | 0.095 | 0.099 |

표준편차별 r (지속성 0.9): 4.074, 3.812, 3.404, 2.882%. 지속성별 r (sd 0.4): 3.956, 3.792, 3.475, 2.882%. 격차/eps² 는 아직 약간 증가 중이며 일부는 격자 효과입니다(2000 점에서 0.096).

**Aiyagari (1994) Table II 대조 (v0.2, `scripts/aiyagari_table2_v0_2.py`).** 원문 Table II 의 σ 는 earnings 의 변동계수(coefficient of variation)로 설명되고, 7상태 Markov 근사의 실제 모멘트가 입력값과 다릅니다(원문 Table I, 각주 38; 검토자 확인). 또 위험회피 μ = 1, 3, 5 열이 있고, 시험 세트의 EIS 0.5 (μ = 2) 와 같은 열은 없습니다. 그래서 두 셀을 μ 를 맞추고, 사슬(Tauchen 1986 의 폭 3 sd; Rouwenhorst)과 σ 해석(로그 sd; 변동계수 → 로그 sd = √log(1+σ²))을 바꿔 계산했습니다. β = 0.96, α = 0.36, δ = 0.08, 차입한도 0, 평균 노동 1 로 정규화 (차입한도 0 이면 가계 문제가 노동 수준에 대해 동차이므로 정규화가 r 을 바꾸지 않음), 2000 점, 상한 400.

| 셀 (σ, ρ, μ) | 원문 값* r / 저축률 (%) | Tauchen, 로그 sd | Tauchen, 변동계수 | Rouwenhorst, 로그 sd | Rouwenhorst, 변동계수 |
|---|---|---|---|---|---|
| 0.2, 0.6, 3 | 3.8767 / 24.25 | **3.8782** / 24.25 | 3.8842 | 3.9127 | 3.9179 |
| 0.4, 0.9, 5 | −0.3456 / 37.63 | −0.0856 / 36.39 | 0.1041 | 0.7265 | 0.9007 |

\* 원문 값은 검토 모형(GPT 6.1 Sol)이 원문 스캔(인쇄 p. 678)에서 읽어 옮긴 것입니다. 우리는 스캔본의 텍스트를 읽을 수 없어 직접 확인하지 못했습니다. 저축률은 s = δα/(r+δ) 로 r 의 함수이므로, 같은 r 에서 맞는 저축률은 독립된 두 번째 검증이 아닙니다.

- **첫 셀**은 Tauchen 사슬로 0.0016%p 이내에서 맞습니다(3.8782 대 3.8767).
- **둘째 셀은 재현되지 않았고 원인은 미해결입니다**(−0.0856 대 −0.3456, 약 0.26%p). 원문 Table I(인쇄 p. 676)은 사슬의 **수준 소득 CV 와 자기상관**을 보고하는데, Tauchen 3 sd 사슬의 값이 이미 원문의 두 자리 반올림 값과 일치합니다.

  | 셀 (σ, ρ) | 원문 Table I: CV / corr (검토 모형이 읽은 값) | Tauchen 3 sd 사슬: 수준 CV / 수준 corr |
  |---|---|---|
  | 0.2, 0.6 | 0.21 / 0.59 | 0.2142 / 0.5935 |
  | 0.4, 0.9 | 0.49 / 0.89 | 0.4885 / 0.8913 |

  따라서 v0.2 에서 쓴 "사슬의 실제 위험 차이로 설명된다"는 해석은 철회합니다. 두 모멘트가 같아도 전이행렬 전체가 같다는 보장은 없으므로 사슬 차이를 배제한 것도 아닙니다. 자산 격자를 바꿔도 차이는 닫히지 않습니다(r: 1000 점·상한 400 에서 −0.0862, 2000·400 에서 −0.0856, 4000·400 에서 −0.0855, 2000·800 에서 −0.0858%; 상단·격자 밖 질량 0). 후보 원인(원문의 모의 시계열 기반 계산, 전이행렬, 기타 구현 차이)은 모두 **[미확인]** 입니다.
- 참고로 Tauchen 폭을 2, 2.5, 3 으로 바꾸면 사슬의 로그 sd 가 0.399, 0.440, 0.468, r 은 0.915, 0.334, −0.086% 로 둘째 셀이 사슬에 매우 민감합니다. 원문 값에 맞추려고 사슬을 조정하지 않았습니다.
- 공개된 Econ-ARK Aiyagari 복제 프로젝트(HARK, DolARK)도 있으며, 그 코드의 정확도를 이번에 감사하지는 않았습니다.

---

## 8. T5: IKC 일관성 시험 (intertemporal Keynesian cross consistency)

**범위.** 같은 블록 야코비안으로 구성한 GE-IKC **일관성 시험**입니다. 독립적인 경제 모형 정확도 인증과 구별합니다. F1·F2 를 잡아냈다는 점에서 실용적입니다.

**명제.** HA 블록이 가계뿐인 모형에서 일반균형 산출은 dY = (I − M)^{-1} [dG + J^{C,Tax} dTax + J^{C,r} dr] 입니다 (M = J^{C,Y}; dTax, dr 은 일반균형 경로). 출처: Auclert, Rognlie, Straub (2024), JPE 132(12), 4068-4121 (NBER WP 25020).

**정확한 조건.** 가계: 차입한도 0, 소득 (Y − Tax) e, Rouwenhorst 7상태(지속성 0.92, sd 0.6), 500 점, σ 0.5. r = 0.02, B = 1 = A (β = 0.9105, 제약 질량 19.8%), G = 0.2, Tax_t = Tax + 0.1 (B_{t−1} − B). 고정 실질금리 모형과 실질금리 준칙 r_t = r + 0.25 (Y_t − 1) + 충격 모형. 미지수 Y, B. 충격 4개, T = 100, 200, 300. 풀이 3종: (a) 옵션을 일치시킨 SSJ, (b) 반복 정제를 더한 같은 알고리즘, (c) 옵션 불일치(F1).

**수치 결과 (상대 오차, T = 100 / 200 / 300).**

| 충격 | (a) SSJ 기본 | (b) 정제 | (c) 옵션 불일치 |
|---|---|---|---|
| G, 고정 금리 | 1.6e-15 / 1.5e-14 / 9.5e-14 | 5.3e-16 / 9.7e-15 / 1.9e-15 | 1.6e-15 / 1.5e-14 / 9.5e-14 |
| 금리, 고정 금리 | 3.9e-15 / 4.8e-14 / 3.1e-13 | 8.5e-16 / 1.7e-14 / 1.4e-13 | **9.9e-5** / 9.9e-5 / 9.9e-5 |
| G, 금리 준칙 | 5.0e-13 / 2.2e-9 / **5.9e-7** | 2.0e-15 / 3.2e-15 / 3.0e-15 | 5.0e-13 / 2.2e-9 / 5.9e-7 |
| 금리 충격, 금리 준칙 | 4.8e-13 / 2.3e-9 / **3.3e-7** | 3.1e-13 / 1.3e-12 / 9.8e-12 | **9.0e-5** / 9.0e-5 / 9.0e-5 |

- 허용치 1e-10. (a) 의 금리 준칙 T ≥ 200 과 (c) 의 금리 충격은 `xfail(strict=True)` 입니다. **이 엄격한 예상 실패는 SSJ 1.0.0 을 고정한 회귀 시험에만 해당합니다.** SSJ 가 고쳐지면 XPASS 가 되어 갱신을 알려 줍니다. 다른 풀이기를 위한 일반 규격에 "실패해야 한다"는 규칙을 옮기지 않습니다.
- F1·F2 의 원인과 최소 재현은 14절.
- 왈라스 법칙 잔차: 자산시장 9e-10 이하(정제 풀이, T=300).

---

## 9. T6: HANK − RANK 차분 항등식과 분해 (differencing identity and decomposition)

**명제.** A = I − M_HA, B = I − M_RA 이면 A^{-1} − B^{-1} = A^{-1}(B − A)B^{-1} (유한 T 리졸벤트 항등식). 같은 정상상태와 같은 정책 경로(dG, dTax, dr)에서 dY_HA − dY_RA = (I−M_HA)^{-1} [ (M_HA − M_RA)(dY_RA − dTax) + (J^r_HA − J^r_RA) dr ].

**범위.** 유한 T 의 리졸벤트 항등식과 같은 경로에 대한 차분 분해를 검증합니다. 무한 시계에서의 유일한 집단별 일반균형 귀속을 검증한 것이 아닙니다(1절 주요 발견 5).

**T6 시험이 각각 확인하는 것 (v0.3 설명).** 행렬 항등식은 가역 행렬이면 언제나 성립하므로 모형이 아니라 선형대수를 확인하는 검사입니다. 오류 주입 시험 v0.1 에서도 아무것도 잡지 않았습니다. 분해 시험은 경로와 야코비안의 일관성을 확인합니다. RA 닫힌 해(Ricardian 승수 1)와 HA 승수 조건(> 1)은 v0.2.1 까지 한 시험이었고 v0.3 부터 별도 시험입니다. cohort 시뮬레이션과 첫 라운드의 T 무관성은 별도 계산 경로지만 같은 가계 코드·정책·분포를 쓰므로 외부 기준이 아닙니다.

**정확한 조건.** T5 의 고정 실질금리 재정 HANK. RA 는 같은 r, σ, C, A(= B) 의 닫힌 해(β_RA = 1/(1+r)). 충격은 부채 조달 G(0.01 × 0.8^t)와 금리 −0.01 × 0.8^t. T = 150, 300.

**수치 결과.**
- 항등식 오차: T=150 3.8e-15, T=300 2.2e-14. 조건수: cond(I − M_RA) 31 / 1131, cond(I − M_HA) 294 / 5933.
- 분해 오차: G 충격 9.8e-12 / 1.9e-10, 금리 충격 1.1e-11 / 2.2e-10.
- RANK 승수는 리카도 동등성 때문에 1 이며, T=150 의 1 + 6.4e-6 은 잘린 부채 경로 꼬리(0.92^T) 때문입니다(T=300 에서 1 + 2.5e-11). HANK 충격 시점 승수 1.979, 현재가치 누적 승수 3.383.

| t | 0 | 1 | 2 | 3 | 4 | 5 |
|---|---|---|---|---|---|---|
| dY_HA − dY_RA | 0.00979 | 0.00932 | 0.00880 | 0.00825 | 0.00770 | 0.00716 |
| 첫 번째 단계 | 0.00316 | 0.00274 | 0.00231 | 0.00190 | 0.00151 | 0.00115 |
| 일반균형 증폭 | 0.00664 | 0.00658 | 0.00648 | 0.00635 | 0.00619 | 0.00601 |

- 부 사분위별 첫 단계 격차(t=0): Q1 1.05e-3, Q2 8.7e-4, Q3 6.5e-4, Q4 5.9e-4. 합계는 야코비안 총계와 3.2e-6 이내, T=150 과 T=300 에서 1e-8 이내로 같습니다.
- 첫 단계 격차의 현재가치: Q1-Q4 에서 −2.09e-3, −8.5e-4, 3.5e-4, 2.58e-3 (합은 0). 집단별로 (I−M_HA)^{-1} 를 적용한 값은 T 에 따라 커집니다(최저 사분위 충격 시점 T=150: −2.0e-3, T=300: −5.8e-3, T=450: −8.7e-3). 그래서 집단별로는 첫 단계까지만 보고하고 일반균형 증폭은 합계로만 보고합니다.

---

## 10. 자신의 HANK 코드를 연결하는 방법 (plugging in your own code)

검사 함수(`hank_limits/checks.py`)는 SSJ 에 의존하지 않고 provider 함수만 받습니다. 규약은 `hank_limits/interface.py` 에 있습니다. 야코비안은 T × T 행렬(행 = 결과 시점 t, 열 = 충격 시점 s), 단위는 수준, dr 은 절대 단위(0.01 = 1%p)이고, **시점 규약은 SSJ 와 같아야 합니다.**

| 검사 | provider 서명 | 반환 |
|---|---|---|
| `check_rank_limit` | `f(eps, T)` | `J={'C': {'Y', 'r'}}`, `ss={'C','A','r','eis'}` |
| `check_tank_limit` | `f(eps, T, dr, chis)` | `M_H, M_S`, `dY={chi: path}`, `lam`, `ss` |
| `check_werning` | `f(b, T, dr, cyclical)` | `dY, r_ss, eis, e_grid, Pi, zeta`, `diag` |
| `check_aiyagari` | `f(sd, rho)` | `r, rho_disc`, 선택 `diag` (잔차) |
| `check_ikc`, `check_differencing` | `f(T, shocks)` | `J`, `ss`, `irfs` |

v0.3.0: 공급자는 선택 필드 `ss_detail`, `J_A`, `J_C`, `ge_target`, `policy`, `irfs[label]` 의 `dC`·`dB`·`drshock` 도 보고할 수 있습니다. `hank_limits/interface.py` 에 공급자별로 검사가 읽는 모든 필드와, 선택 필드마다 쓰이는 V 수량을 적었습니다. V 검사는 이 값들로 정상상태 잔차, 소득 정규화, 명세 고정을 계산해 `spec=` 로 받은 값과 비교하고(이 세트는 `hank_limits/spec.py` 를 넘깁니다. 다른 보정이면 그 보정의 값을 넘기십시오), 공급자가 실제로 푼 외생 경로를 요청한 경로와 대조합니다. V 수량은 자기 식에 필요한 필드로만 계산하므로, 필드가 없으면 그 수량만 빠지고 이유는 `v0_3` → `skipped` 에 남습니다. V 시험은 유한하지 않거나 음수이거나 임계값을 넘는 값이 하나라도 있으면 실패하고, 실패 없이 값이 빠져 있으면 incomplete 로 skip 하며, 공급자가 그 시험의 수량을 하나도 보고하지 않으면 이유와 함께 skip 합니다(`checks.v03_verdict`). T3 에서는 풀이가 청산을 부과한 시장을 `ge_target` 으로 보고하십시오(`'asset_mkt'` 또는 `'goods_mkt'`). 각 검사는 요청한 경로와 실험의 사본을 공급자에게 넘기고, 비교는 따로 보존한 자기 사본으로 합니다.

권장 사항: T1·T4 는 단일 유형 가계 블록 단위 시험으로, T2 는 유형 수를 일반화한 형태로, T5·T6 은 자신의 일반균형 풀이기의 회귀 시험으로 쓰십시오. 숫자 허용치는 이 보정의 값이므로 자신의 보정에 맞춰 다시 정하십시오. 시점 규약이나 단위가 다르면 대부분 eps = 0 정확도 시험(T1, 1e-9 수준)에서 바로 드러납니다.

---

## 11. pytest 출력 요약 (actual pytest output)

아래는 기록한 그대로 둔 v0.2.1 의 출력(41 passed, 2 xfailed)이고, 전체 출력은 `out/pytest_output_v0_2_1.txt` 에 있습니다. T1–T6, D 시험과 그 결과는 그 뒤로 같습니다. 현재 판 v0.3.0 의 출력(114 passed, 4 xfailed: 이 시험들, V 검사, 검사기 회귀 시험)은 `out/pytest_output_v0_3_0.txt` 에 있습니다.

```
tests/test_t5_ikc.py::test_ikc_identity_ssj_solver_rule_long_horizon XFAIL
tests/test_t5_ikc.py::test_ikc_identity_with_mismatched_options XFAIL
tests/test_t5_ikc.py::test_walras_law PASSED
...
tests/test_diagnostics_v0_2.py::test_inverse_iteration_report[audit_t4] PASSED
tests/test_diagnostics_v0_2.py::test_inverse_iteration_report[audit_t1_mid] PASSED
tests/test_diagnostics_v0_2.py::test_ssj_forward_iteration_converges_to_it[audit_t4] PASSED
tests/test_diagnostics_v0_2.py::test_ssj_forward_iteration_converges_to_it[audit_t1_mid] PASSED
tests/test_diagnostics_v0_2.py::test_ssj_stopping_rule_is_not_an_accuracy_bound PASSED
...
=========================== short test summary info ============================
XFAIL tests/test_t5_ikc.py::test_ikc_identity_ssj_solver_rule_long_horizon - Finding F2 (SSJ 1.0.0): partial-pivoting LU on H_U packed in variable-block order shows element growth ~1e7 at T=300 (cond ~2e4); date-interleaved ordering or iterative refinement removes it.
XFAIL tests/test_t5_ikc.py::test_ikc_identity_with_mismatched_options - Finding F1 (SSJ 1.0.0): HetBlock linear impulses recompute Jacobians despite supplied Js (own options, default one-sided h=1e-4); with two-sided Js the r-shock IKC holds only to ~1e-4.
=================== 41 passed, 2 xfailed in 65.36s (0:01:05) ===================
```

---

## 12. FastHetBlock 분포 감사 (distribution audit, v0.2)

`scripts/diagnostics_v0_2.py` A, `out/diagnostics_v0_2.json`. 기준은 기본 설정(τ = 1e-10, 최대 20회, 다듬기 5회).

| 모형 (혼합 속도) | 역반복 횟수 | 음의 질량 | L1 잔차 | τ 1e-12 ~ 1e-6 에서 평균 자산 변화 | SSJ 전진 반복 tol 1e-10: 횟수 / 평균 자산 차이 | tol 1e-12 |
|---|---|---|---|---|---|---|
| T1 PIH eps = 0.1 (느림) | 3 | 0 | 1.2e-16 | ≤ 1.1e-14 | 37,910 / **1.7e-3** | 54,590 / 1.7e-5 |
| T1 PIH eps = 0.0125 (매우 느림) | 4 | 0 | 1.0e-16 | ≤ 1.1e-12 | 600,000 회 안에 정지 안 함 | — |
| T4 sd = 0.2 (빠름) | 3 | 0 | 1.8e-16 | ≤ 1.1e-15 | 3,010 / 1.1e-6 | 3,730 / 1.1e-8 |
| T5 재정 HANK (빠름) | 3 | 0 | 2.0e-16 | ≤ 2.2e-16 | 210 / 4.9e-9 | 250 / 7.5e-11 |

- 역반복 1회만 하고 다듬기를 생략하면 T1 eps = 0.0125 에서 평균 자산이 6.8e-4 다르고, 2회면 2.0e-8, 3회면 5.4e-13 입니다. 기본 설정은 종료 기준(1e-15)을 만족할 때까지 반복합니다.
- SSJ 전진 반복은 허용치를 조일수록 역반복 해로 수렴합니다. SSJ 의 정지 규칙(연속 반복의 차이)은 느린 혼합에서 정확도 한계가 아니므로, 근사 RA 보정에서 자산 목표로 β 를 맞출 때 주의가 필요합니다. 이것은 SSJ 결함 보고가 아니라 사용상 관찰입니다.
- SSJ 정지 규칙을 그대로 따라 반복 횟수를 세는 복제 루프가 SSJ 결과와 원소 단위로 같음을 확인했습니다.

---

## 13. 범위 밖이지만 함께 둔 것 (related material, not part of the test suite)

- `scripts/a1_markov_checks_v0_2.py`: 수반(adjoint) 오차 추정 설계 노트(A1)의 첫 점검. 작은 사슬의 전송항 포함 항등식, 실제 가계 사슬의 보존적 섭동에서 1차 예측 remainder 차수(2.00), 거친→세밀 정상 목표 오차의 항등식과 근사 수반 effectivity. **적응 격자 알고리즘을 구현하거나 검증한 것이 아닙니다.** 결과는 별도 설계 노트에 정리합니다.
- Dynare 의 이질적 가계 기능과 대조할 계획이면 **Dynare 7.1 이상**을 씁니다. 7.1 릴리스(2026-05-28)는 7.0 의 `model(heterogeneity=…)` 블록에서 complementarity 조건 부호가 반대로 해석되던 오류와, 이질적 블록에만 나오는 집계 변수를 인식하지 못하던 오류를 고쳤습니다. 7.0 을 기준으로 고정하면 차입제약 시험을 오염시킬 수 있습니다.

---

## 14. F1·F2·F3 최소 재현 (minimal reproductions)

**F1** (`scripts/repro_F1_v0_2.py`, SSJ 코드만): `sequence_jacobian.hetblocks.hh_sim.hh` 한 블록, T = 300.

| 비교 | 최대 상대 차이 |
|---|---:|
| `impulse_linear(..., Js={'hh': J2})` vs `J2 @ dr` (J2 양측) | 9.3e-5 |
| 같은 호출 vs `J1 @ dr` (J1 기본 단측) | 0 |
| `Js={'hh': 2 × J2}` 를 넣은 호출 vs `Js={'hh': J2}` 호출 | **0** |
| `impulse_linear(..., Js={'hh': J2}, twosided=True)` vs `J2 @ dr` | 0 |

**F2** (`scripts/repro_F2_v0_2.py`, SSJ 코드만, SSJ 자신의 전진 반복 정상분포 사용). `np.linalg.solve` 는 `Block.solve_impulse_linear` 와 같은 경로입니다.

| T | cond₂(H_U) | LU 원소 성장률 | 정규화 역오차 η (SSJ 경로) | 정제 2회 뒤 η |
|---:|---:|---:|---:|---:|
| 100 | 3.6e2 | 1.9e1 | 3.7e-14 | 5.7e-17 |
| 200 | 2.9e3 | 3.7e4 | 1.6e-10 | 7.6e-17 |
| 300 | 2.1e4 | 1.1e7 | 3.1e-8 | 1.0e-16 |

η = ‖H_U x − b‖_∞ / (‖H_U‖_∞ ‖x‖_∞ + ‖b‖_∞). T = 300 에서 날짜별 재배열 시 성장률 1.02, η 4.3e-16; 블록 순서 반전 시 성장률 1.00, η 9.3e-17. H_U 와 b 는 `out/F2_HU_T300_v0_2.npz` 에 있습니다.

**F3** (v0.3, SSJ 코드만). `scripts/repro_F3_v0_1.py` 는 y = e Y^(1+ζ) 에 대해 hetinput 미분을 `HetBlock` 과 같은 방식(`CombinedExtendedFunction([income]).differentiable(ss, h, True).diff({'Y': 1})`)으로 부릅니다. 실제로 쓰이는 미분은 전방 차분과 정확히 같고, 해석 미분과의 오차는 O(h) 입니다. `twosided=None` 으로 부르면 저장된 설정이 쓰여 오차가 O(h²) 가 됩니다.

```
h=0.01: |used-forward| 0.0e+00  |used-central| 2.4e-03  |used-analytic| 2.4e-03  |twosided=None - analytic| 8.4e-06
h=0.001: |used-forward| 0.0e+00  |used-central| 2.4e-04  |used-analytic| 2.4e-04  |twosided=None - analytic| 8.4e-08
h=0.0001: |used-forward| 0.0e+00  |used-central| 2.4e-05  |used-analytic| 2.4e-05  |twosided=None - analytic| 8.4e-10
```

`scripts/repro_F3_hetblock_v0_1.py` 는 SSJ 예제 블록 `hh_sim` 의 `HetBlock.jacobian(..., twosided=True)` 에서 효과를 끝까지 보여 줍니다. (H) 비선형 hetinput y = w^(1+ζ_e) e 에서 w 에 대한 A 의 야코비안, (O) 선형 hetinput 과 hetoutput c² 에서 r 에 대한 C2 의 야코비안입니다. "stored default" 는 `diff` 를 `twosided=None` 으로 부른 경우입니다. 두 경우 모두 차이가 O(h) 입니다.

```
H hetinput  max|J_original - J_stored_default| / max|J_stored_default|: h=1e-03 5.06e-05  h=1e-04 5.05e-06  h=1e-05 5.05e-07
O hetoutput max|J_original - J_stored_default| / max|J_stored_default|: h=1e-03 4.22e-04  h=1e-04 4.22e-05  h=1e-05 4.22e-06
```

코드 경로(sequence-jacobian 1.0.0 wheel): `blocks/het_block.py` 400–401행(hetinput 을 `True` 로 미분 가능하게 만듦)과 405행(hetoutput 은 호출자의 `twosided`), `utilities/function.py` 233행(`diff(..., twosided=False)`), `het_block.py` 302행(hetinput: `.diff({input_shocked: 1})`)과 378행(hetoutput: `.diff(..., outputs=...)`)은 모두 `twosided` 없이 부릅니다. `stage_block.py` 422–423행과 287행도 같습니다. GitHub master `eb4ac4c` 도 같은 코드라는 것은 연구 세션의 확인이며, v0.3 환경에서 다시 확인하지는 않았습니다. SSJ 에 [shade-econ/sequence-jacobian#41](https://github.com/shade-econ/sequence-jacobian/issues/41) 로 보고했습니다.

SSJ 관리자에게 보낸 보고: F1 [shade-econ/sequence-jacobian#39](https://github.com/shade-econ/sequence-jacobian/issues/39)·F2 [shade-econ/sequence-jacobian#40](https://github.com/shade-econ/sequence-jacobian/issues/40) (2026-10-20), F3 [shade-econ/sequence-jacobian#41](https://github.com/shade-econ/sequence-jacobian/issues/41) (2026-10-21 UTC). 제목은 관측된 동작으로 씁니다 ("HetBlock linear impulses recompute Jacobians despite supplied Js").

---

## 15. 그림과 표 (figures and tables)

- 그림: `out/t1_rank_limit.png` … `out/t6_decomposition.png` (로그 축, 극한이 오른쪽).
- 기계 판독 표: `out/tables_v0_2/suite_results_long_v0_2.csv` (열: test, limit_param, limit_value, n_a, a_max, fd_h, fd_twosided, window, horizon_T, solver, metric, value), `out/convergence_v0_2.json`, `out/diagnostics_v0_2.json`, `out/aiyagari_table2_v0_2.json`, `out/repro_F1_v0_2.json`, `out/repro_F2_v0_2.json`.

## 16. 인용 문헌 (references)

- Aiyagari, S. R. (1994). Uninsured Idiosyncratic Risk and Aggregate Saving. QJE 109(3), 659-684. https://ideas.repec.org/a/oup/qjecon/v109y1994i3p659-684..html
- Auclert, A., Bardóczy, B., Rognlie, M., Straub, L. (2021). Using the Sequence-Space Jacobian to Solve and Estimate Heterogeneous-Agent Models. Econometrica 89(5), 2375-2408. https://doi.org/10.3982/ECTA17434 (Definition 1 과 식 (23)이 기대 벡터를 명시하고, 논문 자체도 직접 계산·자동미분과의 야코비안 정확도 대조를 보고합니다. 이 세트는 그 검증을 경제적 극한으로 보완합니다.)
- Auclert, A., Rognlie, M., Straub, L. (2024). The Intertemporal Keynesian Cross. JPE 132(12), 4068-4121. https://doi.org/10.1086/732531 (NBER WP 25020)
- Bayer, C., Luetticke, R., Weiss, M., Winkelmann, Y. (2026). An Endogenous Gridpoint Method for Distributional Dynamics. JME 158, 103895. https://doi.org/10.1016/j.jmoneco.2026.103895 (복권 방식과 다른 분포 갱신의 비교 대상; 본문 대조는 하지 않음)
- Bilbiie, F. O. (2008). JET 140(1), 162-196. https://ideas.repec.org/a/eee/jetheo/v140y2008i1p162-196.html
- Bilbiie, F. O. (2020). The New Keynesian cross. JME 114, 90-108. https://doi.org/10.1016/j.jmoneco.2019.03.003
- Dynare Team (2026). Release of Dynare 7.1. https://forum.dynare.org/t/release-of-dynare-7-1/28376
- Huggett, M. (1993). JEDC 17(5-6), 953-969. https://ideas.repec.org/a/eee/dyncon/v17y1993i5-6p953-969.html
- Santos, M. S. (2000). Accuracy of Numerical Solutions Using the Euler Equation Residuals. Econometrica 68(6), 1377-1402. https://doi.org/10.1111/1468-0262.00165 (T4 Euler 잔차의 배경; 제약 가계·분포·GE 목표 전체에 대한 보장으로 인용하지 않음)
- Tauchen, G. (1986). Finite State Markov-Chain Approximations to Univariate and Vector Autoregressions. Economics Letters 20(2), 177-181.
- Werning, I. (2015). Incomplete Markets and Aggregate Demand. NBER WP 21448. https://www.nber.org/papers/w21448
- Young, E. R. (2010). Solving the Incomplete Markets Model with Aggregate Uncertainty Using the Krusell-Smith Algorithm and Non-Stochastic Simulations. JEDC 34(1), 36-41. https://doi.org/10.1016/j.jedc.2008.11.010
- Econ-ARK, Aiyagari (1994) replication. https://econ-ark.org/materials/aiyagari1994qje/
- sequence-jacobian: https://github.com/shade-econ/sequence-jacobian , https://pypi.org/project/sequence-jacobian/1.0.0/

## 17. AI 사용 공개 (AI usage disclosure)

이 프로젝트에는 생성형 AI 를 폭넓게 사용했습니다. Claude (Anthropic; Claude Opus 5.5, 2026년 10월)가 저자의 지시에 따라 코드·시험·수치 실험·문서 대부분을 작성하고 닫힌 해 기준값을 유도·점검했습니다. GPT 6.1 Sol (OpenAI)이 코드와 문서를 두 차례 검토했고, 그에 따른 변경은 이 README 앞부분에 정리되어 있습니다(이 문서의 "검토자"는 이 모형입니다). 저자는 연구 질문과 포함할 시험을 정하고, 작업을 지휘·검토했으며, 어떤 결과를 어떤 표현으로 보고할지 결정했고, 제출물 전체의 정확성·독창성·라이선스에 책임을 집니다.

## 18. 라이선스 (licence)

MIT. `LICENSE` 파일을 참조하십시오.

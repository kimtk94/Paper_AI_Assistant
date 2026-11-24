# IPF 세브란스 공동연구 — 의사결정 및 실행 기준 (2026-11-21)

**상태:** 내부 초안 완료 / 발송·SFDC 저장하지 않음 / PR 미병합

## 핵심 가설

기저 혈장 FCN3·GPX3와 이후 12개월 FVC 변화의 연관성을 한국인 IPF에서 확인한다. 혈장 PXD010965의 질병 상태 차이와 두 독립 폐 단일세포 데이터에서 RNA 감소를 관찰했지만, FVC 예후 연관성은 아직 **미검증**이다. 임상 검체 자료 없이 신규성이 입증됐다고 주장하지 않는다.

## 기존 정량 근거

- 세브란스 한국인 IPF 107명, 가족성 19명: 유전형 및 FVC/DLCO 종단 분석. https://doi.org/10.1186/s12931-026-03504-w
- PXD010965: IPF 17 / 대조군 20명, 공개 혈장 단백체 167개 재분석. 정규화 민감도가 커서 QC 방식을 원논문과 직접 동일시하지 않는다.
- GSE135893: IPF 12 / 정상 10 donor, Vanderbilt 동일 기관에서 FCN3 endothelial RNA 감소 방향.
- GSE136831: IPF 32 / 정상 28 donor, FCN3 endothelial log₂ ratio −1.644 (IPF 29 / 정상 20), GPX3 fibroblast −1.901 (24 / 12). 3개 *사전 선택 비교*에만 BH를 적용; 임상·배치 공변량 미조정.
- IPF-PRO 2025 MS cohort 299 IPF / 99 controls, 예후 단백체 연구가 이미 존재: https://doi.org/10.1186/s12931-025-03377-5

## 설계

1. **Primary:** (FVC %pred follow-up − baseline) × 365.25 / 일수; 사전 정의 9–15개월 구간. FVC mL/year 보조.
2. **Exposure:** Baseline log₂ plasma FCN3 and GPX3 (2 primary biomarkers, both two-sided), clinical adjustment including baseline FVC, age, sex, anti-fibrotic and smoking; baseline DLCO when available. Avoid covariate proliferation if N small.
3. **Discovery:** DIA global data with FDR and missingness sensitivity; no predictive classification promise.
4. **Technical verification:** initial 6–10 specimens for peptide detectability; then 60–80 analyzable exploratory patients if feasible; PRM/ELISA confirm candidate detectability.
5. **Sample size:** Fisher z approximate N=85 for |r|=.30, power .80, α=.05 under simple unadjusted correlation. Covariate adjustment, 2 candidates, missingness and prediction require more; 80 is feasibility not final clinical validation. At baseline coverage .80 and followup .80, 80 complete cases require 125 screened patients. The 107-person archival cohort may not suffice.
6. **Variables:** subject pseudonym, IPF diagnosis/adjudication, date-plasma and date-FVC, baseline/follow-up FVC% and mL, baseline DLCO, age, sex, smoking, treatment timing, HRCT, liver/kidney function, hemolysis, processing batch, ethics/consent. Keep patient-level data off GitHub.

## Gate

- GO_PILOT only after IRB/consent, same-person sample linkage, assay detectability, usable n≥60 (prefer 80), 9–15 month clinical data and missingness review.
- HOLD if sample/clinical linkage or coverage unknown, n<60, or assay nondetect.
- Screening QC: age must be finite numeric; required covariates cannot be missing/unknown placeholders. Explicit `none` for smoking or antifibrotic treatment remains a recorded value. Hemolysis must be explicitly negative (`no`, `n`, `false`, `0`); positive or unconfirmed status requires review. This adds no age eligibility cutoff or assay-detection claim.
- STOP/REDESIGN if only prevalent IPF vs healthy controls without progression data (already crowded research area).

## Why this is different from prior work

Not a newly discovered diagnostic signature. We prospectively specify FCN3 & GPX3 based on cross-cohort signal, and ask whether they relate to FVC *change* in a Korean IPF clinical cohort. Unlike the existing 2025 IPF-PRO prognostic study, the distinctive asset must be the Korean clinical cohort and independent replication; discovery alone is insufficient.

## Related internal Drive and GitHub

- IPF archival PXD010965: https://drive.google.com/drive/folders/17anMxeR0bbx1vwWf_rM6svRZpVatMmP4
- GSE136831 independent validation: https://drive.google.com/drive/folders/1l646hcA7vmVzIAUiIFr0GTmsz_fBLEfs
- Existing code PR: https://github.com/kimtk94/Paper_AI_Assistant/pull/23

**Commercial/territory safety:** Account/Contact Owner and prior contact history must be confirmed. No automatic outreach, no actual patient uploads to CI, no pricing or authorship guarantees.

## 2026-11-24 update

[CLINICAL_DESIGN_20261124.md](CLINICAL_DESIGN_20261124.md) supersedes earlier planning/approval language and distinguishes preliminary screening, assay qualification and outcome adjudication. It incorporates the #23 evidence replay, recent comparator literature, two-marker multiplicity planning, informative loss/death/transplant, and the v2 local-only readiness screen. **All actual cohort/clinical gates remain HOLD.** The old single-marker N=85 benchmark is not the two-marker design; at assumed |r|=.30 with conservative family-wise planning, 103 complete records are needed before covariate/attrition adjustments.

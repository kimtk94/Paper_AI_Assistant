# Paper AI Assistant 기획서 (논문 검색 + 후속 주제 발굴)

## 0) 비용 전략 (중요)

초기에는 **추가 요금 없이** 검증하는 것을 목표로 한다.

- 1순위: `codex login` 기반 ChatGPT OAuth 인증으로 시작
- API 키를 직접 코드/문서에 넣지 않고 인증 연동으로 운영
- 무료/기본 제공 한도 내에서 MVP 실험 후, 필요 시 유료 확장

- ⚠️ API Key 없이도 시작 가능하지만, 사용량 한도/정책은 플랜에 따라 달라질 수 있으므로 PoC 단계에서 모니터링

## 1) 목표

사용자가 키워드/질문을 입력하면:
1. 관련 논문을 검색하고,
2. 핵심 내용을 요약한 뒤,
3. 그 논문을 기반으로 **후속 연구 주제(다음에 파볼 만한 질문)**를 제안하는 AI를 만든다.

---

## 2) 핵심 사용자 시나리오

### 시나리오 A: 빠른 탐색
- 입력: "멀티모달 RAG 최신 논문 찾아줘"
- 출력:
  - 관련 논문 목록(제목, 저자, 연도, 링크)
  - 각 논문 3~5줄 요약
  - 공통 트렌드 정리
  - 후속 주제 Top 5

### 시나리오 B: 특정 논문 확장
- 입력: "이 논문 기반으로 다음 연구 아이디어 알려줘"
- 출력:
  - 논문의 가정/한계 분석
  - 검증 가능한 가설 후보
  - 데이터셋/실험 설계 제안
  - 난이도/임팩트 기준 우선순위

---

## 3) 기능 요구사항

### 3.1 논문 검색
- 키워드 기반 검색
- 연도, 도메인, 인용수(가능 시) 필터
- 중복 제거 및 관련도 순 정렬

### 3.2 논문 이해/요약
- 초록 기반 1차 요약
- 가능하면 본문/핵심 섹션 기반 2차 요약
- 방법론, 성능, 한계점 분리 요약

### 3.3 후속 주제 생성 (핵심)
- 논문의 한계/미해결 문제를 구조화
- 후속 연구 아이디어를 다음 포맷으로 생성:
  - 주제명
  - 왜 중요한가
  - 검증 가설
  - 실험 방법
  - 성공/실패 기준

### 3.4 결과 품질 제어
- 근거 없는 주장 최소화(출처 표시)
- 중복 아이디어 제거
- 실현 가능성 점수(1~5) 제공

---

## 4) 시스템 구성(초안)

1. **Query Planner**
   - 사용자 질문을 검색 쿼리로 변환
2. **Paper Retriever**
   - 논문 메타데이터/초록 수집
3. **Paper Analyzer (LLM: OpenAI / Codex OAuth 연동 우선)**
   - 핵심 기여/한계/공백 추출
4. **Topic Generator (LLM: OpenAI / Codex OAuth 연동 우선)**
   - 후속 주제 후보 생성
5. **Ranker**
   - 새로움/실험 가능성/영향도 기반 랭킹
6. **Report Builder**
   - Markdown 보고서 생성

---

## 5) 데이터 모델(간단)

## Paper
- id
- title
- authors
- year
- venue
- abstract
- url
- key_findings[]
- limitations[]

## FollowUpTopic
- topic_id
- title
- motivation
- hypothesis
- experiment_plan
- risk
- feasibility_score (1~5)
- impact_score (1~5)
- based_on_papers[]

---

## 6) 출력 포맷 (Markdown)

```md
# 검색 결과 요약

## 관련 논문
1. 제목 (연도) - 링크
   - 핵심 요약

## 공통 인사이트
- ...

## 후속 연구 주제 제안
### 1) 주제명
- 중요성:
- 가설:
- 실험 계획:
- 성공 기준:
- 리스크:
```

---

## 7) 개발 단계 제안

### Phase 1 (MVP)
- 키워드 검색 + 초록 요약 + 후속 주제 3개 생성
- 결과를 Markdown으로 저장

### Phase 2
- 논문 PDF 본문 일부 반영
- 주제 평가 점수 및 랭킹
- 실험 설계 템플릿 자동 생성

### Phase 3
- 사용자 피드백 기반 추천 개선
- 연구 로드맵(3개월/6개월) 자동 제안

---

## 8) MVP 완료 기준 (Definition of Done)

- 동일 질문에서 일관된 형식의 보고서 생성
- 최소 5개 논문 검색 및 요약
- 후속 주제 3개 이상 제시
- 각 주제에 가설 + 실험계획 + 리스크 포함
- 결과를 `.md` 파일로 저장 가능

---

## 9) 다음 액션

1. 검색 대상 소스(예: arXiv, Semantic Scholar 등) 우선순위 확정
2. OpenAI 프롬프트 템플릿 정의 (요약용/주제생성용)
3. MVP CLI 명령 설계 (`search`, `summarize`, `propose-topics`)
4. 샘플 질의 10개로 품질 평가 루프 구축

---

## 10) Multi-omics 구조화 스키마/파서 산출물

아래 파일을 추가해 multi-omics 논문/아이디어 구조화 데이터를 바로 검증할 수 있게 구성했다.

- `schemas/paper_record.schema.json`: Paper Record JSON Schema
- `schemas/idea_record.schema.json`: Idea Record JSON Schema
- `examples/multiomics_records.example.json`: 샘플 입력 템플릿
- `src/multiomics_models.py`: Pydantic 모델 + OpenClaw 입출력 파서(`parse_bundle`, `load_bundle`, `dump_bundle`)

빠른 확인:

```bash
python src/multiomics_models.py
```

`pydantic`이 설치되어 있으면 샘플 JSON 검증 후 레코드 개수를 출력한다.

---

## 11) 3-Assistant 실행 구조 (T2D + AD + Korea affiliation 체크)

`src/paper_assistants.py` 기준으로 아래 3개 Assistant를 한 파이프라인으로 분리했다.

1. **RetrievalAssistant**
   - OpenAlex API에서 `type 2 diabetes`, `Alzheimer's disease` 관련 논문 메타데이터/초록 수집
   - 저자 소속 문자열(`raw_affiliation_strings`) 및 기관 국가코드(`country_code=KR`)를 확인해
     `korea_affiliation_present` 필드를 계산
2. **SummarizerAssistant**
   - 초록 기반으로 핵심 요약/방법론 힌트/한계 코멘트 생성
3. **RelevanceReviewerAssistant**
   - 질병 키워드 매칭으로 주제 적합성(`is_relevant`) 점검

### 실행 예시

```bash
python src/paper_assistants.py \
  --max-results 10 \
  --output outputs/paper_assistant_report.json \
  --mailto you@example.com
```

출력 JSON의 핵심 필드:
- `disease_hits`: `type_2_diabetes`, `alzheimers_disease` 매칭 결과
- `korea_affiliation_present`: 저자 주소/소속 내 Korea 여부
- `korea_affiliation_evidence`: Korea 판단 근거 텍스트
- `summary`, `review`: 요약과 주제 적합성 평가 결과

### 권장 운영 플로우

- 매일/매주 배치로 `paper_assistants.py` 실행
- 결과 JSON을 `src/multiomics_ingest.py`로 append-only 적재
- 적재 데이터 기반으로 후속 아이디어 생성/랭킹(기존 `multiomics_models.py` 스키마 확장)

---

## 12) 논문 검토 효율화 템플릿 (Method / Data DB / Dataset 중심)

논문을 빠르게 검토할 때는 3단계로 진행한다.

1. **1차 스크리닝(5~10분)**: title/abstract/conclusion 위주로 문제-방법-데이터셋 매칭 확인  
2. **2차 구조 파악(20~30분)**: methods / data / experiment split / metric 확인  
3. **3차 심화(선별 논문만)**: ablation, leakage risk, external validation 확인

### 표준 정리 항목

- Method (핵심 접근법)
- Data DB (데이터 출처/저장소)
- Dataset (이름/버전/샘플 수/분할 여부)
- Relevance (타깃 질병 주제 적합성)

### 기존 `paper_assistant_report.json` 리스트 재검토 실행

이미 생성된 JSON 리스트를 대상으로 Method/Data/Dataset 항목을 다시 정리하려면:

```bash
python src/paper_assistants.py \
  --from-report outputs/paper_assistant_report.json \
  --output outputs/paper_assistant_report_reviewed.json \
  --review-md-output outputs/paper_assistant_review.md
```

산출물:
- `outputs/paper_assistant_report_reviewed.json`: `method_data_review` 필드가 추가된 리뷰 결과
- `outputs/paper_assistant_review.md`: Method / Data DB / Dataset / Relevance 표 형식 리포트

---

## 13) 지금까지의 작업 업데이트 (2026-04-11)

현재 저장소 기준으로 완료된 작업을 운영 관점에서 요약하면 아래와 같다.

### 완료된 산출물

- **기획/운영 문서화**
  - `PAPER_AI_ASSISTANT_SPEC.md`: 전체 제품 기획, 3-Assistant 파이프라인, multi-omics 확장 전략 정리
  - `OPENCLAW_SETUP.md`: OpenClaw + Codex CLI + systemd까지 원클릭/단계별 설치 가이드
  - `MULTIOMICS_OPERATIONS.md`: 단건 적재/백필/주기 실행 운영 절차

- **스키마/모델 기반 데이터 구조화**
  - `schemas/paper_record.schema.json`, `schemas/idea_record.schema.json`으로 레코드 구조 고정
  - `src/multiomics_models.py`에 Pydantic 모델 및 번들 입출력 파서 구현
  - `examples/multiomics_records.example.json` 샘플 입력 제공

- **적재 파이프라인 자동화**
  - `src/multiomics_ingest.py` 기반 append-only 적재 로직 구성
  - `ops/07_ingest_multiomics_once.sh`(1회 적재), `ops/08_backfill_multiomics.sh`(백필) 제공
  - `ops/multiomics_ingest.service`, `ops/multiomics_ingest.timer`, `ops/09_setup_multiomics_timer.sh`로 주기 실행 구성

- **논문 수집/요약/검토 파이프라인 구현**
  - `src/paper_assistants.py`에 Retrieval/Summarizer/RelevanceReviewer 3-Assistant 흐름 구현
  - T2D + Alzheimer's 질의 및 Korea affiliation 체크(`korea_affiliation_present`) 포함
  - 기존 리포트 재검토 모드(`--from-report`)와 Markdown 리뷰 출력(`--review-md-output`) 지원

### 현재 운영 가능한 흐름

1. OpenClaw 환경 부트스트랩 (`ops/00_bootstrap_openclaw.sh`)
2. 논문 수집/요약 실행 (`src/paper_assistants.py`)
3. 결과를 multi-omics JSONL로 적재 (`src/multiomics_ingest.py` + ops 스크립트)
4. 필요 시 systemd timer로 배치 자동화 (`ops/09_setup_multiomics_timer.sh`)

### 다음 우선순위 제안

- **중복 제어 강화**: `paper_id`/`idea_id` 기준 dedup 단계 추가
- **품질 메트릭 도입**: 요약 품질/주제 적합성 점수의 추세 모니터링
- **아이디어 랭킹 고도화**: novelty/feasibility/impact 가중치 튜닝
- **리포트 표준화**: 주간 운영 리포트 템플릿(성공/실패 케이스 포함) 정착

## 12) PDF 본문 텍스트 추출 + 연차/질병 분류 스크립트

다운로드한 PDF를 실제로 "연차별 + 질병별"로 정리할 수 있도록 `src/pdf_disease_year_organizer.py`를 추가했다.

### 제공 기능
- PDF 텍스트 추출: `pypdf` 우선, 없으면 `pdfplumber` → `pymupdf(fitz)` → `pdftotext` 순 fallback
- 질병 분류: 텍스트 내 키워드 기반 multi-label 매칭
- 연도 추정: 파일명 우선 + 본문 텍스트 내 연도 후보 보정
- 산출물 생성:
  - 정리된 PDF 폴더: `<output>/<disease>/<year>/...pdf`
  - 추출 텍스트 파일: `<text_dir>/*.txt`
  - JSON 리포트: 파일별 추출/분류 상세
  - Markdown 요약: 질병별/연도별 분포와 실패 목록

### 실행 예시

```bash
python src/pdf_disease_year_organizer.py \
  --pdf-dir data/pdfs \
  --output-dir data/organized_pdfs \
  --text-dir data/extracted_text \
  --report-json data/pdf_organize_report.json \
  --summary-md data/pdf_organize_summary.md
```

### 질병 키워드 커스터마이징

기본 키워드 대신 JSON 파일로 커스텀 매핑을 넣을 수 있다.

```json
{
  "type_2_diabetes": ["type 2 diabetes", "t2d", "insulin resistance"],
  "alzheimers_disease": ["alzheimer", "amyloid", "dementia"]
}
```

```bash
python src/pdf_disease_year_organizer.py \
  --pdf-dir data/pdfs \
  --disease-keywords-json config/disease_keywords.json
```

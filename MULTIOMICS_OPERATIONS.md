# Multi-omics 데이터 적재/주기 실행 가이드

질문하신 3가지를 기준으로 바로 운영 가능한 형태로 정리했습니다.

## 1) 세팅을 어떻게 하면 좋을까?

### 권장 디렉터리
- 입력(최신 파일): `data/incoming/latest_bundle.json`
- 입력(백필 여러 파일): `examples/backfill_inputs/*.json`
- 출력 누적 데이터: `data/multiomics/`
  - `paper_records.jsonl`
  - `idea_records.jsonl`
  - `snapshots/bundle_*.json`
  - `manifests/manifest_*.json`

### 의존성
`pydantic` 기반 엄격 검증이 필요하면:

```bash
pip install -r requirements-multiomics.txt
```

## 2) 데이터를 한번에 쌓을 수 있나?

가능합니다.

### 단일 파일 적재
```bash
sh ops/07_ingest_multiomics_once.sh data/incoming/latest_bundle.json
```

### 여러 파일 백필
```bash
sh ops/08_backfill_multiomics.sh examples/backfill_inputs data/multiomics
```

이 방식은 append-only(JSONL)라서, 과거 데이터 누적/감사 추적에 유리합니다.

## 3) 주기적으로 돌아가도록 할 수 있나?

가능합니다(systemd timer).

### 타이머 설치
```bash
sh ops/09_setup_multiomics_timer.sh
```

기본 주기: 6시간마다 실행(`ops/multiomics_ingest.timer`).

### 동작 방식
- 타이머가 `multiomics_ingest.service`를 호출
- 서비스는 기본적으로 아래 입력을 읽음:
  - `data/incoming/latest_bundle.json`
- 엄격 검증 모드(`--strict-validation`)로 적재

## 운영 팁
- 업스트림 수집기(예: PubMed/EuropePMC fetcher)는 최신 번들을 항상 `data/incoming/latest_bundle.json`로 원자적으로 교체(`mv`)하도록 구성하세요.
- 중복 방지가 필요하면 다음 단계에서 `paper_id`/`idea_id` 기준 dedup 잡을 추가하면 됩니다.

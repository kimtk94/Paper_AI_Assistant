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
`pydantic`을 설치하면 기본 실행에서도 레코드 스키마를 검증합니다:

```bash
pip install -r requirements-multiomics.txt
```

입력은 `papers`와 `ideas` 배열을 담은 번들 JSON입니다. 둘 중 하나만 포함할 수 있으며, 각 배열의 항목은 레코드 객체여야 합니다. `examples/multiomics_records.example.json`을 참고하세요. `paper_assistants.py`의 `results` 보고서는 구조와 필수 필드가 다르므로 직접 적재할 수 없습니다. 확인된 논문 정보를 번들 스키마에 맞게 작성한 뒤 적재해야 합니다.

잘못된 번들 구조와 설치된 Pydantic의 레코드 검증 오류는 저장 전에 실패합니다. Pydantic이 없는 기본 실행만 경고와 함께 구조 검증으로 진행합니다. 운영 시에는 아래와 같이 `--strict-validation`을 사용해 Pydantic이 없을 때도 중단하도록 하세요.

```bash
python3 src/multiomics_ingest.py --input data/incoming/latest_bundle.json --data-dir data/multiomics --strict-validation
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

각 실행에는 UTC 타임스탬프와 UUID를 합친 `run_id`가 부여됩니다. 같은 초에 여러 파일을 백필해도 `bundle_<run_id>.json`과 `manifest_<run_id>.json`을 각각 보존하며 기존 감사 파일을 덮어쓰지 않습니다. 레코드의 `ingested_at`은 기존 UTC 타임스탬프 형식을 유지합니다.

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
- JSONL 두 파일과 감사 파일 전체가 하나의 트랜잭션으로 저장되는 구조는 아닙니다. 쓰기 실패 후 재시도 전에는 스냅샷·매니페스트·누적 행을 확인하세요. 동일 번들을 다시 적재하면 행이 추가됩니다.

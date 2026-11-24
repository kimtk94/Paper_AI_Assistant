# Paper AI Assistant

논문 메타데이터 수집, 연구자·주제 연결, PDF 정리, 검토한 multi-omics 기록 적재를 위한 Python 도구입니다. 현재 요약·분류에는 초록과 키워드에 기반한 규칙이 포함됩니다. 연구 가설과 통계 결과의 해석은 원문 및 데이터 출처를 함께 검토해야 합니다.

## 구성

| 영역 | 진입점 | 역할 |
|---|---|---|
| 논문 검색·검토 | `src/paper_assistants.py` | OpenAlex 등에서 메타데이터 수집, 질병·방법론 단서 정리, PDF 확보 상태 기록 |
| SKKU 문헌·연구자 연결 | `src/skku_pubmed_*.py` | PubMed 소속 검색, ORCID 추적, 논문·연구자·커뮤니티 연결 |
| PDF 분류 | `src/pdf_disease_year_organizer.py` | 질병·연도별 분류 |
| 구조화 기록 적재 | `src/multiomics_models.py`, `src/multiomics_ingest.py` | `papers`/`ideas` 스키마 검증, JSONL 및 실행별 감사 기록 |
| 서버 운영 예시 | `ops/`, `OPENCLAW_SETUP.md` | 설치·서비스·배치 실행 절차 |
| 연구 초안 | [PR #21–25](https://github.com/kimtk94/Paper_AI_Assistant/pulls) | 기본 브랜치와 별도로 검토 중인 공개자료 분석·임상 타당성 프로토콜 |

## 설치 및 로컬 검증

Python 3.11 또는 3.12 기준입니다.

```bash
python -m venv .venv
. .venv/bin/activate
python -m pip install -r requirements.txt
python -m unittest discover -s tests -p 'test_*.py' -v
```

테스트는 합성 자료·임시 디렉터리·모의 응답을 사용합니다. 실제 PubMed/GEO 분석이나 서버 설치를 실행하지 않습니다. 의존성 설치에는 패키지 다운로드가 필요합니다.

## 사용 문서

- [SKKU PubMed 실행 및 해석](SKKU_PUBMED_LINEAGE.md)
- [Multi-omics 적재 및 운영](MULTIOMICS_OPERATIONS.md)
- [기획서와 논문 수집 CLI 예시](PAPER_AI_ASSISTANT_SPEC.md)
- [2026-11-24 저장소·연구 PR 검토 기록](docs/review/2026-11-24.md)

`paper_assistants.py`의 `results` 보고서와 `multiomics_ingest.py`의 `papers`/`ideas` 번들은 구조가 다릅니다. 보고서를 검토하고 [예시 번들](examples/multiomics_records.example.json)의 필수 정보를 작성한 후 적재하세요.

## CI와 외부 자료 검증

`SKKU PubMed CI`는 Python 3.11/3.12에서 모든 기본 코드 테스트, 문법, 노트북 구조를 확인합니다. 외부 PubMed 통합 검증은 Actions의 수동 실행에서 `run_live_smoke`를 선택할 때 실행합니다. PR 검증은 외부 검색 결과나 당일 서비스 상태에 의존하지 않습니다.

PubMed 응답에 오류·누락·중복이 있으면 정상적인 빈 검색 결과로 처리하지 않습니다. XML 캐시는 요청한 논문 ID와 다시 대조하며, 완료되지 않은 검색 구간은 완료 체크포인트에 저장하지 않습니다. 과거 version 1 체크포인트는 검색을 다시 수행하고, 검증에 통과한 XML 캐시는 재사용합니다. 현재 XML 파서는 `PubmedArticle`을 지원하며 지원하지 않는 레코드는 조용히 누락하지 않고 오류를 표시합니다.

연구 PR의 Actions 성공은 해당 작업이 종료됐다는 뜻입니다. 데이터 미확보에 따른 `HOLD`, 분석 코드의 기술적 오류, 실제 연구 가설의 지지 여부는 [검토 기록](docs/review/2026-11-24.md)에서 구분합니다.

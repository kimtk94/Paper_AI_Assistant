# OpenClaw 서버 세팅 가이드

실행은 모두 `sh` 스크립트로 제공하고, 이 문서에는 **사용 방법만** 정리합니다.

## 1) 전체 초기 세팅 (권장)

```bash
sh ops/00_bootstrap_openclaw.sh
```

- 포함 작업: `git pull` → OpenClaw 설치 확인/설치 → 인증 선택(OAuth 또는 API Key) → systemd 등록/기동

## 2) 단계별 수동 실행

### 설치
```bash
sh ops/01_install_openclaw.sh
```

### Codex CLI 설치
```bash
sh ops/02_install_codex_cli.sh
```

### ChatGPT OAuth 로그인 (API Key 없이)
```bash
sh ops/03_login_codex_oauth.sh
```

### API Key 설정
```bash
sh ops/04_set_openai_key.sh '<YOUR_OPENAI_API_KEY>'
```

### systemd 등록/기동
```bash
sh ops/05_setup_systemd_openclaw.sh
```

### 업데이트 + 재시작
```bash
sh ops/06_update_and_restart_openclaw.sh
```

## FAQ

### Q) API Key 없이도 가능한가?
가능합니다. 아래 스크립트로 로그인하면 됩니다.

```bash
sh ops/02_install_codex_cli.sh
sh ops/03_login_codex_oauth.sh
```

단, 플랜/시점에 따라 사용량·모델 제한이 달라질 수 있습니다.

### 설치 에러가 `set: Illegal option -o pipefail` 로 나는 경우
OpenClaw 설치 스크립트는 내부적으로 bash를 요구할 수 있습니다.
최신 스크립트(`sh ops/01_install_openclaw.sh`)는 자동으로 `bash` 파이프 실행을 사용합니다.

## 3) Multi-omics 적재 운영

### 1회 적재
```bash
sh ops/07_ingest_multiomics_once.sh data/incoming/latest_bundle.json
```

### 백필(여러 JSON)
```bash
sh ops/08_backfill_multiomics.sh examples/backfill_inputs data/multiomics
```

### 주기 실행(systemd timer)
```bash
sh ops/09_setup_multiomics_timer.sh
```

상세 운영 설명은 `MULTIOMICS_OPERATIONS.md`를 참고하세요.

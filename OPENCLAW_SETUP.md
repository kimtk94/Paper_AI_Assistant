# OpenClaw 서버 세팅 가이드

실행은 모두 `sh` 스크립트로 제공하고, 이 문서에는 **사용 방법만** 정리합니다.

## 1) 전체 초기 세팅 (권장)

```bash
sh ops/bootstrap_openclaw.sh
```

- 포함 작업: `git pull` → OpenClaw 설치 확인/설치 → 인증 선택(OAuth 또는 API Key) → systemd 등록/기동

## 2) 단계별 수동 실행

### 설치
```bash
sh ops/install_openclaw.sh
```

### ChatGPT OAuth 로그인 (API Key 없이)
```bash
sh ops/login_codex_oauth.sh
```

### API Key 설정
```bash
sh ops/set_openai_key.sh '<YOUR_OPENAI_API_KEY>'
```

### systemd 등록/기동
```bash
sh ops/setup_systemd_openclaw.sh
```

### 업데이트 + 재시작
```bash
sh ops/update_and_restart_openclaw.sh
```

## FAQ

### Q) API Key 없이도 가능한가?
가능합니다. 아래 스크립트로 로그인하면 됩니다.

```bash
sh ops/login_codex_oauth.sh
```

단, 플랜/시점에 따라 사용량·모델 제한이 달라질 수 있습니다.

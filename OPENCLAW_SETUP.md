# OpenClaw 서버 세팅 가이드 (Codex OAuth / ChatGPT 로그인 우선)

추가 요금 없이 먼저 시도할 수 있도록, **API 키 직접 입력 대신 `codex --login` (ChatGPT OAuth)** 경로를 기본으로 정리했습니다.

## 1) 서버에서 코드 최신화

```bash
git pull --ff-only
```

## 2) OpenClaw 설치

> 보안상 URL을 꼭 확인한 뒤 실행하세요.

```bash
curl -fsSL https://openclaw.ai/install.sh | bash
```

설치 후:

```bash
openclaw --version
```

## 3) 인증 방식 선택

### 권장: ChatGPT OAuth 로그인 (키 수동 입력 없음)

```bash
codex --login
```

- 브라우저에서 **Sign in with ChatGPT** 진행
- CLI가 API 자격정보를 자동 연결/생성
- 키를 직접 복붙하지 않아도 됨

### 대안: OpenAI API 키 직접 사용

```bash
cp .env.openclaw.example .env.openclaw
bash ops/set_openai_key.sh '<YOUR_OPENAI_API_KEY>'
```

## 4) 부팅 후 자동 기동(systemd)

```bash
sudo cp ops/openclaw.service /etc/systemd/system/openclaw.service
sudo systemctl daemon-reload
sudo systemctl enable --now openclaw
sudo systemctl status openclaw
```

## 5) 업데이트 루틴

```bash
bash ops/update_and_restart_openclaw.sh
```

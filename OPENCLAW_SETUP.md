# OpenClaw 서버 세팅 가이드 (OpenAI 전용)

이 저장소를 서버에서 `git pull` 한 뒤 바로 진행할 수 있도록, 최소 실행 순서로 정리했습니다.

## 1) 서버에서 코드 최신화

```bash
git pull --ff-only
```

## 2) OpenClaw 설치(공식 설치 스크립트 방식)

> 보안상 **URL/도메인 오타**를 꼭 확인한 뒤 실행하세요.

```bash
curl -fsSL https://openclaw.ai/install.sh | bash
```

설치 후 셸을 다시 열고 버전 확인:

```bash
openclaw --version
```

## 3) 환경변수 파일 준비 (OpenAI만 사용)

```bash
cp .env.openclaw.example .env.openclaw
```

### OpenAI 키 적용 방법

키를 직접 파일에 넣어도 되고, 아래 스크립트로 적용해도 됩니다.

```bash
bash ops/set_openai_key.sh '<YOUR_OPENAI_API_KEY>'
```

> `.env.openclaw`는 `.gitignore`에 포함되어 있어 git에 커밋되지 않습니다.

## 4) 부팅 후 자동 기동(systemd)

원하면 `ops/openclaw.service`를 사용해 systemd 서비스로 등록하세요.

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

위 스크립트는 저장소 업데이트 후 OpenClaw 서비스를 재시작합니다.

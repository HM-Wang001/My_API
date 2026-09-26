
import json
import time
import requests
from fastapi import FastAPI, HTTPException
from pydantic import BaseModel, Field

# 导入你现有的登录模块
from zepp登录接口 import zepp_login

app = FastAPI(title="Zepp 步数 API")

class StepRequest(BaseModel):
    account: str = Field(..., description="Zepp 账号")
    password: str = Field(..., description="密码")
    steps: int = Field(..., ge=1000, le=98880, description="目标步数")

def submit_steps(user_id, app_token, steps):
    timestamp = int(time.time())
    url = f"https://api-mifit-cn.huami.com/v1/data/band_data.json?&t={timestamp}"

    headers = {
        "apptoken": app_token,
        "User-Agent": "MiFit/4.6.0 (iPhone; iOS 14.0.1; Scale/2.00)"
    }

    today = time.strftime("%Y-%m-%d")
    data_json = json.dumps([{
        "date": today,
        "data_hr": "",
        "data": [{
            "stop": 1439,
            "value": "",
            "did": "DA932FFFFE8816E7",
            "tz": 32,
            "src": 17,
            "start": 0
        }],
        "summary": json.dumps({
            "stp": {"ttl": steps, "dis": 144, "cal": 6, "wk": 5},
            "v": 5,
            "goal": 8000
        })
    }])

    post_data = {
        "data_json": data_json,
        "userid": user_id,
        "device_type": "0",
        "last_sync_data_time": "1589917081",
        "last_deviceid": "DA932FFFFE8816E7"
    }

    try:
        response = requests.post(url, headers=headers, data=post_data, timeout=10)
        print(f"HTTP 状态码: {response.status_code}")
        print(f"服务器原始返回: {response.text}")

        if response.status_code == 200:
            resp_json = response.json()
            if resp_json.get("code") == 1:
                return True, f"步数更新成功: {steps}"
            else:
                return False, resp_json.get("message", "未知错误")
        else:
            return False, f"HTTP {response.status_code}"
    except requests.exceptions.Timeout:
        return False, "请求超时"
    except requests.exceptions.RequestException as e:
        return False, f"请求异常: {e}"

@app.post("/api/update-steps")
async def update_steps(req: StepRequest):
    login_result = zepp_login(req.account, req.password)
    if not login_result or not login_result.app_token:
        raise HTTPException(status_code=401, detail="登录失败")

    ok, msg = submit_steps(login_result.user_id, login_result.app_token, req.steps)
    if ok:
        return {"success": True, "message": msg}
    else:
        raise HTTPException(status_code=500, detail=msg)

@app.get("/")
async def root():
    return {"status": "ok", "message": "Zepp 步数 API 运行中"}

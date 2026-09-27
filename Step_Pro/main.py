import os
from flask import Flask, request, jsonify

app = Flask(__name__)

# 从环境变量读取配置（你之前用的 CONFIG JSON 格式）
import json

CONFIG = json.loads(os.environ.get("CONFIG", "{}"))

# 导入 mimotion 的核心逻辑（你下载的仓库文件）
from zepp_helper import (
    login_access_token,
    grant_login_tokens,
    get_user_device_id,
    post_fake_brand_data,
)

API_KEY = "0120"  # 改成你自己的密钥


@app.route("/api/update-steps", methods=["GET"])
def update_steps():
    # 1. 密钥校验
    key = request.args.get("key")
    if key != API_KEY:
        return jsonify({"success": False, "message": "密钥错误"}), 403

    # 2. 获取参数
    account = request.args.get("account") or CONFIG.get("USER")
    password = request.args.get("password") or CONFIG.get("PWD")
    steps = request.args.get("steps", type=int) or 20000

    if not account or not password:
        return jsonify({"success": False, "message": "缺少账号或密码"}), 400

    try:
        # 3. 登录获取 access_token
        import uuid
        device_id = str(uuid.uuid4())
        access_token, err = login_access_token(account, password)
        if err:
            return jsonify({"success": False, "message": f"登录失败: {err}"}), 401

        # 4. 换取 login_token、app_token、user_id
        login_token, app_token, userid, err = grant_login_tokens(
            access_token, device_id, is_phone=False
        )
        if err:
            return jsonify({"success": False, "message": f"换取 token 失败: {err}"}), 401

        # 5. 查询真实设备 ID（解决新账号同步问题的关键）
        real_device_id = get_user_device_id(app_token, userid)

        # 6. 提交步数
        ok, msg = post_fake_brand_data(str(steps), app_token, userid, real_device_id)
        if ok:
            return jsonify({
                "success": True,
                "message": f"步数更新成功: {steps}",
                "user_id": userid,
                "device_id": real_device_id or "默认值"
            })
        else:
            return jsonify({"success": False, "message": f"提交失败: {msg}"}), 500

    except Exception as e:
        import traceback
        return jsonify({"success": False, "message": f"执行异常: {str(e)}", "trace": traceback.format_exc()}), 500


@app.route("/")
def root():
    return jsonify({"status": "ok", "message": "Zepp 步数 API (Flask) 运行中"})
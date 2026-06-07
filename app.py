"""
BioMultiAgent Flask Web App
============================
Run: python app.py
Open: http://localhost:5001
"""

from __future__ import annotations

import logging
import os
import uuid

from flask import Flask, jsonify, render_template, request, session

from bioagent.supervisor import run_bio_agent

logging.basicConfig(level=logging.INFO, format="%(asctime)s  %(message)s")
log = logging.getLogger("biomultiagent")

app = Flask(__name__)
app.secret_key = os.environ.get("SECRET_KEY", os.urandom(24).hex())


@app.route("/")
def home():
    return render_template("index.html")


@app.route("/analyze", methods=["POST"])
def analyze():
    data  = request.get_json(force=True)
    query = (data.get("query") or "").strip()
    if not query:
        return jsonify({"error": "query field is required"}), 400

    session_id = session.get("session_id") or str(uuid.uuid4())
    session["session_id"] = session_id
    log.info("Session=%s  Query: %s", session_id[:8], query[:80])

    try:
        state = run_bio_agent(query, session_id=session_id)
    except Exception as exc:
        log.error("Agent error: %s", exc)
        return jsonify({"error": str(exc)}), 500

    return jsonify({
        "result":    state["final_response"],
        "intent":    state["intent"],
        "sub_tasks": state["sub_tasks"],
        "citations": state["citations"],
        "details":   state["agent_results"],
    })


@app.route("/health")
def health():
    return jsonify({"status": "ok"})


if __name__ == "__main__":
    port = int(os.environ.get("PORT", 5001))
    app.run(host="0.0.0.0", port=port, debug=True)

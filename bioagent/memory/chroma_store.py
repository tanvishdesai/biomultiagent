"""ChromaDB session memory for BioMultiAgent."""

from __future__ import annotations

import hashlib
import time
from pathlib import Path
from typing import Any, Dict, List, Optional


class SessionMemory:
    """
    Persists Q&A interactions per session using ChromaDB PersistentClient.
    Falls back gracefully to an in-memory list when ChromaDB is unavailable.
    """

    def __init__(self, collection_name: str = "bio_session_memory") -> None:
        self._fallback: List[Dict[str, Any]] = []
        self._chroma_ok = False
        persist_dir = Path(__file__).resolve().parents[2] / "data" / "chroma"
        try:
            import chromadb
            self._client = chromadb.PersistentClient(path=str(persist_dir))
            self._collection = self._client.get_or_create_collection(collection_name)
            self._chroma_ok = True
        except ImportError:
            pass
        except Exception:
            pass

    def store(self, session_id: str, query: str, result: str) -> None:
        doc_id = hashlib.md5(f"{session_id}_{time.time()}".encode()).hexdigest()
        doc    = f"Q: {query}\nA: {result}"

        if self._chroma_ok:
            try:
                self._collection.add(
                    documents=[doc],
                    ids=[doc_id],
                    metadatas=[{"session_id": session_id, "ts": int(time.time())}],
                )
                return
            except Exception:
                pass

        self._fallback.append({
            "id":         doc_id,
            "session_id": session_id,
            "document":   doc,
            "ts":         int(time.time()),
        })

    def retrieve(self, session_id: str, query: str, n: int = 3) -> List[str]:
        if self._chroma_ok:
            try:
                res = self._collection.query(
                    query_texts=[query],
                    n_results=n,
                    where={"session_id": session_id},
                )
                return res["documents"][0] if res["documents"] else []
            except Exception:
                pass

        session_docs = [
            r["document"]
            for r in self._fallback
            if r["session_id"] == session_id
        ]
        return session_docs[-n:]

    def clear_session(self, session_id: str) -> None:
        if self._chroma_ok:
            try:
                self._collection.delete(where={"session_id": session_id})
            except Exception:
                pass
        self._fallback = [
            r for r in self._fallback if r["session_id"] != session_id
        ]


_memory: Optional[SessionMemory] = None


def get_memory() -> SessionMemory:
    global _memory
    if _memory is None:
        _memory = SessionMemory()
    return _memory

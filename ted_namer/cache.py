"""Lightweight SQLite-backed JSON cache for external biological database API calls."""

from __future__ import annotations

import json
import sqlite3
import time
from pathlib import Path
from typing import Any, Optional


class ApiCache:
    """Thread/process-safe SQLite key-value cache for API responses."""

    def __init__(
        self,
        db_path: str | Path = ".ted_namer_cache.sqlite",
        ttl_seconds: int = 7 * 24 * 3600,
        enabled: bool = True,
    ) -> None:
        self.db_path = Path(db_path)
        self.ttl_seconds = ttl_seconds
        self.enabled = enabled
        if self.enabled:
            self._init_db()

    def _init_db(self) -> None:
        with sqlite3.connect(self.db_path) as conn:
            conn.execute(
                """
                CREATE TABLE IF NOT EXISTS api_cache (
                    namespace TEXT NOT NULL,
                    cache_key TEXT NOT NULL,
                    payload TEXT NOT NULL,
                    created_at REAL NOT NULL,
                    PRIMARY KEY (namespace, cache_key)
                )
                """
            )
            conn.commit()

    def get(self, namespace: str, key: str) -> Optional[Any]:
        if not self.enabled:
            return None
        try:
            with sqlite3.connect(self.db_path) as conn:
                row = conn.execute(
                    "SELECT payload, created_at FROM api_cache WHERE namespace = ? AND cache_key = ?",
                    (namespace, key),
                ).fetchone()
                if not row:
                    return None
                payload_str, created_at = row
                if self.ttl_seconds > 0 and (time.time() - created_at) > self.ttl_seconds:
                    return None
                return json.loads(payload_str)
        except Exception:
            return None

    def set(self, namespace: str, key: str, value: Any) -> None:
        if not self.enabled:
            return
        try:
            payload_str = json.dumps(value)
            with sqlite3.connect(self.db_path) as conn:
                conn.execute(
                    """
                    INSERT OR REPLACE INTO api_cache (namespace, cache_key, payload, created_at)
                    VALUES (?, ?, ?, ?)
                    """,
                    (namespace, key, payload_str, time.time()),
                )
                conn.commit()
        except Exception:
            pass

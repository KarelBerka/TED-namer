"""Lightweight SQLite-backed JSON cache for external biological database API calls."""

from __future__ import annotations

import json
import sqlite3
import time
from pathlib import Path
from typing import Any, Optional

_HOT_MEM_NAMESPACES = {"cath_h_ipr", "cath_tree_node", "cath_resolution_v3"}


class ApiCache:
    """Fast WAL-backed SQLite key-value cache with in-memory hot-key acceleration."""

    def __init__(
        self,
        db_path: str | Path = ".ted_namer_cache.sqlite",
        ttl_seconds: int = 7 * 24 * 3600,
        enabled: bool = True,
    ) -> None:
        self.db_path = Path(db_path)
        self.ttl_seconds = ttl_seconds
        self.enabled = enabled
        self._conn: Optional[sqlite3.Connection] = None
        self._mem: dict[tuple[str, str], Any] = {}
        self._pending_writes = 0
        if self.enabled:
            self._init_db()

    def _get_conn(self) -> sqlite3.Connection:
        if self._conn is None:
            self._conn = sqlite3.connect(self.db_path, check_same_thread=False, timeout=30.0)
            self._conn.execute("PRAGMA journal_mode=WAL")
            self._conn.execute("PRAGMA synchronous=NORMAL")
            self._conn.execute("PRAGMA temp_store=MEMORY")
        return self._conn

    def _init_db(self) -> None:
        conn = self._get_conn()
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
        mk = (namespace, key)
        if mk in self._mem:
            return self._mem[mk]
        try:
            conn = self._get_conn()
            row = conn.execute(
                "SELECT payload, created_at FROM api_cache WHERE namespace = ? AND cache_key = ?",
                (namespace, key),
            ).fetchone()
            if not row:
                return None
            payload_str, created_at = row
            if self.ttl_seconds > 0 and (time.time() - created_at) > self.ttl_seconds:
                return None
            val = json.loads(payload_str)
            if namespace in _HOT_MEM_NAMESPACES:
                self._mem[mk] = val
            return val
        except Exception:
            return None

    def set(self, namespace: str, key: str, value: Any) -> None:
        if not self.enabled:
            return
        if namespace in _HOT_MEM_NAMESPACES:
            self._mem[(namespace, key)] = value
        try:
            payload_str = json.dumps(value)
            conn = self._get_conn()
            conn.execute(
                """
                INSERT OR REPLACE INTO api_cache (namespace, cache_key, payload, created_at)
                VALUES (?, ?, ?, ?)
                """,
                (namespace, key, payload_str, time.time()),
            )
            self._pending_writes += 1
            if self._pending_writes >= 25:
                conn.commit()
                self._pending_writes = 0
        except Exception:
            pass

    def flush(self) -> None:
        if self._conn is not None and self._pending_writes > 0:
            try:
                self._conn.commit()
                self._pending_writes = 0
            except Exception:
                pass

"""Markdown と Notebook に書かれた相対リンクの参照先が存在することを確かめる。

ディレクトリ再編でリンクが切れるのを防ぐための検査。外部URLと、
ページ内アンカー（``#...``）の見出しの有無は検査しない。

Run from the repository root: python -m unittest discover -s tests -v
"""

import json
import re
from pathlib import Path
import unittest

ROOT = Path(__file__).resolve().parents[1]
LINK = re.compile(r"\]\(([^)\s]+)\)")
SKIP_DIRS = {".git", ".venv", "node_modules"}


def documents():
    for path in sorted(ROOT.rglob("*")):
        if path.suffix in {".md", ".ipynb"} and not SKIP_DIRS & set(path.relative_to(ROOT).parts):
            yield path


def markdown_text(path):
    """Markdownはファイル全体、Notebookはマークダウンセルだけを返す。"""
    text = path.read_text(encoding="utf-8")
    if path.suffix != ".ipynb":
        return text
    cells = json.loads(text)["cells"]
    return "\n".join("".join(c["source"]) for c in cells if c["cell_type"] == "markdown")


def broken_links(path):
    broken = []
    for target in LINK.findall(markdown_text(path)):
        if re.match(r"[a-z][a-z0-9+.-]*:", target, re.IGNORECASE) or target.startswith("#"):
            continue  # 外部URL・mailto・ページ内アンカー
        if not (path.parent / target.split("#", 1)[0]).exists():
            broken.append(target)
    return broken


class DocumentLinkTests(unittest.TestCase):
    def test_relative_links_resolve(self):
        paths = list(documents())
        self.assertTrue(paths)
        for path in paths:
            with self.subTest(path=path.relative_to(ROOT).as_posix()):
                self.assertEqual(broken_links(path), [], "リンク先が存在しない")


if __name__ == "__main__":
    unittest.main()

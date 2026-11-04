"""Convert percent-format sources (notebook_sources/*.py) into Jupyter notebooks."""
from pathlib import Path
import nbformat as nbf

OUT = Path(__file__).resolve().parents[1]
SRC, NBD = OUT / "notebook_sources", OUT / "notebooks"
NBD.mkdir(exist_ok=True)
BOOT = '''import sys
from pathlib import Path
for _p in [Path.cwd(), *Path.cwd().parents]:
    if (_p / "scripts" / "amr_common.py").exists():
        sys.path.insert(0, str(_p / "scripts")); break
    if (_p / "Analysis_Results" / "scripts" / "amr_common.py").exists():
        sys.path.insert(0, str(_p / "Analysis_Results" / "scripts")); break
from amr_common import *'''

def convert(src: Path) -> Path:
    cells, kind, buf = [], None, []
    def flush():
        if kind is None: return
        text = "\n".join(buf).strip("\n")
        if not text: return
        if kind == "md":
            lines = [l[2:] if l.startswith("# ") else ("" if l.strip() == "#" else l) for l in text.split("\n")]
            cells.append(nbf.v4.new_markdown_cell("\n".join(lines)))
        else:
            cells.append(nbf.v4.new_code_cell(text.replace("#BOOTSTRAP", BOOT)))
    for line in src.read_text(encoding="utf-8").split("\n"):
        if line.startswith("# %%"):
            flush(); buf = []
            kind = "md" if "[markdown]" in line else "code"
            continue
        buf.append(line)
    flush()
    nb = nbf.v4.new_notebook(cells=cells)
    nb.metadata["kernelspec"] = {"name": "python3", "display_name": "Python 3 (ipykernel)", "language": "python"}
    nb.metadata["language_info"] = {"name": "python"}
    dst = NBD / (src.stem + ".ipynb")
    nbf.write(nb, dst)
    return dst

if __name__ == "__main__":
    import sys
    NBD.mkdir(exist_ok=True)
    targets = sys.argv[1:] or [p.stem for p in sorted(SRC.glob("*.py"))]
    for t in targets:
        print("built", convert(SRC / f"{t}.py").name)

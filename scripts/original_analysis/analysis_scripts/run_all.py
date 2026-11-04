"""Build and execute every analysis notebook from a fresh kernel, in order.
Usage: python run_all.py [notebook_stem ...]"""
import sys, time, json, re, datetime as dt, traceback
from pathlib import Path
import nbformat as nbf
from nbclient import NotebookClient
from nbconvert import ScriptExporter
sys.path.insert(0, str(Path(__file__).resolve().parent))
import build_notebooks as B

OUT = Path(__file__).resolve().parents[1]
EXP = OUT / "scripts" / "nb_exports"; EXP.mkdir(exist_ok=True)
LOGF = OUT / "logs" / "notebook_run_log.jsonl"

def run(stem):
    nb_path = B.convert(B.SRC / f"{stem}.py")
    nb = nbf.read(nb_path, as_version=4)
    t0 = time.time(); status, err = "success", ""
    try:
        NotebookClient(nb, timeout=3600, kernel_name="python3",
                       resources={"metadata": {"path": str(nb_path.parent)}}).execute()
    except Exception as e:
        status, err = "error", re.sub(r"\x1b\[[0-9;]*m", "", "".join(traceback.format_exception_only(type(e), e)))[-2500:]
    nbf.write(nb, nb_path)
    body, _ = ScriptExporter().from_notebook_node(nb)
    (EXP / f"{stem}.py").write_text(body, encoding="utf-8")
    rec = dict(notebook=stem, status=status, seconds=round(time.time() - t0, 1),
               finished=dt.datetime.now().isoformat(timespec="seconds"), error=err)
    with open(LOGF, "a", encoding="utf-8") as f: f.write(json.dumps(rec) + "\n")
    print(json.dumps(rec)[:3500]); sys.stdout.flush()
    return status

if __name__ == "__main__":
    stems = sys.argv[1:] or [p.stem for p in sorted(B.SRC.glob("*.py"))]
    for s in stems:
        if run(s) != "success":
            sys.exit(1)

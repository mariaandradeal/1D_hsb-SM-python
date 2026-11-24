"""Compile the technical manual with pdfLaTeX/BibTeX, without latexmk/Perl."""
from pathlib import Path
import shutil
import subprocess

directory = Path(__file__).resolve().parents[1]/"docs"/"manual"
for program in ("pdflatex", "bibtex"):
    if shutil.which(program) is None:
        raise SystemExit(f"Install a TeX distribution providing {program}, then retry.")
commands = [["pdflatex", "-interaction=nonstopmode", "-halt-on-error", "main.tex"],
            ["bibtex", "main"],
            ["pdflatex", "-interaction=nonstopmode", "-halt-on-error", "main.tex"],
            ["pdflatex", "-interaction=nonstopmode", "-halt-on-error", "main.tex"]]
for command in commands:
    completed = subprocess.run(command, cwd=directory, capture_output=True, text=True, errors="replace")
    if completed.returncode:
        print(completed.stdout[-6000:])
        print(completed.stderr)
        raise SystemExit(completed.returncode)
print(f"Built {directory/'main.pdf'}")

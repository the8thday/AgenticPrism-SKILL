#!/usr/bin/env python3
"""Install the AgenticPrism runtime next to its Skills.

Creates <repo>/.venv, installs this package (editable, so `git pull` keeps Skills
and code on the same version) with dependency versions constrained by
requirements-lock.txt, then runs a self-check. Uses uv when available (it can
fetch the validated Python 3.13 itself); otherwise the interpreter running this
script, which must be Python 3.12 or newer. Standard library only.

    python3 install.py                 # install or update the runtime
    python3 install.py --dev           # also install test/validation tools
    python3 install.py --link-skills ~/.claude/skills   # symlink each Skill folder there
"""
import argparse
import json
import os
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parent
VENV = ROOT / ".venv"
WINDOWS = os.name == "nt"
BIN = VENV / ("Scripts" if WINDOWS else "bin")
EXE = BIN / ("agentic-prism.exe" if WINDOWS else "agentic-prism")
PY = BIN / ("python.exe" if WINDOWS else "python")
VALIDATED_PYTHON = "3.13"
MIN_PYTHON = (3, 12)


def say(message):
    print(f"[agentic-prism] {message}", flush=True)


def run(args, dry):
    say("$ " + " ".join(str(a) for a in args))
    if not dry:
        subprocess.run([str(a) for a in args], check=True, cwd=ROOT)


def version():
    for line in (ROOT / "pyproject.toml").read_text().splitlines():
        if line.startswith("version"):
            return line.split("=", 1)[1].strip().strip('"')
    raise SystemExit("pyproject.toml has no version")


def create_environment(args):
    uv = None if args.no_uv else shutil.which("uv")
    lock = ROOT / "requirements-lock.txt"
    if VENV.exists() and args.force:
        say(f"removing existing {VENV}")
        if not args.dry_run:
            shutil.rmtree(VENV)
    if uv:
        if not VENV.exists():
            run([uv, "venv", VENV, "--python", args.python or VALIDATED_PYTHON], args.dry_run)
        run([uv, "pip", "install", "--python", PY, "-c", lock, "-e", ROOT], args.dry_run)
        if args.dev:
            run([uv, "pip", "install", "--python", PY, "-r", lock], args.dry_run)
        return
    if args.python:
        raise SystemExit("--python needs uv; without uv, run install.py with the interpreter you want")
    if sys.version_info[:2] < MIN_PYTHON:
        raise SystemExit(f"Python {MIN_PYTHON[0]}.{MIN_PYTHON[1]}+ is required (pinned numpy/scipy); "
                         f"this is {sys.version.split()[0]}. Install uv or a newer Python.")
    if f"{sys.version_info[0]}.{sys.version_info[1]}" != VALIDATED_PYTHON:
        say(f"note: releases are validated on Python {VALIDATED_PYTHON}; using {sys.version.split()[0]}")
    if not VENV.exists():
        run([sys.executable, "-m", "venv", VENV], args.dry_run)
    run([PY, "-m", "pip", "install", "--upgrade", "pip"], args.dry_run)
    run([PY, "-m", "pip", "install", "-c", lock, "-e", ROOT], args.dry_run)
    if args.dev:
        run([PY, "-m", "pip", "install", "-r", lock], args.dry_run)


def self_check(dry):
    say("self-check")
    if dry:
        return {"status": "skipped (dry run)"}
    doctor = subprocess.run([str(EXE), "doctor", "--collection", str(ROOT)], capture_output=True, text=True)
    print(doctor.stdout, end="")
    if doctor.returncode:
        print(doctor.stderr, end="", file=sys.stderr)
        raise SystemExit("doctor reported a problem; see above")
    with tempfile.TemporaryDirectory(prefix="agentic-prism-selfcheck-") as tmp:
        out = Path(tmp) / "run"
        env = {**os.environ, "MPLBACKEND": "Agg"}
        subprocess.run([str(EXE), "analyze", "--config", str(ROOT / "fixtures/groups_synthetic/paired_config.json"),
                        "--output", str(out)], check=True, capture_output=True, env=env)
        subprocess.run([str(EXE), "verify", "--run", str(out)], check=True, capture_output=True)
        if not (out / "report.html").stat().st_size > 1000:
            raise SystemExit("self-check report was not written")
    say("self-check passed: analyze, render and verify work")
    return json.loads(doctor.stdout)


def link_skills(target, dry):
    target = Path(target).expanduser()
    say(f"linking Skill folders into {target}")
    if not dry:
        target.mkdir(parents=True, exist_ok=True)
    for skill in sorted(p for p in (ROOT / "skills").iterdir() if (p / "SKILL.md").exists()):
        link = target / skill.name
        if link.is_symlink() and link.resolve() == skill.resolve():
            say(f"  {link} already points here")
            continue
        if link.exists() or link.is_symlink():
            say(f"  skipped {link}: something else already exists there (not changed)")
            continue
        say(f"  {link} -> {skill}")
        if not dry:
            try:
                link.symlink_to(skill, target_is_directory=True)
            except OSError as exc:
                raise SystemExit(f"Could not create a symlink ({exc}). Do not copy the folders instead: a copied "
                                 f"Skill cannot find its runtime. Point your agent at {ROOT / 'skills'} directly.")


def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--dev", action="store_true", help="also install pytest and validation-only tools from the lock file")
    parser.add_argument("--force", action="store_true", help="delete and recreate .venv")
    parser.add_argument("--python", help=f"Python version for uv to use (default {VALIDATED_PYTHON})")
    parser.add_argument("--no-uv", action="store_true", help="use venv + pip even if uv is installed")
    parser.add_argument("--link-skills", metavar="DIR", help="symlink each Skill folder into DIR (for example ~/.claude/skills)")
    parser.add_argument("--dry-run", action="store_true", help="print the steps without changing anything")
    parser.add_argument("--skip-self-check", action="store_true", help=argparse.SUPPRESS)
    args = parser.parse_args()
    say(f"AgenticPrism {version()} at {ROOT}")
    create_environment(args)
    report = {} if args.skip_self_check else self_check(args.dry_run)
    if args.link_skills:
        link_skills(args.link_skills, args.dry_run)
    say(f"runtime: {EXE}")
    say(f"Skills:  {ROOT / 'skills'} (start with {ROOT / 'skills/agentic-prism/SKILL.md'})")
    if report.get("status") == "warning":
        say("installed with warnings; see the doctor output above")


if __name__ == "__main__":
    main()

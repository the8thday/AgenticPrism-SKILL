"""Runtime self-report: version, interpreter, dependency pins and collection match."""
from pathlib import Path
import importlib.metadata
import platform
import re
import sys
from . import __version__

VALIDATED_PYTHON = "3.13"


def _lock(root):
    path = Path(root) / "requirements-lock.txt"
    if not path.exists():
        return None
    pins = {}
    for line in path.read_text().splitlines():
        match = re.match(r"^([A-Za-z0-9_.-]+)==([^\s;]+)", line.strip())
        if match:
            pins[match.group(1).lower().replace("_", "-")] = match.group(2)
    return pins


def doctor(collection=None):
    package = Path(__file__).resolve().parent
    source_root = package.parents[1] if package.parent.name == "src" else None
    root = Path(collection).resolve() if collection else source_root
    report = {"package_version": __version__, "package_location": str(package), "python": platform.python_version(),
              "python_executable": sys.executable, "platform": platform.platform(), "dependencies": {},
              "notes": [], "problems": [], "status": "ok"}
    for name in ("numpy", "scipy", "pandas", "matplotlib"):
        report["dependencies"][name] = importlib.metadata.version(name)
    if sys.version_info[:2] < (3, 12):
        report["problems"].append("Python 3.12+ is required by the pinned numerical dependencies")
    if platform.python_version().rsplit(".", 1)[0] != VALIDATED_PYTHON:
        report["notes"].append(f"releases are validated on Python {VALIDATED_PYTHON}")
    pins = _lock(root) if root else None
    if pins is None:
        report["notes"].append("requirements-lock.txt not found; dependency pins not checked")
    else:
        mismatched = {}
        for dist in importlib.metadata.distributions():
            name = dist.metadata["Name"].lower().replace("_", "-")
            if name in pins and dist.version != pins[name]:
                mismatched[name] = {"installed": dist.version, "locked": pins[name]}
        report["dependency_pin_mismatches"] = mismatched
        if mismatched:
            report["notes"].append("installed dependency versions differ from requirements-lock.txt; "
                                   "rerun install.py to restore the validated versions")
    if collection:
        pyproject = root / "pyproject.toml"
        version = None
        if pyproject.exists():
            found = re.search(r'^version\s*=\s*"([^"]+)"', pyproject.read_text(), re.M)
            version = found.group(1) if found else None
        report["collection_root"] = str(root)
        report["collection_version"] = version
        if not (root / "skills").is_dir():
            report["problems"].append(f"{root} does not contain a skills/ folder")
        if version != __version__:
            report["problems"].append(f"runtime {__version__} does not match collection {version}; "
                                      "run install.py in the collection root")
        if source_root is not None and source_root != root:
            report["notes"].append("this runtime is installed from a different source tree than the collection")
    if report["problems"]:
        report["status"] = "error"
    elif report.get("dependency_pin_mismatches"):
        report["status"] = "warning"
    return report

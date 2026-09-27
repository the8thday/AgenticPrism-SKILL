import argparse
import json
import sys
from pathlib import Path
from . import __version__
from .workflow import analyze, verify_run


def main():
    p = argparse.ArgumentParser(description="AgenticPrism: validated module entry points")
    p.add_argument("--version", action="version", version=__version__)
    sub = p.add_subparsers(dest="command", required=True)
    a = sub.add_parser("analyze", help="Run a supported analysis (see the module registry for analysis_type values)")
    a.add_argument("--config", required=True)
    a.add_argument("--output", required=True)
    a.add_argument("--no-render", action="store_true")
    r = sub.add_parser("render", help="Render saved results, without fitting")
    r.add_argument("--run", required=True)
    r.add_argument("--style", choices=["standard", "prism_like"])
    v = sub.add_parser("verify", help="Verify analysis artifact hashes")
    v.add_argument("--run", required=True)
    imp = sub.add_parser("import-octet", help="Import verified Octet Results.txt layout with explicit metadata")
    imp.add_argument("--manifest", required=True)
    imp.add_argument("--output", required=True)
    plate = sub.add_parser("import-plate", help="Convert plate-reader grid CSVs and plate-map layers to a long table")
    plate.add_argument("--manifest", required=True)
    plate.add_argument("--output", required=True)
    doc = sub.add_parser("doctor", help="Report runtime version, Python, dependency pins and collection match")
    doc.add_argument("--collection", help="Skill collection root to compare against (the folder holding skills/)")
    args = p.parse_args()
    try:
        if args.command == "analyze":
            out = analyze(args.config, args.output, not args.no_render)
            result = json.loads((out / "results.json").read_text())
            if result.get("analysis_type") == "variance_components":
                print(json.dumps({"output": str(out), "status": result["status"],
                                  "interpretation": "Precision estimates do not establish SOP acceptance"}))
                return 0
            if result.get("analysis_type") in ("method_validation", "ada_sensitivity", "ada_drug_tolerance"):
                facts = json.loads((out / "interpretation_facts.json").read_text())
                print(json.dumps({"output": str(out), "analysis_type": result["analysis_type"], "primary": facts["primary"],
                                  "must_mention": facts["must_mention"],
                                  "interpretation": "Meeting declared criteria is part of a validation, not regulatory acceptance"}, ensure_ascii=False))
                return 0
            failed = sum(f["status"] == "failed" for f in result["fits"])
            counts = ({"n_fit_groups": len(result["fits"]),
                       "n_sensorgrams": sum(f["n_curves"] for f in result["fits"])}
                      if result.get("analysis_type") == "binding_kinetics" else
                      {"n_plates": len(result["fits"]), "n_unknown_wells": len(result["wells"])}
                      if result.get("analysis_type") == "elisa_quantification" else
                      {"n_comparisons": 1} if result.get("analysis_type") == "group_comparison" else
                      {"n_groups": result["fits"][0]["n_groups"], "n_contrasts": len(result["contrasts"])}
                      if result.get("analysis_type") in ("multi_group_comparison", "nonparametric", "mmrm", "repeated_measures", "time_to_event", "tumor_growth") else
                      {"n_curves": len(result["fits"])})
            unreportable = sum(not f.get("reportable", f["status"] != "failed") for f in result["fits"])
            withheld = sum(w["status"] != "quantified" for w in result.get("wells", []))
            print(json.dumps({"run": str(out), **counts, "n_failed": failed,
                              "n_unreportable_fits": unreportable, "n_withheld_wells": withheld,
                              "interpretation": "Read reportable/status and diagnostics; exit zero means execution completed, not scientific validation"}, ensure_ascii=False))
            return 3 if failed else 0
        if args.command == "render":
            cfg = json.loads((Path(args.run)/"config.resolved.json").read_text())
            if cfg["analysis_type"] == "variance_components":
                from .precision_workflow import render_precision
                print(render_precision(Path(args.run), args.style))
            elif cfg["analysis_type"] == "ada_cut_point":
                from .ada_workflow import render_ada
                print(render_ada(Path(args.run), args.style))
            elif cfg["analysis_type"] == "method_validation":
                from .method_validation_workflow import render_method_validation
                print(render_method_validation(Path(args.run), args.style))
            elif cfg["analysis_type"] in ("ada_sensitivity", "ada_drug_tolerance"):
                from .ada_performance import render_ada_performance
                print(render_ada_performance(Path(args.run), args.style))
            elif cfg["analysis_type"] == "binding_kinetics":
                from .kinetics_report import render_kinetics
                print(render_kinetics(Path(args.run), args.style))
            elif cfg["analysis_type"] == "dose_response_4pl":
                from .dose_report import render_dose
                print(render_dose(Path(args.run), args.style))
            elif cfg["analysis_type"] in ("elisa_quantification", "group_comparison", "multi_group_comparison", "nonparametric", "mmrm", "repeated_measures", "time_to_event", "tumor_growth"):
                from .simple_report import render_simple
                print(render_simple(Path(args.run), args.style))
            else:
                from .report import render_report
                print(render_report(Path(args.run), args.style))
        elif args.command == "doctor":
            from .doctor import doctor
            report = doctor(args.collection)
            print(json.dumps(report, ensure_ascii=False, indent=2))
            return 1 if report["status"] == "error" else 0
        elif args.command == "import-plate":
            from .plate_import import import_plates
            print(import_plates(args.manifest, args.output))
        elif args.command == "import-octet":
            from .octet_import import import_octet
            print(import_octet(args.manifest,args.output))
        else:
            verify_run(args.run)
            print("Scientific artifact hashes verified")
        return 0
    except Exception as exc:
        print(f"{type(exc).__name__}: {exc}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    sys.exit(main())

"""Command-line entry points used verbatim in the README and manual."""
from pathlib import Path
import argparse
import sys
from .config import load_config
from .model import run
from .io import read_observations, write_json
from .evaluation import evaluate, signatures


def main(argv=None):
    parser = argparse.ArgumentParser(prog="hsbsm", description="Run the hsB-SM hydrological model")
    parser.add_argument("--version", action="version", version="hsbsm 0.1.0")
    sub = parser.add_subparsers(dest="command", required=True)
    for name in ("run", "calibrate", "sensitivity"):
        p = sub.add_parser(name)
        p.add_argument("config", type=Path)
        p.add_argument("--output", type=Path, default=Path("outputs")/name)
        p.add_argument("--steps", type=int, help="Use only the first N forcing rows")
        p.add_argument("--no-plots", action="store_true")
        if name in ("run", "calibrate"):
            p.add_argument("--observations", type=Path, required=name == "calibrate")
            p.add_argument("--warmup", type=int, default=0, help="Steps excluded from metrics")
        if name in ("calibrate", "sensitivity"):
            p.add_argument("--active", nargs="+", default=["Aq_k", "Aq_f"])
        if name == "calibrate":
            p.add_argument("--maxiter", type=int, default=10)
            p.add_argument("--popsize", type=int, default=5)
            p.add_argument("--seed", type=int, default=42)
            p.add_argument("--split", type=float, default=0.8)
    args = parser.parse_args(argv)
    try:
        cfg = load_config(args.config)
        if args.steps is not None:
            if args.steps <= 0:
                raise ValueError("--steps must be positive")
            cfg.forcing = cfg.forcing.iloc[:args.steps].copy()
        if args.command == "sensitivity":
            from .calibration import sensitivity
            table = sensitivity(cfg, args.active)
            args.output.mkdir(parents=True, exist_ok=True)
            table.to_csv(args.output/"sensitivity.csv", index=False)
            print(f"Sensitivity saved to {args.output.resolve()}")
            return 0
        observations = read_observations(args.observations) if args.observations else None
        if args.command == "calibrate":
            from .calibration import calibrate
            result, report = calibrate(cfg, observations, active=args.active, split=args.split,
                warmup_steps=args.warmup, maxiter=args.maxiter, popsize=args.popsize, seed=args.seed,
                progress=lambda n, f: print(f"Evaluation {n}: 1-KGE = {f:.5g}", flush=True))
        else:
            result = run(cfg, progress=lambda n, total: print(f"Simulation {n}/{total}", flush=True))
        result.save(args.output, plots=not args.no_plots)
        if args.command == "calibrate":
            write_json(args.output/"calibration.json", report)
        if observations is not None:
            summary, annual, comparison = evaluate(result, observations, args.warmup)
            write_json(args.output/"evaluation.json", summary)
            annual.to_csv(args.output/"evaluation_by_year.csv", index=False)
            comparison.to_csv(args.output/"comparison.csv", index=False)
            if not args.no_plots:
                from .plotting import plot_result
                plot_result(result, args.output, comparison)
        write_json(args.output/"signatures.json", signatures(result))
        print(f"Outputs: {args.output.resolve()}")
        print(f"Relative water-balance error: {result.diagnostics['relative_water_balance_error']:.3e}")
        return 0
    except (ValueError, ArithmeticError, OSError, KeyError) as exc:
        print(f"hsbsm: {exc}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())

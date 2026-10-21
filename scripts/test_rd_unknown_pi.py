#!/usr/bin/env python3
"""Compare production RD with the previous RC/normal-RD simulated fits."""
import argparse
import math
from pathlib import Path
import statistics as st
import sys

sys.dont_write_bytecode = True
from test_rd_rc_normal import trace_values, ROOT, run, table


def stats(chains):
    batches, halves = [], []
    for chain in chains:
        size, half = len(chain) // 20, len(chain) // 2
        batches.extend(st.mean(chain[j*size:(j+1)*size]) for j in range(20))
        halves.extend([chain[:half], chain[half:2*half]])
    within = st.mean(st.variance(x) for x in halves)
    between = half * st.variance(st.mean(x) for x in halves)
    rhat = math.sqrt(((half-1)/half*within + between/half)/within) if within else 1.
    return st.mean(batches), math.sqrt(st.variance(batches)/len(batches)), rhat


def summarize(folder, baseline, replicates):
    report = ["# Production RD with learned annotation mixture weights", "",
              "Same simulated data and seeds as the normal-prior sanity check. "
              "Production RD uses fixed annotation variances 0, 0.0001, 0.001, 0.01, 0.1 "
              "with Dirichlet(1,1,1,1,1) weights; annotation pi is 1 minus the zero weight. "
              "Its prior differs from the learned single-normal prior as well as allowing exclusion.", "",
              "| Quantity | RC mean | Normal RD mean | Default RD mean | Default RD MCSE | Default RD split R-hat |",
              "|---|---:|---:|---:|---:|---:|"]
    results, samples = {}, {}
    for model in ["RC", "NormalRD", "DefaultRD"]:
        src, stem = (baseline, "RD") if model == "NormalRD" else (folder, model)
        if model == "RC" and not (folder / "RC_0.parRes").exists():
            src = baseline
        core = [table(src / f"{stem}_{r}.mcmcsamples/CoreParameters.mcmcsamples.txt")
                for r in range(replicates)]
        annotations = [trace_values(src / f"{stem}_{r}.mcmcsamples/AnnoEffects_p2.mcmcsamples.txt")
                       for r in range(replicates)]
        results[model] = {name: stats([[float(row[name]) for row in chain] for chain in core])
                          for name in ["hsq", "NnzSnp", "SigmaSq"]}
        for col, name in enumerate(["Annotation intercept", "Annotation Enriched"]):
            results[model][name] = stats([[row[col] for row in chain] for chain in annotations])
        samples[model] = [table(src / f"{stem}_{r}.snpRes") for r in range(replicates)]
        if model == "DefaultRD":
            for name in core[0][0]:
                if name.startswith("AnnoPi"):
                    results[model][name] = stats([[float(row[name]) for row in chain] for chain in core])
    for name in results["RC"]:
        rd = results["DefaultRD"][name]
        report.append(f"| {name} | {results['RC'][name][0]:.6g} | "
                      f"{results['NormalRD'][name][0]:.6g} | {rd[0]:.6g} | {rd[1]:.3g} | {rd[2]:.4f} |")
    report.extend(["", "Learned annotation inclusion masses (not annotation PIPs):"])
    for name, values in results["DefaultRD"].items():
        if name.startswith("AnnoPi"):
            report.append(f"- {name}: {values[0]:.6g}; MCSE {values[1]:.3g}; split R-hat {values[2]:.4f}.")
    # Annotation coefficient assignment PIP and slab weights, averaged over chains.
    annotation_results = [table(folder / f"DefaultRD_{r}.parSetRes") for r in range(replicates)]
    report.extend(["", "Production RD annotation outputs:"])
    keys = [(row["Parameter"], row["Annotation"]) for row in annotation_results[0]
            if row["Parameter"].startswith(("AnnoPIP_", "AnnoMix_"))
            and row["Annotation"] != "Intercept"]
    for label, annotation in keys:
        mean = st.mean(float(next(row["Mean"] for row in chain
                                  if row["Parameter"] == label and row["Annotation"] == annotation))
                       for chain in annotation_results)
        report.append(f"- {label} {annotation}: {mean:.6g}.")
    for column in ["A1Effect", "PIP"]:
        pooled = {model: [st.mean(float(chain[j][column]) for chain in samples[model])
                          for j in range(len(samples[model][0]))] for model in samples}
        for other in ["RC", "NormalRD"]:
            a, b = pooled[other], pooled["DefaultRD"]
            invalid = sum(not math.isfinite(x) for x in a+b)
            if invalid:
                report.append(f"\n{column}, {other} versus default RD: {invalid} nonfinite pooled values.")
                continue
            ma, mb = st.mean(a), st.mean(b)
            corr = sum((x-ma)*(y-mb) for x,y in zip(a,b)) / math.sqrt(
                sum((x-ma)**2 for x in a)*sum((y-mb)**2 for y in b))
            report.append(f"\n{column}, {other} versus default RD: correlation {corr:.8f}, "
                          f"maximum absolute difference {max(abs(x-y) for x,y in zip(a,b)):.6g}, "
                          f"mean absolute difference {st.mean(abs(x-y) for x,y in zip(a,b)):.6g}.")
    text = "\n".join(report) + "\n"
    (folder / "comparison.md").write_text(text)
    print(text, flush=True)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("baseline", type=Path)
    parser.add_argument("output", type=Path)
    parser.add_argument("--iterations", type=int, default=30000)
    parser.add_argument("--replicates", type=int, default=4)
    parser.add_argument("--two-components", action="store_true")
    parser.add_argument("--run-rc", action="store_true")
    args = parser.parse_args()
    baseline, folder = args.baseline.resolve(), args.output.resolve()
    folder.mkdir(parents=True, exist_ok=True)
    for replicate in range(args.replicates):
        for model in (["RC", "DefaultRD"] if args.run_rc else ["DefaultRD"]):
            bayes = "RD" if model == "DefaultRD" else "RC"
            prefix = folder / f"{model}_{replicate}"
            print(f"Running {model}, replicate {replicate+1}/{args.replicates}", flush=True)
            run([ROOT / "scr/gctb", "--sbayes", bayes,
                 "--gwas-summary", baseline / "sim.ma", "--ldm", baseline / "sim_ld.ldm.sparse",
                 "--annot", baseline / "sim.annot", "--chain-length", args.iterations,
                 "--burn-in", args.iterations // 5, "--thin", "5", "--out-freq", "1000",
                 "--thread", "1", "--num-chains", "1", "--alg", "Gibbs", "--write-mcmc-txt",
                 *(["--sbayes", "RC", "--gamma", "0,1", "--pi", "0.85,0.15",
                    "--sbayes", bayes] if args.two_components else []),
                 "--seed", 123 + replicate * 1009, "--out", prefix], str(prefix) + ".stdout")
    summarize(folder, baseline, args.replicates)


if __name__ == "__main__":
    main()

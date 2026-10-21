#!/usr/bin/env python3
"""Full simulated-GWAS sanity check using an isolated normal-prior RD build.

Run after make -C scr, with EIGEN3_INCLUDE_DIR and BOOST_LIB set.
Production source/defaults are not modified. Outputs persist in the given folder.
"""
import argparse
import math
import os
from pathlib import Path
import random
import statistics
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[1]


def run(args, log):
    with open(log, "w") as out:
        subprocess.run([str(x) for x in args], cwd=ROOT, stdout=out,
                       stderr=subprocess.STDOUT, check=True)
    if "Error:" in Path(log).read_text():
        raise RuntimeError(f"GCTB reported an error; see {log}")


def simulate(folder):
    rng = random.Random(20261010)
    n, m = 1000, 200
    annotation = [int(j % 4 == 0) for j in range(m)]
    genotypes, effects = [], []
    for j in range(m):
        frequency = rng.uniform(.1, .5)
        genotypes.append([int(rng.random() < frequency) +
                          int(rng.random() < frequency) for _ in range(n)])
        effects.append(rng.gauss(0, 1) if rng.random() <
                       (.35 if annotation[j] else .08) else 0.)
    genetic = [sum(effects[j] * genotypes[j][i] for j in range(m)) for i in range(n)]
    mean = sum(genetic) / n
    vg = sum((x - mean)**2 for x in genetic) / (n - 1)
    phenotype = [x + rng.gauss(0, math.sqrt(vg)) for x in genetic]
    ym = sum(phenotype) / n
    phenotype = [x - ym for x in phenotype]
    yy = sum(x*x for x in phenotype)
    prefix = folder / "sim"
    with open(str(prefix) + ".bed", "wb") as bed:
        bed.write(bytes([0x6c, 0x1b, 0x01]))
        for genotype in genotypes:
            codes = [3, 2, 0]  # dosage of BIM A1: 0, 1, 2
            bed.write(bytes(sum(codes[genotype[i+k]] << (2*k) for k in range(4))
                            for i in range(0, n, 4)))
    with open(str(prefix) + ".bim", "w") as bim:
        for j in range(m):
            bim.write(f"1 rs{j+1} 0 {(j+1)*1000} A G\n")
    with open(str(prefix) + ".fam", "w") as fam:
        for i in range(n):
            fam.write(f"F{i+1} I{i+1} 0 0 1 -9\n")
    with open(folder / "sim.annot", "w") as annot, open(folder / "sim.ma", "w") as ma:
        annot.write("SNP Intercept Enriched\n")
        ma.write("SNP A1 A2 freq b se p N\n")
        for j, genotype in enumerate(genotypes):
            xm = sum(genotype) / n
            xx = sum((x-xm)**2 for x in genotype)
            xy = sum((genotype[i]-xm)*phenotype[i] for i in range(n))
            b = xy / xx
            se = math.sqrt((yy - xy*xy/xx) / (n-2) / xx)
            p = math.erfc(abs(b/se) / math.sqrt(2))
            annot.write(f"rs{j+1} 1 {annotation[j]}\n")
            ma.write(f"rs{j+1} A G {xm/2:.9g} {b:.9g} {se:.9g} {p:.9g} {n}\n")
    print(f"Simulated {n} individuals, {m} SNPs, {sum(x != 0 for x in effects)} causal SNPs, target h2=0.5", flush=True)


def normal_build(folder):
    source = (ROOT / "scr/model.cpp").read_text()
    start = source.index("void ApproxBayesRD::sampleUnknowns(const unsigned iter){")
    source = source[:start] + source[start:].replace(
        "void ApproxBayesRD::sampleUnknowns(const unsigned iter){",
        "void ApproxBayesRD::sampleUnknowns(const unsigned iter){\n"
        "    annoEffects.bayesR = false;\n"
        "    piAnno.values.setOnes(); piAnno.value = 1.0f;\n", 1)
    start = source.index("        annoEffects.updateMixtureWeights(sampleAnnoEffectsIndep);", start)
    end = source.index("        addProfileTime(profileTimePiUpdate);", start)
    source = source[:start] + (
        "        sigmaSqAnno.sampleFromFC(annoEffects.ssq,\n"
        "            VectorXf::Constant(annoEffects.numComp, annoEffects.numAnno - 1));\n"
    ) + source[end:]
    (folder / "model_normal.cpp").write_text(source)
    compiler = os.environ.get("CXX", "/opt/homebrew/opt/llvm/bin/clang++" if sys.platform == "darwin" else "g++")
    flags = ["-std=c++11", "-O2", "-fopenmp", "-DUNIX", "-I", ROOT / "scr",
             "-I", os.environ["EIGEN3_INCLUDE_DIR"], "-I", os.environ["BOOST_LIB"]]
    if sys.platform == "darwin":
        flags += ["-stdlib=libc++", "-Dfopen64=fopen"]
    run([compiler, *flags, "-c", folder / "model_normal.cpp", "-o", folder / "model_normal.o"],
        folder / "compile.log")
    objects = sorted(p for p in (ROOT / "scr").glob("*.o") if p.name != "model.o")
    run([compiler, *flags, *objects, folder / "model_normal.o", "-lz", "-o", folder / "gctb_normal"],
        folder / "link.log")


def table(path):
    rows = [line.split() for line in path.read_text().splitlines() if line.strip()]
    return [dict(zip(rows[0], row)) for row in rows[1:]]


def trace_values(path):
    """Read parameter values from either legacy or named MCMC text traces."""
    lines = Path(path).read_text().splitlines()
    named = bool(lines and lines[0].split()[:2] == ['Iteration', 'Chain'])
    return [[float(x) for x in line.split()[2 if named else 0:]]
            for line in lines[1 if named else 0:]]


def summarize(folder, replicates):
    """Batch-mean MCSEs and ordinary split R-hat, plus pooled SNP agreement."""
    def stats(chains):
        batches = []
        halves = []
        for chain in chains:
            size = len(chain) // 20
            batches.extend(statistics.mean(chain[j*size:(j+1)*size]) for j in range(20))
            half = len(chain) // 2
            halves.extend([chain[:half], chain[half:2*half]])
        mean = statistics.mean(batches)
        mcse = math.sqrt(statistics.variance(batches) / len(batches))
        within = statistics.mean(statistics.variance(x) for x in halves)
        between = half * statistics.variance(statistics.mean(x) for x in halves)
        rhat = math.sqrt(((half-1)/half*within + between/half) / within) if within else 1.
        return mean, mcse, rhat

    summaries = {}
    snps = {}
    for model in ["RC", "RD"]:
        core = [table(folder / f"{model}_{r}.mcmcsamples/CoreParameters.mcmcsamples.txt")
                for r in range(replicates)]
        annotations = [trace_values(folder / f"{model}_{r}.mcmcsamples/AnnoEffects_p2.mcmcsamples.txt")
                       for r in range(replicates)]
        summaries[model] = {name: stats([[float(row[name]) for row in chain] for chain in core])
                            for name in ["hsq", "SigmaSq", "NnzSnp"]}
        for col, name in enumerate(["Annotation intercept", "Annotation Enriched"]):
            summaries[model][name] = stats([[row[col] for row in chain] for chain in annotations])
        snps[model] = [table(folder / f"{model}_{r}.snpRes") for r in range(replicates)]

    report = ["# RD normal-prior versus RC sanity check", "",
              "RD uses its retained independent Gaussian-slab sampler with annotation pi=1, "
              "and the same learned annotation-variance update as RC. Production sources are unchanged.",
              "", "Simulation: 1,000 individuals, 200 SNPs, 27 causal variants, target h2=0.5; "
              "one binary enrichment annotation plus a flat-prior intercept. "
              "LD is computed from the same genotypes, with all correlations retained.", "",
              f"{replicates} independent chains per model. MCSE uses 20 contiguous batches per chain; "
              "R-hat is the ordinary split-chain statistic (not rank normalized).", "",
              "| Quantity | RC mean | RD mean | Difference / combined MCSE | RC R-hat | RD R-hat |",
              "|---|---:|---:|---:|---:|---:|"]
    for name in summaries["RC"]:
        rc, rd = summaries["RC"][name], summaries["RD"][name]
        z = abs(rc[0]-rd[0]) / math.hypot(rc[1], rd[1])
        report.append(f"| {name} | {rc[0]:.6g} | {rd[0]:.6g} | {z:.3f} | {rc[2]:.4f} | {rd[2]:.4f} |")
    report.append("")
    for column in ["A1Effect", "PIP"]:
        pooled = {model: [statistics.mean(float(chain[j][column]) for chain in snps[model])
                          for j in range(len(snps[model][0]))] for model in snps}
        a, b = pooled["RC"], pooled["RD"]
        invalid = sum(not math.isfinite(x) for x in a+b)
        if invalid:
            report.append(f"{column}: {invalid} nonfinite pooled values; correlation is unavailable.")
            continue
        ma, mb = statistics.mean(a), statistics.mean(b)
        corr = sum((x-ma)*(y-mb) for x,y in zip(a,b)) / math.sqrt(
            sum((x-ma)**2 for x in a)*sum((y-mb)**2 for y in b))
        report.append(f"Pooled SNP {column}: correlation {corr:.8f}; "
                      f"maximum absolute difference {max(abs(x-y) for x,y in zip(a,b)):.6g}.")
    report.append("\nThe generated RD AnnoMix and AnnoPi output labels belong to its production "
                  "constructor and are stale in this isolated test build; they are excluded from this comparison. "
                  "The inclusion probabilities passed to the annotation sampler are fixed at 1.")
    text = "\n".join(report) + "\n"
    (folder / "comparison.md").write_text(text)
    print(text, flush=True)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("output", type=Path)
    parser.add_argument("--iterations", type=int, default=10000)
    parser.add_argument("--replicates", type=int, default=4)
    parser.add_argument("--reuse-build", action="store_true")
    parser.add_argument("--two-components", action="store_true",
                        help="Use SNP mixture 0,1 to avoid weakly identified intermediate slabs")
    args = parser.parse_args()
    folder = args.output.resolve()
    folder.mkdir(parents=True, exist_ok=True)
    simulate(folder)
    if not args.reuse_build:
        print("Building isolated RD normal-prior binary...", flush=True)
        normal_build(folder)
    elif not (folder / "gctb_normal").exists():
        raise FileNotFoundError("--reuse-build requires gctb_normal in the output directory")
    run([ROOT / "scr/gctb", "--bfile", folder / "sim", "--make-sparse-ldm",
         "--chisq", "0", "--thread", "1", "--out", folder / "sim_ld"], folder / "ld.stdout")
    for replicate in range(args.replicates):
        for model in ["RC", "RD"]:
            prefix = folder / f"{model}_{replicate}"
            print(f"Running {model}, replicate {replicate+1}/{args.replicates}...", flush=True)
            run([folder / "gctb_normal", "--sbayes", model,
                 "--gwas-summary", folder / "sim.ma", "--ldm", folder / "sim_ld.ldm.sparse",
                 "--annot", folder / "sim.annot", "--chain-length", args.iterations,
                 "--burn-in", args.iterations // 5, "--thin", "5", "--out-freq", "1000",
                 "--thread", "1", "--num-chains", "1", "--alg", "Gibbs",
                 "--write-mcmc-txt",
                 # Legacy --pi parsing recognizes RC but omits RD; set the
                 # mixture while RC is selected, then restore the requested model.
                 *(["--sbayes", "RC", "--gamma", "0,1", "--pi", "0.85,0.15",
                    "--sbayes", model] if args.two_components else []),
                 "--seed", 123 + replicate * 1009, "--out", prefix], str(prefix) + ".stdout")
    print(f"Full fit outputs: {folder}", flush=True)
    summarize(folder, args.replicates)


if __name__ == "__main__":
    main()

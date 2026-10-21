#!/usr/bin/env python3
"""Simulate RC or RD annotation effects and an independent Gaussian GWAS."""
import argparse
import json
import math
from pathlib import Path
import random
import statistics as st
import struct
import sys

sys.dont_write_bytecode = True
from test_rd_unknown_pi import ROOT, run, stats, table
from test_rd_rc_normal import trace_values

VARIANCES = [0., .0001, .001, .01, .1]
GAMMA = [0., .001, .01, .1, 1.]


def categorical(rng, weights):
    u, cumulative = rng.random(), 0.
    for g, weight in enumerate(weights):
        cumulative += weight
        if u < cumulative:
            return g
    return len(weights)-1


def simulate(folder, snps, annotations, sample_size, seed, annotation_model="RD",
             rc_variances=None, annotation_reference=None, append_null_annotations=0):
    rng = random.Random(seed)
    names = [f"A{j+1:02d}" for j in range(annotations)]
    weights, groups, effects, intercepts = [], [], [], []
    coverage = [rng.uniform(.15, .4) for _ in names]
    reference_rows=None
    if annotation_reference is not None:
        reference=Path(annotation_reference).read_text().splitlines()
        if reference[0].split() != ['SNP', 'Intercept']+names:
            raise ValueError('Reference annotation header does not match requested annotation count')
        reference_rows=[[int(x) for x in line.split()[2:]] for line in reference[1:]]
        if len(reference_rows) != snps or any(len(row)!=annotations for row in reference_rows):
            raise ValueError('Reference annotation dimensions do not match the simulation')
        coverage=[st.mean(row[j] for row in reference_rows) for j in range(annotations)]
    for c in range(4):
        if annotation_model == "RC":
            groups.append([-1 for _ in names])
            effects.append([rng.gauss(0., math.sqrt(rc_variances[c])) for _ in names])
        else:
            draws = [rng.gammavariate(1., 1.) for _ in VARIANCES]
            weights.append([x/sum(draws) for x in draws])
            groups.append([categorical(rng, weights[c]) for _ in names])
            effects.append([rng.gauss(0., math.sqrt(VARIANCES[g])) if g else 0.
                            for g in groups[c]])
        intercepts.append(st.NormalDist().inv_cdf(.4 if c == 0 else .65)
                          - sum(p*a for p,a in zip(coverage, effects[c])))
    designs, memberships, probabilities = [], [], []
    for i in range(snps):
        row = reference_rows[i] if reference_rows is not None else [int(rng.random() < p) for p in coverage]
        conditional = [st.NormalDist().cdf(intercepts[c] + sum(a*x for a,x in zip(effects[c], row)))
                       for c in range(4)]
        mass, prefix = [], 1.
        for p in conditional:
            mass.append(prefix*(1-p))
            prefix *= p
        mass.append(prefix)
        designs.append(row)
        probabilities.append(mass)
        memberships.append(categorical(rng, mass))
    # Set a fixed SNP base variance using the design-conditional expected Vg.
    # Do not rescale sampled effects: their Gaussian component distributions
    # must remain exactly those specified by the generating model.
    expected_gamma = sum(sum(p*g for p,g in zip(row, GAMMA)) for row in probabilities)
    sigma_sq = .5/expected_gamma
    beta = [rng.gauss(0., math.sqrt(sigma_sq*GAMMA[g])) if g else 0. for g in memberships]
    vg, vare = sum(b*b for b in beta), .5
    varp = vg + vare
    # The sparse GCTB GWAS likelihood uses total phenotype variance / N.
    observed = [b+rng.gauss(0., math.sqrt(varp/sample_size)) for b in beta]
    active_annotations=annotations
    active_weights=weights
    forced_null_names=[f'Null{j+1:02d}' for j in range(append_null_annotations)]
    if append_null_annotations:
        # Append known-negative predictors after generating the GWAS, using
        # an independent RNG, so the original data and causal effects survive.
        null_rng=random.Random(seed ^ 0x4E554C4C)
        null_coverage=[null_rng.uniform(.15,.4) for _ in forced_null_names]
        for row in designs:
            row.extend(int(null_rng.random()<p) for p in null_coverage)
        names.extend(forced_null_names)
        for c in range(4):
            groups[c].extend([0]*append_null_annotations)
            effects[c].extend([0.]*append_null_annotations)
        annotations += append_null_annotations
        if annotation_model=='RD':
            weights=[[(active_annotations*w + (append_null_annotations if g==0 else 0))/annotations
                      for g,w in enumerate(level)] for level in active_weights]
    with open(folder / "sim.annot", "w") as out:
        out.write("SNP Intercept " + " ".join(names) + "\n")
        for i,row in enumerate(designs):
            out.write(f"rs{i+1} 1 " + " ".join(map(str,row)) + "\n")
    with open(folder / "sim.ma", "w") as out:
        out.write("SNP A1 A2 freq b se p N\n")
        for i,b in enumerate(observed):
            # Allele-scale GWAS with 2pq=.5; fit uses standardized SNP effects.
            se = math.sqrt(varp/sample_size/.5)
            raw = b/math.sqrt(.5)
            p = max(math.erfc(abs(raw/se)/math.sqrt(2)), 1e-300)
            out.write(f"rs{i+1} A G .5 {raw:.12g} {se:.12g} {p:.12g} {sample_size}\n")
    with open(folder / "identity.ldm.sparse.info", "w") as info, open(folder / "identity.ldm.sparse.bin", "wb") as binary:
        info.write("Chrom ID GenPos PhysPos A1 A2 A1Freq Index WindStart WindEnd WindSize WindWidth N SamplVar LDsum\n")
        for i in range(snps):
            info.write(f"1 rs{i+1} 0 {(i+1)*1000} A G .5 {i} {i} {i} 1 0 {sample_size} 0 1\n")
            binary.write(struct.pack("=If", i, 1.))
    with open(folder / "annotation_truth.tsv", "w") as out:
        out.write("Annotation Level Component Variance Effect\n")
        for c in range(4):
            out.write(f"Intercept p{c+2} flat NA {intercepts[c]:.12g}\n")
            for j,name in enumerate(names):
                component="normal" if annotation_model=="RC" and groups[c][j]==-1 else groups[c][j]
                variance=rc_variances[c] if component=="normal" else VARIANCES[groups[c][j]]
                out.write(f"{name} p{c+2} {component} {variance} {effects[c][j]:.12g}\n")
    with open(folder / "snp_truth.tsv", "w") as out:
        out.write("SNP Component StandardizedEffect\n")
        for i,(g,b) in enumerate(zip(memberships,beta)):
            out.write(f"rs{i+1} {g} {b:.12g}\n")
    truth = dict(seed=seed, snps=snps, annotations=annotations, sample_size=sample_size,
                 weights=weights, variances=VARIANCES, gamma=GAMMA, sigma_sq=sigma_sq,
                 genetic_variance=vg, residual_variance=vare, hsq=vg/varp,
                 summary_noise_variance=varp/sample_size,
                 intercepts=intercepts, names=names, groups=groups, effects=effects,
                 component_counts=[memberships.count(g) for g in range(5)])
    truth.update(annotation_model=annotation_model, rc_variances=rc_variances if annotation_model=="RC" else None,
                 annotation_reference=str(annotation_reference) if annotation_reference else None,
                 active_annotations=active_annotations, forced_null_names=forced_null_names,
                 active_weights=active_weights if annotation_model=='RD' else None)
    (folder / "truth.json").write_text(json.dumps(truth, indent=2)+"\n")
    print(f"Simulated {snps} SNPs, {annotations} annotations, counts={truth['component_counts']}, h2={truth['hsq']:.4f}", flush=True)
    print("Annotation variances:" if annotation_model=="RC" else "Annotation weights:",
          rc_variances if annotation_model=="RC" else weights, flush=True)
    return truth


def recovery(folder, truth, models, replicates):
    rows = ["Model Level Annotation TrueEffect Mean SD PIP SplitRhat"]
    generating_model=truth.get('annotation_model', 'RD')
    prior_description=("Annotation effects are independent N(0,tau_c^2) draws with one shared variance "
                       "per level, and no exact-zero spike. True variances p2–p5 are "
                       +', '.join(map(str,truth['rc_variances']))
                       +". Variances are fixed generating parameters, not draws from RC's inverse-chi-square hyperprior; "
                       "RC learns its variances using its production hyperprior. "
                       if generating_model=='RC' else
                       "Annotation weights are drawn independently from Dirichlet(1,1,1,1,1), "
                       "then annotation assignments and effects from RD's exact five-component prior. ")
    if truth.get('forced_null_names'):
        prior_description=(f"The original {truth['active_annotations']} annotations follow the "
                           f"{generating_model} generating model; {len(truth['forced_null_names'])} "
                           "additional annotations have effects fixed at zero at all four levels. "
                           "This shared-null stress test introduces cross-level sparsity beyond default RD's "
                           "independent per-level assignments. The original GWAS and causal effects are unchanged. "
                           +prior_description)
    report = [f"# Multiple annotations generated under {generating_model}", "",
              f"{truth['snps']} independent SNPs; {truth['annotations']} overlapping binary annotations "
              f"plus intercept; Gaussian summary-data sample size {truth['sample_size']:,}; "
              f"realized h2={truth['hsq']:.5f}.", "",
              prior_description +
              "SNP membership follows the sequential probit model shared by RC and RD; SNP effects follow the "
              "0,0.001,0.01,0.1,1 variance mixture. SNP base variance is fixed for generation. "
              "Fits learn it with --hsq-percentage-model false. LD is exact identity; "
              "no genotypes or LD reference sample are generated. Summary noise variance "
              "is total phenotype variance / N, matching the sparse GWAS likelihood.", "",
              "| Model | Level | Effect correlation | Effect RMSE | 95% mean±1.96 SD coverage | Max coefficient split R-hat |",
              "|---|---|---:|---:|---:|---:|"]
    diagnostics=[]
    for model in models:
        summaries = [table(folder / f"{model}_{r}.parSetRes") for r in range(replicates)]
        core = [table(folder / f"{model}_{r}.mcmcsamples/CoreParameters.mcmcsamples.txt") for r in range(replicates)]
        for c in range(4):
            label=f"AnnoEffects_p{c+2}"
            chains=[trace_values(folder / f"{model}_{r}.mcmcsamples/{label}.mcmcsamples.txt") for r in range(replicates)]
            intercept_result=stats([[row[0] for row in chain] for chain in chains])
            diagnostics.append(f"{model} p{c+2} intercept split R-hat: {intercept_result[2]:.4f}.")
            means, rhats, covered = [], [], 0
            for j,name in enumerate(truth['names']):
                samples=[[row[j+1] for row in chain] for chain in chains]
                mean, mcse, rhat = stats(samples)
                sd=st.stdev(x for chain in samples for x in chain)
                true=truth['effects'][c][j]
                pip=st.mean(float(next(row['Mean'] for row in chain if row['Parameter']==f'AnnoPIP_p{c+2}' and row['Annotation']==name)) for chain in summaries) if model=='RD' else float('nan')
                means.append(mean); rhats.append(rhat)
                covered += abs(mean-true) <= 1.96*sd
                rows.append(f"{model} p{c+2} {name} {true:.9g} {mean:.9g} {sd:.9g} {pip:.9g} {rhat:.6g}")
            true=truth['effects'][c]; ma,mb=st.mean(true),st.mean(means)
            corr=sum((x-ma)*(y-mb) for x,y in zip(true,means))/math.sqrt(sum((x-ma)**2 for x in true)*sum((y-mb)**2 for y in means))
            rmse=math.sqrt(st.mean((x-y)**2 for x,y in zip(true,means)))
            report.append(f"| {model} | p{c+2} | {corr:.4f} | {rmse:.4f} | {covered}/{len(true)} | {max(rhats):.4f} |")
        diagnostics.append("")
        for name in ['hsq','NnzSnp','SigmaSq']:
            result=stats([[float(row[name]) for row in chain] for chain in core])
            diagnostics.append(f"{model} {name}: mean {result[0]:.6g}, MCSE {result[1]:.3g}, split R-hat {result[2]:.4f}.")
        invalid=sum(not math.isfinite(float(row['PIP'])) for r in range(replicates) for row in table(folder/f'{model}_{r}.snpRes'))
        diagnostics.append(f"{model}: {invalid} nonfinite SNP PIP entries across saved chain summaries.")
    report.extend(diagnostics)
    report.extend(["", "## RD mixture weights", "",
                   "| Level | Generating weights | Realized annotation assignment fractions | Posterior mean weights | Max weight split R-hat |",
                   "|---|---|---|---|---:|"])
    summaries=[table(folder/f'RD_{r}.parSetRes') for r in range(replicates)]
    for c in range(4):
        values=[]
        for component in ['zero','var0.0001','var0.001','var0.01','var0.1']:
            values.append(st.mean(float(next(row['Mean'] for row in chain if row['Parameter']==f'AnnoMix_p{c+2}' and row['Annotation']==component)) for chain in summaries))
        fractions=[truth['groups'][c].count(g)/truth['annotations'] for g in range(5)]
        draws=[trace_values(folder/f'RD_{r}.mcmcsamples/AnnoMix_p{c+2}.mcmcsamples.txt') for r in range(replicates)]
        rhat=max(stats([[row[g] for row in chain] for chain in draws])[2] for g in range(5))
        fmt=lambda x: ', '.join(f'{v:.3f}' for v in x)
        generating=fmt(truth['weights'][c]) if generating_model=='RD' else f"single normal, variance {truth['rc_variances'][c]}"
        assignment=fmt(fractions) if generating_model=='RD' else "not applicable"
        report.append(f"| p{c+2} | {generating} | {assignment} | {fmt(values)} | {rhat:.4f} |")
    report.extend(["", "Mixture components are ordered by variance: 0, 0.0001, 0.001, 0.01, 0.1. "
                   "Nearby tiny-variance components may be hard to distinguish even when coefficient recovery is good. "
                   "Coverage above uses a Gaussian approximation, not posterior quantiles; RD includes an exact-zero mass."])
    metadata=json.loads((folder/'fit_settings.json').read_text())
    report.extend(['', f"Fit settings: {replicates} independently seeded chains per model; "
                   f"{metadata['iterations']:,} iterations, {metadata['burn_in']:,} burn-in, thin 5. "
                   "Split R-hat is the ordinary, non-rank-normalized statistic."])
    (folder/'annotation_recovery.tsv').write_text('\n'.join(rows)+'\n')
    text='\n'.join(report)+'\n'
    (folder/'recovery.md').write_text(text)
    print(text, flush=True)


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('output',type=Path)
    parser.add_argument('--snps',type=int,default=3000)
    parser.add_argument('--annotations',type=int,default=40)
    parser.add_argument('--sample-size',type=int,default=5000000)
    parser.add_argument('--seed',type=int,default=20261011)
    parser.add_argument('--iterations',type=int,default=10000)
    parser.add_argument('--replicates',type=int,default=4)
    parser.add_argument('--include-rc',action='store_true')
    parser.add_argument('--annotation-model',choices=['RD','RC'],default='RD')
    parser.add_argument('--rc-variances',default='0.01,0.03,0.06,0.1',
                        help='Four fixed generating Gaussian variances for RC, p2 through p5')
    parser.add_argument('--annotation-reference',type=Path,
                        help='Reuse an existing sim.annot matrix with matching dimensions')
    parser.add_argument('--append-null-annotations',type=int,default=0,
                        help='Append predictors with exactly zero effects at all four levels')
    parser.add_argument('--summarize-only',action='store_true')
    args=parser.parse_args()
    rc_variances=[float(x) for x in args.rc_variances.split(',')]
    if len(rc_variances)!=4 or any(not math.isfinite(x) or x<=0 for x in rc_variances):
        parser.error('--rc-variances requires four positive finite variances')
    if args.append_null_annotations<0:
        parser.error('--append-null-annotations must be nonnegative')
    folder=args.output.resolve(); folder.mkdir(parents=True,exist_ok=True)
    models=['RD','RC'] if args.include_rc else ['RD']
    if args.summarize_only:
        truth=json.loads((folder/'truth.json').read_text())
    else:
        truth=simulate(folder,args.snps,args.annotations,args.sample_size,args.seed,
                       args.annotation_model,rc_variances,args.annotation_reference,args.append_null_annotations)
        (folder/'fit_settings.json').write_text(json.dumps(dict(iterations=args.iterations,
            burn_in=args.iterations//5, replicates=args.replicates, models=models),indent=2)+'\n')
        for model in models:
            for r in range(args.replicates):
                prefix=folder/f'{model}_{r}'
                print(f'Fitting {model}, chain {r+1}/{args.replicates}',flush=True)
                run([ROOT/'scr/gctb','--sbayes',model,'--gwas-summary',folder/'sim.ma',
                     '--ldm',folder/'identity.ldm.sparse','--annot',folder/'sim.annot',
                     '--hsq-percentage-model','false','--chain-length',args.iterations,
                     '--burn-in',args.iterations//5,'--thin','5','--out-freq','1000',
                     '--thread','1','--num-chains','1','--write-mcmc-txt',
                     '--sbayes','RC','--pi','0.6,0.1,0.1,0.1,0.1','--sbayes',model,
                     '--seed',123+r*1009,'--out',prefix],str(prefix)+'.stdout')
    recovery(folder,truth,models,args.replicates)


if __name__=='__main__':
    main()

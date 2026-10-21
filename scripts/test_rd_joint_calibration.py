#!/usr/bin/env python3
"""Fit shared-indicator RD to existing simulations and compare inclusion calibration."""
import argparse
import json
import math
from pathlib import Path
import statistics as st
import sys
sys.dont_write_bytecode = True
from test_rd_unknown_pi import ROOT, run, stats, table
from test_rd_rc_normal import trace_values


def analyze(folder):
    truth = json.loads((folder/'truth.json').read_text())
    settings = json.loads((folder/'fit_settings.json').read_text())
    names = truth['names']
    causal = [any(truth['effects'][c][j] != 0 for c in range(4)) for j in range(len(names))]
    summaries = [table(folder/f'RD_Joint_{r}.parSetRes') for r in range(settings['replicates'])]
    pips = [st.mean(float(next(row['Mean'] for row in s if row['Parameter']=='AnnoPIP' and row['Annotation']==name)) for s in summaries) for name in names]
    weights = [st.mean(float(next(row['Mean'] for row in s if row['Parameter']=='AnnoMix' and row['Annotation']==name)) for s in summaries) for name in ['zero','var0.0001','var0.001','var0.01','var0.1']]
    report = ['# Shared-indicator RD calibration', '', f"Generating model: {truth['annotation_model']} plus shared nulls. Same GWAS data and fit settings as independent RD.", '', '| Model | Mean PIP, causal | Mean PIP, null | Brier score |', '|---|---:|---:|---:|']
    independent = table(folder/'annotation_pip_calibration.tsv')
    oldp = [float(row['RDAnyPIP']) for row in independent]
    for label, values in [('RD independent (any level)',oldp),('RD_Joint',pips)]:
        report.append(f"| {label} | {st.mean(p for p,y in zip(values,causal) if y):.4f} | {st.mean(p for p,y in zip(values,causal) if not y):.4f} | {st.mean((p-y)**2 for p,y in zip(values,causal)):.5f} |")
    report += ['', '| Joint PIP threshold | Calls | True positives | False positives | Observed false discovery proportion | Predicted FDR |', '|---|---:|---:|---:|---:|---:|']
    for cutoff in [.5,.9,.95,.99]:
        selected = [j for j,p in enumerate(pips) if p>=cutoff]
        tp = sum(causal[j] for j in selected); fp = len(selected)-tp
        observed = f'{fp/len(selected):.4f}' if selected else 'NA'
        predicted = f'{st.mean(1-pips[j] for j in selected):.4f}' if selected else 'NA'
        report.append(f'| {cutoff} | {len(selected)} | {tp} | {fp} | {observed} | {predicted} |')
    report += ['', 'Posterior mean weights (zero, 0.0001, 0.001, 0.01, 0.1): '+', '.join(f'{w:.4f}' for w in weights)+'.', '', '| Level | Effect correlation | Effect RMSE | Maximum coefficient split R-hat |', '|---|---:|---:|---:|']
    traces = []
    for c in range(4):
        chains = [trace_values(folder/f'RD_Joint_{r}.mcmcsamples/AnnoEffects_p{c+2}.mcmcsamples.txt') for r in range(settings['replicates'])]
        traces.append(chains)
        estimates = [stats([[row[j+1] for row in chain] for chain in chains]) for j in range(len(names))]
        true = truth['effects'][c]; means = [x[0] for x in estimates]
        mt,me = st.mean(true),st.mean(means)
        corr = sum((t-mt)*(e-me) for t,e in zip(true,means))/math.sqrt(sum((t-mt)**2 for t in true)*sum((e-me)**2 for e in means))
        rmse = math.sqrt(st.mean((t-e)**2 for t,e in zip(true,means)))
        report.append(f'| p{c+2} | {corr:.4f} | {rmse:.4f} | {max(x[2] for x in estimates):.4f} |')
    # The shared assignment must switch all four non-intercept coefficients together.
    for r in range(settings['replicates']):
        assert len({len(traces[c][r]) for c in range(4)})==1
        for i in range(len(traces[0][r])):
            for j in range(1,len(names)+1):
                assert len({traces[c][r][i][j]==0 for c in range(4)})==1, 'Shared inclusion violated'
        assert all(math.isfinite(float(row['PIP'])) for row in table(folder/f'RD_Joint_{r}.snpRes'))
    report += ['', 'All saved draws pass the shared-inclusion check; all SNP PIPs are finite. Joint PIPs use conditional inclusion probabilities averaged over posterior samples. One simulated dataset per generating scenario; four chains are computational replicates, not independent datasets.']
    (folder/'rd_joint_calibration.tsv').write_text('Annotation TrueAnyNonzero JointPIP IndependentAnyPIP\n'+'\n'.join(f'{name} {int(y)} {p:.9g} {q:.9g}' for name,y,p,q in zip(names,causal,pips,oldp))+'\n')
    result = '\n'.join(report)+'\n'
    (folder/'rd_joint_calibration.md').write_text(result)
    print(result,flush=True)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('folders',type=Path,nargs='+')
    parser.add_argument('--summarize-only',action='store_true')
    args = parser.parse_args()
    for folder in args.folders:
        folder = folder.resolve()
        settings = json.loads((folder/'fit_settings.json').read_text())
        if not args.summarize_only:
            for r in range(settings['replicates']):
                prefix = folder/f'RD_Joint_{r}'
                print(f'Fitting {folder.name}: RD_Joint chain {r+1}',flush=True)
                run([ROOT/'scr/gctb','--sbayes','RC','--pi','0.6,0.1,0.1,0.1,0.1','--sbayes','RD_Joint',
                     '--gwas-summary',folder/'sim.ma','--ldm',folder/'identity.ldm.sparse','--annot',folder/'sim.annot',
                     '--hsq-percentage-model','false','--chain-length',settings['iterations'],'--burn-in',settings['burn_in'],
                     '--thin','5','--out-freq','1000','--thread','1','--num-chains','1','--write-mcmc-txt',
                     '--seed',123+r*1009,'--out',prefix],str(prefix)+'.stdout')
        analyze(folder)


if __name__=='__main__':
    main()

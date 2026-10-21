#!/usr/bin/env python3
"""Check RD inclusion PIPs using known annotation truth and synchronized draws."""
import argparse
import json
import math
from pathlib import Path
import statistics as st
import sys

sys.dont_write_bytecode=True
from test_rd_rc_normal import table, trace_values


def check(folder):
    truth=json.loads((folder/'truth.json').read_text())
    settings=json.loads((folder/'fit_settings.json').read_text())
    n, replicates=len(truth['names']), settings['replicates']
    causal=[int(any(truth['effects'][c][j]!=0 for c in range(4))) for j in range(n)]
    level_causal=[[int(a!=0) for a in truth['effects'][c]] for c in range(4)]
    summaries=[table(folder/f'RD_{r}.parSetRes') for r in range(replicates)]
    level_pips=[[st.mean(float(next(row['Mean'] for row in summary
                                   if row['Parameter']==f'AnnoPIP_p{c+2}' and row['Annotation']==name))
                        for summary in summaries) for name in truth['names']] for c in range(4)]
    totals=[0]*n
    draws=0
    for r in range(replicates):
        levels=[trace_values(folder/f'RD_{r}.mcmcsamples/AnnoEffects_p{c+2}.mcmcsamples.txt')
                for c in range(4)]
        if len({len(level) for level in levels})!=1:
            raise ValueError('Annotation coefficient traces are not synchronized')
        for i in range(len(levels[0])):
            for j in range(n):
                totals[j]+=any(levels[c][i][j+1]!=0 for c in range(4))
            draws+=1
    any_pips=[x/draws for x in totals]
    null=[j for j,y in enumerate(causal) if not y]
    active=[j for j,y in enumerate(causal) if y]
    rd_brier=st.mean((p-y)**2 for p,y in zip(any_pips,causal))
    rc_brier=st.mean((1-y)**2 for y in causal)
    def subgroup(values, indices):
        return f'{st.mean(values[j] for j in indices):.4f}' if indices else 'NA'
    report=['# Annotation inclusion calibration with known nulls','',
            f"Generating effects: {truth.get('annotation_model','RD')}; {n} annotations, "
            f"{len(active)} causal at some level and {len(null)} null at all levels.", '',
            'RD annotation-level PIP is the posterior probability of an exactly nonzero coefficient at any level, '
            'computed from synchronized draws, not from products of marginal PIPs. RC has no zero-effect spike; '
            'its implied nonzero probability is 1. RC sign probabilities are not inclusion PIPs.', '',
            '| Fitted model | Mean PIP, causal annotations | Mean PIP, null annotations | Any-level Brier score |',
            '|---|---:|---:|---:|',
            f'| RD | {subgroup(any_pips,active)} | {subgroup(any_pips,null)} | {rd_brier:.5f} |',
            f'| RC (implied PIP=1) | {"1.0000" if active else "NA"} | {"1.0000" if null else "NA"} | {rc_brier:.5f} |', '',
            'Brier score measures binary probability prediction error, combining calibration and resolution.', '',
            '| RD PIP threshold | Calls | True positives | False positives | Observed FDR | PIP-predicted FDR |',
            '|---|---:|---:|---:|---:|---:|']
    for cutoff in [.5,.9,.95,.99]:
        selected=[j for j,p in enumerate(any_pips) if p>=cutoff]
        tp=sum(causal[j] for j in selected); fp=len(selected)-tp
        fdr=f'{fp/len(selected):.4f}' if selected else 'NA'
        predicted=f'{st.mean(1-any_pips[j] for j in selected):.4f}' if selected else 'NA'
        report.append(f'| {cutoff:.2f} | {len(selected)} | {tp} | {fp} | {fdr} | {predicted} |')
    report.extend(['', '| Probability level | True nonzero fraction | Mean RD PIP | Mean RD PIP on wholly null annotations | Level Brier score |',
                   '|---|---:|---:|---:|---:|---:|'])
    for c in range(4):
        score=st.mean((p-y)**2 for p,y in zip(level_pips[c],level_causal[c]))
        report.append(f'| p{c+2} | {st.mean(level_causal[c]):.4f} | {st.mean(level_pips[c]):.4f} | {subgroup(level_pips[c],null)} | {score:.5f} |')
    report.extend(['','| RD any-level PIP bin | Count | Mean PIP | Observed causal fraction |','|---|---:|---:|---:|'])
    for lo,hi in [(0,.5),(.5,.8),(.8,.9),(.9,.95),(.95,.99),(.99,1.00001)]:
        indices=[j for j,p in enumerate(any_pips) if lo<=p<hi]
        if indices:
            report.append(f'| [{lo:.2f}, {min(hi,1):.2f}{"]" if hi>1 else ")"} | {len(indices)} | '
                          f'{st.mean(any_pips[j] for j in indices):.4f} | {st.mean(causal[j] for j in indices):.4f} |')
    report.extend(['','## Separate RC sign-evidence screen','',
                   'For comparison of screening behavior only, call an annotation if any level has '
                   'posterior sign probability ≥0.975 or ≤0.025. This is an unadjusted sign screen across four levels, '
                   'not a calibrated inclusion probability or an FDR-controlled procedure.'])
    rc_summaries=[table(folder/f'RC_{r}.parSetRes') for r in range(replicates)]
    sign_evidence=[]
    for name in truth['names']:
        certainty=[]
        for c in range(4):
            pos=st.mean(float(next(row['PostProbAboveZero'] for row in summary
                                   if row['Parameter']==f'AnnoEffects_p{c+2}' and row['Annotation']==name))
                        for summary in rc_summaries)
            certainty.append(max(pos,1-pos))
        sign_evidence.append(max(certainty))
    selected=[j for j,p in enumerate(sign_evidence) if p>=.975]
    tp=sum(causal[j] for j in selected);fp=len(selected)-tp
    report.append(f'RC sign screen: {len(selected)} calls, {tp} true positives, {fp} false positives.')
    report.extend(['','These results describe one realization per generating scenario. Shared null annotations '
                   'introduce cross-level sparsity not represented by default RD’s independent inclusion assignments. '
                   'Read inclusion metrics alongside chain diagnostics in recovery.md.'])
    rows=['Annotation TrueAnyNonzero RDAnyPIP RCImpliedPIP RCSignEvidence']
    for j,name in enumerate(truth['names']):
        rows.append(f'{name} {causal[j]} {any_pips[j]:.9g} 1 {sign_evidence[j]:.9g}')
    (folder/'annotation_pip_calibration.tsv').write_text('\n'.join(rows)+'\n')
    text='\n'.join(report)+'\n'
    (folder/'pip_calibration.md').write_text(text)
    print(text,flush=True)


if __name__=='__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('folders',type=Path,nargs='+')
    args=parser.parse_args()
    for folder in args.folders:
        check(folder.resolve())

#!/usr/bin/env python3
"""Integration check for sampled annotation indicators and trace alignment."""
import argparse
from pathlib import Path
import sys
import struct
sys.dont_write_bytecode = True
from test_rd_rc_normal import ROOT, run, trace_values


def matrix(path):
    return trace_values(path)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('output',type=Path)
    args = parser.parse_args()
    folder = args.output.resolve()
    folder.mkdir(parents=True,exist_ok=True)
    data = ROOT/'test/rd_null_annotations'
    names = (data/'sim.annot').read_text().splitlines()[0].split()[1:]
    for model, chains in [('RD_Joint',1),('RD_Joint',2),('RD',1),('RD',2)]:
        prefix = folder/f'{model}_{chains}'
        run([ROOT/'scr/gctb','--sbayes','RC','--pi','0.6,0.1,0.1,0.1,0.1','--sbayes',model,
             '--gwas-summary',data/'sim.ma','--ldm',data/'identity.ldm.sparse','--annot',data/'sim.annot',
             '--hsq-percentage-model','false','--chain-length','300','--burn-in','23','--thin','5',
             '--thread',chains,'--num-chains',chains,'--write-mcmc-bin','--write-mcmc-txt','--seed','123','--out',prefix],str(prefix)+'.stdout')
        samples = Path(str(prefix)+'.mcmcsamples')
        effects = [matrix(samples/f'AnnoEffects_p{c+2}.mcmcsamples.txt') for c in range(4)]
        labels = ['AnnoDelta'] if model=='RD_Joint' else [f'AnnoDelta_p{c+2}' for c in range(4)]
        expected = [[iteration,chain] for iteration in range(25,301,5) for chain in range(1,chains+1)]
        for c,label in enumerate(labels):
            trace = samples/f'{label}.mcmcsamples.txt'
            lines = trace.read_text().splitlines()
            assert lines[0].split()==['Iteration','Chain']+names
            records = [[int(x) for x in line.split()] for line in lines[1:]]
            delta = [row[2:] for row in records]
            assert [row[:2] for row in records]==expected
            assert not Path(str(trace)+'.columns.txt').exists()
            assert not Path(str(trace)+'.rows.txt').exists()
            assert len(delta)==len(expected)
            assert all(len(row)==len(names) for row in delta)
            assert all(x in (0,1) for row in delta for x in row)
            assert all(row[0]==1 for row in delta), 'Intercept must always be included'
            assert {row[j] for row in delta for j in range(1,len(names))}=={0,1}
            for level in (range(4) if model=='RD_Joint' else [c]):
                assert len(effects[level])==len(delta)
                assert all(bool(d[j])==(a[j]!=0) for d,a in zip(delta,effects[level]) for j in range(1,len(names)))
        for trace in samples.glob('*.mcmcsamples.txt'):
            lines = trace.read_text().splitlines()
            header = lines[0].split()
            assert header[:2]==['Iteration','Chain'], trace
            rows = [line.split() for line in lines[1:]]
            assert [[int(x) for x in row[:2]] for row in rows]==expected, trace
            assert all(len(row)==len(header) for row in rows), trace
            assert all(len(row)==len(header)-2 for row in trace_values(trace)), trace
        binary = samples/'SnpEffects.mcmcsamples.bin'
        raw = binary.read_bytes()
        nrow,ncol,nnz = struct.unpack_from('III',raw)
        assert nrow==len(expected) and ncol==3000
        assert len(raw)==12+12*nnz
        rowlines = Path(str(binary)+'.rows.txt').read_text().splitlines()
        assert rowlines[0].split()==['Row','Iteration','Chain']
        assert [[int(x) for x in line.split()] for line in rowlines[1:]]==[[j+1]+r for j,r in enumerate(expected)]
        collines = Path(str(binary)+'.columns.txt').read_text().splitlines()
        assert collines[0].split()==['Column','Name']
        assert [line.split() for line in collines[1:]]==[[str(j+1),f'rs{j+1}'] for j in range(3000)]
        assert all(r<nrow and c<ncol for r,c,v in struct.iter_unpack('IIf',raw[12:]))
        print(f'{model}, {chains} chain(s): named text traces, binary metadata and delta alignment passed',flush=True)
    print('Annotation delta output checks passed',flush=True)


if __name__=='__main__':
    main()

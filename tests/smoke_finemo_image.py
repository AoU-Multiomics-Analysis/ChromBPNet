"""Run real Fi-NeMo on planted allele motifs and check gains, losses, and empty rows."""
import csv
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import zipfile

import h5py
import numpy as np

SCRIPTS = Path('/opt/finemo_pipeline/scripts')
sys.path.insert(0, str(SCRIPTS))
from finemo_io import sha256


def run(script, *args):
    subprocess.run([sys.executable, str(SCRIPTS/script), *map(str,args)], check=True)


def table(path):
    with path.open() as stream:
        return list(csv.DictReader(stream,delimiter='\t'))


def main():
    os.environ['OMP_NUM_THREADS'] = '2'
    os.environ['NUMBA_NUM_THREADS'] = '2'
    with tempfile.TemporaryDirectory() as tmp:
        root=Path(tmp)
        motif='ACGTCAGTGCAT'
        width=64
        cwm=np.zeros((len(motif),4),dtype=np.float32)
        for i,b in enumerate(motif):
            cwm[i,'ACGT'.index(b)]=1
        motifs=root/"motifs ' $(touch NEVER).h5"
        with h5py.File(motifs,'w') as handle:
            pattern=handle.create_group('pos_patterns/pattern_0')
            pattern.create_dataset('sequence',data=cwm)
            pattern.create_dataset('contrib_scores',data=cwm)
            pattern.create_dataset('hypothetical_contribs',data=cwm)
        # Interleaved REF/ALT: lost, gained, retained with weaker ALT, no calls.
        raw=np.zeros((8,4,width),dtype=np.int8)
        raw[:,0,:]=1
        hyp=np.zeros_like(raw,dtype=np.float32)
        for row,scale in [(0,2),(3,2),(4,2),(5,1)]:
            raw[row,:,26:38]=cwm.T
            hyp[row,:,26:38]=cwm.T*scale
        rows=root/'rows.tsv'
        lines=[]
        for i in range(8):
            lines.append([i,['lost','gained','retained','empty'][i//2], 'REF' if i%2==0 else 'ALT',
                          'chr1',101+i//2*100,'A','T',32,33,68+i//2*100,100+i//2*100])
        with rows.open('w') as stream:
            csv.writer(stream,delimiter='\t',lineterminator='\n').writerows(lines)
        coordinates=root/'coordinates.h5'
        with h5py.File(coordinates,'w') as handle:
            handle.create_dataset('raw/seq',data=raw)
            handle.create_dataset('reference_positions',data=np.stack([np.arange(68+i//2*100,132+i//2*100) for i in range(8)]))
            handle.attrs['region_sha256']=sha256(rows)
        average=root/'mean.h5'
        with h5py.File(average,'w') as handle:
            handle.create_dataset('raw/seq',data=raw)
            handle.create_dataset('shap/seq',data=hyp)
            handle.attrs.update(cell_type='synthetic_CD4',head='counts',fold_count=5,averaging='arithmetic_mean_signed',region_sha256=sha256(rows))
        provenance=root/'motif_provenance.json'
        provenance.write_text(json.dumps(dict(cell_type='synthetic_CD4',head='counts')))
        matches=root/'matches.tsv'
        matches.write_text('pattern\ttarget_motif\tcandidate_tf\tq_value\n'
                           'pos_patterns.pattern_0\tm1\tSynthetic_TF\t0.001\n'
                           'pos_patterns.pattern_0\tm2\tRelated_TF\t0.002\n')
        prepared=root/'prepared'
        run('prepare_finemo.py','--contributions',average,'--regions',rows,'--motifs',motifs,
            '--motif-provenance',provenance,'--coordinates',coordinates,'--kind','variants',
            '--cell-type','synthetic_CD4','--output-dir',prepared)
        common=['--regions',prepared/'regions.npz','--motifs',motifs,'--preparation',prepared/'metadata.json',
                '--threads','2','--batch-size','4']
        denied=subprocess.run([sys.executable,str(SCRIPTS/'call_finemo.py'),*map(str,common),
                               '--output-dir',str(root/'denied')],capture_output=True,text=True)
        assert denied.returncode!=0 and 'GPU is required' in denied.stderr
        calls=root/'calls'
        run('call_finemo.py',*common,'--allow-cpu','--output-dir',calls)
        summary=root/'summary'
        run('summarize_finemo.py','--hits',calls/'hits.tsv','--qc',calls/'peaks_qc.tsv',
            '--regions',rows,'--coordinates',coordinates,'--preparation',prepared/'metadata.json',
            '--call-metadata',calls/'metadata.json','--tf-matches',matches,'--output-dir',summary)
        results={r['variant_id']:r for r in table(summary/'variant_summary.tsv')}
        assert set(results)=={'lost','gained','retained','empty'},results
        assert int(results['lost']['n_lost'])>=1,results
        assert int(results['gained']['n_gained'])>=1,results
        assert int(results['retained']['n_retained'])>=1,results
        assert int(results['empty']['n_ref_hits'])==int(results['empty']['n_alt_hits'])==0
        retained=[r for r in table(summary/'motif_changes.tsv') if r['variant_id']=='retained' and r['status']=='retained']
        assert retained and float(retained[0]['delta_hit_coefficient_global'])<0,retained
        assert 'Synthetic_TF' in retained[0]['candidate_tfs'] and 'Related_TF' in retained[0]['candidate_tfs']
        assert not (root/'NEVER').exists()
        for bundle in [calls/'raw_calls.zip',summary/'report.zip']:
            with zipfile.ZipFile(bundle) as archive:
                assert archive.testzip() is None
                assert not any(name.endswith('.zip') for name in archive.namelist())
        # Exercise absolute genomic coordinates with the existing discovery inputs.
        peaks=root/'peaks.bed'
        peaks.write_text(''.join(f'chr1\t{80+i*100}\t{120+i*100}\tp{i}\t0\t.\t0\t0\t0\t20\n' for i in range(8)))
        with h5py.File(average,'r+') as handle:
            handle.attrs['region_sha256']=sha256(peaks)
        provenance.write_text(json.dumps(dict(cell_type='synthetic_CD4',head='counts',prepared_peak_sha256=sha256(peaks))))
        peak_prepared=root/'peak_prepared'
        run('prepare_finemo.py','--contributions',average,'--regions',peaks,'--motifs',motifs,
            '--motif-provenance',provenance,'--cell-type','synthetic_CD4','--output-dir',peak_prepared)
        peak_calls=root/'peak_calls'
        run('call_finemo.py','--regions',peak_prepared/'regions.npz','--motifs',motifs,
            '--preparation',peak_prepared/'metadata.json','--threads','2','--allow-cpu','--output-dir',peak_calls)
        peak_hits=table(peak_calls/'hits_unique.tsv')
        assert {int(row['peak_id']) for row in peak_hits} == {0,3,4,5},peak_hits
        print('PASS: real Fi-NeMo calls, REF/ALT gains/losses/strength, empty variants, TF annotation, genomic coordinates, report, and GPU requirement')


if __name__=='__main__':
    main()

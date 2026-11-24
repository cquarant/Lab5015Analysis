#!/usr/bin/env python3
"""Check selected production coverage, REF7 baselines and report assets."""
import argparse
from collections import Counter
import csv
from html.parser import HTMLParser
import json
from pathlib import Path
import re
from urllib.parse import urlsplit

RUNS=(4194,4257,4281,4291)
REFS=(4,5,7,8,9,10,11,12)
THRESHOLDS=(3,5,7,10,15,20)
SAFE_COUNTS={4194:848605,4257:1040888,4281:1176158,4291:776485}


class LocalLinks(HTMLParser):
    def __init__(self):
        super().__init__()
        self.links=[]

    def handle_starttag(self, tag, attrs):
        self.links.extend(urlsplit(value).path for key,value in attrs if key in ('href','src')
                          and value and not value.startswith(('#','../','http')))


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('work',type=Path)
    parser.add_argument('--safe-baseline',type=Path,required=True)
    parser.add_argument('--spill-manifest',type=Path,default=Path(
        '/eos/home-c/cquarant/TB_Sept2025_Studies4DN/timing_closure_withref_TOFHIR/input_spills.csv'))
    args=parser.parse_args()
    import ROOT

    def file_uuid(path):
        root=ROOT.TFile.Open(str(path),'READ')
        assert root and not root.IsZombie(),path
        result=root.GetUUID().AsString()
        root.Close()
        return result

    work=args.work
    published=json.loads((work/'published_validation/validation.json').read_text())
    with args.spill_manifest.open() as stream:
        spills=list(csv.DictReader(stream))
    assert {r['run'] for r in published if r['passed']}==set(RUNS)
    subsets=[]
    aggregate_groups=0
    max_residual_covariance_identity=0.
    max_catr_identity=0.
    for run in RUNS:
        expected=[row for row in spills if int(row['run'])==run and float(row['vov_first_entry'])==3.]
        expected_filenames=Counter(Path(row['reco_path']).name for row in expected)
        expected_raw_entries=sum(int(row['entries']) for row in expected)
        for ref in REFS:
            folder=work/'production_final'/str(run)/f'ref{ref:02d}'
            receipt=json.loads((folder/'complete.json').read_text())
            assert all((folder/path).stat().st_size==size for path,size in receipt['sizes'].items())
            assert receipt['run']==run and receipt['ref']==ref
            replay=json.loads((folder/'closure/summary.json').read_text())
            actual={key:file_uuid(folder/name) for key,name in
                    [('reference','refCharacterization.root'),('corrections','step2.root')]}
            assert actual==replay['input_uuids'],(run,ref)
            flat=json.loads((work/'flat_final'/str(run)/f'ref{ref:02d}/complete.json').read_text())
            assert flat['reference_uuid']==actual['reference'],(run,ref)
            assert flat['closure_uuid']==file_uuid(folder/'closure/timing_closure.root'),(run,ref)
            step1_log=(folder/'moduleCharacterization_step1.log').read_text()
            assert Counter(re.findall(r'Adding file: (\S+)',step1_log))==expected_filenames,(run,ref)
            raw_entries=re.search(r'external bar loop: reading entry 0 / (\d+)',step1_log)
            assert raw_entries and int(raw_entries[1])==expected_raw_entries,(run,ref)
            if ref==7:
                assert receipt['groups']==96 and receipt['accepted']==SAFE_COUNTS[run],receipt
                current=replay['groups']
                safe=json.loads((args.safe_baseline/str(run)/'summary.json').read_text())['groups']
                assert len(current)==len(safe)==96
                for new,old in zip(current,safe):
                    for key in ('label','counts','source_entries_sha256'):
                        assert new[key]==old[key],(run,new['label'],key)
                    for stage in ('tw','full'):
                        for key in ('moments','gaussian','fits'):
                            assert new['stages'][stage][key]==old['stages'][stage][key],(run,new['label'],stage,key)
            subsets.append({key:receipt[key] for key in ('run','ref','groups','skipped','accepted')})
        for bar in range(16):
            for threshold in THRESHOLDS:
                label=f'bar{bar:02d}_Vov3.00_th{threshold:02d}_energyBin01'
                group=json.loads((work/'aggregate'/str(run)/label/'summary.json').read_text())
                aggregate_groups+=1
                assert group['run']==run and group['bar']==bar and group['threshold']==threshold
                assert set(group['reference_bars'])|{int(r) for r in group['excluded']}==set(REFS)
                assert sum(group['counts_by_ref'].values())==group['pairs']
                assert group['pairs']-group['unique_events']==group['duplicated_pairs']
                for stage in ('tw','full'):
                    for sample in ('union','centered','single_ref_per_event'):
                        result=group['stages'][stage][sample]
                        cov=result['covariance']
                        max_residual_covariance_identity=max(max_residual_covariance_identity,
                            abs(cov['residual_covariance_ps2']-cov['from_average_minus_matched_ps2']))
                        max_catr_identity=max(max_catr_identity,abs(result['moments']['closure_ps2']))
    assert max_residual_covariance_identity<1e-6 and max_catr_identity<1e-6
    report=work/'report'
    for page in ('index.html','results.html','selection.html','aggregation.html','diagnostics.html','methodology.html','commands.html'):
        parser=LocalLinks()
        parser.feed((report/page).read_text())
        for path in parser.links:
            assert (report/path).is_file(),path
    for run in RUNS:
        for threshold in THRESHOLDS:
            for estimator in ('moments','gaussian'):
                for suffix in ('png','pdf'):
                    for name in ('delta', 'methods_centered_vs_ref7'):
                        assert (report/'plots'/f'{name}_{run}_{threshold}_{estimator}.{suffix}').stat().st_size>0
        for sample in ['union','centered','single_ref_per_event','published']+[f'ref{r:02d}' for r in REFS]:
            for threshold in THRESHOLDS:
                for estimator in ('moments','gaussian'):
                    stem=f'timing_{estimator}_Vov3.00_th{threshold:02d}'
                    for suffix in ('png','pdf'):
                        assert (report/'plots'/str(run)/sample/f'{stem}.{suffix}').stat().st_size>0
    with (report/'data/all_results.csv').open() as stream:
        rows=list(csv.DictReader(stream))
    keys={(int(r['run']),r['sample'],int(r['bar']),int(r['threshold']),r['stage'],r['estimator']) for r in rows}
    assert len(keys)==len(rows)
    assert len(rows)==4*(sum(r['groups'] for r in subsets)+4*aggregate_groups)
    validation=dict(published_REF7=published,safe_REF7_counts=SAFE_COUNTS,selected_REF_bars=REFS,
                    subsets=subsets,aggregate_groups=aggregate_groups,result_rows=len(rows),
                    max_residual_covariance_identity_ps2=max_residual_covariance_identity,
                    max_catr_identity_ps2=max_catr_identity)
    (report/'data/validation.json').write_text(json.dumps(validation,indent=2)+'\n')
    print(f'Validated {len(subsets)} run–REF samples, {aggregate_groups} aggregate groups, {len(rows)} table rows and all plot assets.')


if __name__=='__main__':
    main()

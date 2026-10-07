"""Verify released bytes and source-to-segment mappings, without refitting data."""
from pathlib import Path
import csv, hashlib, json, math, sys

ROOT=Path(__file__).resolve().parent

def rows(relative):
    with (ROOT/relative).open(encoding='utf-8-sig',newline='') as stream:
        return list(csv.DictReader(stream))

def close(a,b,tol=1e-12):
    a,b=float(a),float(b)
    return (math.isnan(a) and math.isnan(b)) or abs(a-b)<=tol

def main():
    manifest=json.loads((ROOT/'metadata/checksums.json').read_text(encoding='utf-8'))
    errors=[]
    for name, expected in manifest['sha256'].items():
        path=(ROOT/name).resolve()
        if not path.is_relative_to(ROOT) or not path.is_file():
            errors.append('Missing or invalid path: '+name)
        elif hashlib.sha256(path.read_bytes()).hexdigest()!=expected:
            errors.append('Hash differs: '+name)
    fits=rows('data/fitting_segments.csv')
    comparisons=rows('data/figure6_comparisons.csv')
    assert len(fits)==17 and len(comparisons)==6
    assert sum(x['channel']=='yaw' for x in fits)==9
    assert all(x['participant_id']=='participant_02' for x in fits)
    assert sum(x['participant_id']=='participant_01' for x in comparisons)==4
    record_paths=[x['source_record_csv'] for x in fits+comparisons]
    assert len(set(record_paths))==23
    total_rows=sum(len(rows(p)) for p in record_paths)
    assert total_rows==36754
    for meta in fits:
        record=rows(meta['source_record_csv'])
        segment=rows(meta['segment_csv'])
        start=int(meta['source_start_sample_1based'])-1
        end=int(meta['source_end_sample_1based'])
        selected=record[start:end]
        assert len(segment)==len(selected)==int(meta['samples'])
        angle='smoothed_horizontal_left_deg' if meta['channel']=='yaw' else 'smoothed_vertical_left_deg'
        sign=-1 if meta['channel']=='yaw' else 1
        origin=float(meta['segment_time_origin_record_s'])
        sdk_origin=float(selected[0]['sdk_time_left_s'])
        for s,r in zip(segment,selected):
            assert close(s['angle_deg'],sign*float(r[angle]),0)
            assert close(s['time_s'],float(r['elapsed_time_s'])-origin)
            assert close(s['left_sdk_time_s'],float(r['sdk_time_left_s'])-sdk_origin)
    for meta in comparisons:
        record=rows(meta['source_record_csv'])
        segment=rows(meta['selected_segment_csv'])
        start=int(meta['source_start_row_1based'])-1
        end=int(meta['source_end_row_1based'])
        selected=record[start:end]
        assert len(segment)==len(selected)==int(meta['selected_samples'])
        for s,r in zip(segment,selected):
            for key in record[0]:
                assert close(s[key],r[key],0)
            assert close(s['segment_time_s'],float(r['elapsed_time_s'])-float(meta['segment_time_origin_record_s']))
    if errors:
        print('\n'.join(errors),file=sys.stderr)
        print('Rerunning analyses can alter output hashes; use a clean versioned checkout for byte-level verification.',file=sys.stderr)
        return 1
    print(json.dumps({'file_hashes':'PASS','checked_files':len(manifest['sha256']),
                      'source_recordings':23,'saved_source_rows':total_rows,
                      'fitting_segments':17,'comparison_windows':6,
                      'source_segment_mappings':'PASS'},indent=2))
    return 0

if __name__=='__main__':
    raise SystemExit(main())

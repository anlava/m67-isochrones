#!/usr/bin/env python3
"""Validate all three source data files."""
import sys
import csv
import io

def validate_members(path):
    """Validate m67_members_hr2023.tsv"""
    print(f"=== Validating {path} ===")
    errors = []
    data_rows = 0
    header_found = False
    
    with open(path, 'r') as f:
        lines = f.readlines()
    
    # Find header line
    header_line_idx = None
    for i, line in enumerate(lines):
        stripped = line.strip()
        if stripped.startswith('Name\t') and 'GaiaDR3' in stripped:
            header_line_idx = i
            header_found = True
            expected_cols = ['Name', 'GaiaDR3', 'Prob', 'inrt', 'Gmag', 'BPmag', 'RPmag']
            cols = stripped.rstrip('\n').split('\t')
            if cols != expected_cols:
                errors.append(f"Line {i+1}: header mismatch: {cols}")
            break
    
    if not header_found:
        errors.append("Header line not found!")
        return errors
    
    print(f"  Header at line {header_line_idx + 1}")
    
    # Parse data rows using csv reader with tab delimiter
    data_start = header_line_idx + 1
    # Skip units line and separator line
    for i in range(header_line_idx + 1, len(lines)):
        stripped = lines[i].strip()
        if stripped == '' or stripped.startswith('---') or stripped == 'mag\tmag\tmag':
            data_start = i + 1
            break
    
    print(f"  Data starts at line {data_start + 1}")
    
    reader = csv.reader(lines[data_start:], delimiter='\t')
    for row_num, row in enumerate(reader, start=data_start + 1):
        # Skip empty rows
        if not row or all(c.strip() == '' for c in row):
            continue
        
        data_rows += 1
        
        if len(row) != 7:
            errors.append(f"Line {row_num}: expected 7 columns, got {len(row)}: '{row}'")
            continue
        
        name, gaia_dr3, prob, inrt, gmag, bpmag, rpmag = row
        
        # Check GaiaDR3 is numeric
        try:
            int(gaia_dr3.strip())
        except ValueError:
            errors.append(f"Line {row_num}: GaiaDR3 not integer: '{gaia_dr3}'")
        
        # Check Prob is float
        try:
            float(prob.strip())
        except ValueError:
            errors.append(f"Line {row_num}: Prob not float: '{prob}'")
        
        # Check inrt is boolean
        inrt_val = inrt.strip()
        if inrt_val not in ('True', 'False'):
            errors.append(f"Line {row_num}: inrt not boolean: '{inrt_val}'")
        
        # Check magnitudes are numeric
        for col_name, val in [('Gmag', gmag), ('BPmag', bpmag), ('RPmag', rpmag)]:
            try:
                float(val.strip())
            except ValueError:
                errors.append(f"Line {row_num}: {col_name} not numeric: '{val}'")
    
    print(f"  Data rows: {data_rows}")
    print(f"  Errors: {len(errors)}")
    for e in errors[:10]:
        print(f"    {e}")
    return errors

def validate_parameters(path):
    """Validate m67_parameters_cavallo2024.tsv"""
    print(f"\n=== Validating {path} ===")
    errors = []
    data_rows = 0
    
    with open(path, 'r') as f:
        lines = f.readlines()
    
    # Find header line
    header_line_idx = None
    for i, line in enumerate(lines):
        stripped = line.strip()
        if stripped.startswith('Cluster\t') and 'RA_ICRS' in stripped:
            header_line_idx = i
            break
    
    if header_line_idx is None:
        errors.append("Header line not found!")
        return errors
    
    print(f"  Header at line {header_line_idx + 1}")
    
    # Find data start (skip units and separator)
    data_start = header_line_idx + 1
    for i in range(header_line_idx + 1, len(lines)):
        stripped = lines[i].strip()
        if stripped == '' or stripped.startswith('---'):
            data_start = i + 1
            break
    
    print(f"  Data starts at line {data_start + 1}")
    
    reader = csv.reader(lines[data_start:], delimiter='\t')
    for row_num, row in enumerate(reader, start=data_start + 1):
        if not row or all(c.strip() == '' for c in row):
            continue
        
        data_rows += 1
        cols = [c.strip() for c in row]
        
        if len(cols) < 22:
            errors.append(f"Line {row_num}: expected >=22 columns, got {len(cols)}")
            continue
        
        cluster = cols[0]
        if cluster != 'NGC_2682':
            errors.append(f"Line {row_num}: cluster name is '{cluster}', expected 'NGC_2682'")
        
        # Check key values
        try:
            logAge50 = float(cols[10])
            fe_h50 = float(cols[12])
            av50 = float(cols[15])
            dmod50 = float(cols[17])
            print(f"  logAge50={logAge50}, [Fe/H]50={fe_h50}, AV50={av50}, dMod50={dmod50}")
        except (ValueError, IndexError) as e:
            errors.append(f"Line {row_num}: error reading key values: {e}")
    
    print(f"  Data rows: {data_rows}")
    print(f"  Errors: {len(errors)}")
    for e in errors:
        print(f"    {e}")
    return errors

def validate_isochrone(path):
    """Validate m67_parsec12s_av018.dat"""
    print(f"\n=== Validating {path} ===")
    errors = []
    data_rows = 0
    termination_point = None
    
    with open(path, 'r') as f:
        lines = f.readlines()
    
    # Find header line (starts with # Zini)
    header_line_idx = None
    for i, line in enumerate(lines):
        if line.startswith('# Zini'):
            header_line_idx = i
            break
    
    if header_line_idx is None:
        errors.append("Header line not found!")
        return errors
    
    # Parse header columns
    header_line = lines[header_line_idx].strip()
    # The header has column names separated by spaces/tabs
    # Let's count expected columns from the header
    # Zini MH logAge Mini int_IMF Mass logL logTe logg label McoreTP C_O period0 period1 period2 period3 period4 pmode Mloss tau1m X Y Xc Xn Xo Cexcess Z mbolmag Gmag G_BPmag G_RPmag
    expected_cols = ['Zini', 'MH', 'logAge', 'Mini', 'int_IMF', 'Mass', 'logL', 'logTe', 'logg', 
                     'label', 'McoreTP', 'C_O', 'period0', 'period1', 'period2', 'period3', 'period4',
                     'pmode', 'Mloss', 'tau1m', 'X', 'Y', 'Xc', 'Xn', 'Xo', 'Cexcess', 'Z', 
                     'mbolmag', 'Gmag', 'G_BPmag', 'G_RPmag']
    print(f"  Expected columns: {len(expected_cols)}")
    
    # Find data start
    data_start = header_line_idx + 1
    for i in range(header_line_idx + 1, len(lines)):
        stripped = lines[i].strip()
        if stripped != '' and not stripped.startswith('#'):
            data_start = i
            break
    
    print(f"  Data starts at line {data_start + 1}")
    
    for i in range(data_start, len(lines)):
        line = lines[i].strip()
        if line == '' or line.startswith('#'):
            continue
        
        data_rows += 1
        cols = line.split()
        
        if len(cols) != len(expected_cols):
            errors.append(f"Line {i+1}: expected {len(expected_cols)} columns, got {len(cols)}")
            continue
        
        try:
            zini = float(cols[0])
            mh = float(cols[1])
            log_age = float(cols[2])
            mini = float(cols[3])
            mass = float(cols[5])
            log_l = float(cols[6])
            label = int(cols[9])
            gmag = float(cols[29])
            g_bp = float(cols[30])
            g_rp = float(cols[31])
        except (ValueError, IndexError) as e:
            errors.append(f"Line {i+1}: error parsing values: {e}")
            continue
        
        # Check for termination point
        if label == 9 and log_l == -9.999:
            termination_point = i + 1
            print(f"  Termination point found at line {termination_point}: label={label}, logL={log_l}")
    
    print(f"  Data rows: {data_rows}")
    print(f"  Termination point: line {termination_point}")
    print(f"  Errors: {len(errors)}")
    for e in errors[:10]:
        print(f"    {e}")
    if len(errors) > 10:
        print(f"    ... and {len(errors)-10} more")
    return errors

if __name__ == '__main__':
    from pathlib import Path
    base = Path(__file__).resolve().parents[1] / "data" / "raw"
    
    all_errors = []
    all_errors.extend(validate_members(base / 'm67_members_hr2023.tsv'))
    all_errors.extend(validate_parameters(base / 'm67_parameters_cavallo2024.tsv'))
    all_errors.extend(validate_isochrone(base / 'm67_parsec12s_av018.dat'))
    
    print(f"\n=== TOTAL ERRORS: {len(all_errors)} ===")
    if all_errors:
        sys.exit(1)
    else:
        print("All files validated OK.")
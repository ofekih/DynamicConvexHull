import json
import math
import sys

def main():
    if len(sys.argv) < 2:
        print("Usage: python analyze_complexity.py <benchmark_json_file>")
        sys.exit(1)

    try:
        with open(sys.argv[1], 'r') as f:
            content = f.read()
            # Fix potentially truncated JSON
            content = content.strip()
            if content.endswith(','):
                content = content[:-1]
            if not content.endswith('}'):
                if not content.endswith(']'):
                     content += ']'
                content += '}'
            data = json.loads(content)
    except Exception as e:
        print(f"Error reading JSON: {e}")
        # Try a more aggressive fix if the standard simple fix failed
        try:
            # Find the last closing brace and cut off everything after
            last_brace = content.rfind('}')
            if last_brace != -1:
                content = content[:last_brace+1] + ']}'
                data = json.loads(content)
            else:
                sys.exit(1)
        except:
            print("Could not recover JSON data.")
            sys.exit(1)

    print(f"{'Benchmark':<30} | {'N':<8} | {'Time (ns)':<12} | {'T/N':<10} | {'T/N^2':<10} | {'T/N^3':<10} | {'T/log^2N':<10}")
    print("-" * 110)

    for bm in data['benchmarks']:
        name = bm['name']
        if "BM_FindStabbingLine" not in name or "_mean" in name or "_median" in name or "_stddev" in name or "BigO" in name or "RMS" in name:
            continue
            
        try:
            n_part = name.split('/')[-1]
            n = int(n_part)
        except ValueError:
            continue
            
        time_ns = bm['real_time']
        
        if n > 1:
            log_n = math.log2(n)
            t_n = time_ns / n
            t_n2 = time_ns / (n**2)
            t_n3 = time_ns / (n**3)
            t_log2n = time_ns / (log_n ** 2)
        else:
            t_n = t_n2 = t_n3 = t_log2n = 0

        print(f"{name:<30} | {n:<8} | {time_ns:<12.2f} | {t_n:<10.2f} | {t_n2:<10.2f} | {t_n3:<10.2f} | {t_log2n:<10.2f}")
        
    print("\nInterpretation:")
    print("If T/f(N) is roughly constant as N increases, then Time ~ O(f(N)).")
    print("Look for the column that stays the most stable across different N values.")

if __name__ == "__main__":
    main()

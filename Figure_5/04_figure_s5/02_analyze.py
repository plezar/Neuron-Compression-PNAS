import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy import stats

def analyze_middle_20_slices(data, use_unnormalized=True):
    # grab middle 20 slices from each series and average them
    series_averages = []
    
    # pick which column to use
    fosb_column = 'mean_fosb_per_pixel_astrocytes' if use_unnormalized else 'fosb_per_pixel_ratio'
    
    for (replicate_id, series_num), series_data in data.groupby(['replicate_id', 'series_number']):
        series_data = series_data.sort_values('slice_index')
        
        # get middle 20 slices
        total_slices = len(series_data)
        mid_point = total_slices // 2
        start_idx = mid_point - 10
        end_idx = start_idx + 20
        
        middle_data = series_data.iloc[start_idx:end_idx]
        avg_ratio = middle_data[fosb_column].mean()
        
        series_averages.append({
            'replicate_id': replicate_id,
            'series_number': series_num,
            'group': 'experimental' if '_e' in replicate_id.lower() else 'control',
            'avg_fosb_ratio': avg_ratio,
            'total_nuclei': middle_data['total_nuclei'].iloc[0],
            'astrocyte_count': middle_data['astrocyte_count'].iloc[0],
            'astrocyte_percentage': 100 * middle_data['astrocyte_count'].iloc[0] / middle_data['total_nuclei'].iloc[0],
            'slice_range': f"{start_idx}-{end_idx}"
        })
    
    return pd.DataFrame(series_averages)

def create_summary_plot(data):
    plt.figure(figsize=(10, 8))
    
    ctrl_data = data[data['group'] == 'control']['avg_fosb_ratio']
    exp_data = data[data['group'] == 'experimental']['avg_fosb_ratio']
    
    # boxplot
    box_data = [ctrl_data, exp_data]
    plt.boxplot(box_data, labels=['control', 'experimental'])
    
    # add points with jitter so you can see them
    x_ctrl = np.random.normal(1, 0.04, size=len(ctrl_data))
    x_exp = np.random.normal(2, 0.04, size=len(exp_data))
    plt.scatter(x_ctrl, ctrl_data, alpha=0.6, color='blue', label='control')
    plt.scatter(x_exp, exp_data, alpha=0.6, color='red', label='experimental')
    
    # stats
    t_stat, p_val = stats.ttest_ind(exp_data, ctrl_data)
    effect_size = exp_data.mean() / ctrl_data.mean()
    
    plt.title('middle 20 slices analysis\n' + 
             f'control (n={len(ctrl_data)}): {ctrl_data.mean():.3f} ± {ctrl_data.std():.3f}\n' +
             f'experimental (n={len(exp_data)}): {exp_data.mean():.3f} ± {exp_data.std():.3f}\n' +
             f'p = {p_val:.2e}, effect size: {effect_size:.2f}x')
    
    plt.ylabel('fosb ratio (per pixel)')
    plt.grid(True, alpha=0.3)
    
    plt.savefig('middle_20_slices_analysis.png', dpi=300, bbox_inches='tight')
    plt.close()
    
    return {
        'control_mean': ctrl_data.mean(),
        'control_std': ctrl_data.std(),
        'control_n': len(ctrl_data),
        'experimental_mean': exp_data.mean(),
        'experimental_std': exp_data.std(),
        'experimental_n': len(exp_data),
        'p_value': p_val,
        'effect_size': effect_size
    }

def main():
    print("loading analysis results...")
    data = pd.read_csv('/users/charlessander/desktop/patzke lab computing/datta lab/actual final/exact_tiff_nuclei_nd2_fosb_slice_data.csv')
    
    print("\nanalyzing middle 20 slices...")
    results = analyze_middle_20_slices(data, use_unnormalized=False)
    
    print("\ncreating plot...")
    stats = create_summary_plot(results)
    
    # save outputs
    results.to_csv('middle_20_slices_data.csv', index=False)
    pd.DataFrame([stats]).to_csv('middle_20_slices_stats.csv', index=False)
    
    print("\ndone! files saved:")
    print("  middle_20_slices_analysis.png")
    print("  middle_20_slices_data.csv")
    print("  middle_20_slices_stats.csv")
    
    print(f"\ncontrol: {stats['control_mean']:.3f} ± {stats['control_std']:.3f} (n={stats['control_n']})")
    print(f"experimental: {stats['experimental_mean']:.3f} ± {stats['experimental_std']:.3f} (n={stats['experimental_n']})")
    print(f"p-value: {stats['p_value']:.2e}")
    print(f"effect size: {stats['effect_size']:.2f}x")

if __name__ == '__main__':
    main()

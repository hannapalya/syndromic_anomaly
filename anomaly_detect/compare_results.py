#!/usr/bin/env python3
"""
Compare LSTM results between two CSV files with the same structure
"""

import pandas as pd
import numpy as np

def load_results(filepath):
    """Load results from CSV file"""
    try:
        df = pd.read_csv(filepath, index_col=0)
        print(f"Loaded {len(df)} signals from {filepath}")
        return df
    except FileNotFoundError:
        print(f"File not found: {filepath}")
        return None
    except Exception as e:
        print(f"Error loading {filepath}: {e}")
        return None

def compare_results(old_df, new_df):
    """Compare two result dataframes"""
    if old_df is None or new_df is None:
        print("Cannot compare - one or both files failed to load")
        return None
    
    # Ensure both have the same columns
    common_columns = list(set(old_df.columns) & set(new_df.columns))
    if not common_columns:
        print("No common columns found between the two files")
        return None
    
    print(f"Comparing on {len(common_columns)} common columns: {common_columns}")
    
    # Get common signals
    common_signals = list(set(old_df.index) & set(new_df.index))
    if not common_signals:
        print("No common signals found between the two files")
        return None
    
    print(f"Comparing {len(common_signals)} common signals")
    
    comparison_data = []
    
    for signal in sorted(common_signals):
        old_row = old_df.loc[signal]
        new_row = new_df.loc[signal]
        
        for col in common_columns:
            if col in old_row and col in new_row and pd.notna(old_row[col]) and pd.notna(new_row[col]):
                old_val = float(old_row[col])
                new_val = float(new_row[col])
                improvement = new_val - old_val
                pct_change = (improvement / old_val * 100) if old_val != 0 else 0
                
                comparison_data.append({
                    'Signal': signal,
                    'Metric': col,
                    'Old_Value': old_val,
                    'New_Value': new_val,
                    'Improvement': improvement,
                    'Pct_Change': pct_change
                })
    
    return pd.DataFrame(comparison_data)

def print_detailed_comparison(comparison_df):
    """Print detailed comparison results"""
    print("\n" + "=" * 100)
    print("DETAILED COMPARISON: OLD vs NEW LSTM RESULTS")
    print("=" * 100)
    
    # Group by signal
    for signal in sorted(comparison_df['Signal'].unique()):
        signal_data = comparison_df[comparison_df['Signal'] == signal]
        print(f"\n--- Signal {signal} ---")
        print(f"{'Metric':<20} {'Old':<10} {'New':<10} {'Change':<10} {'% Change':<10}")
        print("-" * 70)
        
        for _, row in signal_data.iterrows():
            print(f"{row['Metric']:<20} {row['Old_Value']:<10.3f} {row['New_Value']:<10.3f} "
                  f"{row['Improvement']:<+10.3f} {row['Pct_Change']:<+10.1f}%")

def print_summary_comparison(comparison_df):
    """Print summary statistics"""
    print("\n" + "=" * 100)
    print("SUMMARY STATISTICS")
    print("=" * 100)
    
    # Group by metric
    for metric in sorted(comparison_df['Metric'].unique()):
        metric_data = comparison_df[comparison_df['Metric'] == metric]
        
        old_vals = metric_data['Old_Value'].values
        new_vals = metric_data['New_Value'].values
        improvements = metric_data['Improvement'].values
        
        print(f"\n{metric.upper()}:")
        print(f"  Old Average: {np.mean(old_vals):.3f} ± {np.std(old_vals):.3f}")
        print(f"  New Average: {np.mean(new_vals):.3f} ± {np.std(new_vals):.3f}")
        print(f"  Average Improvement: {np.mean(improvements):+.3f} ± {np.std(improvements):.3f}")
        print(f"  Signals Improved: {sum(1 for imp in improvements if imp > 0)}/{len(improvements)}")
        print(f"  Best Improvement: {np.max(improvements):+.3f}")
        print(f"  Worst Change: {np.min(improvements):+.3f}")

def print_improvement_analysis(comparison_df):
    """Analyze which signals improved the most"""
    print("\n" + "=" * 100)
    print("IMPROVEMENT ANALYSIS")
    print("=" * 100)
    
    # Calculate overall improvement score for each signal
    signal_scores = {}
    
    for signal in comparison_df['Signal'].unique():
        signal_data = comparison_df[comparison_df['Signal'] == signal]
        
        # Weight different metrics differently based on your priorities
        # Sensitivity, POD, Timeliness are equally important
        # Specificity is a constraint (should be >= 0.98)
        priority_metrics = ['Sensitivity_All', 'POD_All', 'Timeliness_All', 
                           'Sensitivity_NonBH', 'POD_NonBH', 'Timeliness_NonBH']
        
        score = 0
        count = 0
        
        for metric in priority_metrics:
            metric_data = signal_data[signal_data['Metric'] == metric]
            if not metric_data.empty:
                score += metric_data['Improvement'].iloc[0]
                count += 1
        
        if count > 0:
            signal_scores[signal] = score / count
        else:
            signal_scores[signal] = 0
    
    # Sort by improvement
    sorted_signals = sorted(signal_scores.items(), key=lambda x: x[1], reverse=True)
    
    print("\nSignals ranked by overall improvement:")
    print(f"{'Signal':<8} {'Improvement Score':<18} {'Key Changes'}")
    print("-" * 60)
    
    for signal, score in sorted_signals:
        signal_data = comparison_df[comparison_df['Signal'] == signal]
        changes = []
        
        for _, row in signal_data.iterrows():
            if row['Improvement'] > 0.05:
                changes.append(f"↑{row['Metric']}")
            elif row['Improvement'] < -0.05:
                changes.append(f"↓{row['Metric']}")
        
        print(f"{signal:<8} {score:<+18.3f} {', '.join(changes[:3]) if changes else 'Minimal'}")

def print_metric_breakdown(comparison_df):
    """Show breakdown by metric type"""
    print("\n" + "=" * 100)
    print("METRIC BREAKDOWN")
    print("=" * 100)
    
    # Group metrics by type
    metric_groups = {
        'Sensitivity': [col for col in comparison_df['Metric'].unique() if 'Sensitivity' in col],
        'Specificity': [col for col in comparison_df['Metric'].unique() if 'Specificity' in col],
        'POD': [col for col in comparison_df['Metric'].unique() if 'POD' in col],
        'Timeliness': [col for col in comparison_df['Metric'].unique() if 'Timeliness' in col],
        'Threshold': [col for col in comparison_df['Metric'].unique() if 'Threshold' in col]
    }
    
    for group_name, metrics in metric_groups.items():
        if metrics:
            print(f"\n{group_name.upper()}:")
            group_data = comparison_df[comparison_df['Metric'].isin(metrics)]
            
            for metric in metrics:
                metric_data = group_data[group_data['Metric'] == metric]
                avg_improvement = metric_data['Improvement'].mean()
                improved_count = sum(1 for imp in metric_data['Improvement'] if imp > 0)
                total_count = len(metric_data)
                
                print(f"  {metric:<25}: {avg_improvement:+.3f} avg improvement, "
                      f"{improved_count}/{total_count} signals improved")

def main():
    print("LSTM Results Comparison Tool")
    print("=" * 50)
    
    # Load the two result files
    print("\nLoading results...")
    old_df = load_results("results/LSTM_first_combined.csv")  # Assuming you have this
    new_df = load_results("LSTM_results_combined.csv")
    
    if old_df is None:
        print("\nTo use this tool:")
        print("1. Create a CSV file with your old results in the same format as LSTM_results_combined.csv")
        print("2. Save it as 'results/LSTM_first_combined.csv'")
        print("3. Run this script again")
        return
    
    # Compare results
    comparison_df = compare_results(old_df, new_df)
    
    if comparison_df is not None and not comparison_df.empty:
        print_detailed_comparison(comparison_df)
        print_summary_comparison(comparison_df)
        print_improvement_analysis(comparison_df)
        print_metric_breakdown(comparison_df)
        
        # Save comparison to CSV
        comparison_df.to_csv("LSTM_comparison_detailed.csv", index=False)
        print(f"\nDetailed comparison saved to: LSTM_comparison_detailed.csv")
        
        # Create a summary pivot table
        pivot_df = comparison_df.pivot(index='Signal', columns='Metric', values='Improvement')
        pivot_df.to_csv("LSTM_comparison_pivot.csv")
        print(f"Pivot table saved to: LSTM_comparison_pivot.csv")
    
    else:
        print("No valid comparison data found")

if __name__ == "__main__":
    main()
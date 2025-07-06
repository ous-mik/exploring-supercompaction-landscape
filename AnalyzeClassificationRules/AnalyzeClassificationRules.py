import re
import csv
import pandas as pd

def analyze_rules(rules_file_path, output_base_name="CompactionPhenotype_feature_importance"):
    '''
    Analyzes FastGentleBoosting rules (exported from CellProfiler Analyst) from a text file to extract feature importance.

    Arguments:
        rules_file_path (str): Path to the text file containing the rules.
        output_base_name (str): Base name for the output CSV/XLSX files.
    '''
    parsed_rules = []
    
    with open(rules_file_path, 'r') as file:
        for line in file:
            line = line.strip()
            if line.startswith("IF"):
                # Regex to capture feature, threshold, and scores for two or three classes
                match = re.match(r'IF \((.+?) > (.+?), \[(.+?)\], \[(.+?)\]\)', line, flags=re.IGNORECASE)
                if match:
                    feature_name, threshold, scores_true_str, scores_false_str = match.groups()
                    threshold = float(threshold)
                    scores_true = list(map(float, scores_true_str.split(',')))
                    scores_false = list(map(float, scores_false_str.split(',')))
                    num_classes = len(scores_true)
                    score_changes = [scores_true[i] - scores_false[i] for i in range(num_classes)]
                    parsed_rules.append({
                        'feature_name': feature_name,
                        'threshold': threshold,
                        'score_changes': score_changes
                    })
                else:
                    print(f"Warning: Could not parse rule: {line}")

    feature_importance = {}
    num_classes = 0
    if parsed_rules:
        num_classes = len(parsed_rules[0]['score_changes'])

    total_abs_importance = [0.0] * num_classes

    for rule in parsed_rules:
        feature_name = rule['feature_name']
        score_changes = rule['score_changes']
        threshold = rule['threshold']

        if feature_name not in feature_importance:
            feature_importance[feature_name] = {'score_changes': [0.0] * num_classes, 'thresholds': []}

        for i in range(len(score_changes)):
            feature_importance[feature_name]['score_changes'][i] += abs(score_changes[i])
            total_abs_importance[i] += abs(score_changes[i])
        
        feature_importance[feature_name]['thresholds'].append(threshold)

    # Prepare plotting data with normalized importance
    plotting_data = [["Feature Name"] + ["Threshold"] + [f"Class {i+1}" for i in range(num_classes)] + [f"Normalized Absolute Importance Class {i+1}" for i in range(num_classes)]]
    for rule in parsed_rules:
        feature = rule['feature_name']
        threshold = rule['threshold']
        score_changes = rule['score_changes']
        normalized_scores = [abs(score) / total_abs_importance[score_changes.index(score)] if total_abs_importance[score_changes.index(score)]!= 0 else 0 for score in score_changes]
        plotting_data.append([feature] + [threshold] + score_changes + normalized_scores)

    # Save plotting data to CSV
    plotting_csv_file = f"{output_base_name}_plotting.csv"
    with open(plotting_csv_file, 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerows(plotting_data)
        print(f"Plotting data saved to {plotting_csv_file}")

    try:
        df_plotting = pd.DataFrame(plotting_data[1:], columns=plotting_data[0])
        plotting_excel_file = f"{output_base_name}_plotting.xlsx"
        df_plotting.to_excel(plotting_excel_file, index=False)
        print(f"Plotting data saved to {plotting_excel_file}")
    except ImportError:
        print("Pandas not installed. Skipping XLSX export for plotting data.")

    # Prepare data for summary table
    summary_data = [["Feature Name"] + [f"Aggregate Importance Class {i+1}" for i in range(num_classes)] + [f"Normalized Absolute Aggregate Importance Class {i+1}" for i in range(num_classes)] + ["Number of Rules", "Min Threshold", "Max Threshold"]]
    for feature, data in feature_importance.items():
        aggregate_scores = data['score_changes']
        thresholds = data['thresholds']
        num_rules = len(thresholds)
        min_threshold = min(thresholds) if thresholds else None
        max_threshold = max(thresholds) if thresholds else None
        normalized_importance = [abs(score) / total_abs_importance[aggregate_scores.index(score)] if total_abs_importance[aggregate_scores.index(score)]!= 0 else 0 for score in aggregate_scores]
        summary_data.append([feature] + aggregate_scores + normalized_importance + [num_rules, min_threshold, max_threshold])

    # Save summary data to CSV
    summary_csv_file = f"{output_base_name}_summary.csv"
    with open(summary_csv_file, 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerows(summary_data)
        print(f"Summary table data saved to {summary_csv_file}")

    try:
        df_summary = pd.DataFrame(summary_data[1:], columns=summary_data[0])
        summary_excel_file = f"{output_base_name}_summary.xlsx"
        df_summary.to_excel(summary_excel_file, index=False)
        print(f"Summary table data saved to {summary_excel_file}")
    except ImportError:
        print("Pandas not installed. Skipping XLSX export for summary data.")

if __name__ == "__main__":
    rules_file = "CompactionPhenotypeFilter_MainAnalysis.txt" # Replace with the actual path to your rules file
    analyze_rules(rules_file)
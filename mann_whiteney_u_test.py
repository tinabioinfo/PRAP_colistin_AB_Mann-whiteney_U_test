import pandas as pd
import scipy.stats as stats
import statsmodels.stats.multitest as smm
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np

# Load the data
file_path = '/home/tinas/download_AMR_data/colistin_analysis/AB/prokka_annotation/all_gff_files/output_roary_407_with_protein_sequences/PRAP_statistical_test/PRAP_for_statistical_test.csv'
df = pd.read_csv(file_path, index_col=0)

# Extract phenotypes and clean the data
phenotypes = df.iloc[-1, :]
df = df.drop(df.tail(1).index)

# Split the data into resistant and susceptible groups
resistant_indices = phenotypes[phenotypes == 1].index
susceptible_indices = phenotypes[phenotypes == 0].index

resistant_df = df[resistant_indices].astype(float)
susceptible_df = df[susceptible_indices].astype(float)

# Initialize lists to store results
p_values = []
effect_sizes = []

# Perform Mann-Whitney U tests and calculate effect sizes for each gene
for gene in df.index:
    resistant_vals = resistant_df.loc[gene].dropna()
    susceptible_vals = susceptible_df.loc[gene].dropna()
    u_stat, p_val = stats.mannwhitneyu(resistant_vals, susceptible_vals, alternative='two-sided')
    p_values.append(p_val)
    # Calculate effect size as mean difference
    effect_size = resistant_vals.mean() - susceptible_vals.mean()
    effect_sizes.append(effect_size)

# Adjust p-values for multiple testing
adjusted_p_values = smm.multipletests(p_values, alpha=0.05, method='fdr_bh')[1]

# Create a DataFrame to store test results
results_df = pd.DataFrame({
    'Gene': df.index,
    'p_value': p_values,
    'adjusted_p_value': adjusted_p_values,
    'effect_size': effect_sizes
})

# Save the results to a CSV file
results_file_path = '/home/tinas/download_AMR_data/colistin_analysis/AB/prokka_annotation/all_gff_files/output_roary_407_with_protein_sequences/PRAP_statistical_test/out_code.csv'
results_df.to_csv(results_file_path, index=False)

# Visualizations
# Histogram of p-values
plt.figure(figsize=(10, 6))
sns.histplot(p_values, bins=30, kde=False)
plt.title('Histogram of p-values')
plt.xlabel('p-value')
plt.ylabel('Frequency')
plt.savefig('/home/tinas/download_AMR_data/colistin_analysis/AB/prokka_annotation/all_gff_files/output_roary_407_with_protein_sequences/PRAP_statistical_test/p_value_histogram.png')

# Preparing data for the volcano plot
results_df['-log10(p_value)'] = -np.log10(results_df['p_value'])

# Volcano plot of statistical results
plt.figure(figsize=(10, 6))
sns.scatterplot(x='effect_size', y='-log10(p_value)', data=results_df, hue=results_df['adjusted_p_value'] < 0.05)
plt.title('Volcano Plot of Gene Differences')
plt.xlabel('Effect Size (Mean Difference)')
plt.ylabel('-log10(p-value)')
plt.axhline(-np.log10(0.05), color='red', linestyle='--')  # Significance threshold line
plt.savefig('/home/tinas/download_AMR_data/colistin_analysis/AB/prokka_annotation/all_gff_files/output_roary_407_with_protein_sequences/PRAP_statistical_test/volcano_plot.png')

# Box plots for a selection of significant genes
significant_genes = results_df[results_df['adjusted_p_value'] < 0.05]['Gene']
for gene in significant_genes.head(5):  # Limit to top 5 significant genes for demonstration
    plt.figure(figsize=(10, 6))
    sns.boxplot(x='variable', y='value', data=pd.melt(df.loc[[gene]]), palette='Set2')
    plt.title(f'Protein Identity Distribution: {gene}')
    plt.xticks([0, 1], ['Resistant', 'Susceptible'])
    plt.ylabel('Protein Identity')
    plt.savefig(f'/home/tinas/download_AMR_data/colistin_analysis/AB/prokka_annotation/all_gff_files/output_roary_407_with_protein_sequences/PRAP_statistical_test/{gene}_boxplot.png')


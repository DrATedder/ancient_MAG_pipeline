mport matplotlib.pyplot as plt
import pandas as pd

# Read data from CSV file
file_path = '~/combined_MAG_data.csv'  # Replace 'your_data_file.csv' with the actual file path
df = pd.read_csv(file_path)

# Sort the DataFrame by 'abundance'
df_sorted = df.sort_values(by='mean_abundance')

# Create a figure with two subplots
fig, axes = plt.subplots(nrows=1, ncols=2, figsize=(12, 4))

# Plot 1: 'presence' based on 'abundance'
colors = ['red' if p == 0 else 'blue' for p in df_sorted['presence']]
axes[0].scatter(df_sorted['mean_abundance'], range(len(df_sorted)), c=colors, marker='o')
axes[0].set_xlabel('Abundance')
axes[0].set_ylabel('Sample Index')
axes[0].set_title('MAG presence (Blue) and absence (Red) based on sorted abundance')

# Plot 2: 'presence' based on 'abundance' and 'genome_size'
colors = ['red' if p == 0 else 'blue' for p in df['presence']]
scatter = axes[1].scatter(df['mean_abundance'], df['genome_size'], c=colors, marker='o')
axes[1].set_xlabel('Abundance')
axes[1].set_ylabel('Genome Size')
axes[1].set_title('Presence (Red) and Absence (Blue) based on Abundance and Genome Size')
axes[1].legend(handles=[plt.Line2D([0], [0], marker='o', color='w', markerfacecolor='red', markersize=10, label='Absence'),
                        plt.Line2D([0], [0], marker='o', color='w', markerfacecolor='blue', markersize=10, label='Presence')])

# Adjust layout
plt.tight_layout()

# Save the figure as a single PNG file
plt.savefig('combined_plot_ancient.png')

# Show the plot (optional)
plt.show()

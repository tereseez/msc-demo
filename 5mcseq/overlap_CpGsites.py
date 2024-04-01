#!/usr/bin/env python
#start overlap_CpGsites.py
import pandas as pd
from matplotlib import pyplot as plt
from matplotlib_venn import venn2
import glob
from scipy.stats import pearsonr

def add_numerical_index2file(file_names, split_str):
    """
    Input: list of files with string separator in filenames
    Output: df of filenames associated with a numerical value (chr)
    """
    a_df = pd.DataFrame(data = file_names, columns = ['file_name'])
    a_df2 = a_df.file_name.str.split('/', expand = True)
    new_a = pd.concat([a_df, a_df2.iloc[:,-1]], axis = 1).copy()
    new_a['numerical'] = new_a.iloc[:,-1].apply(
                            lambda x: x.split(split_str)[1]
                            .replace('chr','').replace('X','21')
                            .replace('Y','22').replace('M','23')
                            .replace('T','')
                            )
    new_a.numerical = new_a.numerical.astype(int)
    return new_a.copy()

run1_files = glob.glob('out_met_count/A_2/CpG_context/*')
run2_files = glob.glob('out_met_count/run2_junbai/\
        A_2/me_counts/CpG_context/*')

#combine filename dataframes into one shared dataframe
a_run1 = add_numerical_index2file(run1_files,'_')
a_run2 = add_numerical_index2file(run2_files,'.')
combined_file_names = pd.merge(a_run1, a_run2, on = 'numerical',
                               how = 'left').copy()

#sort combined dataframe by numerical column (chr)
sorted_a = combined_file_names.sort_values(by = ['numerical'])
sorted_a.reset_index(inplace = True, drop = True)

#Empty sets for venn diagram
run1_CpG = set()
run2_CpG = set()

#Empty dataframe to track all chr
all_chromosomes = pd.DataFrame()

#compare CpG sites chromosome wise
for idx, row in sorted_a.iterrows():
    run1_file = row['file_name_x']
    run2_file = row['file_name_y']
    run1_df = pd.read_csv(run1_file, sep = '\t', header = None,
                          names = ["chrom", "start", "end", "percentage",
                                   "total_counts", "strand"])
    run2_df = pd.read_csv(run2_file, sep = '\t', header = None,
                          names = ["chrom", "start", "end", "percentage",
                                   "total_counts", "strand"])
    run1_CpG.update(set(run1_df.iloc[:,1]))
    run2_CpG.update(set(run2_df.iloc[:,1]))
    combined_datasets = pd.merge(run1_df, run2_df, on = 'start',
                                 how = 'inner')
    all_chromosomes = pd.concat([all_chromosomes, combined_datasets],
                                axis = 0)

all_chromosomes.reset_index(inplace = True, drop = True)

#Make venn diagram
total = len(run1_CpG)
v = venn2([run1_CpG, run2_CpG], set_labels=('Run 1 ', 'Run 2'),
      set_colors=('#ed7d31', 'darksalmon'), alpha = 0.5,
      subset_label_formatter = lambda x: str(x) +
      "\n(" + f"{(x/total):1.0%}" + ")")
v.get_label_by_id('01').set_text('')
v.get_label_by_id('10').set_text('')
plt.title('A2 comparison: CpG sites')
plt.show()

#Calculate correlation coefficient
r, p = pearsonr(all_chromosomes['total_counts_x'],
                all_chromosomes['total_counts_y'])

#scatter plot of total counts per CpG site
f, ax = plt.subplots(figsize=(7,5))
plt.scatter(all_chromosomes['total_counts_x'],
            all_chromosomes['total_counts_y'],
            color="#ed7d31", alpha=0.5)
plt.xlabel('Run 1')
plt.ylabel('Run 2')
plt.title('A2 comparison: total count per CpG site')
ax.plot([0, 1], [0, 1], transform=ax.transAxes, linestyle='--',
        color="darksalmon")
plt.show()
#end overlap_CpGsites.py

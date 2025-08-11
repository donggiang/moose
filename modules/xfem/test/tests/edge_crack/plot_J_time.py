import pandas as pd
import matplotlib.pyplot as plt


# plot

plt.figure(figsize=(10,6))


# # csv file name
# csv_path = 'edge_crack_2d_en_ref.csv'
# # read csv file
# df = pd.read_csv(csv_path, skiprows=[1])
# x2= df['time']#
# y_columns = ['J_1_1']
# for y_col in y_columns:
#     plt.plot(x2, df[y_col], linestyle='-',marker='o', linewidth=2, label=f'enrichment', color='b')


# csv file name
csv_path = 'ADedge_crack_2d_en_inc_ri3_ro5.csv'
# read csv file
df = pd.read_csv(csv_path, skiprows=[1])
x2= df['time']#
y_columns = ['J_1_1']
for y_col in y_columns:
    plt.plot(x2, df[y_col], linestyle='-',marker='o', linewidth=2, label=f'enrichment', color='r')



# csv file name
csv_path = 'edge_crack_2d_en_inc_ri3_ro5.csv'
# read csv file
df = pd.read_csv(csv_path, skiprows=[1])
x2= df['time']#
y_columns = ['J_1_1']
for y_col in y_columns:
    plt.plot(x2, df[y_col], linestyle='-',marker='none', linewidth=2, label=f'enrichment', color='b')


# ad labels and title
plt.xlabel('Time', fontsize=14)
plt.ylabel('J', fontsize=14)
plt.legend(fontsize=14, frameon=False)
#plt.xlim(0, 1200)
#plt.ylim(0, 3)
#tick customization
# Customize the outer frame (spines)
ax = plt.gca()  # Get current axis
spine_width = 1.5  # Specify the spine width
ax.spines['top'].set_linewidth(spine_width)
ax.spines['bottom'].set_linewidth(spine_width)
ax.spines['left'].set_linewidth(spine_width)
ax.spines['right'].set_linewidth(spine_width)

# Ticks customization
ax.tick_params(axis='both', which='major', labelsize=12)

#save the plot
out_path = 'J_vs_time.png'
plt.savefig(out_path)
plt.show()

print('Plotting is completed successfully!!!')




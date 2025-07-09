import pandas as pd
import matplotlib.pyplot as plt

def compute_crack_length(c, n, B, time):
    length = [0]  # Initialize the length array with the first element as 0
    print(time)
    for i in range(0, len(time)-1):
        if i>0:
            dt = time[i] - time[i-1]
        else:
            dt = time[i]
        if B[i]>0:
           #print('c value: ')
           #print (B[i])
           #if i==0:
            #B[i] = B[i]/4
           if(i==0):
            B[i] = B[i]/100000
           increment = c * pow(B[i], n) * dt

        else:
            increment = 0
       # print('increment')
       # print(increment)

      #  if i>0:
        length.append(length[-1] + increment)
        # else :
        #     length.append(0)
        print(length[i])
    return length


plt.figure(figsize=(10,6))

y_columns = ['C_1_1']


c= 0.6e-2
e =0.9 #1.1 #0.87 #0.867#factor/(factor+1)
############
# csv_path = 'CCG_en_cr_cl_h0p2_dt2_rd0p59_ri2_ro6p5_0p18em3_e1p35_c6em24_n6p8.csv'
# # read csv file
# df3 = pd.read_csv(csv_path)
# x3= df3['time']#
# for y_col in y_columns:
#     plt.plot(x3, df3[y_col], linestyle='-',marker='o', linewidth=3, label=f'h=0.2', color='r')


# csv_path = 'CCG_en_cr_cl_h0p16_dt2_rd0p59_ri2_ro6p5_0p18em3_e1p35_c6em24_n6p8.csv'
# # read csv file
# df2 = pd.read_csv(csv_path)
# x2= df2['time']#
# for y_col in y_columns:
#     plt.plot(x2, df2[y_col], linestyle='-',marker='o', linewidth=3, label=f'h=0.16', color='g')

# csv_path = 'CCG_en_cr_cl_h0p125_dt5_rd0p59_ri2_ro6p5_0p18em3_e1p35_c6em24_n6p8.csv'
# # read csv file
# df = pd.read_csv(csv_path)
# x= df['time']#
# for y_col in y_columns:
#     plt.plot(x, df[y_col], linestyle='-',marker='o', linewidth=3, label=f'h=0.125', color='b')
############
# #########
csv_path = 'CCG_en_cr_cl_h0p5_d0p1_m0p7_n0p9_c5em24_e6p3.csv'
# read csv file
df4 = pd.read_csv(csv_path)
x4= df4['time']#
for y_col in y_columns:
    plt.plot(x4, df4[y_col], linestyle='-',marker='none', linewidth=4, label=f'h=0.16', color='purple')



#print (df[y_col])
# ad labels and title
plt.xlabel('Time', fontsize=14)
plt.ylabel('C(t)', fontsize=14)
plt.legend(fontsize=14, frameon=False)
#plt.xlim(0, 140)
plt.ylim(-5, 5)
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
out_path = 'Ct_vs_time.png'
plt.savefig(out_path)
plt.show()

print('Plotting C(t) is completed successfully!!!')


plt.figure(figsize=(10,6))



csv_path = 'vm_experiment.csv'
df_exp = pd.read_csv(csv_path, header=None)
plt.plot(df_exp[0], df_exp[1], linestyle='--',   linewidth=4, label='WM-Zhao\'s experiment', color='orange')



length = compute_crack_length(c, e, df4[y_col], x4)
print (length)
print('complete print length')
for y_col in y_columns:
    plt.plot(x4, length, linestyle='-',marker='none', linewidth=3, label=f'h=0.16', color='purple')



# ad labels and title
plt.xlabel('Time', fontsize=14)
plt.ylabel('crack length)', fontsize=14)
plt.legend(fontsize=14, frameon=False)
#plt.xlim(0, 1200)
#plt.ylim(0, 5)
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
out_path = 'cracklength_vs_time.png'
plt.savefig(out_path)
plt.show()

print('Plotting is completed successfully!!!')



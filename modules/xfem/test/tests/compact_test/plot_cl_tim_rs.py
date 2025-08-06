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
          # if(i==0):
           # B[i] = B[i]/100000
  
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


c= 0.15e-3
e =0.87 #1.1 #0.87 #0.867#factor/(factor+1)




############
# #########
csv_path = 'nonAD_CCG_en_cr_cl_h0p16_d1p6_m0p25em3_n0p87_c2em23_e6p83_n1_ri3p2_ro6p4.csv'
# read csv file
df4 = pd.read_csv(csv_path)

# Loop over each row, starting from the second row
for i in range(1, len(df4["C_1_1"])):
    current_value = df4["C_1_1"][i]
    previous_value = df4["C_1_1"][i-1]

    if i>0:
        # Check if current value is 10 times greater than the previous value or negative
        if current_value > 2 * previous_value or current_value < 0:
            # Modify the current value as needed
            # For example, set it to the previous value or some other logic
            df4["C_1_1"][i] = previous_value

x4= df4['time']#
for y_col in y_columns:
    plt.plot(x4, df4[y_col], linestyle='-',marker='none', linewidth=4, label=f'h=0.1', color='purple')



############
# #########
csv_path = 'nonAD_CCG_en_cr_cl_h0p125_d1p6_m0p15em3_n0p87_c2em23_e6p83_n_ri2p5_ro5.csv'
# read csv file
df2 = pd.read_csv(csv_path)

# Loop over each row, starting from the second row
for i in range(1, len(df2["C_1_1"])):
    current_value = df2["C_1_1"][i]
    previous_value = df2["C_1_1"][i-1]

    if i>0:
        # Check if current value is 10 times greater than the previous value or negative
        if current_value > 2 * previous_value or current_value < 0:
            # Modify the current value as needed
            # For example, set it to the previous value or some other logic
            df2["C_1_1"][i] = previous_value

x2= df2['time']#
for y_col in y_columns:
    plt.plot(x2, df2[y_col], linestyle='-',marker='none', linewidth=4, label=f'h=0.1', color='yellow')

csv_path = 'nonAD_CCG_en_cr_cl_h0p1_d0p8_m0p15em3_n0p87_c2em23_e6p83_n_ri2_ro4.csv'
# read csv file
df3 = pd.read_csv(csv_path)

# Loop over each row, starting from the second row
for i in range(1, len(df3["C_1_1"])):
    current_value = df3["C_1_1"][i]
    previous_value = df3["C_1_1"][i-1]

    if i>0:
        # Check if current value is 10 times greater than the previous value or negative
        if current_value > 2 * previous_value or current_value < 0:
            # Modify the current value as needed
            # For example, set it to the previous value or some other logic
            df3["C_1_1"][i] = previous_value

x3= df3['time']#
for y_col in y_columns:
    plt.plot(x3, df3[y_col], linestyle='-',marker='none', linewidth=4, label=f'h=0.16', color='red')


#print (df[y_col])
# ad labels and title
plt.xlabel('Time', fontsize=14)
plt.ylabel('C(t)', fontsize=14)
plt.legend(fontsize=14, frameon=False)
#plt.xlim(0, 140)
#plt.ylim(-5, 5)
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


length = compute_crack_length(c, e, df2[y_col], x2)
print (length)
print('complete print length')
for y_col in y_columns:
    plt.plot(x2, length, linestyle='-',marker='none', linewidth=3, label=f'h=0.125', color='yellow')

length = compute_crack_length(c, e, df3[y_col], x3)
print (length)
print('complete print length')
for y_col in y_columns:
    plt.plot(x3, length, linestyle='-',marker='none', linewidth=3, label=f'h=0.1', color='red')

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



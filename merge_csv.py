import os

path = './rta_test/csv/IRF_North_z60_0.5h/'
# merge files ---!
csv_file = os.listdir(path=path)
csv_file = sorted(csv_file)
print(len(csv_file))
for idx, file in enumerate(csv_file):
    if '.log' in file:
        csv_file.pop(idx)
print(len(csv_file))

csv_merged = 'sensitivity.csv'
fout = open(path + csv_merged, 'w+')
# first file ---!
for line in open(path + csv_file[0]):
    fout.write(line)
# remaining files ---!    
for idx, file in enumerate(csv_file):
    f = open(path + file)
    next(f)  # skip the header ---!
    for line in f:
        fout.write(line)
    f.close()
fout.close()

print(csv_merged)
print('done')
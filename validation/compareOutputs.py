import ROOT as r
import sys
import os

print(">>> STARTING THE VALIDATION...")
print("\n")

if len(sys.argv) < 3:
    print("> Not enough inputs to do the comparison")
    print("> They should be at least 2")
    print("> Please, run 'python compareOutputs.py output1 output2'")
    print("> Exiting...")
    print("\n")
    exit()


## Read the inputs

input1 = sys.argv[1]
if input1[-1] != '/':
    input1 += '/'
input2 = sys.argv[2]
if input2[-1] != '/':
    input2 += '/'
print(">> Comparing 2 inputs:")
print(">>   - " + input1)
print(">>   - " + input2)
print("\n")



## Loop over files:

files1 = []
for f in os.listdir(input1):
    if '.root' not in f:
        continue
    files1.append(f)

files2 = []
for f in os.listdir(input2):
    if '.root' not in f:
        continue
    files2.append(f)

if files1 != files2:
    print("> The provided inputs do not contain the same ROOT files")
    print("> Exiting...")
    print("\n")
    exit()

checks = [] # fail: 0 ; success: 1; cannot compare: -1

for filename in files1:
    print(">> Validating " + filename + '...')
    tfile1 = r.TFile(input1 + filename)
    tfile2 = r.TFile(input2 + filename)

    histos1 = [th1f.GetName() for th1f in tfile1.GetListOfKeys()]
    histos2 = [th1f.GetName() for th1f in tfile2.GetListOfKeys()]

    if histos1 != histos2:
        print(" > The histograms are not the same")
        print(" > Exiting...")
        exit()

    for histo in histos1:
        h1 = tfile1.Get(histo)
        h2 = tfile2.Get(histo)
        if h1.GetEntries() < 1 and h2.GetEntries() < 1:
            checks.append(-1)
            continue
        if h2.GetMean() < 1e-15:
            checks.append(-1)
            continue
        if h2.GetRMS() < 1e-15:
            checks.append(-1)
            continue

        if h1.GetEntries() != h2.GetEntries() or h1.GetMean() != h2.GetMean() or h1.GetRMS() != h2.GetRMS():
            checks.append(0)
        else:
            checks.append(1)

## Print results:

print(">>> VALIDATION RESULTS:")
print(">>  Compared %d distributions from which" % len(checks))
print(">>    - %d are indetermined" % checks.count(-1))
print(">>    - %d successed" % checks.count(1))
print(">>    - %d failed" % checks.count(0))




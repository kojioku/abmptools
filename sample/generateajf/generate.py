import sys

argvs = sys.argv
inpdb = argvs[1]

# atom section
'''
    --mp2d,  --nbo,  --esp
    --mp2d, --nbo
    --mp2d, --esp
    --mp2d
    --hf, --nbo --esp
    --hf, --nbo
    --hf, --esp
    --hf
'''

istrs = []

methods = ['MP2D', 'HF']
esps = ['--resp', '']
nbos = ['--nonbo', '']
for method in methods:
    for esp in esps:
        for nbo in nbos:
            istr = '--method ' + method + ' ' + esp + ' ' + nbo
            istrs.append(istr)

#  monmer section
'''
    --mp2
    --hf
    --mp3
    --lrd
'''

methods = ['HF', 'MP2', 'MP3', 'HF+D']
for method in methods:
    istr = '--method ' + method
    istrs.append(istr)


piedas = ['--nopieda', '']
solvs = ['-pb', '']
bsses = ['-bsse', '']
methods = ['HF', 'MP2', 'MP3', 'HF+D']

for pieda in piedas:
    for solv in solvs:
        for bsse in bsses:
            for method in methods:
                istr = pieda + ' ' + solv + ' ' + bsse + ' --method ' + method
                istrs.append(istr)

# print(istrs)


f = open('generate.sh', 'w')
for aaa in istrs:
    print('python -m abmptools.generateajf.py -i', inpdb, aaa, file=f)

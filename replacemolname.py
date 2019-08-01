import sys
import os

## --- user input ---
replacelists = [['001', 'LIG'], ['002', 'CYC'], ['003', 'NA'], ['004', 'WAT']]

# replacelists = [['*', '1', 'BEX']]
# infos = [ ['*', '1', 'LIG'], ['*', '2', 'CYC'], ['*', '3', 'AAA']]


## --- end ---

argvs = sys.argv
print(argvs)

for i in range(len(argvs)):
    if i == 0:
        continue
    infile = argvs[i]
    head, ext = os.path.splitext(infile)
    print(head, ext)

    if ext != '.pdb':
        out = head.split('.pdb')[0] + ext + '-sed.pdb'
    else:
        out = head + '-sed.pdb'
    print(out)

    lines = open(infile, 'r').readlines()
    # print(lines)


    outf = open(out, 'w')
    # for repinfo in replacelists:
        # print('{0:>3}'.format(repinfo[0]) + '     ' + repinfo[1] + ' ')
    for line in lines:
        for repinfo in replacelists:
            if len(repinfo) == 3:
                line = line.replace(' {0:>3}'.format(repinfo[0]) + '     ' + repinfo[1] + ' ', '{0:>3}'.format(repinfo[2]) + '     ' + repinfo[1] + ' ')
            if len(repinfo) == 2:
                line = line.replace(' {0:>3}'.format(repinfo[0]) + ' ', ' {0:>3}'.format(repinfo[1]) + ' ')

        print(line[:-1], file=outf)

    print(out, 'was created.')
# sed 's/  \*     1 /LIG     1 /g' $1 |sed 's/  \*     2 /CYC     2 /g' > $out
# echo $out was created.


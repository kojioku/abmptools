#!/bin/bash
pdb=$1
basis=${2:-STO-3G}
head=${pdb%.*}

python -m abmptools.generateajf -i $pdb -basis "$basis" --method MP2D --resp --nonbo
python -m abmptools.generateajf -i $pdb -basis "$basis" --method MP2D --resp
python -m abmptools.generateajf -i $pdb -basis "$basis" --method MP2D  --nonbo
python -m abmptools.generateajf -i $pdb -basis "$basis" --method MP2D
python -m abmptools.generateajf -i $pdb -basis "$basis" --method HF --resp --nonbo
python -m abmptools.generateajf -i $pdb -basis "$basis" --method HF --resp
python -m abmptools.generateajf -i $pdb -basis "$basis" --method HF  --nonbo
python -m abmptools.generateajf -i $pdb -basis "$basis" --method HF
python -m abmptools.generateajf -i $pdb -basis "$basis" --method HF
python -m abmptools.generateajf -i $pdb -basis "$basis" --method MP2
python -m abmptools.generateajf -i $pdb -basis "$basis" --method MP3
python -m abmptools.generateajf -i $pdb -basis "$basis" --method HF+D
python -m abmptools.generateajf -i $pdb -basis "$basis" --nopieda -pb -bsse --method HF
python -m abmptools.generateajf -i $pdb -basis "$basis" --nopieda -pb -bsse --method MP2
python -m abmptools.generateajf -i $pdb -basis "$basis" --nopieda -pb -bsse --method MP3
python -m abmptools.generateajf -i $pdb -basis "$basis" --nopieda -pb -bsse --method HF+D
python -m abmptools.generateajf -i $pdb -basis "$basis" --nopieda -pb  --method HF
python -m abmptools.generateajf -i $pdb -basis "$basis" --nopieda -pb  --method MP2
python -m abmptools.generateajf -i $pdb -basis "$basis" --nopieda -pb  --method MP3
python -m abmptools.generateajf -i $pdb -basis "$basis" --nopieda -pb  --method HF+D
python -m abmptools.generateajf -i $pdb -basis "$basis" --nopieda  -bsse --method HF
python -m abmptools.generateajf -i $pdb -basis "$basis" --nopieda  -bsse --method MP2
python -m abmptools.generateajf -i $pdb -basis "$basis" --nopieda  -bsse --method MP3
python -m abmptools.generateajf -i $pdb -basis "$basis" --nopieda  -bsse --method HF+D
python -m abmptools.generateajf -i $pdb -basis "$basis" --nopieda   --method HF
python -m abmptools.generateajf -i $pdb -basis "$basis" --nopieda   --method MP2
python -m abmptools.generateajf -i $pdb -basis "$basis" --nopieda   --method MP3
python -m abmptools.generateajf -i $pdb -basis "$basis" --nopieda   --method HF+D
python -m abmptools.generateajf -i $pdb -basis "$basis"  -pb -bsse --method HF
python -m abmptools.generateajf -i $pdb -basis "$basis"  -pb -bsse --method MP2
python -m abmptools.generateajf -i $pdb -basis "$basis"  -pb -bsse --method MP3
python -m abmptools.generateajf -i $pdb -basis "$basis"  -pb -bsse --method HF+D
python -m abmptools.generateajf -i $pdb -basis "$basis"  -pb  --method HF
python -m abmptools.generateajf -i $pdb -basis "$basis"  -pb  --method MP2
python -m abmptools.generateajf -i $pdb -basis "$basis"  -pb  --method MP3
python -m abmptools.generateajf -i $pdb -basis "$basis"  -pb  --method HF+D
python -m abmptools.generateajf -i $pdb -basis "$basis"   -bsse --method HF
python -m abmptools.generateajf -i $pdb -basis "$basis"   -bsse --method MP2
python -m abmptools.generateajf -i $pdb -basis "$basis"   -bsse --method MP3
python -m abmptools.generateajf -i $pdb -basis "$basis"   -bsse --method HF+D
python -m abmptools.generateajf -i $pdb -basis "$basis"    --method HF
python -m abmptools.generateajf -i $pdb -basis "$basis"    --method MP2
python -m abmptools.generateajf -i $pdb -basis "$basis"    --method MP3
python -m abmptools.generateajf -i $pdb -basis "$basis"    --method HF+D

mkdir ${head} 2>/dev/null
mv *${head}* $head
cp $head/$pdb .

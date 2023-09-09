#!/bin/bash
pdb=$1
head=${pdb%.*}

python -m abmptools.generateajf -i $pdb --method MP2D --resp --nonbo
python -m abmptools.generateajf -i $pdb --method MP2D --resp 
python -m abmptools.generateajf -i $pdb --method MP2D  --nonbo
python -m abmptools.generateajf -i $pdb --method MP2D  
python -m abmptools.generateajf -i $pdb --method HF --resp --nonbo
python -m abmptools.generateajf -i $pdb --method HF --resp 
python -m abmptools.generateajf -i $pdb --method HF  --nonbo
python -m abmptools.generateajf -i $pdb --method HF  
python -m abmptools.generateajf -i $pdb --method HF
python -m abmptools.generateajf -i $pdb --method MP2
python -m abmptools.generateajf -i $pdb --method MP3
python -m abmptools.generateajf -i $pdb --method HF+D
python -m abmptools.generateajf -i $pdb --nopieda -pb -bsse --method HF
python -m abmptools.generateajf -i $pdb --nopieda -pb -bsse --method MP2
python -m abmptools.generateajf -i $pdb --nopieda -pb -bsse --method MP3
python -m abmptools.generateajf -i $pdb --nopieda -pb -bsse --method HF+D
python -m abmptools.generateajf -i $pdb --nopieda -pb  --method HF
python -m abmptools.generateajf -i $pdb --nopieda -pb  --method MP2
python -m abmptools.generateajf -i $pdb --nopieda -pb  --method MP3
python -m abmptools.generateajf -i $pdb --nopieda -pb  --method HF+D
python -m abmptools.generateajf -i $pdb --nopieda  -bsse --method HF
python -m abmptools.generateajf -i $pdb --nopieda  -bsse --method MP2
python -m abmptools.generateajf -i $pdb --nopieda  -bsse --method MP3
python -m abmptools.generateajf -i $pdb --nopieda  -bsse --method HF+D
python -m abmptools.generateajf -i $pdb --nopieda   --method HF
python -m abmptools.generateajf -i $pdb --nopieda   --method MP2
python -m abmptools.generateajf -i $pdb --nopieda   --method MP3
python -m abmptools.generateajf -i $pdb --nopieda   --method HF+D
python -m abmptools.generateajf -i $pdb  -pb -bsse --method HF
python -m abmptools.generateajf -i $pdb  -pb -bsse --method MP2
python -m abmptools.generateajf -i $pdb  -pb -bsse --method MP3
python -m abmptools.generateajf -i $pdb  -pb -bsse --method HF+D
python -m abmptools.generateajf -i $pdb  -pb  --method HF
python -m abmptools.generateajf -i $pdb  -pb  --method MP2
python -m abmptools.generateajf -i $pdb  -pb  --method MP3
python -m abmptools.generateajf -i $pdb  -pb  --method HF+D
python -m abmptools.generateajf -i $pdb   -bsse --method HF
python -m abmptools.generateajf -i $pdb   -bsse --method MP2
python -m abmptools.generateajf -i $pdb   -bsse --method MP3
python -m abmptools.generateajf -i $pdb   -bsse --method HF+D
python -m abmptools.generateajf -i $pdb    --method HF
python -m abmptools.generateajf -i $pdb    --method MP2
python -m abmptools.generateajf -i $pdb    --method MP3
python -m abmptools.generateajf -i $pdb    --method HF+D

mkdir ${head} 2>/dev/null
mv *${head}* $head
cp $head/$pdb .

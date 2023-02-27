#!/bin/bash 
DIR="$(echo $PWD)"

rm input.ogo leap.ls packmol.in 2> /dev/null
echo -e "###OGOLEM###\n<GEOMETRY>\nNumberOfParticles=ZZZZ" | tee -a input.ogo > /dev/null
echo -e "tolerance 2.0\nfiletype pdb\noutput pack.pdb" | tee -a packmol.in > /dev/null
echo -e "source leaprc.gaff" | tee -a leap.ls > /dev/null



for i in $(seq 1 2 $#)
  do
    arg=${!i}
    echo $arg
    c=$(($i+1))
    FRC=false                # is true if frcmod file is not empty
    base=$(name $arg)
    echo $base
    babel -ixyz $base.xyz -opdb $base.pdb 2> /dev/null
    
    antechamber -i $base.pdb -fi pdb -o $base.mol2 -fo mol2 -c bcc -nc ${!c} -rn $base > /dev/null
    rm ANTECHAMBER* ATOMTYPE.INF sqm.*
    
    if [ ! -e $base.mol2 ]
      then
        echo >&2 "Can't create mol2 files"
        exit 1
    fi
    
    parmchk2 -i $base.mol2 -f mol2 -o $base.frcmod > /dev/null
    
    if [ $(wc -l < $base.frcmod) -le 15 ]
    then
      rm $base.frcmod
    else
      FRC=true
    fi
    
    
    babel -imol2 $base.mol2 -opdb $base.pdb 2> /dev/null
    
    
    while read -r line
    do
      if [ $(echo $line | awk '{print $1}') == "ATOM" ]
      then
        element=$(echo $line | awk '{print $3}')
        element=$(printf '%s' "$element" | sed 's/[0-9]//g')
        
        echo "${line:0:76}$element" | tee -a tmp.pdb > /dev/null
      else
        echo $line | tee -a tmp.pdb > /dev/null
      fi
    done < "$base.pdb"
    
    mv tmp.pdb $base.pdb
    
    if [ $FRC == true ]
    then
      echo -e "source leaprc.gaff\nloadamberparams $base.frcmod\n$base = loadmol2 $base.mol2\ncheck $base\nsaveoff $base $base.off\nquit" > tmp.ls
    else
      echo -e "source leaprc.gaff\n$base = loadmol2 $base.mol2\ncheck $base\nsaveoff $base $base.off\nquit" > tmp.ls
    fi
    
    tleap -f tmp.ls > /dev/null
    grep -i "Warning" leap.log
#    if [[ $(grep "Warning" leap.log) != '' ]]
#      then
#        cp leap.log $base.log
#    fi
    rm leap.log tmp.ls
    
    echo -e "<MOLECULE>\nMoleculeRepetitions=XXX$i\nMoleculePath=$base.xyz\n</MOLECULE>" | tee -a input.ogo > /dev/null
    echo -e "\nstructure $base.pdb\n  number XXX$i\n  resnumbers 2\n  inside box 0. 0. 0. 20. 20. 20.\nend structure" | tee -a packmol.in > /dev/null
    if [ $FRC == true ]
    then
      echo -e "loadamberparams $base.frcmod" | tee -a leap.ls > /dev/null
    fi
    echo "loadoff $base.off" | tee -a leap.ls > /dev/null
    
    i=$((i+1))
  done

echo -e "</GEOMETRY>\nNumberOfGlobIterations=AAAA\nPoolSize=BBBB\nLocOptAlgo=amber" | tee -a input.ogo > /dev/null
echo -e "pack = loadpdb pack.pdb\ncheck pack\nsaveamberparm pack input.prmtop inpcrd\nquit" | tee -a leap.ls > /dev/null


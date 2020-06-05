#!/bin/zsh
cd ..
mkdir -p output/matrices
mkdir build & cd build
cmake ../
make
for x in 1 2 3 4 5;
do
    for geo in "emboite" "non_emboite"
    do
        for type in 0 1 3 4 5;
        do
            echo ./src/Mtf_$geo $x $type;
            ./src/Mtf_$geo $x $type;
        done
    done
done

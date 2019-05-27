#!/bin/zsh
cd ..
mkdir meshes
mkdir build & cd build
cmake ../
make
for x in 1 2 3 4 5;
do
    echo ./src/Gen_maillage_emboite $x;
    ./src/Gen_maillage_emboite $x;
    echo ./src/Gen_non_maillage_emboite $x;
    ./src/Gen_maillage_non_emboite $x;
done

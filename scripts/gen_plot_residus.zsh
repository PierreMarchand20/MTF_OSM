#!/bin/zsh
mkdir -p ../output/figures
for x in 3;
do
    for geo in "emboite" "non_emboite"
    do
        for type in 0 1 3;
        do
            echo python3 plot_gmres_residus.py --geo $geo  --ni $x --type $type --save ../output/figures/ --show 0
            python3 plot_gmres_residus.py --geo $geo  --ni $x --type $type --save ../output/figures/ --show 0
        done
    done
done

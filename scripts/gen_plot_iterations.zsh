#!/bin/zsh
mkdir -p ../output/figures
for x in 5;
do
    for geo in "emboite" "non_emboite"
    do
        for type in 0 1 3;
        do
            echo python3 plot_gmres_iterations.py --geo $geo  --nbr $x --type $type --save ../output/figures/ --show 0
            python3 plot_gmres_iterations.py --geo $geo  --nbr $x --type $type --save ../output/figures/ --show 0
        done
    done
done

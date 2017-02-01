# Script to test the robustness of the robust formulation
Omega=$(awk 'BEGIN{for(i=0;i<=3.1;i+=0.1)print i}')
for n in $Omega
#for Omega in $(seq 0.1 0.1 0.2)
do
    for sigma in 1
    do
        ./bin/mmkp -f data/I01.txt -o $n -s $sigma
        cat robust.txt >> summaryRobust.csv
    done
done

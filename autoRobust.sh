# Script to test the robustness of the robust formulation
# Omega=$(awk 'BEGIN{for(i=0;i<=5.1;i+=0.2)print i}')
Omega=$(awk 'BEGIN{for(i=0.2;i<=5.1;i+=0.2)print i}')
for n in $Omega
#for Omega in $(seq 0.1 0.1 0.2)
do
    for sigma in 1
    do
        ./bin/mmkp -f ../data/I07.txt -o $n -s $sigma -c 0.25 -n 14 -z 0.8 -u 0.05
        cat robust.txt >> summaryRobust.csv
    done
done

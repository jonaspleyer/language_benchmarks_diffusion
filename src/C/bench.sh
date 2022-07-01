#!/usr/bin/zsh

N_warmup=3
N_runs=7
report_file="report.txt"

rm $report_file

############ GCC ############
# Do the warmup run for gcc
for i in $(seq 1 $N_warmup);
do
    /usr/bin/time -f "%C %E" ./diffusion_gcc;
done

# Now do the actual run
for i in $(seq 1 $N_runs);
do
    /usr/bin/time -f "%C %E" -ao $report_file ./diffusion_gcc;
done

############ CLANG ############
# Do the warmup run fog clang
for i in $(seq 1 $N_warmup);
do
    /usr/bin/time -f "%C %E" ./diffusion_clang;
done

# Now do the actual run
for i in $(seq 1 $N_runs);
do
    /usr/bin/time -f "%C %E" -ao $report_file ./diffusion_clang;
done


for i in $bay_is; do x_id=$(grep $i ../Simons_pilot_sample_info | awk '{print $7"-"$5"-"$6}'); num_g=$(cat
    ${i}*.bed.txt | awk '{if($2 < 15) print $0}' | wc -l); echo -e "${x_id}\t$num_g" >> baylor_num_genes_lt15.txt ; done

for i in $nyg_id; do x_id=$(grep ${i}, ../BloodVsSaliva_Descr.csv | awk '{split($0,a,","); print a[2]"-"a[4]}');
    num_g=$(cat ${i}-D1*.bed.txt | awk '{if($2 < 15) print $0}' | wc -l); echo -e "${x_id}\t$num_g" >>
    nyg_num_genes_lt15.txt ; done


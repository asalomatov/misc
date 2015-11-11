
cat *blood*.txt | awk '{print $1}' | sort | uniq -c | sort -nr | head -50 > nyg_blood_freq_genes_lt15.txt


for file in *.fastq.gz; do
    new_name=$(echo "$file" | sed 's/_[^_]*_/_S1_/')
    mv "$file" "$new_name"
done

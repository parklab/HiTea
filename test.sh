hitea -i 'examples/test.bam' -e 'MboI' -g 'hg38' -w test_out -o test

size=$(($(wc -c < "test_out/test.candidate.insertions.bed") +0))
if [[ $size < 10 ]]; then
  echo "ERROR : The run was unsuccessful"
  exit 1  
else
  mv test_out/test.candidate.insertions.bed examples/test.candidate.insertions.bed
  mv test_out/HiTEA_Report.html examples/HiTEA_Report_test.html
  rm -r -f test_out
  exit 0
fi
   
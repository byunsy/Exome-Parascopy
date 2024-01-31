# Exome-Parascopy



### Reference set
Building the reference set is an essential step for our main algorithm. This step might take some time depending on the number of input samples and number of exons.
```
python build_referencesets.py all.counts.EUR.exons.tsv > EUR.betafit.out
```

### True copy numbers 
We also need to get the true copy numbers of samples that will serve as ground-truth values when calculating accuracy. Here, we use `res.samples.bed` from Parascopy output.

By default, copy number estimate at the mid-point of the gene is selected as the single point estimate of the true copy number of the gene.
```
python get_true_cn.py ground-truths/res.samples.EUR.SMN1.bed > EUR.SMN1.truecn.out
```

If instead we wish to select a specific index of copy number intervals as the point estimate, we can also run 
```
python get_true_cn.py ground-truths/res.samples.EUR.SMN1.bed -d 1 > EUR.SMN1.truecn.out
```
where `-d` represents the 0-based index.

### Run main algorithm 

Now, we can run the main algorithm and also obtain the general accuracy statistics ouput using the following command:
```
python main.py \
-i EUR.SMN1.counts.exons.tsv \
-o EUR.acc.out \
-p EUR.betafit.out \
--truth EUR.SMN1.truecn.out \
-c 4
```
where `-c` represents the normal reference aggregate copy number. 

This will produce an output file like `EUR.acc.out` which will contain comprehensive summaries of `connected components`, `fractionalCN`, `integerCN`, `trueCN`, etc. 

### Benchmark: compute precision and recall

To compute precision and recall of a single population or a single output file, run
```
python compute_prec_rec.py -i EUR.SMN1.out -c 4
```

To compute precision and recall across multiple populations or multiple files altogether, run
```
python compute_prec_rec.py -I EUR.SMN1.out AFR.SMN1.out AMR.SMN1.out EAS.SMN1.out SAS.SMN1.out -c 4
```

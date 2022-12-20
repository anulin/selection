The tool searches for selection on polygenic traits. It consists of two scripts: scores.dll and getKEGGscores.py

Scores.dll

  Input:
    Two snp counts files in the followind format:
    SNP_ID  total_count countofallele_1 2 3 4 Nameofallele_1 2 3
    1_534247	96	192	0	0	0	C	T	.

    Two numbers of sample sizes

  Usage:
  
    dotnet scoresvcf.dll snpcountsfile1 snpcountsfile2 samplesize1 samplesize2 [--optional flags]
    --help - outputs help
    --maf - minor allele frequency filter (default 0.05)
    --mi - minimal fraction of individuals available for snp (default 1)
    
  Output:
  
    Scores file for snps with the following format:
      SNP_ID  SNP_score  Nameofallele_1 2 3
      1_754063 0.618055382764516 G T .
      
getKEGGscores.py
  Input:
  
    -Scores file for snps produced by Scores.dll
    -(optional) snp to gene annotation in the following format:
      Comment(can be alllele nucleotide)  SNP_ID  Gene
      .	1_871269	SAMD11
    by default annotated_exomeCeu.txt is used
    
    -(optional) gene to pathway annotation file in the following format:
      Gene  Pathway_ID  Human-readable_pathway_name
      AKR1A1	hsa00010	Glycolysis / Gluconeogenesis	
      ADH1A	hsa00010	Glycolysis / Gluconeogenesis	
    by default KEGG to Genes Upd_ShortNams.txt is used

  Usage:
  
    getKEGGscores scoresfile --a annotationfile --n networksfile --p 15
    scoresfile is produced by scores.dll 
    --a (optional) - custom file with snp to gene annotation
    --n (optional) - custom file with gene to pathway annotation
    --p (optional) - percentile of snp scores for each pathway that will be calculated as a pathway score. Default is 15%
    --h - help
    
  Output:
  
    Pathway scores and their p-values (obtained with permutation test) without adjustments.
    
    Format:
    
      Name  ID  P-value score snps genes
      Epstein-Barr virus infection hsa05169 0.000326349	0.107253818	413	202	

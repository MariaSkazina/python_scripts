## fastq_barbershop.py

fastq_barbershop.py is a Python 3 script for .fastq files. The script was created by Maria Skazina and [Olesya Dogonasheva](https://github.com/xomaiya) fastq_barbershop.py works in two different modes. Filtering mode allows to select reads of proper GC content and length, while trimming mode enable to cut nucleotides with low quality or cut the specified number of nucleotides from the start or from the end of read.


**Launching the script in filtering mode:**

```
python3 fastq_barbershop.py filter --gc_bounds 40 60 --min_length 50 your.fastq
```

Here your.fastq is a positional argument. Also we recommend to use --gc_bounds right after filter, for proper work of script. This command will give you an output file your__passed.fastq with reads having GC content rate between 40 and 60% and with length no more than 50 bp.

If you would like to keep the reads which were filtered out, you can use the flag --keep_filtered, like this:

```
python3 fastq_barbershop.py filter --gc_bounds 40 60 --min_length 50 --keep_filtered your.fastq
```

In this case you will have two files in output - with passed (your__passed.fastq) and with failed (your__failed.fastq) reads.


The additional help about this mode is available by command:

```
python3 fastq_barbershop.py filter -h
```


**Launching the script in the trimming mode:**

```
 python3 fastq_barbershop.py trimmer CROP --length 5 test_data.fastq 
```

This command will crop 5 nucleotides from the your reads and save it in the file named 'output'. 

```
 python3 fastq_barbershop.py trimmer SLIDINGWINDOW --threshold 32 --window 3 test_data.fastq 
```

This command crops your reads depending on the mean quality. It screens your read from the start and when the mean quality drops lower the threshold (32) it crops everithing on the right part of the read and the nucleotides in the window with size 3, where the dramatic drop of quality happened.

The additional help about this mode is available by command:

```
python3 fastq_barbershop.py trimmer -h
```





  

  

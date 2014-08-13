## Bioinformatics Models and Algorithms

In the fall of 2013 I took [this][0] course in bioinformatics at the University of California, Santa Cruz.  I took this course to gage my interest in the field of bioinformatics and also to sharpen my programing skills.  Here I present the programs I wrote for this course.

### Course Description

Covers bioinformatics models and algorithms: the use of computational techniques to convert the masses of information from biochemical experiments (DNA sequencing, DNA chips, and other high-throughput experimental methods) into useful information. Emphasis is on DNA and protein sequence alignment and analysis.

### Program Descriptions

Each program has its own README with more details (including the original assignment description).  Below is a list of the assignments with a brief description.

[HW1: Word Count][1] -  Break input text into words and count the frequency of each word.

[HW2: Parsing Sequences][2] - Parsers for reading some of the most common sequence formats.

[HW4: Simple Markov Chains][4] -  Given many protein sequences, I built stochastic models of these sequences as zero-order and first-order Markov chains, and measured the information gain of the first-order model over the zero-order model.

[HW5: Palindromes in Genomic Sequences][5] - DNA palindromes, which are reverse-complement palindromes (that is, the word is the same on both strands of the DNA—reversing the word and complementing it gives you back the word), seem to have biological importance.  Do these occur more often or less often in a given genome? A statistical model is created to predict how often we would expect to see a certain DNA palindromes to occur on chance, then this is compared to the actual frequency.

[HW6: Null Models][6] - An exercise in the importance of selecting a good null model.  This project examines four null models.  A researcher at UCSC had discovered an ORF (open reading frame) 388-codons long.  He wondered if this ORF coded for a protein, or if it just occurred by chance.  How unexpected is the 388-codon-long ORF? Using different null models, we arrive at different answers (though some null models are clearly better choices).

[HW7: Protein Information][7] - The Purpose of this assignment was to research a protein given a gene which coded for this protein.  We were given a DNA sequence of a PCR product that contained the gene for the protein.  From there I used various on-line tools to analyze the protein.  The code I wrote to find the longest ORF in the DNA is included here.

[HW8: Affine-gap Alignment][8] - In this assignment I wrote a program to align proteins.

[HW9 : Degenerate Codons][9] - This assignment is intended to deepen your understanding of cassette mutagenesis, the genetic code, and libraries of mutant genes.

### Future work

The individual README files for each assignment were created form the [html][10] assignment (converted to markdown with pandoc).  I am in the process of editing these - looking for errors in the automatic conversion with pandoc and also adding information about the project and the usage of the program.

Review program documentation.  A good exercise since I wrote these programs almost a year ago (if the documentation is good it should still make sense).

[0]: http://users.soe.ucsc.edu/~karplus/bme205/f13/index.html
[1]: https://github.com/lebailly/BME205/tree/master/HW1
[2]: https://github.com/lebailly/BME205/tree/master/HW2
[4]: https://github.com/lebailly/BME205/tree/master/HW4
[5]: https://github.com/lebailly/BME205/tree/master/HW5
[6]: https://github.com/lebailly/BME205/tree/master/HW6
[7]: https://github.com/lebailly/BME205/tree/master/HW7
[8]: https://github.com/lebailly/BME205/tree/master/HW8
[9]: https://github.com/lebailly/BME205/tree/master/HW9
[10]: http://users.soe.ucsc.edu/~karplus/bme205/f13/index.html#hw
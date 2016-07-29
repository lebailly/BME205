## Bioinformatics Models and Algorithms

In the fall of 2013 I took [this][0] course in bioinformatics at the University of California, Santa Cruz.  I took this course to gage my interest in the field of bioinformatics and also to sharpen my programing skills.  While it was a very enjoyable course which taught me lots about programing and modeling, I decided pursuing other avenues of data science was a better fit for me (instead of going for a second graduate degree in bioinformatics).  The programs I wrote for this class give a sample of my coding and modeling abilities, which is why I put them in a public repo.

### Course Description

This course covered bioinformatics models and algorithms.  Specifically, it covered the use of computational techniques to convert the masses of information from biochemical experiments (DNA sequencing, DNA chips, and other high-throughput experimental methods) into useful information. Emphasis was on DNA and protein sequence alignment and analysis.

### Program Descriptions

Each program has its own README with more details.  These were taken from the [assignment descriptions][10] and translated to markdown using [pandoc][11].  Below is a list of the assignments with a brief description.

[HW1: Word Count][1] - A warm up exercise where we break input text into words and count the frequency of each word.

[HW2: Parsing Sequences][2] - Building parsers for reading some of the most common sequence formats.

[HW3: Fellowship Application][3] - The purpose of this assignment was to practice writing a fellowship application.  I have included the link to the assignment, but since it since no code or modeling was involved with this project I did not include anything else in the repo for this assignment.

[HW4: Simple Markov Chains][4] -  Given many protein sequences, I built stochastic models of these sequences as zero-order and first-order Markov chains, and measured the information gain of the first-order model over the zero-order model.

[HW5: Palindromes in Genomic Sequences][5] - DNA palindromes, which are reverse-complement palindromes (that is, the word is the same on both strands of the DNA—reversing the word and complementing it gives you back the word), seem to have biological importance.  Do these occur more often or less often in a given genome? A statistical model is created to predict how often we would expect to see a certain DNA palindromes to occur on chance, then this is compared to the actual frequency.

[HW6: Null Models][6] - An exercise in the importance of selecting a good null model.  This project examines four null models.  A researcher at UCSC had discovered an ORF (open reading frame) 388-codons long.  He wondered if this ORF coded for a protein, or if it just occurred by chance.  How unexpected is the 388-codon-long ORF? Using different null models, we arrive at different answers (though some null models are clearly better choices).

[HW7: Protein Information][7] - The purpose of this assignment was to research a protein given a gene which coded for this protein.  I was given a DNA sequence of a PCR product that contained the gene for the protein.  From there I used various on-line tools to analyze the protein.  The code I wrote to find the longest ORF in the DNA is included here.

[HW8: Affine-gap Alignment][8] - In this assignment I wrote a program to align proteins.  Given a [matrix][13] which scores how biologically likely one amino acid can be substituted for another, we can score how well two different proteins align.  This program uses [dynamic programing][12] to determine the optimal alignment and the best possible alignment score.  This was one of the more challenging assignments, but also my favorite.

[HW9: Degenerate Codons][9] - In this assignment I studied cassette mutagenesis (the production of mutants within a region).  I wrote a program to produce a table containing all possible sets of amino acids that can be encoded by a single degenerate codon.

[0]: http://users.soe.ucsc.edu/~karplus/bme205/f13/index.html
[1]: https://github.com/lebailly/BME205/tree/master/HW1
[2]: https://github.com/lebailly/BME205/tree/master/HW2
[3]: https://users.soe.ucsc.edu/%7Ekarplus/bme205/f13/Fellowship.html
[4]: https://github.com/lebailly/BME205/tree/master/HW4
[5]: https://github.com/lebailly/BME205/tree/master/HW5
[6]: https://github.com/lebailly/BME205/tree/master/HW6
[7]: https://github.com/lebailly/BME205/tree/master/HW7
[8]: https://github.com/lebailly/BME205/tree/master/HW8
[9]: https://github.com/lebailly/BME205/tree/master/HW9
[10]: http://users.soe.ucsc.edu/~karplus/bme205/f13/index.html#hw
[11]: http://pandoc.org/
[12]: https://en.wikipedia.org/wiki/Dynamic_programming
[13]: https://github.com/lebailly/BME205/blob/master/HW8/BLOSUM62.txt
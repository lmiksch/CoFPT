# Cotranscriptional Folding Path designer

This script was created during the course software project in the bioinformatics master at the university of vienna. 

CotFPD is a python based script which consits of two steps: 

* domain_seq_generator: Given a folding path as an input, it creates a domain based sequence which should fold like the given input.  
* ir_domain_translator: With the domain based sequence, this part now translates the domain based sequence into a nucleotide sequence using the [Infrared Package](https://www.lix.polytechnique.fr/~will/Software/Infrared/Doc/index.html). 


## Prerequisites

Following packages need to be installed to use ir_domain_translator: 

[ViennaRNA Package](https://github.com/ViennaRNA/ViennaRNA)

[Infrared Package](https://www.lix.polytechnique.fr/~will/Software/Infrared/Doc/index.html)


## Getting Started

To use the scripts you need to clone the repository using git: 

```bash
# Clone this repository
$ git clone https://github.com/lmiksch/Software_Project

#Change your directory 
$ cd Software_Project
```

## Using domain_seq_generator

Takes an abstract folding path as input and calculates a domain based sequence which should have the same general co transcriptional folding pathway as the input. 

In an additonal step it calculates the cotranscriptional folding pathway using the Nussinov Algorithm. It then compares the input to the output of the nussinov algorithm. 

Input is given as a txt file where each line break indicates a new transcription step using a dot-bracket annotation.

#### Example input: 

```
.
()
.()
()()
```
Example command:

```bash
$ python3 domain_seq_generator.py -i input_test.txt

```

Output: 

A file will be generated which consists of the domain based sequence and additional informations. The output file can then be used for the second part. 

Limitations:

Due to certain restrictions of our algorithm and domain design, some cotranscritpional pathways are not possible to design using our algorithm and domain design. 
This will be validated after the domain sequence design step. 
In the event, that the calculated path from the nussinov algorithm does not match the input path. The domain based sequence is still put out. Be aware that this sequence has a similar folding path as the input path but will ultimatley fold in such a way as calculated in the nussinov algorithm. 


## domain_translator

Using the [Infrared Package](https://www.lix.polytechnique.fr/~will/Software/Infrared/Doc/index.html) the domain level sequence gets translated into a nucleotide level sequence. 

Input: 

* Using the output file of domain_seq_generator.py 
	```bash
	$ pyhton3 ir_domain_translator -i domain_seq_out.txt
	```
* Using the command line to input a domain level sequence
	```bash 
	$ python3 ir_domain_translator
	>Input here a domain level sequence: 
	```

Ouput: 

A textfile will be created in where a nucleotide sequence is displayed aswell as the score of the corresponding sequence. 





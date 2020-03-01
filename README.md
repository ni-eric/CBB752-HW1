# CBB752-HW1
Homework Assignment 1 for Biomedical Data Science 2020 ([CBB752](http://cbb752b20.gersteinlab.org/))

Implementation of the Smith-Waterman local alignment algorithm.

## Prerequisites:
- Python 2.7.x

## Installation:
Install git and then:
```
git clone git://github.com/ni-eric/CBB752-HW1.git
```

## Input:
2 Sequences for comparison, which contain no spaces and are separated by a newline
### Example Input:
```
ERICNI
EEERICCCNIIIIII
```

### Usage:      
```
python hw1.py -i <input file> -s <score file> -o <open gap penalty> -e <extend gap panalty> -O <output file name>
```
The program to run is 'hw1.py'. It takes 5 parameters:
1. Input Filename: 2 sequences for performing aligment on (Required)
2. Score Matrix Filename: Substitution matrix for scoring alignments (Required)
3. Open Gap Penalty: Penalty for opening a gap when creating alignment (Optional, default = -2)
4. Extend Gap Penalty: Penalty for extending an already open gap (Optional, default = -1)
5. Output Filename: File name for output (Optional, default = 'output.txt')

### Example Usages:  
```
python hw1.py -i input.txt -s blosum62.txt
python hw1.py -i input.txt -s blosum62.txt -o -3 -e 0 -O 'output.txt'
python hw1.py -i input.txt -s blosum62.txt -O 'output.txt'
```

## Output:
The program will output back your original input sequences, followed by the calculated scoring matrix according to the given subtitution matrix and open/gap penalty scores. Finally, there is also the best local alignment score, followed by a human-readable alignment.
### Example Output:
```
************************************************************
*    Input Sequences                                       *
************************************************************

Sequence 1	ERICNI
Sequence 2	EEERICCCNIIIIII

************************************************************
*    Score Matrix                                          *
************************************************************

		E	R	I	C	N	I	
	0	0	0	0	0	0	0	
E	0	5	3	2	1	0	0	
E	0	5	5	3	2	1	0	
E	0	5	5	3	2	2	0	
R	0	3	10	8	7	6	5	
I	0	2	8	14	12	11	10	
C	0	1	7	12	23	21	20	
C	0	0	6	11	21	20	20	
C	0	0	5	10	20	18	19	
N	0	0	4	9	19	26	24	
I	0	0	3	8	18	24	30	
I	0	0	2	7	17	23	28	
I	0	0	1	6	16	22	27	
I	0	0	0	5	15	21	26	
I	0	0	0	4	14	20	25	
I	0	0	0	4	13	19	24	

************************************************************
*    Best Local Alignment                                  *
************************************************************

Maximum Alignment Score:[30]
  [ERI--CNI]
   |||  ||| 
EE[ERICCCNI]IIIII

```
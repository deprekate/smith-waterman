Introduction
------------

The included perl script is a standalone script for aligning two sequences using the Smith-Waterman algorithm.

Smith-Waterman Example
--------------

To run you supply two sequences to the script as command line arguements:
```sh
./smith-waterman.pl ATCGGG ATCCGG
```
Output is the simple text alignment of the two sequences followed by the percent identity.
```sh
ATCGGG
ATCCGG
0.833333333333333
```

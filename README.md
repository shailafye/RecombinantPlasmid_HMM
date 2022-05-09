# RecombinantPlasmid_HMM
By: Shaila Fye, Stella Park, Roshni Merugu

05/09/2022

Final group project for Computational Genomics. 
HMM to detect foreign insert/gene in recombinant plasmid based on GC content and restriction site verification.

-----------

###### **Files in RecombinantPlasmid_HMM Directory:**

**detect_insert_HMM.py** 
* This file has the code to instantiate and run the InsertDetector class.
* InsertDetector takes in a sequence and returns the top 1 or 2 detected intervals with confidence scoring
* To test, run the main method at the bottom of the file. 
* The test is dependent on recombinant_seqs.txt file in this directory
* This class can also be tested by inputting a recombinant plasmid sequence
* The output is printed to the terminal
* In the main method, there is also an option to save to csv called output_intervals.csv 

**recombinant_seqs.txt**
* This file is the data we used to test and run InsertDetector. It includes 30 sequences of recombinant plasmids
* We simulated and generated this data by using 3 different plasmids and 4 different genes
* We inserted the gene at different locations in the plasmid
* Each sequence is separated by a new line character

**output_intervals.csv**
* This is where output can be saved to in the main method of detect_insert_HMM.py file
* Saved and pushed results to github as an example of InsertDetector output


-----------

Everything is run on Python 3.7 in PyCharm IDE.



"""
This is the code for the User to input a sequence and an HMM will detect where the
insert DNA is detected

Assumption is that plasmid is longer than insert and plasmid GC content is greater than insert GC content
"""

class InsertDetectHMM:

    def __init__(self, sequence):
        self.sequence = sequence

        # run first HMM on this sequence
        self.detected_intervals = self.first_HMM()

        # run verification steps


        # run second HMM on new intervals for length of insert and GC of insert


        # check confidence scoring


    def first_HMM(self):
        # Returns detected intervals

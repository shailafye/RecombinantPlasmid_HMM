"""
This is the code for the User to input a sequence and an HMM will detect where the
insert DNA is detected

Assumption is that plasmid is longer than insert and plasmid GC content is greater than insert GC content
"""

import re


class InsertDetectHMM:

    def __init__(self, sequence):
        self.sequence = sequence.upper()
        self.inter_dict_score = {}  # all interval scores
        self.top_dict = {}  # top 2 or 3 intervals

        # run first HMM on this sequence
        self.detected_intervals = self.first_HMM()

        # run verification steps

        # run second HMM on new intervals for length of insert and GC of insert

        # check confidence scoring

    def first_HMM(self):
        # Returns detected intervals in [[a,b],[c,d]...] form
        detected_intervals = []
        ## HMM BELOW ###

        #####
        self.detected_intervals = detected_intervals
        return detected_intervals

    # Function to get top three intervals
    def filter_best_intervals(self):
        # first check to see if each one has at least 1 restriction site detected on either side
        # then return top intervals

        # if invalid = length 0 then all good, if not, remove
        invalid_interval, valid_interval = self.check_restriction_site_valid()
        if len(invalid_interval) != 0:
            print("REMOVE", invalid_interval)
            for i in invalid_interval:
                if i[0] in self.detected_intervals:
                    self.detected_intervals.remove(i[0])
                if i[1] in self.detected_intervals:
                    self.detected_intervals.remove(i[1])

        # code to return top 3 intervals with scoring (higher the better)
        return self.best_intervals(self.detected_intervals)

    # three helper functions below to get best intervals and restriction sites
    def best_intervals(self, insert_intervals):
        # insert_intervals = [[0, 569], [1605, 1964], [7216, 7389]]
        # tests combinations of 2 and 3 fragments to get best combo of intervals
        # higher score = better
        interval_lengths = []
        prev_index = 0
        for i in insert_intervals:
            interval_lengths.append(i[0] - prev_index)
            interval_lengths.append(i[1] - i[0])
            prev_index = i[1]
        interval_lengths.append(len(self.sequence) - prev_index + interval_lengths[0])
        interval_lengths = interval_lengths[1:]

        naming_intervals = [item for sublist in insert_intervals for item in sublist]
        # all intervals and scores = self.inter_dict_score

        if len(insert_intervals) > 1:
            for i in range(1, len(interval_lengths), 2):
                a = i - 1
                b = i
                c = (i + 1) % len(interval_lengths)
                res = (interval_lengths[b] / (interval_lengths[a] + interval_lengths[c] + interval_lengths[b]))
                self.inter_dict_score[1 - res] = [naming_intervals[a],
                                                  naming_intervals[(a + 3) % len(interval_lengths)]]
        if len(insert_intervals) > 2:
            for i in range(1, len(interval_lengths), 2):
                a = i - 1
                b = i
                c = (i + 1) % len(interval_lengths)
                d = (i + 2) % len(interval_lengths)
                e = (i + 3) % len(interval_lengths)
                res = (interval_lengths[b] + interval_lengths[d]) / (
                        interval_lengths[a] + interval_lengths[b] + interval_lengths[c] + interval_lengths[d] +
                        interval_lengths[e])
                self.inter_dict_score[1 - res] = [naming_intervals[a],
                                                  naming_intervals[(a + 5) % len(interval_lengths)]]
        # need to add if there are more than 3 intervals
        if len(insert_intervals) > 3:
            print("MORE THAN 3 Intervals")
            # CHECK!
        # top intervals
        keys = list(self.inter_dict_score.keys())
        keys.sort(reverse=True)
        top_keys = keys[0:3]
        for k in top_keys:
            self.top_dict[k] = self.inter_dict_score[k]
        return self.top_dict

    def find_rest_sites(self, offset=0):
        # offset is the index value your sequence starts at
        rest_dict = {
            "sbfi": "TGCA", "bglii": "GATC", "eari": "TGG", "bpuei": "CC", "bspdi": "CG", "maubi": "CGCG",
            "bsrgi": "GTAC", "nsii": "TGCA",
            "pcii": "CATG", "bspei": "CCGG", "tati": "GTAC", "bsssi": "TCGT", "pfoi": "CCCGG", "pspomi": "GGCC"
        }

        sites_found = {}
        for site in rest_dict:
            res = [m.start() for m in re.finditer(rest_dict[site], self.sequence)]
            res = [r + offset for r in res]
            if len(res) > 0:
                sites_found[site] = res
        return sites_found

    def check_restriction_site_valid(self):
        # if two intervals share no restriction sites --> return to look like: [0, 569], [1605, 1964]
        invalid_interval = []
        valid_interval = []
        test_seq_double = self.sequence + self.sequence
        # returns list of invalid intervals --> if list is empty then use orignal detected_intervals list

        for i in range(len(self.detected_intervals)):
            interval = self.detected_intervals[i]
            int_1 = [interval[0] - 20, interval[0] + 50]
            # if interval is negative because  circular ex)[-20,50]
            if int_1[0] < 0:
                int_1 = [len(self.sequence) - abs(int_1[0]), len(self.sequence) + abs(int_1[1])]

            # find all possible restriction sites for the first interval
            sites_1 = self.find_rest_sites(test_seq_double[int_1[0]:int_1[1]], int_1[0])

            for j in range(len(self.detected_intervals)):
                interval_end = self.detected_intervals[j]
                int_2 = [interval_end[1] - 20, interval_end[1] + 50]
                # print("int2",int_2)
                # print(interval, interval_end)
                sites_2 = self.find_rest_sites(test_seq_double[int_2[0]:int_2[1]], int_2[0])
                common_sites = list(set(sites_1.keys()) & set(sites_2.keys()))
                if len(common_sites) < 1:
                    # print("NO RESTRICTION SITE COMMON")
                    invalid_interval.append([interval, interval_end])
                else:
                    valid_interval.append([interval, interval_end])

        return invalid_interval, valid_interval
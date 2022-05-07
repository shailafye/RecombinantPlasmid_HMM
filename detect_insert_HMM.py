"""
This is the code for the User to input a sequence and an HMM will detect where the
insert DNA is detected

Assumption is that plasmid is longer than insert and plasmid GC content is greater than insert GC content
"""

import re
from collections import Counter
import pandas as pd
import csv
import numpy as np


class DetectInsertHMM:

    def __init__(self, sequence):
        self.sequence = sequence.upper()
        self.inter_dict_score = {}  # all interval scores
        self.filtered_intervals = {} # top 2 or 3 intervals for the sequence to run again on HMM2
        self.detected_intervals = [] # first detected intervals HMM
        #self.second_detected_intervals = [] # second detected intervals HMM

        # run first HMM on this sequence
        self.detected_intervals = self.first_HMM()

        # run verification steps
        if len(self.detected_intervals) > 0:
            self.filtered_intervals = self.filter_best_intervals()
        else:
            print("No INSERT detected!")
            return

        # run second HMM on new intervals for length of insert and GC of insert
        self.final_intervals_score = {}
        for i in self.filtered_intervals:
            print(i, self.filtered_intervals[i])
            intervals_output = self.second_HMM(self.filtered_intervals[i])
            print(intervals_output)
            # check confidence scoring on new intervals and multiply
            #sc = self.confidence_score_secondHMM(ground=self.filtered_intervals[i], predict=intervals_output)
            #print(sc)



    def first_HMM(self):
        # Returns detected intervals in [[a,b],[c,d]...] form

        ### Calculate GC content ###

        # recomb dna length and gc content
        gc_recomb = gc_content(self.sequence)
        length_recomb = len(self.sequence)

        # gc content of insert and plasmid
        gc_insert = 1 - gc_recomb  # under assumption that plasmid gc content is always higher # actinobacteria
        gc_plasmid = gc_recomb

        # at content of insert and plasmid
        at_insert = 1-gc_insert
        at_plasmid = 1-gc_plasmid

        # length of insert and plasmid
        temp_length_insert = 0.3 * length_recomb # double check average -> 30%?
        temp_length_plasmid = length_recomb
        temp_length_total = temp_length_insert + temp_length_plasmid
        length_insert = temp_length_insert
        length_plasmid = temp_length_plasmid

        # length portion of insert in recombinant DNA
        fraction_insert = temp_length_insert / temp_length_total
        fraction_plasmid = 1 - fraction_insert

        ### Matricies ###

        observations = self.sequence
        observations = observations.replace('\r', '')
        states = ("insert", "vector")
        start_p = {
            "insert": fraction_insert,
            "vector": fraction_plasmid
        }
        trans_p = {
            "insert": {"insert": 1 - (1 / length_insert), "vector": 1 / length_insert},
            "vector": {"insert": 1 / length_plasmid, "vector": 1 - (1 / length_plasmid)},
        }
        emit_p = {
            "insert": {"A": at_insert / 2, "T": at_insert / 2, "C": gc_insert / 2,
                       "G": gc_insert / 2},
            "vector": {"A": at_plasmid / 2, "T": at_plasmid / 2, "C": gc_plasmid / 2,
                       "G": gc_plasmid / 2},
        }

        ### RUN HMM BELOW ###
        detected_intervals = viterbi_algorithm(observations, states, start_p, trans_p, emit_p)

        #####
        self.detected_intervals = detected_intervals
        return detected_intervals

    def second_HMM(self, candidate):
        # feedback interval sequence --> what first HMM predicts and scored
        test_insert_interval = ""

        if candidate[0] <= candidate[1]:
            test_insert_interval += self.sequence[candidate[0]:candidate[1]]
        else:
            test_insert_interval += self.sequence[candidate[0]:]
            test_insert_interval += self.sequence[:candidate[1]]

        # length of insert and plasmid
        length_plasmid = len(self.sequence) - len(test_insert_interval)
        length_insert = len(test_insert_interval)

        # gc content of insert and plasmid
        gc_insert = gc_content(test_insert_interval)
        gc_plasmid = (gc_content(self.sequence) * len(self.sequence) - gc_insert * length_insert) / length_plasmid

        # at content of insert and plasmid
        at_insert = 1 - gc_insert
        at_plasmid = 1 - gc_plasmid

        # length portion of insert in recombinant DNA
        fraction_insert = length_insert / (length_insert + length_plasmid)
        fraction_plasmid = 1 - fraction_insert

        observations = self.sequence
        observations = observations.replace('\r', '')
        states = ("insert", "vector")
        start_p_fb = {
            "insert": fraction_insert,
            "vector": fraction_plasmid
        }
        trans_p_fb = {
            "insert": {"insert": 1 - (1 / length_insert), "vector": 1 / length_insert},
            "vector": {"insert": 1 / length_plasmid, "vector": 1 - (1 / length_plasmid)},
        }
        emit_p_fb = {
            "insert": {"A": at_insert / 2, "T": at_insert / 2, "C": gc_insert / 2, "G": gc_insert / 2},
            "vector": {"A": at_plasmid / 2, "T": at_plasmid / 2, "C": gc_plasmid / 2, "G": gc_plasmid / 2},
        }

        ### RUN SECOND HMM BELOW ###
        second_detected_intervals = viterbi_algorithm(observations, states, start_p_fb, trans_p_fb, emit_p_fb)

        #####
        #self.second_detected_intervals = second_detected_intervals
        return second_detected_intervals

    def get_intersection(self, ground, predict):
        ground_ranges = []
        for g in ground:
            ground_ranges.append(list(range(g[0], g[1] + 1)))
        predict_ranges = []
        for p in predict:
            predict_ranges.append(list(range(p[0], p[1] + 1)))

        intersect = sorted(set(ground_ranges).intersection(predict_ranges))
        return len(intersect)

    def get_union(self, ground, predict):
        ground_ranges = []
        for g in ground:
            ground_ranges.append(list(range(g[0], g[1] + 1)))
        predict_ranges = []
        for p in predict:
            predict_ranges.append(list(range(p[0], p[1] + 1)))

        intersect = sorted(set(ground_ranges).union(predict_ranges))
        return len(intersect)

    def confidence_score_secondHMM(self, ground, predict):
        # pass in intervals: truth and observed
        # calculate the score (higher > more confident)
        # return score as float
        score = self.get_intersection(ground, predict) / self.get_union(ground, predict)
        return score

    # Function to get top three intervals
    def filter_best_intervals(self):
        # if there is only one interval, just return it with score == 1
        if len(self.detected_intervals) == 1:
            return {1: self.detected_intervals}

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
        top_dict = {}  # top 2 or 3 intervals
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
            top_dict[k] = self.inter_dict_score[k]
        return top_dict

    def find_rest_sites(self, seq, offset=0):
        # offset is the index value your sequence starts at
        rest_dict = {
            "sbfi": "TGCA", "bglii": "GATC", "eari": "TGG", "bpuei": "CC", "bspdi": "CG", "maubi": "CGCG",
            "bsrgi": "GTAC", "nsii": "TGCA",
            "pcii": "CATG", "bspei": "CCGG", "tati": "GTAC", "bsssi": "TCGT", "pfoi": "CCCGG", "pspomi": "GGCC"
        }

        sites_found = {}
        for site in rest_dict:
            res = [m.start() for m in re.finditer(rest_dict[site], seq)]
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


# calculate GC content
def gc_content(sequence):
    res = Counter(sequence)
    gc_sum = res["g"] + res["c"] + res["G"] + res["C"]
    # at_sum = res["a"]+res["t"]+res["A"]+res["T"]
    total_len = len(sequence)
    # assert total_len == at_sum+ gc_sum
    gc_perc = gc_sum / total_len
    return gc_perc


# viterbi algorithm
def viterbi_algorithm(observations, states, start_p, trans_p, emit_p):
    V = [{}]
    for st in states:
        V[0][st] = {"log_prob": np.log(start_p[st] * emit_p[st][observations[0]]), "prev": None}  # storing start probability of first position
    for t in range(1, len(observations)):
        V.append({})
        for st in states:
            log_max_tr_prob = V[t - 1][states[0]]["log_prob"] + np.log(trans_p[states[0]][st])  # maximum transition probability for site 1
            prev_st_selected = states[0]
            for prev_st in states[1:]:  # maximum transition probability for the rest
                log_tr_prob = V[t - 1][prev_st]["log_prob"] + np.log(trans_p[prev_st][st])
                if log_tr_prob > log_max_tr_prob:
                    log_max_tr_prob = log_tr_prob
                    prev_st_selected = prev_st

            log_max_prob = log_max_tr_prob + np.log(emit_p[st][observations[t]])
            V[t][st] = {"log_prob": log_max_prob, "prev": prev_st_selected}

    opt = []
    log_max_prob = -np.inf
    best_st = None

    for st, data in V[-1].items():  # selecting the best prediction
        if data["log_prob"] > log_max_prob:
            log_max_prob = data["log_prob"]
            best_st = st

    opt.append(best_st)  # best state sequence
    previous = best_st

    insert_indices = []

    dict_vec_insert = {"vector": 0, "insert": 1}

    total = 0
    for t in range(len(V) - 2, -1, -1):
        opt.append(dict_vec_insert[V[t + 1][previous]["prev"]])
        if dict_vec_insert[V[t + 1][previous]["prev"]] == 1:
            total += 1
            insert_indices.append(t)
        previous = V[t + 1][previous]["prev"]

    def interval_extract(lists):
        lists = sorted(set(lists))
        range_start = previous_number = lists[0]

        for number in lists[1:]:
            if number == previous_number + 1:
                previous_number = number
            else:
                yield [range_start, previous_number]
                range_start = previous_number = number
        yield [range_start, previous_number]
    insert_result_intervals = list(interval_extract(insert_indices))

    return insert_result_intervals




if __name__ == '__main__':
    all_recomb_seqs = []
    with open("recombinant_seqs.txt") as file:
        for line in file:
            all_recomb_seqs.append(line.rstrip())

    test_num = 15
    test_seq1 = all_recomb_seqs[test_num]
    print(len(test_seq1))
    obj1 = DetectInsertHMM(test_seq1)
    #print(obj1.detected_intervals)
    #print(obj1.filtered_intervals)



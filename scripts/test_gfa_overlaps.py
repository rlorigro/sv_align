import argparse
from collections import defaultdict


class CigarStats:
    def __init__(self, n_matches, n_mismatches, n_inserts, n_deletes):
        self.n_matches = n_matches
        self.n_mismatches = n_mismatches
        self.n_inserts = n_inserts
        self.n_deletes = n_deletes

    def __str__(self):
        return str(self.n_matches) + '\t' + str(self.n_mismatches )+ '\t' + str(self.n_inserts )+ '\t' + str(self.n_deletes)


class Cigar:
    def __init__(self, cigar_string):
        self.tuples = list()
        self.parse_as_tuples(cigar_string)

    def parse_as_tuples(self, cigar_string):
        count_string = ""

        for c in cigar_string:
            if c.isalpha():
                operation = c

                self.tuples.append((operation, int(count_string)))
                count_string = ""

            if c.isdigit():
                count_string += c

    def __str__(self):
        return str(self.tuples)


class Link:
    def __init__(self, reversal_a, reversal_b, cigar):
        self.reversal_a = self.parse_reversal_char_as_bool(reversal_a)
        self.reversal_b = self.parse_reversal_char_as_bool(reversal_b)
        self.cigar = cigar

        self.reversal_bool_to_char = ["+","-"]

    def __str__(self):
        s = self.reversal_bool_to_char[self.reversal_a] + " " + self.reversal_bool_to_char[self.reversal_b] + " " + str(self.cigar)
        return s

    def parse_reversal_char_as_bool(self, reversal_char):
        if reversal_char == '+':
            return False

        elif reversal_char == "-":
            return True

        else:
            exit("ERROR: invalid reversal character cannot be parsed: " + reversal_char)



def parse_gfa(gfa_path):
    sequences = dict()
    links = defaultdict(dict)

    with open(gfa_path, 'r') as file:
        for line in file:
            data = line.strip().split('\t')

            # Sequence data
            if line[0] == 'S':
                node = data[1]
                sequence = data[2]

                sequences[node] = sequence

            # Overlap data
            elif line[0] == 'L':
                node_a = data[1]
                reversal_a = data[2]
                node_b = data[3]
                reversal_b = data[4]

                cigar = Cigar(data[5])

                link = Link(reversal_a, reversal_b, cigar)
                links[node_a][node_b] = link

    return sequences, links


def count_alignment_lengths(cigar):
    length_a = 0
    length_b = 0

    for cigar_tuple in cigar.tuples:
        if cigar_tuple[0] == "M":
            length_a += cigar_tuple[1]
            length_b += cigar_tuple[1]

        elif cigar_tuple[0] == "I":
            length_b += cigar_tuple[1]

        elif cigar_tuple[0] == "D":
            length_a += cigar_tuple[1]

        else:
            exit("ERROR: invalid cigar operation: " + cigar_tuple[0])

    return length_a, length_b


def count_true_alignment_stats(sequence_a, sequence_b, link):
    alignment_length_a, alignment_length_b = count_alignment_lengths(link.cigar)
    alignment_string_a = ""
    alignment_string_b = ""
    alignment_string_ab = ""
    warning = ""

    true_cigar_stats = CigarStats(0,0,0,0)
    gfa_cigar_stats = CigarStats(0,0,0,0)

    invalid_alignment_length = False
    if alignment_length_a > len(sequence_a):
        invalid_alignment_length = True
        warning += "ERROR: invalid alignment length " + str(alignment_length_a) + " for sequence length " + str(len(sequence_a)) + '\n'
    if alignment_length_b > len(sequence_b):
        invalid_alignment_length = True
        warning += "ERROR: invalid alignment length " + str(alignment_length_b) + " for sequence length " + str(len(sequence_b)) + '\n'

    if link.reversal_a:
        index_a = alignment_length_a
        offset_a = -1
    else:
        index_a = len(sequence_a) - alignment_length_a
        offset_a = 1

    if link.reversal_b:
        index_b = len(sequence_b) - 1
        offset_b = -1
    else:
        index_b = 0
        offset_b = 1

    if not invalid_alignment_length:
        for cigar_tuple in link.cigar.tuples:
            for i in range(cigar_tuple[1]):
                if cigar_tuple[0] == "M":
                    gfa_cigar_stats.n_matches += 1

                    if sequence_a[index_a] == sequence_b[index_b]:
                        true_cigar_stats.n_matches += 1
                        alignment_string_ab += "|"
                        alignment_string_a += sequence_a[index_a]
                        alignment_string_b += sequence_b[index_b]
                    else:
                        true_cigar_stats.n_mismatches += 1
                        alignment_string_ab += " "
                        alignment_string_a += sequence_a[index_a]
                        alignment_string_b += sequence_b[index_b]

                    index_a += offset_a
                    index_b += offset_b


                elif cigar_tuple[0] == "I":
                    gfa_cigar_stats.n_inserts += 1
                    true_cigar_stats.n_inserts += 1
                    alignment_string_ab += " "
                    alignment_string_a += "-"
                    alignment_string_b += sequence_b[index_b]
                    index_b += offset_b

                elif cigar_tuple[0] == "D":
                    true_cigar_stats.n_deletes += 1
                    gfa_cigar_stats.n_deletes += 1
                    alignment_string_ab += " "
                    alignment_string_a += sequence_a[index_a]
                    alignment_string_b += "-"
                    index_a += offset_a

                else:
                    exit("ERROR: invalid cigar operation: " + cigar_tuple[0])

    alignment_string = alignment_string_a + '\n' + alignment_string_ab + '\n' + alignment_string_b

    return true_cigar_stats, gfa_cigar_stats, alignment_string, warning


def main(gfa_path):
    sequences, links = parse_gfa(gfa_path)

    for node_a in links:
        for node_b in links[node_a]:
            link = links[node_a][node_b]

            true_cigar_stats, gfa_cigar_stats, alignment_string, warning = \
                count_true_alignment_stats(sequences[node_a], sequences[node_b], link)

            if true_cigar_stats.n_mismatches > 0 or len(warning) > 0:
                warning += "ERROR: mismatch found in alignment"

                print(node_a + '\t' + node_b)
                print(link)
                print(true_cigar_stats)
                print(gfa_cigar_stats)
                print(alignment_string)
                print(warning)
                print()


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--input",
        type=str,
        required=True,
        help="path of file containing quadrant bounds"
    )

    args = parser.parse_args()

    main(gfa_path=args.input)
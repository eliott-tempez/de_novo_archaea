import re
import numpy as np
import psa
from my_functions.genomic_functions import extract_denovo_info, get_sequence_from_loci


from my_functions.paths import GENOMES_LIST


def get_extended_matched_seq(genome, contig, start, end, strand, missing_cter, missing_nter):
    """
    Get the sequence (nucleotides) of the match, extended on each relevant side.
    """
    # Extend sequence on both sides
    if strand == "+":
        add_right = (missing_cter * 3) * 2 if missing_cter > 4 else 0
        add_left = (missing_nter * 3) * 2 if missing_nter > 4 else 0
    elif strand == "-":
        add_left = (missing_cter * 3) * 2 if missing_cter > 4 else 0
        add_right = (missing_nter * 3) * 2 if missing_nter > 4 else 0
    match_start = start - add_left
    match_end = end + add_right
    # Extract extended sequence
    seq = get_sequence_from_loci(genome, contig, match_start, match_end, strand)
    # Get the limits of the actual match
    limit_start = add_left if strand == "+" else add_right
    limit_end = (len(seq) - add_right) if strand == "+" else (len(seq) - add_left)
    return seq, limit_start, limit_end


def in_frame(pos, start_or_end):
    """Get the closest number that is a 3 multiple"""
    if pos % 3 == 0:
        return pos
    is_start = start_or_end == "start"
    if is_start:
        if (pos - 1) % 3 == 0:
            return pos - 1
        return pos - 2
    if (pos + 1) % 3 == 0:
        return pos + 1
    return pos + 2


def smith_waterman(ref_seq, subject_seq):
    aln_dict = {}
    # Remove stops
    ref_seq = re.sub(r"\*", "", str(ref_seq))
    subject_seq = re.sub(r"\*", "", str(subject_seq))
    # Align
    aln = psa.water(moltype = "prot", qseq = ref_seq, sseq = subject_seq)
    # Store alignment in dict
    aln_dict["pval"] = aln.pvalue()
    aln_dict["psim"] = aln.psimilarity
    aln_dict["length"] = aln.length
    aln_dict["qstart"], aln_dict["qend"] = aln.qstart, aln.qend
    aln_dict["sstart"], aln_dict["send"] = aln.sstart, aln.send
    # Uncomment to show the raw alignment output
    #aln_dict["raw"] = aln.raw
    return aln_dict


def is_significant(alignment):
    if not alignment:
        return False
    aln_len = alignment["length"]
    aln_pval = alignment["pval"]
    if aln_len >= 5 and aln_pval <= 10e-3:
        return True
    return False


def get_absolute_match(aln, start_pos_query, start_pos_subject):
    # For the qstart
    relative_qstart = aln["qstart"]
    absolute_qstart = relative_qstart + start_pos_query
    aln["qstart"] = absolute_qstart
    # For the sstart
    frame = aln["frame"]
    relative_sstart = aln["sstart"]
    absolute_sstart = relative_sstart * 3 + start_pos_subject - frame
    aln["sstart"] = absolute_sstart
    # For the qend
    relative_qend = aln["qend"]
    absolute_qend = relative_qend + start_pos_query
    aln["qend"] = absolute_qend
    # For the send
    relative_send = aln["send"]
    absolute_send = relative_send * 3 + start_pos_subject + frame
    aln["send"] = absolute_send
    return aln


def get_segments_from_set(indexes, threshold):
    """
    Take a set of continuous numbers and extract the segments of consecutive numbers.
    """
    continuous_segments = []
    unmatched_indexes = sorted(list(indexes))
    if len(indexes) < threshold:
        return []
    segment_start = unmatched_indexes[0]
    for i in range(len(unmatched_indexes) - 1):
        if unmatched_indexes[i+1] - unmatched_indexes[i] > 1:
            segment_end = unmatched_indexes[i]
            if segment_end - segment_start >= threshold - 1:
                continuous_segments.append((segment_start, segment_end))
            segment_start = unmatched_indexes[i+1]
    if unmatched_indexes[i+1] - unmatched_indexes[i] == 1:
        if unmatched_indexes[i+1] - segment_start >= threshold - 1:
            continuous_segments.append((segment_start, unmatched_indexes[i+1]))
    return continuous_segments


def get_uncovered_segments(matches, query_len, subject_len):
    """
    Get all the segments > 4 aa in the query and subject sequences that haven't found a match.
    """
    covered_segments_query = []
    covered_segments_subject = []
    for match in matches:
        qstart = match["qstart"]
        qend = match["qend"]
        covered_segments_query += list(range(qstart, qend))
        sstart = match["sstart"]
        send = match["send"]
        covered_segments_subject += list(range(sstart, send))
    uncovered_pos_query = set(range(query_len)) - set(covered_segments_query)
    uncovered_pos_subject = set(range(subject_len)) - set(covered_segments_subject)
    uncovered_segments_query = get_segments_from_set(uncovered_pos_query, 5)
    uncovered_segments_subject = get_segments_from_set(uncovered_pos_subject, 15)
    # Get the start and end points
    qstarts, qends, sstarts, sends = [], [], [], []
    for segment in uncovered_segments_query:
        qstarts.append(segment[0])
        qends.append(segment[1])
    for segment in uncovered_segments_subject:
        sstarts.append(segment[0])
        sends.append(segment[1])
    return qstarts, qends, sstarts, sends


def recursively_align(query_seq, subject_seq_nu, start_pos_query_l, end_pos_query_l, start_pos_subject_l, end_pos_subject_l, matches):
    best_aln = None
    # For all segment combinations
    for i in range(len(start_pos_query_l)):
        start_pos_query = in_frame(start_pos_query_l[i], "start")
        end_pos_query = in_frame(end_pos_query_l[i], "end")
        for j in range(len(start_pos_subject_l)):
            start_pos_subject = in_frame(start_pos_subject_l[j], "start")
            end_pos_subject = in_frame(end_pos_subject_l[j], "end")
            # Get the cut sequences
            query_seq_segment = query_seq[start_pos_query:end_pos_query]
            subject_seq_nu_segment = subject_seq_nu[start_pos_subject:end_pos_subject]
            # Align in all 3 frames and keep the best one
            best_pval = 1
            best_psim = 0
            for frame in [0, 1, 2]:
                subject_seq_aa = subject_seq_nu_segment[frame:(frame-3)].translate(table=11)
                aln = smith_waterman(query_seq_segment, subject_seq_aa)
                pval = aln["pval"]
                psim = aln["psim"]
                if pval < best_pval:
                    best_pval = pval
                    best_psim = psim
                    best_aln = aln
                    best_aln.update({"frame": frame})
                    best_start_pos_query = start_pos_query
                    best_start_pos_subject = start_pos_subject
                # If we have the same pval, chose the match with the highest similarity
                elif pval == best_pval:
                    if psim > best_psim:
                        best_psim = psim
                        best_aln = aln
                        best_aln.update({"frame": frame})
                        best_start_pos_query = start_pos_query
                        best_start_pos_subject = start_pos_subject

    # Check the best alignment is qualitative
    if is_significant(best_aln):
        # Replace the relative positions by absolute ones
        best_aln = get_absolute_match(best_aln, best_start_pos_query, best_start_pos_subject)
        # Save the alignment
        matches.append(best_aln)
        # Get the segments that aren't covered by a match
        qstarts, qends, sstarts, sends = get_uncovered_segments(matches, len(query_seq), len(subject_seq_nu))
        recursively_align(query_seq, subject_seq_nu, qstarts, qends, sstarts, sends, matches)


def order_matches(matches):
    qstarts = np.array([m["qstart"] for m in matches])
    sorted_index = np.argsort(qstarts)
    sorted_matches = []
    for i in sorted_index:
        sorted_matches.append(matches[i])
    return sorted_matches


def look_for_frameshifts(denovo_seq, denovo_start, denovo_end, extended_match_seq, extended_start, extended_end):
    extend_left = extended_start != 0
    extend_right = extended_end != len(extended_match_seq)
    matches_left, matches_right = [], []
    # Start with the left
    if extend_left:
        left_denovo = denovo_seq[:denovo_start]
        extended_left_seq = extended_match_seq[:extended_start]
        # Get the matches recursively
        recursively_align(left_denovo, extended_left_seq, [0], [len(left_denovo)], [0], [len(extended_left_seq)], matches_left)
        # Order the matches
        matches_left = order_matches(matches_left)
    # Do the right
    if extend_right:
        right_denovo = denovo_seq[denovo_end:]
        extended_right_seq = extended_match_seq[extended_end:]
        # Get the matches recursively
        recursively_align(right_denovo, extended_right_seq, [0], [len(right_denovo)], [0], [len(extended_right_seq)], matches_right)
        matches_right = [get_absolute_match(r, denovo_end, extended_end) for r in matches_right]
        # Order the matches
        matches_right = order_matches(matches_right)
    return matches_left, matches_right


def print_results(denovo, all_matches, qcov):
    # Get the frame covers in string form
    strings = {0: ["-"] * 100, 1: ["-"] * 100, 2: ["-"] * 100}
    min_start = min([dic["sstart"] for dic in all_matches])
    max_end = max([dic["send"] for dic in all_matches])
    factor = 100 / (max_end - min_start)
    n_nucl_subject = 0
    for match in all_matches:
        frame = match["frame"]
        str_start = round((match["sstart"] - min_start) * factor)
        str_end = round((match["send"] - min_start) * factor)
        for i in range(str_start, str_end):
            strings[frame][i] = "*"
    print(f"{''.join(strings[0])}\t{denovo}\n")
    print(f"{''.join(strings[1])}\tqcov = {qcov}%\n")
    print(f"{''.join(strings[2])}\n\n")




if __name__ == "__main__":
    with open(GENOMES_LIST, "r") as f:
        genomes = f.readline().split()
    genomes = [re.sub('"', '', g) for g in genomes]

    # Get the denovo info
    denovo_dict = {}
    for genome in genomes:
        new_denovo = extract_denovo_info(genome)

        # For each denovo gene
        for denovo in new_denovo:
            """#--------------------------------------------------
            if "IOIKJFFK_00291_gene_mRNA" not in denovo:
                continue
            #--------------------------------------------------"""
            new_denovo[denovo]["genome"] = genome

            # Get the match sequence
            gene_len = len(new_denovo[denovo]["sequence"])
            loci = new_denovo[denovo]["loci"]
            contig = loci[0]
            start = loci[1]
            end = loci[2]
            strand = loci[3]
            outgroup = new_denovo[denovo]["ancestor_sp"]
            # Look if the Cter has been entirely matched
            missing_cter = gene_len - new_denovo[denovo]["qend"]
            # Look if the Nter has been entirely matched
            missing_nter = new_denovo[denovo]["qstart"]
            # Get the extended match sequence in outgroup
            match_seq_extended, extended_start, extended_end = get_extended_matched_seq(outgroup, contig, start, end, strand, missing_cter, missing_nter)
            # Add the info to the dict
            new_denovo[denovo]["extended_match_seq"] = match_seq_extended
            new_denovo[denovo]["extended_start"] = extended_start
            new_denovo[denovo]["extended_end"] = extended_end


        # Add to global dict
        denovo_dict.update(new_denovo)


    """#--------------------------------------------------------------
    # Keep only gene of interest
    dict_interest = {}
    dict_interest["IOIKJFFK_00291_gene_mRNA"] = denovo_dict["IOIKJFFK_00291_gene_mRNA"]
    denovo_dict = dict_interest
    print(denovo_dict)
    #--------------------------------------------------------------"""

    
    for denovo in denovo_dict:
        extended_match_seq = denovo_dict[denovo]["extended_match_seq"]
        extended_start = denovo_dict[denovo]["extended_start"]
        extended_end = denovo_dict[denovo]["extended_end"]
        denovo_seq = denovo_dict[denovo]["sequence"]
        denovo_start = denovo_dict[denovo]["qstart"]
        denovo_end = denovo_dict[denovo]["qend"]

        # Look for frameshift on both sides
        frameshifts_left, frameshifts_right = look_for_frameshifts(denovo_seq, denovo_start, denovo_end, extended_match_seq, extended_start, extended_end)
        # Keep gene only if there are matches
        if frameshifts_left == [] and frameshifts_right == []:
            continue
        # Get list of all matches for all frames
        origin_match = {"qstart": denovo_start, "qend": denovo_end, "sstart": extended_start, "send": extended_end, "frame": 0}
        all_matches = frameshifts_left + [origin_match] + frameshifts_right
        # Get the total qcov
        bases_covered = []
        for dic in all_matches:
            bases_covered += list(range(dic["qstart"], dic["qend"]))
        total_qcov = round((len(set(bases_covered)) / len(denovo_seq) * 100), 1)
        # Print the results
        print_results(denovo, all_matches, total_qcov)
        print("\n")
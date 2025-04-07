import re
import psa
from my_functions.genomic_functions import extract_denovo_info, get_sequence_from_loci


from my_functions.paths import GENOMES_LIST


def get_extended_matched_seq(genome, contig, start, end, strand, missing_right, missing_left):
    """
    Get the sequence (nucleotides) of the match, extended on each relevant side.
    """
    # Extend sequence on bith sides
    add_right = (missing_right * 3) * 2 if missing_right > 4 else 0
    add_left = (missing_left * 3) * 2 if missing_left > 4 else 0
    match_start = start - add_left
    match_end = end + add_right
    # Extract extended sequence
    seq = get_sequence_from_loci(genome, contig, match_start, match_end, strand)
    # Get the limits of the actual match
    limit_start = add_left
    limit_end = len(seq) - add_right
    return seq, limit_start, limit_end


def align_extensions(ref_seq, seq_f0, seq_f1, seq_f2):
    """
    Try aligning the match extensions in all 3 frames with the de novo sequence that isn't aligned with the match.
    """
    matches = []
    # Remove stops
    ref_seq = re.sub(r"\*", "", str(ref_seq))
    # Get optimal local alignment for each seq couple
    subjects = [seq_f0, seq_f1, seq_f2]
    for i in range(3):
        subject_seq = re.sub(r"\*", "", str(subjects[i]))
        # Do a smith waterman alignment
        aln = psa.water(moltype = "prot", qseq = ref_seq, sseq = subject_seq)
        # Keep only alignments with pval <10-3 and length >= 5aa
        pval = aln.pvalue()
        if pval < 10e-3 and aln.length >= 5:
            qstart, qend = aln.qstart, aln.qend
            sstart, send = aln.sstart, aln.send
            frame = i
            matches.append([frame, qstart, qend, sstart, send, pval])
    return matches


def get_uncovered_segments(extension_sequence, local_aln_matches):
    """Get the fragments longer than 5aa in the query extension that weren't aligned"""
    extension_len = len(extension_sequence)
    matches_segments = []
    # Get all the aa indexes covered by the alignments
    for match in local_aln_matches:
        qstart = match[1]
        qend = match[2]
        matches_segments += list(range((qstart - 1), qend))
    matches_segments = set(matches_segments)
    # Get all the indexes uncovered by the alignments
    unmatched_indexes = set(list(range(extension_len))) - matches_segments
    if unmatched_indexes < 5:
        return []
    # Transform indexes to continuous segments
    continuous_segments = []
    unmatched_indexes = sorted(list(unmatched_indexes))
    segment_start = unmatched_indexes[0]
    for i in range(len(unmatched_indexes - 1)):
        if unmatched_indexes[i+1] - unmatched_indexes[i] > 1:
            segment_end = unmatched_indexes[i]
            if segment_end - segment_start >= 5:
                continuous_segments.append((segment_start, segment_end))
            segment_start = unmatched_indexes[i+1]
    if unmatched_indexes[i+1] - unmatched_indexes[i] == 1:
        if unmatched_indexes[i+1] - segment_start >= 5:
            continuous_segments.append((segment_start, unmatched_indexes[i+1]))
    return continuous_segments


def recursively_align(matches, unaligned_segments, ref_seq, seq_f0, seq_f1, seq_f2):
    if unaligned_segments == []:
        return []
    unaligned_segments = get_uncovered_segments(ref_seq, matches)



def look_for_frameshifts(denovo_seq, denovo_start, denovo_end, extended_match_seq, extended_start, extended_end):
    extend_left = extended_start != 0
    extend_right = extended_end != len(extended_match_seq)
    frameshifts = {}
    # Start with the left
    if extend_left:
    # Translate in the 3 frames
        left_denovo = denovo_seq[:denovo_start]
        left_f0 = extended_match_seq[:extended_start].translate(table=11)
        left_f1 = extended_match_seq[1:extended_start+1].translate(table=11)
        left_f2 = extended_match_seq[2:extended_start+2].translate(table=11)
        # Get the matches
        matches = recursively_align(matches, left_denovo, left_f0, left_f1, left_f2)
        matches = align_extensions(left_denovo, left_f0, left_f1, left_f2)
        # See if the qcov is completely covered
        segments_left = get_uncovered_segments(left_denovo, matches)
        while segments_left != []:
            for segment in segments_left:
                left_denovo_segment = left_denovo[segment[0]:segment[1]]








    # Do the right
    if extend_right:
        right_denovo = denovo_seq[denovo_end:]
        right_f0 = extended_match_seq[extended_end:].translate(table=11)
        right_f1 = extended_match_seq[extended_end-1:-1].translate(table=11)
        right_f2 = extended_match_seq[extended_end-2:-2].translate(table=11)
        # Get the matches
        matches = align_extensions(right_denovo, right_f0, right_f1, right_f2)
        # See if the qcov is completely covered
        segments_left = get_uncovered_segments(right_denovo, matches)






if __name__ == "__main__":
    with open(GENOMES_LIST, "r") as f:
        genomes = f.readline().split()
    genomes = [re.sub('"', '', g) for g in genomes]

    # Get the denovo info
    denovo_dict = {}
    for genome in genomes:
        new_denovo = extract_denovo_info(genome)

        for denovo in new_denovo:
            if denovo == "HPMEPLIM_01176_gene_mRNA":
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
                missing_right = gene_len - new_denovo[denovo]["qend"] - 1
                # Look if the Nter has been entirely matched
                missing_left = new_denovo[denovo]["qstart"]
                # Get the extended match sequence in outgroup
                match_seq_extended, extended_start, extended_end = get_extended_matched_seq(outgroup, contig, start, end, strand, missing_right, missing_left)
                new_denovo[denovo]["extended_match_seq"] = match_seq_extended
                new_denovo[denovo]["extended_start"] = extended_start
                new_denovo[denovo]["extended_end"] = extended_end

        # Add to global dict
        denovo_dict.update(new_denovo)


    #--------------------------------------------------------------
    # Keep only gene of interest
    dict_interest = {}
    dict_interest["HPMEPLIM_01176_gene_mRNA"] = denovo_dict["HPMEPLIM_01176_gene_mRNA"]
    denovo_dict = dict_interest
    print(denovo_dict)
    #--------------------------------------------------------------

    
    for denovo in denovo_dict:
        extended_match_seq = denovo_dict[denovo]["extended_match_seq"]
        extended_start = denovo_dict[denovo]["extended_start"]
        extended_end = denovo_dict[denovo]["extended_end"]
        denovo_seq = denovo_dict[denovo]["sequence"]
        denovo_start = denovo_dict[denovo]["qstart"]
        denovo_end = denovo_dict[denovo]["qend"]

        # Look for frameshift on both sides
        frameshifts = look_for_frameshifts(denovo_seq, denovo_start, denovo_end, extended_match_seq, extended_start, extended_end)
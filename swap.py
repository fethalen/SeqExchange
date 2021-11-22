'Module for swapping bases and sequences among a set of alignments.'

import random

def pick_otus_randomly(msa, count, exclude=None):
    '''Takes a MultipleSequenceAlignment object, an integer, and a list of OTUs
    as an input. Randomly choose and return a set of OTUs from the provided
    alignment. Choose as many OTUs as the provided integer. OTUs in the exclude
    list cannot be chosen. Return an empty set if no OTU has been found.
    '''
    otus_in_msa = msa.otus()

    for otu in exclude:
        if otu in otus_in_msa:
            otus_in_msa.remove(otu)

    if len(otus_in_msa) < count:
        return set()

    return random.sample(otus_in_msa, count)

def replace_receiver_seqs(msa, receivers):
    '''Takes a MultipleSequenceAlignment object and a set of OTUs as an input.
    If a OTU in the list has more than one sequence within this alignment, then
    randomly choose another OTU's sequence and assign it to that OTU.
    '''
    receiver_count = receivers_in_msa(msa, receivers)
    affected_receivers = set()
    swap_count = 0

    for otu in receiver_count:
        if receiver_count[otu] <= 1:
            continue

        swap_count += receiver_count[otu] - 1

        # remove old sequences
        receiver_seqs = get_seqs_from_otus(msa, otu)
        for seq in random.sample(receiver_seqs, swap_count):
            msa.sequences.remove(seq)

        # assign new sequences
        affected_receivers.add(otu)
        donors = pick_otus_randomly(msa, swap_count, receivers)
        seqs = get_seqs_from_otus(msa, donors)
        remove_seqs_from_otu(msa, donors)
        msa = randomly_assign_seqs(msa, seqs, set(otu))

    print('-> swapping {} sequences in the alignment {}'.format(
        swap_count, msa))
    print('-> randomly assigning the {} sequences to {} taxa'.format(
        swap_count, len(affected_receivers)))

    return msa

def get_seqs_from_otus(msa, otus):
    '''Takes a MultipleSequenceAlignment object and a set of OTUs as an input.
    Returns a set of all the sequences that belong to that OTU.
    '''
    seqs_from_otus = set()
    for seq in msa.sequences:
        if seq.otu in otus:
            seqs_from_otus.add(seq)
    return seqs_from_otus

def remove_seqs_from_otu(msa, otus):
    '''Takes a MultipleSequenceAlignment object and a set of OTUs as an input.
    Returns a MultipleSequenceAlignment where all of the sequences that belong
    to any OTU in the provided list is removed.
    '''
    for seq in msa.sequences:
        if seq.otu in otus:
            msa.remove_sequence(seq)

def randomly_assign_seqs(msa, sequences, otus):
    '''Takes a MultipleSequenceAlignment object, a set of sequences, and set of
    OTUs as an input. The Sequence objects within the provided list are
    randomly assigned as sequences to the provided OTUs. The text string
    '_from_' and the old OTU to which the sequence originally belonged to is
    also added to the end of the sequence identifier.
    '''
    for seq in sequences:
        otu = random.sample(otus, 1)[0]
        new_description = '{}@{}_from_{}'.format(
            otu, seq.identifier, seq.otu)
        seq.description = new_description
        msa.add_sequence(seq)
    return msa

def receivers_in_msa(msa, receivers):
    '''Takes a MultipleSequenceAlignment object and a set of OTUs as an input.
    Returns a dictionary where each receiver in the provided set of OTUs act as
    keys and the respective value is how many sequences that OTU have in the
    provided alignment.
    '''
    receiver_count = {otu : 0 for otu in receivers}

    for seq in msa.sequences:
        if seq.otu in receivers:
            receiver_count[seq.otu] += 1

    return receiver_count

def cross_contaminate(msa, receivers, likelihood, max_substitutions):
    '''Takes a MultipleSequenceAlignment object, a set of OTUs, a floating
    point indicating likelihood and an intiger that corresponds to the
    maximally allowed number of substitutions. For each receiver that is
    present within the alignment, there is <likelihood> chance that one of
    the receiver's sequences becomes a cross-contaminant: First, a copy of the
    sequence is created, then up to <max_substitutions> bases are rearranged
    within the sequence and the sequence is then assigned to another OTU that
    is already present within the alignment. A sequence from that OTU is then
    selected by random and deleted. Returns the MultipleSequenceAlignment,
    containing one, or more, newly generated cross-contaminant sequences.
    '''
    receiver_count = receivers_in_msa(msa, receivers)
    contamination_count = 0

    for otu in receiver_count:
        if receiver_count[otu] == 0:
            continue

        if bool(random.randrange(100) > (likelihood * 100)):
            continue

        contamination_count += 1

        # randomly select one of the OTU's sequences and make a copy
        otu_seqs = list(get_seqs_from_otus(msa, otu))
        original_seq = random.sample(otu_seqs, 1)[0]
        contaminated_seq, substitutions = add_noise(
            original_seq.sequence_data, max_substitutions)
        contaminant = random.sample(msa.sequences, 1)[0]

        # select a contaminant, remove one of it's sequences and use the OTU
        # for the contamination
        contaminant_identifier = \
            original_seq.identifier + '_contaminated_' + otu + '_' + \
            str(substitutions) + '_subs'
        contaminant_description = contaminant.otu + '@' + contaminant_identifier
        msa.add_sequence(
            None, contaminant_description, contaminated_seq)
        msa.remove_sequence(contaminant)

    print('-> contaminated {} sequence{} in {}'.format(
        contamination_count,
        '' if contamination_count == 1 else 's',
        msa))
    return msa

def add_noise(sequence, max_substitutions):
    '''Takes a string of sequence data (nucleotides or amino acids) and an
    integer as an input. Randomly pick up to <max_substitutions> positions that
    receives another positions' nucleotide or amino acid (or whatever
    character happens to be at that position). Returns the newly generated
    sequence and the number of subtitutions that were added.
    '''
    substitutions = random.randint(1, max_substitutions)
    positions = random.sample(range(1, len(sequence)), substitutions)
    replacements = random.sample(range(1, len(sequence)), substitutions)

    sequence = list(sequence)
    for index, position in enumerate(positions):
        sequence[position] = sequence[replacements[index]]

    return ''.join([position for position in sequence]), substitutions

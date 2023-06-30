import pandas as pd


def is_int(string):
    try:
        int(string)
        return True
    except ValueError:
        return False


def breaks_a_number(prefix, string):
    if len(string) == len(prefix):
        return False
    else:
        boundary = string[len(prefix) - 1:len(prefix) + 1]
        return is_int(boundary)


def acceptable_prefix(prefix, a, b):
    # Short prefixes are no good
    if len(prefix) < 3:
        return False

    # Check if prefix splits a number, reject it then
    if breaks_a_number(prefix, a) or breaks_a_number(prefix, b):
        return False

    return True


def longest_common_prefix(a, b):
    if len(b) < len(a):
        a, b = b, a

    for i in reversed(range(len(a) + 1)):
        prefix = a[:i]
        if b.startswith(prefix) and acceptable_prefix(prefix, a, b):
            return prefix

    return ''


def longest_common_prefixes_sequence(sequence):
    sequence = sorted(sequence)
    if len(sequence) < 2:
        return sequence

    prefixes = set()

    for a, b in zip(sequence, sequence[1:]):
        prefix = longest_common_prefix(a, b)
        if prefix:
            prefixes.add(prefix)

    for item in sequence:

        has_prefix = False
        for prefix in prefixes:
            if item.startswith(prefix):
                has_prefix = True
                break

        if not has_prefix:
            prefixes.add(item)

    if sorted(prefixes) == sequence:
        return prefixes
    else:
        return longest_common_prefixes_sequence(prefixes)


def concatenate_int_suffixes(suffixes):
    try:
        int_suffixes = [int(s) for s in suffixes]
    except ValueError:
        return suffixes

    int_suffixes = sorted(int_suffixes)

    previous = int_suffixes[0]
    min_ = int_suffixes[0]

    concatenated_suffixes = []
    for s in int_suffixes[1:]:
        if s == previous + 1:
            previous = s
        else:
            if previous == min_:
                concatenated_suffixes.append(str(previous))
            elif previous == min_ + 1:
                concatenated_suffixes.extend([str(min_), str(previous)])
            else:
                concatenated_suffixes.append('{}-{}'.format(min_, previous))

            previous = s
            min_ = s

    if previous == min_:
        concatenated_suffixes.append(str(previous))
    elif previous == min_ + 1:
        concatenated_suffixes.extend([str(min_), str(previous)])
    else:
        concatenated_suffixes.append('{}-{}'.format(min_, previous))

    return concatenated_suffixes


def create_label(gene_names):

    # Used to shorten certain troublesome gene heatmap_header_labels
    MANUAL_ALIAS_MAPPINGS = {
        'HIST[1H2AG,1H2AH,1H2AI,1H2AJ,1H2AK,1H2AL,1H2AM,2H2AA3,2H2AA4,2H2AC]/H2AFJ': 'HIST1H2A[G-M]/HIST2H2A[A3,A4,C]/H2AFJ',
        'HIST[1H4A,1H4B,1H4C,1H4E,1H4F,1H4I,1H4J,1H4L,2H4A,4H4]': 'HIST1H4[A-C,E,F,I,J,L]/HIST2H4A/HIST4H4',
        'HIST[2H3A,2H3C,2H3D,2H3PS2,3H3]': 'HIST2H3[A,C,D,PS2]/HIST3H3',
        'HIST[1H2BH,1H2BL,2H2BF]': 'HIST1H2B[H,L]/HIST2H2BF',
        'HIST[1H2AB,1H2AE,3H2A]': 'HIST1H2A[B,E]/HIST3H2A',
        'HIST1H1D/HIST1H1E': 'HIST1H1[D,E]',
        'HIST1H3[A,B,C,D,E,F,G,H,I,J]': 'HIST1H3[A-J]'
    }

    prefixes = longest_common_prefixes_sequence(gene_names)

    # Sort prefixes by the first occurrence, then alphabetically by matching gene name:
    def prefix_sortkey(p):
        for i, gn in enumerate(gene_names):
            if gn.startswith(p):
                return i, gn

        raise Exception('Prefix not found in gene names list')

    prefixes = sorted(prefixes, key=prefix_sortkey)

    components = []
    for prefix in prefixes:
        matching_gene_names = [gn for gn in gene_names if gn.startswith(prefix)]
        suffixes = [gn[len(prefix):] for gn in matching_gene_names]

        suffix_threshold = 3
        if '' in suffixes:
            components.append(prefix)
            # Remove the empty suffix from suffix and corresponding matching gene name
            suffixes = list(filter(lambda x: x != '', suffixes))
            suffix_threshold -= 1 # If we had the empty protein suffix, subtract one from threshold
        suffixes = sorted(suffixes)
        if not suffixes:
            continue
        elif len(suffixes) < suffix_threshold:
            components.extend(['{}{}'.format(prefix, s) for s in suffixes])
        else:
            components.append('{}[{}]'.format(prefix,
                                              ','.join(concatenate_int_suffixes(suffixes))))

    alias = '/'.join(components)
    alias = MANUAL_ALIAS_MAPPINGS.get(alias, alias)
    return alias


def compile_protein_id_map(series, sep=';'):

    groups = []
    seen = set()
    for ix in series:
        for protein in ix.split(sep):
            if protein in seen:
                raise Exception('{} occurs in series more than once'.format(protein))

            seen.add(protein)
            groups.append([protein, ix])

    groups = pd.DataFrame(groups, columns=['Protein ID', series.name])
    return groups


def assign_protein_labels(protein_meta):

    names = protein_meta['Gene names']

    names_split = names.fillna('').apply(lambda x: x.split(';'))

    labels = names_split.apply(create_label)

    # Make empty heatmap_header_labels to be named after the index
    labels.loc[labels == ''] = protein_meta.loc[labels == ''].index.str.replace(';', '/')

    # Get all duplicated heatmap_header_labels
    duplicated = labels.apply(lambda x: x != '' and (labels == x).sum() > 1)

    duplicated_labels = labels[duplicated].unique()

    for duplicate in duplicated_labels:
        subproteins = protein_meta[labels == duplicate]

        # Assign unique identifiers for heatmap_header_labels, the highest score protein getting the first number
        subproteins = subproteins.sort_values(by='Score', ascending=False)

        new_labels = ['{} ({})'.format(duplicate, i) for i in range(1, len(subproteins)+1)]
        labels.loc[subproteins.index] = new_labels

    labels.name = 'Gene label'

    return labels
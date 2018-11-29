import sys
import operator


class Node:
    def __init__(self, content):
        self.content = content
        self.next = None


class LinkedList:
    def __init__(self):
        self.head = None
        self.tail = None


pep_mass = [57, 71, 87, 97, 99, 101, 103, 113, 114, 115, 128, 129, 131, 137, 147, 156, 163, 186]


def create_cycle(peptide):
    l_list = LinkedList()
    first = True
    current_node = None
    for acid in peptide:
        if first:
            first = False
            current_node = Node(acid)
            l_list.head = current_node
        else:
            current_node.next = Node(acid)
            current_node = current_node.next
    current_node.next = l_list.head
    return l_list


def get_cycle_by_length(acid, length):
    total = 0
    current_node = acid
    for i in range(length):
        total += current_node.content
        current_node = current_node.next
    return total


def find_cyclospectrum(peptide):
    cyclospectrum = [0]
    pep_length = len(peptide)
    peptide_list = create_cycle(peptide)
    cyclospectrum.append(get_cycle_by_length(peptide_list.head, pep_length))
    current_node = peptide_list.head
    for j in range(pep_length):
        i = pep_length - 1
        while i != 0:
            cyclospectrum.append(get_cycle_by_length(current_node, i))
            i -= 1
        current_node = current_node.next
    cyclospectrum.sort()
    return cyclospectrum


def create_linear_cycle(peptide):
    l_list = LinkedList()
    first = True
    current_node = None
    for acid in peptide:
        if first:
            first = False
            current_node = Node(acid)
            l_list.head = current_node
        else:
            current_node.next = Node(acid)
            current_node = current_node.next
    return l_list


def find_linear_mass(acid, length):
    total = 0
    current_node = acid
    for i in range(length):
        if current_node is None:
            return None
        else:
            total += current_node.content
            current_node = current_node.next
    return total


def linear_spectrum(peptide):
    l_specturm = [0]
    pep_length = len(peptide)
    linear_list = create_linear_cycle(peptide)
    l_specturm.append(find_linear_mass(linear_list.head, pep_length))
    current_node = linear_list.head
    while current_node is not None:
        i = pep_length - 1
        while i != 0:
            total = find_linear_mass(current_node, i)
            if total is not None:
                l_specturm.append(total)
            i -= 1
        current_node = current_node.next
    l_specturm.sort()
    return l_specturm


def parent_mass(specturm):
    top_val = max(specturm)
    return top_val


def mass(peptide):
    total = 0
    for x in peptide:
        total += x
    return total


def expand(peptides, mass_list):
    new_list = []
    if len(peptides) == 0:
        for mass in mass_list:
            new_list.append([mass])
        return new_list
    for pep in peptides:
        for x in mass_list:
            pos_set = pep.copy()
            pos_set.append(x)
            new_list.append(pos_set)
    return new_list


def trim(leaderboard, spectrum, n):
    if len(leaderboard) == 0:
        return []
    scores = []
    i = 0
    for pep in leaderboard:
        scores.append((i,score(pep, spectrum)))
        i += 1
    scores.sort(key=lambda tup: tup[1], reverse=True)
    new_list = []
    tie_score = None
    sc_count = 0
    for sc_tup in scores:
        index = sc_tup[0]
        sc = sc_tup[1]
        if sc_count < n:
            new_list.append(leaderboard[index])
            sc_count += 1
            tie_score = sc
        elif sc_count >= n:
            if tie_score == sc:
                new_list.append(leaderboard[index])
                sc_count += 1
            else:
                return new_list
    return new_list


def leaderboard_cyclopetide_sequencing(spectrum, mass_list, n):
    leaderboard = []
    leader_peptide = None
    empty_list = False
    while not empty_list:
        next_list = []
        leaderboard = expand(leaderboard, mass_list)
        for peptide in leaderboard:
            if peptide == [99, 71, 137, 57, 72, 57]:
                print('stop here')
            if mass(peptide) == parent_mass(spectrum):
                if score(peptide, spectrum) > score(leader_peptide, spectrum):
                    leader_peptide = peptide
            elif mass(peptide) < parent_mass(spectrum):
                next_list.append(peptide)
        leaderboard = trim(next_list, spectrum, n)
        if len(leaderboard) == 0:
            empty_list = True
    return leader_peptide


def score(peptide, exp_spec):
    if peptide is None:
        return 0
    exp_copy = exp_spec.copy()
    pep_spec = find_cyclospectrum(peptide)
    score = 0
    for acid in pep_spec:
        if acid in exp_copy:
            score += 1
            exp_copy.remove(acid)
    return score


def get_data(filename):
    with open(filename, 'r') as file:
        spec = file.readline().strip().split(' ')
    spectrum = []
    for x in spec:
        spectrum.append(int(x))
    return spectrum


def spectrum_convolution(spectrum):
    count_dict = {}
    for i in range(1, len(spectrum)):
        for j in range(i):
            if not spectrum[i] == spectrum[j]:
                difference = abs(spectrum[i] - spectrum[j])
                if difference in count_dict.keys():
                    count_dict[difference] += 1
                else:
                    count_dict[difference] = 1
    sorted_dict = sorted(count_dict.items(), key=operator.itemgetter(1))
    return sorted_dict[::-1]


def filter_convolution(convolution, m):
    filter = []
    count = 0
    tie_count = None
    for convo in convolution:
        val = convo[0]
        val_count = convo[1]
        if count < m:
            if 57 <= val <= 200:
                filter.append(val)
                count += 1
                tie_count = val_count
        elif count >= m:
            if val_count == tie_count:
                if 57 <= val <= 200:
                    filter.append(val)
                    count += 1
            else:
                return filter
    return filter


def convolution_cyclopeptide_sequencing(spectrum, m, n):
    spectrum.sort()
    spec_conv = spectrum_convolution(spectrum)
    filter_conv = filter_convolution(spec_conv, m)
    results = leaderboard_cyclopetide_sequencing(spectrum, filter_conv, n)
    return results


def write_out(results):
    with open('output.txt', 'w') as file:
        for res in results:
            file.write(str(res) + '-')
    print(results)


if __name__ == '__main__':
    spectrum = get_data(sys.argv[1])
    m = int(sys.argv[2])
    n = int(sys.argv[3])
    result = convolution_cyclopeptide_sequencing(spectrum, m, n)
    write_out(result)


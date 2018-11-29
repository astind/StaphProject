import sys
import random


mass_table = {
    'G': 57,
    'A': 71,
    'S': 87,
    'P': 97,
    'V': 99,
    'T': 101,
    'C': 103,
    'I': 113,
    'L': 113,
    'N': 114,
    'D': 115,
    'K': 128,
    'Q': 128,
    'E': 129,
    'M': 131,
    'H': 137,
    'F': 147,
    'R': 156,
    'Y': 163,
    'W': 186
}


m_to_p = {
    57: 'G',
    71: 'A',
    87: 'S',
    97: 'P',
    99: 'V',
    101: 'T',
    103: 'C',
    114: 'N',
    115: 'D',
    129: 'E',
    131: 'M',
    137: 'H',
    147: 'F',
    156: 'R',
    163: 'Y',
    186: 'W'
}


def mass_to_pep(masses):
    pep = ''
    for mass in masses:
        if mass == 113:
            if random.random() < 0.5:
                pep += 'I'
            else:
                pep += 'L'
        elif mass == 128:
            if random.random() < 0.5:
                pep += 'K'
            else:
                pep += 'Q'
        else:
            pep += m_to_p[mass]
    return pep


pep_mass = [57, 71, 87, 97, 99, 101, 103, 113, 114, 115, 128, 129, 131, 137, 147, 156, 163, 186]


class Node:
    def __init__(self, content):
        self.content = content
        self.next = None


class LinkedList:
    def __init__(self):
        self.head = None
        self.tail = None


def get_data(filename):
    with open(filename, 'r') as file:
        spectrum = file.readline().strip()
    spectrum = spectrum.split(' ')
    spec_int = []
    for x in spectrum:
        spec_int.append(int(x))
    return spec_int


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


def expand(peptides):
    new_list = []
    if len(peptides) == 0:
        for mass in pep_mass:
            new_list.append([mass])
        return new_list
    for pep in peptides:
        for x in pep_mass:
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
    range_length = min(n, len(leaderboard))
    last_score, last_index = None, None
    for j in range(range_length):
        new_list.append(leaderboard[scores[j][0]])
        if j == (range_length - 1):
            last_score = scores[j][1]
            last_index = j
    if last_index < (len(leaderboard) - 1):
        if scores[last_index + 1][1] == last_score:
            last_index += 1
            while scores[last_index][1] == last_score:
                new_list.append(leaderboard[scores[last_index][0]])
                last_index += 1
                if last_index == len(leaderboard):
                    break
    return new_list


def leaderboard_cyclopetide_sequencing(spectrum, n):
    leaderboard = []
    leader_peptide = None
    empty_list = False
    while not empty_list:
        next_list = []
        leaderboard = expand(leaderboard)
        for peptide in leaderboard:
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
    pep_spec = linear_spectrum(peptide)
    score = 0
    for acid in pep_spec:
        if acid in exp_copy:
            score += 1
            exp_copy.remove(acid)
    return score


def write_output(answers):
    with open('output.txt', 'w') as file:
        for ans in answers:
            file.write(str(ans) + '-')
    print(answers)


if __name__ == '__main__':
    # N, spectrum = get_data(sys.argv[1])
    spectrum = get_data(sys.argv[1])
    n = int(sys.argv[2])
    ans = leaderboard_cyclopetide_sequencing(spectrum, n)
    # write_output(ans)
    print(mass_to_pep(ans))


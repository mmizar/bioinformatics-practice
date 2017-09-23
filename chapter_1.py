# Implementation of algorithms desribed in the book
# Bioinformatics Algorithms: An Active Learning Approach
# by Phillip Compeau & Pavel Pevzner

import sys
import math

# Count the frequency of a pattern within a sequence
def pattern_count(text, pattern):
    count = 0
    last_index = 0
    for i in range (0, len(text) - len(pattern) + 1):
        if text[i:i + len(pattern)] == pattern: # and i >= last_index + len(pattern):
            count += 1
            last_index = i
    return count

# Find the most frequent k-mers in a string
def frequent_words(text, k):
    frequent_patterns = set()
    count = []
    for i in range (0, len(text) - k):
        pattern = text[i:i + k]
        count.insert(i, pattern_count(text, pattern))
    max_count = max(count)
    for i in range (0, len(text) - k):
        if count[i] == max_count:
            frequent_patterns.add(text[i:i + k])
    return frequent_patterns

# Find the reverse compliment of a DNA string
def reverse_compliment(pattern):
    reverse = ""
    for nucleotide in pattern:
        reverse += get_compliment(nucleotide)
    return reverse[::-1]

def get_compliment(nucleotide):
    if nucleotide == "C":
        return "G"
    if nucleotide == "G":
        return "C"
    if nucleotide == "A":
        return "T"
    if nucleotide == "T":
        return "A"
    else: return ""

# Find all the locations of a pattern in a genome
def pattern_matching(pattern, genome):
    locations = []
    for i in range (0, len(genome) - len(pattern) + 1):
        if genome[i:i + len(pattern)] == pattern:
            locations.append(i)
    return locations

def pattern_frequencies(text, k):
    count = {}
    for i in range (0, len(text) - k):
        pattern = text[i:i + k]
        count[pattern] = pattern_count(text, pattern)
    return count

def computing_frequencies(text, k):
    count = []
    for i in range (0, int(math.pow(4, k))):
        count.insert(i, 0)
    for i in range (0, len(text) - k):
        pattern = text[i:i + k]
        count[pattern_to_number(pattern)] = pattern_count(text, pattern)
    return count

# Find patterns forming clumps in a string
def clump_finding(genome, k, l, t):
    frequent_patterns = set()
    clump = []
    for i in range (0, int(math.pow(4, k))):
        clump.insert(i, 0)
    frequencies = {}
    text = genome[0:l]
    frequencies = pattern_frequencies(text, k)
    for i in range (0, int(math.pow(4, k))):
        if frequencies.get(number_to_pattern(i, k)) >= t:
            clump[i] = 1
    for i in range (1, len(genome) - l + 1):
        first = genome[i - 1:i - 1 + k]
        last = genome[i + l - k:i + l]
        if first in frequencies:
            frequencies[first] -= 1
        if last in frequencies:
            frequencies[last] += 1
        else:
            frequencies[last] = 1
        if frequencies.get(last) >= t:
            clump[pattern_to_number(last)] = 1
    for i in range (0, int(math.pow(4, k))):
        if clump[i] == 1:
            frequent_patterns.add(number_to_pattern(i, k))
    print(number_to_pattern(pattern_to_number("TAGCCTAACAAG"), 12))
    print(pattern_to_number("TAGCCTAACAAG"))
    print(pattern_to_number("TAGCCTAACAAA"))
    return frequent_patterns

def pattern_to_number(pattern):
    total = int(math.pow(4, len(pattern)))
    chunks = total / 4
    number = 0
    for i in range (0, len(pattern)):
        if pattern[i] == 'A':
            number += chunks * 0
        elif pattern[i] == 'C':
            number += chunks * 1
        elif pattern[i] == 'G':
            number += chunks * 2
        elif pattern[i] == 'T':
            number += chunks * 3
        chunks /= 4
    return number

def number_to_pattern(number, k):
    pattern = ""
    for i in range (0, k):
        remainder = number % 4
        if remainder == 0:
            pattern = "A" + pattern
        elif remainder == 1:
            pattern = "C" + pattern
        elif remainder == 2:
            pattern = "G" + pattern
        elif remainder == 3:
            pattern = "T" + pattern
        number /= 4
    return pattern

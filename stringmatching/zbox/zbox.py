#!/usr/bin/env python

from argparse import ArgumentParser
import copy

def main(pattern, text):
    r = 0
    l = 0

    # Here is the basic algorithm:
    # march through the text one letter at a time
    # for each position we need to find the value "Z",
    # which is the number of letters of prefix in "pattern"
    # that matches "text" at the given position.
    # If Z is len(pattern), it is an exact match.
    # At any given time, r is the rightmost known endpoint
    # of any Z box, and l is the corresponding
    # starting point
    i = 0
    Z = [None] * len(text)
    matches = set()
    while i < len(text):
        # Search for a new pattern start
        if text[i] == pattern[0]:
            matches.add(i)
        else:
            Z[i] = 0
        for match in copy.copy(matches):
            if pattern[i - match] != text[i]:
                Z[match] = i - match
                matches.remove(match)
            if i - match + 1 >= len(pattern) or i + 1 >= len(text):
                Z[match] = i - match + 1
                matches.remove(match)
        i += 1
    print Z

if __name__ == '__main__':
    parser = ArgumentParser('zbox')
    parser.add_argument('pattern')
    parser.add_argument('text')
    args = parser.parse_args()
    main(args.pattern, args.text)

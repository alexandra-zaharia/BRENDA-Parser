# -*- coding: utf-8 -*-


"""
=============================
BRENDA Enzyme Database Parser
=============================

:Author:
    Moritz Emanuel Beber
    Alexandra Zaharia
:Date:
    2011-01-27
    2018-12-08
:Copyright:
    Copyright |c| 2011, Jacobs University Bremen gGmbH, all rights reserved.
    Copyright |c| 2018, Alexandra Zaharia, all rights reserved.
:File:
    parser.py

.. |c| unicode:: U+A9
"""

__all__ = ["BRENDAParser"]

import re
import codecs

from collections import defaultdict
from recordclass import recordclass

from brenda.utils import ArgumentError, ProgressMeter, is_ec_number, has_ec_number, \
    replace_abnormal_comment, init_tags, find_parentheses_indexes


class Enzyme:
    """
    An object that encompasses all information about a kind of enzyme uniquely
    identified by its EC number.
    """

    def __init__(self, ec_number, comment):
        """Initializes an Enzyme instance."""
        self.ec_number = ec_number
        self.comment = comment
        self.proteins = dict()
        self.references = dict()
        self.entries = dict()

    def __str__(self):
        return self.ec_number

    def __repr__(self):
        return self.ec_number


Current = recordclass(
    'Current', ['proteins', 'comment', 'information', 'references', 'ec_number', 'line_number'])


class EntryComment:
    """Encapsulates a comment to an entry in a BRENDA information field."""

    def __init__(self, message, proteins=None, references=None):
        """Initializes an EntryComment instance."""
        self.msg = message
        self.proteins = proteins
        self.references = references

    def __str__(self):
        return self.msg

    def __repr__(self):
        return self.msg


class Entry(EntryComment):
    """Encapsulates an entry in a BRENDA information field."""

    def __init__(self, message, current):
        """Initializes an Entry instance."""
        if not isinstance(current, Current):
            raise ArgumentError('Entry: expected type Current')
        super().__init__(message=message, proteins=current.proteins, references=current.references)
        self.information = current.information
        self.comment = current.comment


class Protein:
    """Encapsulates an entry in a BRENDA PROTEIN (PR) field."""

    _counter = 1

    def __init__(self, organism, current):
        """Initializes a Protein instance."""
        if not isinstance(current, Current):
            raise ArgumentError('Protein: expected type Current')
        self._index = self._counter
        self.__class__._counter += 1
        self.organism = organism
        self.identifiers = current.proteins
        self.references = current.references
        self.information = current.information
        self.comment = current.comment

    def __str__(self):
        return self.organism

    def __repr__(self):
        return '<%s.%s, %d>' % (self.__module__, self.__class__.__name__, id(self))


class BRENDAParser:
    """Encapsulates the parsing of a BRENDA database plain text file."""

    _sections = {
        'ACTIVATING_COMPOUND': 'AC',
        'APPLICATION': 'AP',
        'COFACTOR': 'CF',
        'CLONED': 'CL',
        'CRYSTALLIZATION': 'CR',
        'ENGINEERING': 'EN',
        'GENERAL_STABILITY': 'GS',
        'IC50_VALUE': 'IC5',
        'INHIBITORS': 'IN',
        'KI_VALUE': 'KI',
        'KM_VALUE': 'KM',
        'LOCALIZATION': 'LO',
        'METALS_IONS': 'ME',
        'MOLECULAR_WEIGHT': 'MW',
        'NATURAL_SUBSTRATE_PRODUCT': 'NSP',
        'OXIDATION_STABILITY': 'OS',
        'ORGANIC_SOLVENT_STABILITY': 'OSS',
        'PH_OPTIMUM': 'PHO',
        'PH_RANGE': 'PHR',
        'PH_STABILITY': 'PHS',
        'PI_VALUE': 'PI',
        'POSTTRANSLATIONAL_MODIFICATION': 'PM',
        'PROTEIN': 'PR',
        'PURIFICATION': 'PU',
        'REACTION': 'RE',
        'REFERENCE': 'RF',
        'RENATURED': 'REN',
        'RECOMMENDED_NAME': 'RN',
        'REACTION_TYPE': 'RT',
        'SPECIFIC_ACTIVITY': 'SA',
        'SYSTEMATIC_NAME': 'SN',
        'SUBSTRATE_PRODUCT': 'SP',
        'STORAGE_STABILITY': 'SS',
        'SOURCE_TISSUE': 'ST',
        'SUBUNITS': 'SU',
        'SYNONYMS': 'SY',
        'TURNOVER_NUMBER': 'TN',
        'TEMPERATURE_OPTIMUM': 'TO',
        'TEMPERATURE_RANGE': 'TR',
        'TEMPERATURE_STABILITY': 'TS'}

    def __init__(self, filename, encoding='utf8'):
        """Initializes a BRENDAParser instance."""
        object.__init__(self)
        self._filename = filename
        self._file_handle = None
        self._encoding = encoding
        self._progress = None

        self._tags = init_tags()  # compiled regex's

        # Current information represents current line number, current EC number
        # (Enzyme instance), and current entry information (proteins,
        # information, comment, and references).
        self._current = Current(None, None, None, None, None, None)

        self._skip = False  # skip to next EC number?
        self.enzymes = None  # dict of EC numbers

    def __enter__(self):
        """Opens file and initializes progress meter."""
        self._file_handle = codecs.open(self._filename, mode='rb', encoding=self._encoding)
        tmp = self._file_handle.readlines()
        self._file_handle.close()
        self._file_handle = tmp
        self._progress = ProgressMeter('Parsing flat file', len(tmp))
        self.enzymes = defaultdict(list)
        self._current.line_number = 0
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        """Closes the file handle."""
        if not isinstance(self._file_handle, list):
            self._file_handle.close()
        return False

    def _reset_parser(self):
        """Resets parser fields that are used for information, comment,
        protein and reference extraction.
        """
        self._current.proteins = None
        self._current.comment = None
        self._current.information = None
        self._current.references = None

    @staticmethod
    def is_section_redundant(section_id):
        """Returns True for redundant section IDs (already present in Enzyme),
        currently PROTEIN and REFERENCE.
        """
        return str(section_id) == 'PROTEIN' or str(section_id) == 'REFERENCE'

    def parse(self):
        """Parses multiple Enzyme sections.

        :return: dict of Enzyme objects
        """
        section_name = ''  # long section identifier, e.g. 'PROTEIN'
        section_contents = list()  # contents of the section identified by section_name
        short_entry = ''  # two- or three- letter section identifier, e.g. 'PR' for 'PROTEIN'
        entry = list()  # contents of an entry identified by short_entry
        parser = self._parse_generic_entry

        for line in self._file_handle:
            if self._current.line_number % 1000 == 0:
                self._progress.update(self._current.line_number)
            self._current.line_number += 1

            line = line.rstrip()
            if not line or line.startswith('*'):
                continue

            content = line.split(None, 1)
            if content[0] == 'ID' and has_ec_number(content[1]):
                self._parse_id(content[1])
            elif content[0] in self._sections.keys():  # handle new section
                # Finish handling previous section
                if content[1:]:  # not a new section, actually
                    entry.append(line.lstrip())
                    continue
                if entry:
                    section_contents.append(parser(' '.join(entry)))
                if section_contents and not self.is_section_redundant(section_name):
                    if self._current.ec_number is None:  # skip to next EC due to invalid ID
                        self._skip = True
                        continue
                    self._current.ec_number.entries[section_name] = section_contents

                # Prepare to process current section
                section_contents = list()
                entry = list()
                section_name = content[0]
                parser = self._determine_parser_from_section_name(section_name)
                short_entry = self._sections.get(section_name, False)
                if not short_entry:
                    raise ArgumentError('Unrecognised entry: \'{}\' @ #%{}'
                                        .format(line, self._current.line_number))
            elif content[0] == short_entry:  # handle previous and current entries
                if entry:
                    section_contents.append(parser(' '.join(entry)))
                entry = content[1:]
            elif content[0] == '///':  # handle end of EC number description
                if self._skip:
                    self._skip = False
                    continue
                # end one enzyme entry
                if entry:
                    section_contents.append(parser(' '.join(entry)))
                if section_contents and not self.is_section_redundant(section_name):
                    self._current.ec_number.entries[section_name] = section_contents
                self._current.ec_number = None
            else:
                entry.append(line.lstrip())
        # convert to normal dictionary again
        res = dict(self.enzymes)
        self._progress.close()
        return res

    def _determine_parser_from_section_name(self, section_name):
        """Returns the appropriate parser depending on the current section.

        :param section_name: name of the current section
        :return: reference to parsing function
        """
        if section_name == 'PROTEIN':
            return self._parse_protein
        if section_name == 'REFERENCE':
            return self._parse_reference
        return self._parse_generic_entry

    def _parse_generic_entry(self, text):
        """Parses an entry of a specific information field.

        :param text: generic entry in BRENDA flat file
        :return: an Entry instance corresponding to the parsed text
        """
        self._reset_parser()
        text = self.extract_information(text)
        text = self.extract_proteins(text)
        text = self.extract_references(text)
        text = self.extract_comment(text)
        return Entry(text.strip(), self._current)

    def extract_information(self, text):
        """Extracts and stores information and returns the text resulting from
        removing the information field.

        The information format is defined by the instance's _tags.information
        attribute. If information tags occur more than once in the text, only
        the last occurrence is considered to be valid information (the previous
        ones might belong to the reaction or to comments). If information is
        present in the text, it is stored in the instance's _current.information
        attribute.

        :param text: text that may contain information
        :return: text resulting from information extraction
        """
        mobj = self._tags.information.search(text)

        if mobj:
            occurrences = [x for x in re.finditer(self._tags.information, text)]
            start, end = (occurrences[-1].start(), occurrences[-1].end()) if len(occurrences) > 1\
                else (mobj.start(), mobj.end())
            self._current.information = text[(start + 1):(end - 1)]
            if not self._current.information.strip():
                self._current.information = None
            text = text[:start] + text[end:]
        else:
            self._current.information = None

        return text.strip()

    @staticmethod
    def has_protein_field_structure(text):
        """Determines whether the given text is a potential proteins field.

        :param text: text that may represent a proteins field
        :return: True if text contains protein references, False otherwise
        """
        pattern = re.compile(r'#(\d)+(,*(\s)*\d+)*#')
        return pattern.match(text)

    def _clean_extra_hash_characters(self, text):
        """If extra hash ('#') characters are present in the given text, remove
        them and return the resulting text.

        Hash characters delimit protein references. If they occur an odd number
        of times, the extra occurrence needs to be identified and suppressed.

        :param text: text that may contain extra hash characters to remove
        :return: text resulting from the suppression of the extra occurrence of
            the hash character
        """
        if '#' in text:
            hashes = [x.start() for x in re.finditer(r'#', text)]
            if len(hashes) == 1:
                text = text.replace('#', '')
            elif len(hashes) % 2:
                for i, _ in enumerate(hashes):
                    if i == len(hashes) - 1:
                        break
                    subtext = text[hashes[i]:(hashes[i+1] + 1)]
                    if not self.has_protein_field_structure(subtext):
                        text = text[:hashes[i+1]] + text[(hashes[i+1] + 1):]
                        break
        return text.strip()

    def extract_proteins(self, text):
        """Extracts and stores protein references and returns the text resulting
        from removing the proteins field.

        The proteins field format is defined by the instance's _tags.protein
        attribute. It is assumed that the proteins field occurs at the beginning
        of the text (since this method is used to extract protein fields from
        reaction-type entries). If a proteins field is present in the text, the
        protein references are stored in the instance's _current.proteins
        attribute.

        :param text: text that may contain protein references (at the beginning)
        :return: text resulting from the extraction of protein references
        """
        if text.startswith('#'):
            text = self._clean_extra_hash_characters(text)
            pobj = self._tags.protein.search(text)
            if not pobj or pobj.start() != 0:
                raise ArgumentError('Protein reference missing: \'{}\' @ #%{}'
                                    .format(text, self._current.line_number))
            _, self._current.proteins = self._extract_numbers(pobj.group(), self._tags.protein)
            if self._current.proteins is not None:
                text = text[pobj.end():]
        return text.strip()

    def extract_references(self, text):
        """Extracts and stores literature references and returns the text
        resulting from removing the references field.

        The references field format is defined by the instance's _tags.reference
        attribute. It is assumed that the references field occurs at the end of
        the text (since this method is used to extract reference fields from
        reaction-type entries). If a references field is present in the text,
        the literature references are stored in the instance's
        _current.references attribute.

        :param text: text that may contain literature references (at the end)
        :return: text resulting from the extraction of literature references
        """
        if text.endswith('>'):
            if self._tags.reference.search(text):
                ref_matches = [x for x in self._tags.reference.finditer(text)]
                if ref_matches and ref_matches[-1].end() == len(text):
                    _, self._current.references = \
                        self._extract_numbers(ref_matches[-1].group(), self._tags.reference)
                if self._current.references is not None:
                    text = text[:ref_matches[-1].start()]
        return text.strip()

    @staticmethod
    def _clean_extra_pipe_characters(text):
        """If extra pipe ('|') characters are present in the given text, remove
        them and return the resulting text.

        Pipe characters delimit (not necessarily additional) comments in BRENDA
        flat files. If these pipe characters are present in odd number, the
        extra occurrences need to be removed prior to comment fusion and
        extraction.

        :param text: text that may contain extra pipe characters to remove
        :return: text resulting from the suppression of the first occurrence of
            a pipe character, if there is an odd number of occurrences
        """
        if '|' in text:
            pipes = [x.start() for x in re.finditer(r'\|', text)]
            if len(pipes) == 1:
                text = text.replace('|', '')
            elif len(pipes) % 2:  # suppress first occurrence
                text = text[:pipes[0]] + text[(pipes[0] + 1):]
        return text.strip()

    def has_comment_structure(self, text):
        """Determines whether the given text is a potential comment.

        :param text: text that may represent a comment
        :return: True if text contains protein or reference fields, False
            otherwise
        """
        return self._tags.protein.search(text) or self._tags.reference.search(text)

    def _fuse_abnormal_comment(self, text):
        """Fuses abnormal comments in the text to normal ones, or replaces
        abnormal comments by normal ones.

        An abnormal comment is delimited by pipe ('|') characters), whereas
        normal comments are delimited by parentheses.

        :param text: text that may contain an abnormal comment
        :return: text resulting from fusing the abnormal to the normal comment
            if a normal comment is already present in the text, or the text
            resulting from transforming the abnormal comment into a normal
            comment if a normal comment is not already present in the text
        """
        aobj = self._tags.abnormal_comment.search(text)
        if not aobj:
            return text

        abnormal_comment = text[(aobj.start() + 1):(aobj.end() - 1)]
        if not abnormal_comment.strip():
            return text.replace('|', '')

        cobj = self._tags.comment.search(text)
        if not cobj:
            return replace_abnormal_comment(text, aobj)

        comment_end = [x.start() for x in re.finditer(r'\)', text) if x.start() < aobj.start()]
        if not comment_end:
            return replace_abnormal_comment(text, aobj)

        comm_end = comment_end[-1]
        comm_start = [x.start() for x in re.finditer(r'\(', text) if x.start() < comm_end][-1]
        comment = text[comm_start:comm_end]
        if self.has_comment_structure(comment):
            return text[:comm_end] + '; ' + aobj.group(1) + ')'

        return replace_abnormal_comment(text, aobj)

    def _extract_numbers(self, text, pattern):
        """Extracts and returns numbers from a text according to a given pattern,
        as well as the text resulting from removing the numbers.

        The pattern may be one defining proteins or references (the instance's
        _tags.protein or _tags.reference attributes).

        :param text: text that may contain numbers according to a pattern
        :param pattern: a re pattern
        :return: text resulting from number extraction
        """
        if not isinstance(pattern, re.Pattern):
            raise ArgumentError('Expected re.Pattern: {}'.format(pattern))

        numbers = None
        mobj = pattern.search(text)
        if mobj:
            numbers = [int(num.group(0)) for num in self._tags.numbers.finditer(mobj.group(1))]
            text = text[:mobj.start()] + text[mobj.end():]

        return text.strip(), numbers

    def _guess_comment_indexes(self, text):
        """Determines the left and right indexes in text corresponding to a
        comment.

        If a comment cannot reliably be determined in text, one or both
        return values are None.

        :param text: text for which the comment indexes are to be determined
        :return: two integers representing the left and right indexes that
            delimit a comment in text, or None for either index or both indexes
            if a comment cannot reliably be determined
        """
        left = None
        right = None
        indexes_left, indexes_right = find_parentheses_indexes(text)

        if len(indexes_left) == 1:
            left = indexes_left[0]
        else:
            for i, _ in enumerate(indexes_left):
                next_idx = None
                if i == len(indexes_left) - 1:
                    next_indexes = [idx for idx in indexes_right if idx > indexes_left[-1]]
                    if next_indexes:
                        next_idx = next_indexes[0]
                else:
                    next_idx = indexes_left[i + 1]
                if next_idx:
                    subtext = text[indexes_left[i]:next_idx]
                    if self.has_comment_structure(subtext):
                        left = indexes_left[i]
                        break

        if len(indexes_right) == 1:
            right = indexes_right[0]
        else:
            for i in range(len(indexes_right) - 1, -1, -1):
                prev_idx = indexes_left[-1] if i == 0 else indexes_right[i - 1]
                subtext = text[prev_idx:indexes_right[i]]
                if self.has_comment_structure(subtext):
                    right = indexes_right[i]
                    break

        return left, right

    def _delimit_comment(self, text, enforce_structural_check=True):
        """Delimits and returns a comment in the given text, as well as the
        text resulting from removing the comment field.

        Preliminary preprocessing is assumed in order to assure that the text
        contains a single comment tag.

        :param text: text that may contain a comment
        :param enforce_structural_check: whether a comment should be checked
            for having protein and/or literature reference fields; True by
            default, should be false for Protein entries
        :return: text resulting from comment extraction, and the extracted
            comment
        """
        comment = None

        if self._tags.comment.search(text):
            left, right = self._guess_comment_indexes(text)
            if left and right:
                comment = text[(left + 1):right]
                if not enforce_structural_check or self.has_comment_structure(comment):
                    text = text[:left] + text[(right+1):]

        return text.strip(), comment

    def extract_comment(self, text, enforce_structural_check=True):
        """Extracts and stores a comment and returns the text resulting from
        removing the comment field.

        The comment format is defined by the instance's _tags.comment attribute.
        If a comment field is present in the text, it is stored in the
        instance's _current.comment attribute.

        :param text: text that may contain a comment
        :param enforce_structural_check: whether a comment should be checked
            for having protein and/or literature reference fields; True by
            default, should be false for Protein entries
        :return: text resulting from comment extraction
        """
        text = self._clean_extra_pipe_characters(text)
        text = self._fuse_abnormal_comment(text)
        text, comment = self._delimit_comment(text, enforce_structural_check)
        if comment is not None:
            self._current.comment = self._parse_comment(comment)
        return text.strip()

    def _parse_comment(self, text):
        """Parses a comment field in the BRENDA flat file.

        :param text: text representing the comment to be parsed
        :return: an EntryComment instance corresponding to the parsed comment
        """
        text = text.strip()
        if not text:
            text = None
        elif text[0] == '(' and text[-1] == ')':
            text = text[1:-1].strip()

        proteins = self._get_numbers_in_comment(text, self._tags.protein)
        references = self._get_numbers_in_comment(text, self._tags.reference)

        return EntryComment(text, proteins, references)

    def _get_numbers_in_comment(self, comment, pattern):
        """Returns a list of numbers present in the given comment according
        to the specified pattern.

        :param comment: text that may contain numbers
        :param pattern: pattern describing the numbers
        :return: the list of numbers in comment
        """
        if not isinstance(pattern, re.Pattern):
            raise ArgumentError('Expected re.Pattern: {}'.format(pattern))

        if comment is None or not comment.strip():
            return None

        numbers = list()

        matches = [match.group() for match in pattern.finditer(comment)]
        for number in matches:
            _, values = self._extract_numbers(number, pattern)
            numbers.extend([v for v in values if v not in numbers])

        return numbers

    def _parse_id(self, text):
        """Parses an EC number present in text.

        If the EC number is correctly parsed, a new Enzyme object is created.
        The text  may contain a comment. In addition, the newly created Enzyme
        instance is stored for partial EC numbers in the enzymes dict.

        :param text: text that may contain an EC number and a comment
        """
        comment = None
        mobj = self._tags.comment.search(text)

        if mobj:
            comment = self._parse_comment(mobj.group(1))
            if comment.msg and not comment.msg.strip():
                comment = None
            text = text[:mobj.start()] + text[mobj.end():]

        text = text.strip()
        if is_ec_number(text):
            self._current.ec_number = Enzyme(text, comment.msg if comment else None)
            ec_num = text.split(".")
            for i in range(1, len(ec_num) + 1):
                self.enzymes[".".join(ec_num[:i])].append(self._current.ec_number)

    def _parse_protein(self, text):
        """Parses a PROTEIN (PR) entry from the BRENDA flat file.

        :param text: text that represents a PROTEIN (PR) entry
        """
        self._reset_parser()
        text = self.extract_information(text)
        text = self.extract_comment(text, enforce_structural_check=False)

        mobj = self._tags.protein.search(text)
        if not mobj:
            raise ArgumentError(
                'Protein reference missing: \'{}\' @ #{}'.format(text, self._current.line_number))

        protein_id = int(mobj.group(1))
        text = text[:mobj.start()] + text[mobj.end():]
        text, accessions = self._extract_accessions(text)

        self._current.proteins = sorted(list(set(accessions)))
        text, self._current.references = self._extract_numbers(text, self._tags.reference)
        self._current.ec_number.proteins[protein_id] = Protein(text.strip(), self._current)

    def _extract_accessions(self, text):
        """Extracts and returns protein accessions from the given text, as well
        as the text resulting from removing the protein accession numbers.

        :param text: text that may contain accession numbers
        :return: text resulting from accession extraction, and the resulting
            accession numbers
        """
        accessions = list()

        accession_occurrences = self._tags.accession.findall(text)
        if accession_occurrences:
            accessions = [acc for acc, data_bank in accession_occurrences]

            first_acc = accession_occurrences[0][0]
            last_acc = accession_occurrences[-1][0]
            last_bank = accession_occurrences[-1][1]

            possible_end = text.index(last_acc) + len(last_acc)
            if last_bank is None:
                end = possible_end
            else:
                tmp_text = text[possible_end:]
                end = possible_end + tmp_text.index(last_bank) + len(last_bank)

            text = text[:text.index(first_acc)] + text[end:]

        return text.strip(), accessions

    def _parse_reference(self, text):
        pass

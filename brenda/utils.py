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
    utils.py

.. |c| unicode:: U+A9
"""

import re
import errno
import sys
from collections import namedtuple


class ArgumentError(Exception):
    """
    Error class that is meant to be raised when the arguments provided to a
    function are incorrect.
    """

    def __init__(self, msg, *args, **kw_args):
        """Creates an ArgumentError instance.

        A variable number of arguments may be passed. They will all be used to
        format msg. So take care that the number and type of additional
        arguments matches the format markers in msg.

        :param msg: unformatted string, i.e. it may contain multiple string
            format markers
        :param args:
        :param kw_args:
        """
        Exception.__init__(self, *args, **kw_args)
        self.args = (msg,) + args
        self.errno = errno.EINVAL
        self.strerror = msg % args

    def __str__(self):
        return self.strerror


def has_ec_number(text):
    return re.search(r'[1-7](\.\d+){2}\.\d+', text) is not None


def is_ec_number(text):
    return re.match(r'[1-7](\.\d+){2}\.\d+$', text) is not None


def init_tags():
    Tags = namedtuple(
        'Tags',
        ['protein', 'comment', 'abnormal_comment', 'information', 'reference',
         'numbers', 'accession'])

    # UniProt accession, see https://www.uniprot.org/help/accession_numbers
    _accession = r'[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9][A-Z][A-Z0-9]{2}[0-9]|' \
                 '[A-NR-Z][0-9][A-Z][A-Z0-9]{2}[0-9][A-Z][A-Z0-9]{2}[0-9]'

    return Tags(
        protein=re.compile(r'#(.+?)#', re.UNICODE),
        comment=re.compile(r' \((.*)\)', re.UNICODE),
        abnormal_comment=re.compile(r'\|(.*?)\|', re.UNICODE),
        information=re.compile(r'{(.*?)\}', re.UNICODE),
        reference=re.compile(r'<(.+?)>', re.UNICODE),
        numbers=re.compile(r'\d+', re.UNICODE),
        accession=re.compile(r'(%s)\s+(uniprot|unipro|swissprot|genbank|trembl|embl)*' %
                             _accession, re.UNICODE | re.I))


def find_parentheses_indexes(text):
    return [x.start() for x in re.finditer(r'\(', text)], \
           [x.start() for x in re.finditer(r'\)', text)]


def replace_abnormal_comment(text, aobj):
    return text[:aobj.start()] + \
           '(' + text[(aobj.start() + 1):(aobj.end() - 1)] + ')' + text[aobj.end():]


class ProgressMeter:
    """Displays a progress meter."""
    def __init__(self, label, end=None, **kw_args):
        super(ProgressMeter, self).__init__(**kw_args)
        if end:
            self.end = float(end)
            self.label = label

    def update(self, current):
        sys.stdout.write("\r{} {:.1%}".format(self.label, current / self.end))
        sys.stdout.flush()

    def close(self):
        sys.stdout.write("\r{} {:.1%}\n".format(self.label, 1.0))
        sys.stdout.flush()


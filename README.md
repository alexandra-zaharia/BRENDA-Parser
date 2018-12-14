# BRENDA Parser

This Python module provides classes that encapsulate information provided by the [BRENDA database](https://www.brenda-enzymes.org) and a parser that generates relevant content from the BRENDA [flat file](https://www.brenda-enzymes.org/download_brenda_without_registration.php) (`brenda_download.txt`) that can be downloaded for free.

The starting point of this BRENDA parser is [a project by Moritz Emanuel Beber](https://github.com/Midnighter/BRENDA-Parser) that I have improved by handling various parsing issues and edge cases in current BRENDA flat files (e.g., v2018.2). The project by Moritz has a new branch in the making, so you might want to check it out when it's finished. For the time being, I needed to have the BRENDA flat file v2018.2 correctly parsed, and this project achieves just that.

## Usage

In order to use this module outside of its directory, you will have to add it to your system path:

```python
import sys
sys.path.insert(0, '/path/to/BRENDA-Parser')
```

Provided you saved the BRENDA [flat file](https://www.brenda-enzymes.org/download_brenda_without_registration.php) under the name `brenda_download.txt`, parse it using BRENDA-Parser as follows:

```python
>>> from brenda.parser import BRENDAParser
>>> with BRENDAParser('brenda_download.txt') as parser:
...     brenda = parser.parse()
```

The method `parse` returns a dictionary with all EC numbers (enzymes) that were parsed. Every key in the dictionary is a string representing an EC number starting form the 6 general classes down to the individual enzymes. The values in the dictionary are always lists, e.g.:

```python
>>> len(brenda['1'])
2049
>>> len(brenda['1.1.1.1'])
1
>>> brenda['1.1.1.1'][0]
1.1.1.1
```

In the following, we will briefly survey how information in the BRENDA flat file is parsed and stored; to these means, BRENDA-Parser provides the classes [`Enzyme`](#ec-number-information-enzyme-class), [`Protein`](#protein-information-protein-class), [`Entry`](#brenda-entries-entry-class), and [`EntryComment`](#brenda-comments-entrycomment-class).

## API

### EC number information (`Enzyme` class)

For every EC number, the following fields may be accessed:

  * `ec_number` -- a string designating the EC number.
  * `comment` -- either `None` or a string representing a comment associated to an EC number.
  * `proteins` -- a dict storing protein information for an EC number. Dict keys are numerical protein identifiers as they appearin the BRENDA flat file between `#...#` tags. Dict values are [`Protein`](#protein-information-protein-class) objects.
  * `references` -- a dict storing literature references for an EC number (see [TODO](#TODO)).
  * `entries` -- a dict storing every entry in the BRENDA flat file for a given EC number, with the exception of protein (see `proteins` above) and literature (see `references` above) references. Dict keys are section names in the BRENDA flat file, such as REACTION, SUBSTRATE_PRODUCT, etc. Dict values are lists of [`Entry`](#brenda-entries-entry-class) objects.

For example, suppose you wish to list all EC numbers having an associated comment in the BRENDA flat file, such as:

```
...
ID      1.1.1.109 (transferred to EC 1.3.1.28)
...
ID      5.4.99.10 (deleted, included in EC 5.4.99.11)
...
```

Here is how to retrieve EC numbers with associated comments:

```python
>>> ec_comment = [ec for ec in brenda if brenda[ec][0].comment is not None]
>>> len(ec_comment)
910
>>> print(brenda[ec_comment[0]][0].ec_number, brenda[ec_comment[0]][0].comment)
1.1.1.109 transferred to EC 1.3.1.28
```

### Protein information (`Protein` class)

For a given protein (recall that proteins are stored in a dict where keys are protein numerical identifiers between `#...#` tags in the BRENDA flat file), the following fields may be accessed:

  * `information` -- either `None` or a string representing information associated to the current protein (`{...}` tags in the BRENDA flat file).
  * `comment` -- [`EntryComment`](#brenda-comments-entrycomment-class) object representing a comment associated to the current protein.
  * `organism` -- string representing the species for which the current protein is reported present.
  * `identifiers` -- list of strings representing UniProt accessions for the current protein.
  * `references` -- list of numerical identifiers representing literature references for the current protein (`<...>` tags in the BRENDA flat file).

For example, suppose you wish to display all protein information associated to EC number 6.6.1.2:

```python
>>> for prot_id in brenda['6.6.1.2'][0].proteins: 
...     print(prot_id) 
...     protein = brenda['6.6.1.2'][0].proteins[prot_id] 
...     print('\tInformation:', protein.information) 
...     print('\tComment:', protein.comment) 
...     print('\tOrganism:', protein.organism) 
...     print('\tIdentifiers:', protein.identifiers) 
...     print('\tReferences:', protein.references) 
1
        Information: None
        Comment: nomen rejiciendum
        Organism: Pseudomonas denitrificans
        Identifiers: []
        References: [1]
2
        Information: None
        Comment: None
        Organism: Salmonella enterica
        Identifiers: []
        References: [4]
3
        Information: None
        Comment: None
        Organism: Desulfovibrio vulgaris
        Identifiers: []
        References: [4]
4
        Information: None
        Comment: None
        Organism: Brucella melitensis
        Identifiers: []
        References: [3]
5
        Information: None
        Comment: nomen rejiciendum
        Organism: Pseudomonas denitrificans
        Identifiers: ['P29929', 'P29933', 'P29934', 'Q9HZQ3']
        References: [2]
```

### BRENDA entries (`Entry` class)

The BRENDA flat file contains several sections for every EC number, such as REACTION, SUBSTRATE_PRODUCT, COFACTOR, etc. This information is encapsulated as a list of `Entry` objects. Each `Entry` has the following fields:

  * `msg` -- the `Entry` content.
  * `information` -- additional information on the `Entry` (`{...}` tags in the BRENDA flat file).
  * `comment` -- [`EntryComment`](#brenda-comments-entrycomment-class) object representing a comment associated to the current `Entry`.
  * `proteins` -- list of list of numerical identifiers representing protein references for the current `Entry` (`#...#` tags in the BRENDA flat file).
  * `references` -- list of list of numerical identifiers representing literature references for the current `Entry` (`<...>` tags in the BRENDA flat file).

For example, suppose you wish to display information for the first SP (SUBSTRATE_PRODUCT) entry for EC number 1.1.1.261:

```python
>>> entry = brenda['1.1.1.261'][0].entries['SUBSTRATE_PRODUCT'][0]
>>> print(entry.msg)
dihydroxyacetone phosphate + NAD(P)H = sn-glycerol-1-phosphate + NAD(P)+
>>> print(entry.information)                                                                                              
ir
>>> print(entry.comment)                                                                                                  
#1,2,3,4,5,7# method specific to glycerol-1-phosphate <2>; #2# significant lower reverse reaction <2>; #3# reverse reaction 1/16 of the forward reaction <2>; #3# formation of glycerol-1-phosphate is the natural direction of the reaction <3,5>; #3# key enzyme in the biosynthesis of the enantiomeric glycerophosphate backbone of ether phospholipids of archaebacteria <6>; #3# involved in the biosynthesis of archeal lipids <7>; #8# key enzyme in the formation of archeal enantiomeric polar lipid structures <8>; #8# reverse reaction only with NAD+ as cofactor <8>
>>>  print(entry.proteins)                                                                                                 
[1, 2, 3, 4, 5, 7, 8, 9]
>>> print(entry.references)                                                                                               
[2, 3, 4, 5, 6, 7, 8, 9, 10, 11]
```

### BRENDA comments (`EntryComment` class)

Comments in the BRENDA flat file may provide protein and literature references. It is possible to access these references through the `proteins` and `references` fields of an `EntryComment` object. Here is how to access them for the comment in the above example:

```python
>>> entry = brenda['1.1.1.261'][0].entries['SUBSTRATE_PRODUCT'][0]
>>> print(entry.comment.proteins)
[1, 2, 3, 4, 5, 7, 8]
>>> print(entry.comment.references)
[2, 3, 5, 6, 7, 8]
```

Note that protein and literature references in comments are stored in the order in which they are encountered.

## Dependencies

BRENDA-Parser needs Python 3, as well as `recordclass`. Optionally, `nose` is needed to run the tests. You may install them with `pip`:

```bash
pip install recordclass nose
```

## Tests

Run the tests in the `tests/` directory as follows:

```bash
cd BRENDA-Parser
nosetests tests/
```

## What is different with respect to the forked project

  * UniProt accessions:
      * They are now [correctly](https://www.uniprot.org/help/accession_numbers) parsed (both 6- and 10-character accession types).
      * They are sometimes present in the BRENDA flat file with different data bank identifiers ('UniProt', even 'Unipro' (!), 'SwissProt', 'GenBank', 'TrEMBL', 'EMBL'), with different combinations of upper- and lowercase letters. This is now taken care of.
      * Sometimes the data bank identifier is not present, e.g. BRENDA contains entries such as `PR	#50# Plasmodium falciparum Q965D6  <17>`. This case is now handled.
  * EC numbers with associated comments are now handled (see [example](#ec-number-information-enzyme-class) above).
  * The BRENDA flat file is notoriously badly formatted.
      * Empty information (`{}`) and comment (`()` or `||`, see below) fields are now ignored (whitespace included).
      * `#` artifacts appearing in entries where they should not appear. Extra `#` characters are removed.
      * Unmatched opening and closing parentheses in compound names led the initial parser to incorrectly delimit comments (typically between `(...)` tags). Unmatched parentheses in compound names are now ignored.
      * Unmatched opening and closing parentheses in comments are now ignored and the comments are correctly retrieved.
      * Some entries in the BRENDA flat file use the undocumented comment tags `|...|`. Sometimes only comments delimited by either `(...)` or `|...|` are present, and sometimes both are present. Such comments are now fused together.
      * To make matters more interesting, there may also be extra `|` characters, seemingly inserted at random by BRENDA maintainers. This parser cleans up such spurious characters, then retrieves the (hopefully) correct comment.
  * The SOAP-parsing part has been removed.
  * The refactor/lexer branch has been removed (check out the [forked repository](https://github.com/Midnighter/BRENDA-Parser) for progress on that branch).

## TODO

* Parse and store references (RF entries).

# -*- coding: utf-8 -*-


"""
=============================
BRENDA Enzyme Database Parser
=============================

:Author:
    Alexandra Zaharia
:Date:
    2018-12-08
:Copyright:
    Copyright |c| 2018, Alexandra Zaharia, all rights reserved.
:File:
    test_parser.py

.. |c| unicode:: U+A9
"""

import unittest
import os
from brenda.parser import BRENDAParser

input_test = os.path.join('resources', 'brenda_test.txt')


class TestBrendaParser(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        super(TestBrendaParser, cls).setUpClass()
        with BRENDAParser(input_test) as cls.bp:
            cls.brenda = cls.bp.parse()
        cls.parser = BRENDAParser(None)

    def setUp(self):
        self.parser._reset_parser()

    def test_parsed_brenda_has_enzymes(self):
        self.assertIn('1.1.1.261', self.brenda.keys())
        self.assertIn('6.6.1.2', self.brenda.keys())
        self.assertIn('1.1.1.666', self.brenda.keys())
        self.assertIn('1.1.1.777', self.brenda.keys())
        self.assertIn('1.1.1.888', self.brenda.keys())
        
    def test_number_of_proteins_1_1_1_261(self):
        entry = self.brenda['1.1.1.261'][0]
        self.assertEqual(len(entry.proteins), 12)

    def test_number_of_proteins_6_6_1_2(self):
        entry = self.brenda['6.6.1.2'][0]
        self.assertEqual(len(entry.proteins), 7)
        
    def test_organisms_1_1_1_261(self):
        entry = self.brenda['1.1.1.261'][0]
        
        organisms = [
            'Methanosarcina barkeri',
            'Pyrococcus sp.',
            'Methanothermobacter thermautotrophicus',
            'Pyrococcus furiosus',
            'Halobacterium salinarum',
            'Thermoplasma acidophilum',
            'Thermoplasma sp.',
            'Aeropyrum pernix',
            'Sulfolobus tokodaii',
            'Pyrobaculum calidifontis',
            'Methanocaldococcus jannaschii',
            'Commonote archaea']

        for i in range(1, 13):
            self.assertEqual(entry.proteins[i].organism, organisms[i-1])
            
    def test_comment_6_6_1_2(self):
        entry = self.brenda['6.6.1.2'][0]
            
        comments = [
            'nomen rejiciendum',
            None,
            None, 
            None,
            'nomen rejiciendum',
            'bogus',
            'don\'t forget P29933 is subunit blaaah'
        ]
        
        for i in range(1, 8):
            if comments[i-1] is None:
                print('{} {}'.format(i, entry.proteins[i].comment))
                self.assertIsNone(entry.proteins[i].comment)
            else:
                self.assertEqual(entry.proteins[i].comment.msg, comments[i-1])
                
    def test_identifiers_6_6_1_2(self):
        entry = self.brenda['6.6.1.2'][0]
            
        identifiers = [
            [], ['A0A022YWF9'], ['Q01234'], ['P12345'], 
            ['P29933', 'P29934', 'Q9HZQ3', 'P29929'],
            ['P29933', 'P29934', 'Q9HZQ3', 'P29929'],
            ['P29933', 'P29934', 'Q9HZQ3', 'P29929']
        ]
        
        for i in range(1, 8):
            print(entry.proteins[i].identifiers)
            self.assertEqual(set(entry.proteins[i].identifiers),
                             set(identifiers[i-1]))

    def test_information_6_6_1_2(self):
        entry = self.brenda['6.6.1.2'][0]

        for i in range(1, 7):
            self.assertIsNone(entry.proteins[i].information)

        self.assertIsNotNone(entry.proteins[7].information)
        self.assertEqual(entry.proteins[7].information, 'some information')

    def test_enzyme_comment(self):
        entry666 = self.brenda['1.1.1.666'][0]
        entry777 = self.brenda['1.1.1.777'][0]
        entry888 = self.brenda['1.1.1.888'][0]

        self.assertIsNone(entry666.comment)
        self.assertIsNone(entry777.comment)
        self.assertEqual(entry888.comment, 'transferred from 1.1.1.999')

    def test_reaction_with_comment(self):
        text = \
            'a long-chain acyl-[acyl-carrier protein] + reduced flavodoxin '\
            '+ O2 = a (7S)-7-hydroxy-long-chain-acyl-[acyl-carrier protein] '\
            '+ oxidized flavodoxin + H2O (#1# 1a <1>)'
        reaction = self.parser._parse_generic_entry(text)

        self.assertIsNone(reaction.information)

        self.assertEqual(
            reaction.msg,
            'a long-chain acyl-[acyl-carrier protein] + reduced flavodoxin '
            '+ O2 = a (7S)-7-hydroxy-long-chain-acyl-[acyl-carrier protein] '
            '+ oxidized flavodoxin + H2O')

        self.assertEqual(reaction.comment.msg, '#1# 1a <1>')
        self.assertEqual(reaction.comment.proteins, [1])
        self.assertEqual(reaction.comment.references, [1])

        self.assertIsNone(reaction.proteins)
        self.assertIsNone(reaction.references)

    def test_reaction_with_abnormal_comment(self):
        text = \
            '#1# hexadecanoyl-[acyl-carrier protein] + reduced flavodoxin '\
            '+ O2 = 11-hydroxyhexadecanoyl--[acyl-carrier protein] + '\
            '12-hydroxyhexadecanoyl-[acyl-carrier protein] + '\
            '13-hydroxytetradecanoxl-[acyl-carrier protein] + '\
            '14-hydroxytetradecanoxl-[acyl-carrier protein] + '\
            '15-hydroxytetradecanoxl-[acyl-carrier protein] + '\
            'oxidized flavodoxin + H2O |#1# the enzyme produces mainly the '\
            '11- to 15-hydroxy C16 fatty acids <6>| {r} <6>'
        reaction = self.parser._parse_generic_entry(text)

        self.assertEqual(
            reaction.msg,
            'hexadecanoyl-[acyl-carrier protein] + reduced flavodoxin + O2 '
            '= 11-hydroxyhexadecanoyl--[acyl-carrier protein] + '
            '12-hydroxyhexadecanoyl-[acyl-carrier protein] + '
            '13-hydroxytetradecanoxl-[acyl-carrier protein] + '
            '14-hydroxytetradecanoxl-[acyl-carrier protein] + '
            '15-hydroxytetradecanoxl-[acyl-carrier protein] + '
            'oxidized flavodoxin + H2O')

        self.assertEqual(reaction.information, 'r')

        self.assertEqual(
            reaction.comment.msg,
            '#1# the enzyme produces mainly the '
            '11- to 15-hydroxy C16 fatty acids <6>'
        )

        self.assertEqual(reaction.comment.proteins, [1])
        self.assertEqual(reaction.comment.references, [6])

        self.assertEqual(reaction.proteins, [1])
        self.assertEqual(reaction.references, [6])

    def test_reaction_with_long_protein_and_reference_list(self):
        text = \
            '#1,2,3,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,'\
            '52,53,54 55,56,57,58,59,60,61,62,63,64,65,66,67,68,69,70,71,72,'\
            '73,74,75,76,77,78 79,80,81,82,83,84,86# '\
            'L-Arg + pyruvate + NADH = N2-(D-1-carboxyethyl)-L-Arg + NAD+ '\
            '+ H2O (#35,37,40,64,72,77# r <3,8,14,15,17,24>; '\
            '#84# 2 enzyme forms. One catalyzes the reverse reaction with '\
            'NAD+, the second with NADP+ <33>) |#34# i.e. octopine <1>| '\
            '{} <1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,'\
            '23,24,25,26 27,28,29,30,31,32,33,34,35,36,37>'
        reaction = self.parser._parse_generic_entry(text)

        self.assertEqual(
            reaction.msg,
            'L-Arg + pyruvate + NADH = N2-(D-1-carboxyethyl)-L-Arg + NAD+ + H2O'
        )

        self.assertIsNone(reaction.information)

        self.assertEqual(
            reaction.comment.msg,
            '#35,37,40,64,72,77# r <3,8,14,15,17,24>; '
            '#84# 2 enzyme forms. One catalyzes the reverse reaction with '
            'NAD+, the second with NADP+ <33>; #34# i.e. octopine <1>')

        self.assertEqual(reaction.comment.proteins,
                         [35, 37, 40, 64, 72, 77, 84, 34])
        self.assertEqual(reaction.comment.references,
                         [3, 8, 14, 15, 17, 24, 33, 1])

        proteins = [1, 2, 3]
        proteins.extend([x for x in range(34, 85)])
        proteins.append(86)
        self.assertEqual(reaction.proteins, proteins)

        self.assertEqual(reaction.references, [x for x in range(1, 38)])

    def test_reaction_with_PDB_id_at_the_beginning_of_continuation_line(self):
        self.assertIn('1.1.1.35', self.brenda)
        self.assertIn('REACTION', self.brenda['1.1.1.35'][0].entries)
        reaction = self.brenda['1.1.1.35'][0].entries['REACTION'][0]
        self.assertEqual(reaction.msg, '(S)-3-hydroxyacyl-CoA + NAD+ = 3-oxoacyl-CoA + NADH + H+')
        self.assertEqual(
            reaction.comment.msg,
            '#4# highly conserved residue Ser137 is essential for activity <29>; #5# molecular '
            'reaction mechanism, a His-Glu pair is essential for actalysis in the active site '
            '<33>; #4# reaction mechanism via enolate intermediate <37>; #19# structure-function '
            'analysis and catalytic mechanism, overview <59>; #25# structure-function analysis '
            'using the crystal structure, PDB ID 4j0f, and catalytic mechanism, overview <59>')

    def test_reaction_with_bogus_pipe_character_among_reactants(self):
        text = '#5# lithocholic acid + NADPH + H+ | = ursodeoxycholic acid + NADP+ <9>'
        reaction = self.parser._parse_generic_entry(text)
        self.assertEqual(reaction.proteins, [5])
        self.assertEqual(reaction.references, [9])
        self.assertIsNone(reaction.comment)
        self.assertEqual(
            reaction.msg, 'lithocholic acid + NADPH + H+  = ursodeoxycholic acid + NADP+')

    def test_reaction_with_bogus_pipe_and_abnormal_comment(self):
        text = "#4# (ribonucleotide)n-2',3'-cyclophosphate + 5'-hydroxy-(ribonucleotide)m + GTP " \
               "+ H2O = (ribonucleotide)n+m + GMP + diphosphate | (#4# substrate mimicks the " \
               "broken tRNAGlu(UUC) anticodon stem-loop generated by Kluyveromyces lactis " \
               "gamma-toxin, leaving 2,3-cyclic phosphate and 5-OH ends <1,7>) |#4# overall " \
               "reaction <1,7>| <1,7>"
        reaction = self.parser._parse_generic_entry(text)
        self.assertEqual(reaction.proteins, [4])
        self.assertEqual(reaction.references, [1, 7])
        self.assertEqual(
            reaction.comment.msg,
            '#4# substrate mimicks the broken tRNAGlu(UUC) anticodon stem-loop generated by '
            'Kluyveromyces lactis gamma-toxin, leaving 2,3-cyclic phosphate and 5-OH ends <1,7>; '
            '#4# overall reaction <1,7>')
        self.assertEqual(reaction.comment.proteins, [4])
        self.assertEqual(reaction.comment.references, [1, 7])
        self.assertEqual(
            reaction.msg,
            "(ribonucleotide)n-2',3'-cyclophosphate + 5'-hydroxy-(ribonucleotide)m + GTP + H2O = "
            "(ribonucleotide)n+m + GMP + diphosphate")

    def test_reaction_with_curly_brackets_and_information(self):
        text = '#7,28,41,42,48,56# cefoperazone + H2O = (2R)-2-[(R)-carboxy{[(2S)-2-[(4-ethyl-2 ' \
               '3-dioxopiperazine-1-carbonyl)amino]-2-(4 hydroxyphenyl)acetyl]amino}methyl]' \
               '-5-{[(1-methyl-1H-tetrazol-5 yl)sulfanyl]methyl}-3,6-dihydro-2H-1,3-thiazine' \
               '-4-carboxylic acid (#41# 11% relative activity to cephaloridine <4>) {} ' \
               '<4,8,9,45,77,153>'
        reaction = self.parser._parse_generic_entry(text)
        self.assertEqual(reaction.proteins, [7, 28, 41, 42, 48, 56])
        self.assertEqual(reaction.references, [4, 8, 9, 45, 77, 153])
        self.assertEqual(reaction.comment.msg, '#41# 11% relative activity to cephaloridine <4>')
        self.assertEqual(reaction.comment.proteins, [41])
        self.assertEqual(reaction.comment.references, [4])
        self.assertIsNone(reaction.information)
        self.assertEqual(
            reaction.msg,
            'cefoperazone + H2O = (2R)-2-[(R)-carboxy{[(2S)-2-[(4-ethyl-2 3-dioxopiperazine'
            '-1-carbonyl)amino]-2-(4 hydroxyphenyl)acetyl]amino}methyl]-5-{[(1-methyl-1H'
            '-tetrazol-5 yl)sulfanyl]methyl}-3,6-dihydro-2H-1,3-thiazine-4-carboxylic acid')

    def test_reaction_with_unmatched_number_of_parentheses(self):
        text = '#1,2# dihydroclavaminate + 2-oxoglutarate + O2 = clavaminate + succinate + ' \
               'CO2 + H2O (#2# stereochemical course of oxygen insertion <5>; #1,2# cyclization ' \
               '<1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17>; #2# syn-elimination of the ' \
               'requisite hydrogens <12>) |#2# 5S) enantiomer <3,4>| ' \
               '<1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17> '
        reaction = self.parser._parse_generic_entry(text)
        self.assertEqual(reaction.proteins, [1, 2])
        self.assertEqual(reaction.references, [x for x in range(1, 18)])
        self.assertEqual(reaction.comment.msg,
                         '#2# stereochemical course of oxygen insertion <5>; #1,2# cyclization '
                         '<1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17>; #2# syn-elimination of '
                         'the requisite hydrogens <12>; #2# 5S) enantiomer <3,4>')
        self.assertEqual(reaction.comment.proteins, [2, 1])
        self.assertEqual(reaction.comment.references,
                         [5, 1, 2, 3, 4, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17])
        self.assertIsNone(reaction.information)
        self.assertEqual(
            reaction.msg,
            'dihydroclavaminate + 2-oxoglutarate + O2 = clavaminate + succinate + CO2 + H2O')

    def test_reaction_with_unmatched_brackets_in_comment_and_matched_brackets_in_equation(self):
        text = \
            '#1,4,5,6,8,9,10,14,15,16,18,19,20,25,28,29,31,32,34# (E,E)-farnesyl ' \
            'diphosphate + isopentenyl diphosphate = geranylgeranyl diphosphate + ' \
            'diphosphate (#6,14# no activity with dimethylallyl diphosphate <8,25>; ' \
            '#29# 66% of activity with geranyl diphosphate <24>; #10# preferred ' \
            'allylic substrate <15>; #25# enzyme regulates taxane biosynthesis <25>; ' \
            '#34# bifunctional enzyme EC2.5.1.29/EC2.5.1.10, the FPP/GGPP product ' \
            'ratio increases with the rise of the reaction temperature. The enzyme ' \
            'contributes to adjust the membrane composition to the cell growth ' \
            'temperature modulating its substrate and product specificities <32>; ' \
            '#34# bifunctional enzyme EC2.5.1.29/EC2.5.1.10, the FPP/GGPP product ' \
            'ratio increases with the rise of the reaction temperature <32>) |#8# ' \
            'exclusive product <19>; #1,5,19# E,E)-geranylgeranyl diphosphate ' \
            '<2,7,10>| {}'
        reaction = self.parser._parse_generic_entry(text)
        self.assertEqual(reaction.proteins,
                         [1, 4, 5, 6, 8, 9, 10, 14, 15, 16, 18, 19, 20, 25, 28, 29, 31, 32, 34])
        self.assertIsNone(reaction.references)
        self.assertEqual(
            reaction.comment.msg,
            '#6,14# no activity with dimethylallyl diphosphate <8,25>; #29# 66% of activity with '
            'geranyl diphosphate <24>; #10# preferred allylic substrate <15>; #25# enzyme '
            'regulates taxane biosynthesis <25>; #34# bifunctional enzyme EC2.5.1.29/EC2.5.1.10, '
            'the FPP/GGPP product ratio increases with the rise of the reaction temperature. '
            'The enzyme contributes to adjust the membrane composition to the cell growth '
            'temperature modulating its substrate and product specificities <32>; #34# '
            'bifunctional enzyme EC2.5.1.29/EC2.5.1.10, the FPP/GGPP product ratio increases '
            'with the rise of the reaction temperature <32>; #8# exclusive product <19>; '
            '#1,5,19# E,E)-geranylgeranyl diphosphate <2,7,10>')
        self.assertEqual(reaction.comment.proteins, [6, 14, 29, 10, 25, 34, 8, 1, 5, 19])
        self.assertEqual(reaction.comment.references, [8, 25, 24, 15, 32, 19, 2, 7, 10])
        self.assertIsNone(reaction.information)
        self.assertEqual(
            reaction.msg,
            '(E,E)-farnesyl diphosphate + isopentenyl diphosphate = '
            'geranylgeranyl diphosphate + diphosphate')

    def test_reaction_with_many_parentheses_in_equation(self):
        text = \
            '#10# Manalpha(1-6)(Manalpha(1-3))Manbeta(1-4)GlcNAcbeta(1-4)GlcNAc + H2O = ' \
            'Manalpha(1-6)Manbeta(1-4)GlcNAcbeta(1-4)GlcNAc + alpha-D-mannose (#10# NaBH4 ' \
            'reduced, cleaves the Manalpha(1-6)Man linkage only after its Manalpha(1-3) ' \
            'residue is removed <4>) |#10# NaBH4 reduced, no product: ' \
            'Manalpha(1-3)Manbeta(1-4)GlcNAcbeta(1-4)GlcNAc <4>| {} <4,7>'
        reaction = self.parser._parse_generic_entry(text)
        self.assertEqual(reaction.proteins, [10])
        self.assertEqual(reaction.references, [4, 7])
        self.assertEqual(
            reaction.comment.msg,
            '#10# NaBH4 reduced, cleaves the Manalpha(1-6)Man linkage only after its '
            'Manalpha(1-3) residue is removed <4>; #10# NaBH4 reduced, no product: '
            'Manalpha(1-3)Manbeta(1-4)GlcNAcbeta(1-4)GlcNAc <4>')
        self.assertEqual(reaction.comment.proteins, [10])
        self.assertEqual(reaction.comment.references, [4])
        self.assertIsNone(reaction.information)
        self.assertEqual(reaction.msg,
                         'Manalpha(1-6)(Manalpha(1-3))Manbeta(1-4)GlcNAcbeta(1-4)GlcNAc + H2O = '
                         'Manalpha(1-6)Manbeta(1-4)GlcNAcbeta(1-4)GlcNAc + alpha-D-mannose')

    def test_reaction_with_extra_hash_character(self):
        text = \
            '#5# colloidal chitin + H2O = N-acetylglucosamine + ' \
            'N,N#-diacetylchitobiose + ? |#5# ChiB produces relatively large amounts ' \
            'of chitin monomer and dimer, but very low amounts of trimer <10>| <10>'
        reaction = self.parser._parse_generic_entry(text)
        self.assertEqual(reaction.proteins, [5])
        self.assertEqual(reaction.references, [10])
        self.assertEqual(
            reaction.comment.msg,
            '#5# ChiB produces relatively large amounts '
            'of chitin monomer and dimer, but very low amounts of trimer <10>')
        self.assertEqual(reaction.comment.proteins, [5])
        self.assertEqual(reaction.comment.references, [10])
        self.assertIsNone(reaction.information)
        self.assertEqual(
            reaction.msg,
            'colloidal chitin + H2O = N-acetylglucosamine + N,N-diacetylchitobiose + ?')


if __name__ == '__main__':
    unittest.main()

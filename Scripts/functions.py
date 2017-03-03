#!/usr/bin/env python
"""
    EffectorP: predicting fungal effector proteins from secretomes using machine learning
    Copyright (C) 2015-2016 Jana Sperschneider	

    This program is free software; you can redistribute it and/or modify  
    it under the terms of the GNU General Public License as published by  
    the Free Software Foundation; either version 3 of the License, or     
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.

    Please also see the CSIRO Disclaimer provided with EffectorP (LICENCE.txt).

    Contact: jana.sperschneider@csiro.au
"""
# -----------------------------------------------------------------------------------------------------------
# -----------------------------------------------------------------------------------------------------------
# -----------------------------------------------------------------------------------------------------------
import os
import sys
import re
import io
import getopt
# -----------------------------------------------------------------------------------------------------------
# -----------------------------------------------------------------------------------------------------------
# -----------------------------------------------------------------------------------------------------------
def usage():
    """ Function: usage()

        Purpose:  Print helpful information for the user.        
        
        Input:    None.
    
        Return:   Print options for running EffectorP.       
    """
    print '''
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# EffectorP :: predicting fungal effector proteins from secretomes using machine learning
# EffectorP 1.0 (July 2015); http://effectorp.csiro.au/
# Copyright (C) 2015-2016 Jana Sperschneider, CSIRO.
# Freely distributed under the GNU General Public License (GPLv3).
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    '''
    print "Usage for EffectorP: ", 
    print "python EffectorP.py [-options] -i <input_file>"
    print 
    print "where basic options are:"
    print "-h : show brief help on version and usage" 
    print 
    print "options for output format:"
    print "-s : short output format that provides predictions for all proteins as one tab-delimited table [default long format]"
    print
    print "options directing output:"
    print "-o <f> : direct output to file <f>, not stdout"
    print "-E <f> : save predicted effectors to FASTA file <f>"        
    print
    print "# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -"
    print
    sys.exit()    

    return
# -----------------------------------------------------------------------------------------------------------
def scan_arguments(commandline):
    """ Function: scan_arguments()

        Purpose:  Scan the input options given to the EffectorP program.        
        
        Input:    Input options given by the user.
    
        Return:   Parsed options.
    """
    try:
        opts, args = getopt.getopt(commandline, "hso:E:i:", ["help"])        
    except getopt.GetoptError as err:
        # print help information and exit:
        print str(err) # will print something like "option -a not recognized"
        usage()
        sys.exit(2)

    FASTA_FILE = None
    short_format = False
    output_file = None
    effector_output = None

    i_count, o_count, E_count, P_count = 0, 0, 0, 0
   
    for opt, arg in opts:
        if opt in ("-o"):
            output_file = arg
            o_count += 1
        elif opt in ("-s"):
            short_format = True
        elif opt in ("-i"):
            FASTA_FILE = arg
            i_count += 1
        elif opt in ("-E"):
            effector_output = arg
            E_count += 1
        elif opt in ("-h", "--help"):
            usage()
        else:
            assert False, "unhandled option"

    if i_count > 1 or o_count > 1 or E_count > 1:
       usage()

    return FASTA_FILE, short_format, output_file, effector_output
# -----------------------------------------------------------------------------------------------------------
# This is the ARFF file header in WEKA format
ARFF_HEADER = '''@RELATION effectors
@ATTRIBUTE Tiny NUMERIC
@ATTRIBUTE Small NUMERIC
@ATTRIBUTE Aliphatic NUMERIC
@ATTRIBUTE Aromatic NUMERIC
@ATTRIBUTE Nonpolar NUMERIC
@ATTRIBUTE Polar NUMERIC
@ATTRIBUTE Charged NUMERIC
@ATTRIBUTE Basic NUMERIC
@ATTRIBUTE Acidic NUMERIC
@ATTRIBUTE A NUMERIC
@ATTRIBUTE C NUMERIC
@ATTRIBUTE D NUMERIC
@ATTRIBUTE E NUMERIC
@ATTRIBUTE F NUMERIC
@ATTRIBUTE G NUMERIC
@ATTRIBUTE H NUMERIC
@ATTRIBUTE I NUMERIC
@ATTRIBUTE K NUMERIC
@ATTRIBUTE L NUMERIC
@ATTRIBUTE M NUMERIC
@ATTRIBUTE N NUMERIC
@ATTRIBUTE P NUMERIC
@ATTRIBUTE Q NUMERIC
@ATTRIBUTE R NUMERIC
@ATTRIBUTE S NUMERIC
@ATTRIBUTE T NUMERIC
@ATTRIBUTE V NUMERIC
@ATTRIBUTE W NUMERIC
@ATTRIBUTE Y NUMERIC
@ATTRIBUTE MW NUMERIC
@ATTRIBUTE Charge NUMERIC
@ATTRIBUTE Length NUMERIC
@ATTRIBUTE class {effector,non-effector}

@DATA
'''
# -----------------------------------------------------------------------------------------------------------
def get_seqs_ids_fasta(FASTA_FILE):
    """ Function: get_seqs_ids_fasta()

        Purpose:  Given a FASTA format file, this function extracts
                  the list of identifiers and the list of sequences 
                  in the order in which they appear in the FASTA file.
              
        Input:    Path to FASTA format file.
    
        Return:   List of identifiers and list of sequences in the order 
                  in which they appear in the FASTA file.
    """ 
    identifiers = []
    sequences = []

    with open(FASTA_FILE) as f: 
        content = f.readlines()

        for position, line in enumerate(content):
            if '>' in line:
                identifiers.append(line)
                seq = []
                following_lines = content[position + 1:]
                for next_line in following_lines:
                    if '>' not in next_line:
                        seq.append(next_line.strip())
                    else:
                        break
                sequence = "".join(seq)
                sequence = sequence.replace('*', '')
                sequences.append(sequence)

    return identifiers, sequences
# -----------------------------------------------------------------------------------------------------------
def write_FASTA_short_ids(f_output, ORIGINAL_IDENTIFIERS, ORIGINAL_SEQUENCES):
    """ Function: write_FASTA_short_ids()

        Purpose:  Given a list of identifiers and the corresponding list 
                  of sequence, write these to a FASTA file using short
                  identifiers such as protein1, protein2, .... This is 
                  done because some programs like pepstats do not like 
                  long identifier names as input.
              
        Input:    Path to desired FASTA format output file, list of 
                  identifiers and list of corresponding sequences.
    
        Return:   List of short identifiers.
    """

    with open(f_output, 'w') as f:
        SHORT_IDENTIFIERS = []
        # Change identifiers to protein1, protein2, ...
        # and write to temporary file
        SET = zip(ORIGINAL_IDENTIFIERS, ORIGINAL_SEQUENCES)
    
        for index,  (identifier, sequence) in enumerate(SET):
            short_id = '>protein' + str(index)
            SHORT_IDENTIFIERS.append(short_id)
            f.writelines(short_id + '\n')
            f.writelines(sequence + '\n')

    return SHORT_IDENTIFIERS
# -----------------------------------------------------------------------------------------------------------
def pepstats(SHORT_IDENTIFIERS, SEQUENCES, pepstats_file):
    """ Function: pepstats()

        Purpose:  Given a set that contains the list of identifiers and 
                  the corresponding list of sequences, scan the given 
                  pepstats result file to extract protein properties.
              
        Input:    Set that contains the list of identifiers and 
                  the corresponding list of sequences and peptstats 
                  result file. 
    
        Return:   Dictionary of protein properties for each protein in the 
                  set.
    """

    pepstats_dic = {}

    with open(pepstats_file) as f: 
        content = f.readlines()
        for start, line in enumerate(content):
            if 'PEPSTATS of ' in line:
                TARGET_ID = line.split('PEPSTATS of ')[1]
                TARGET_ID = TARGET_ID.split('from 1 to')[0]
                TARGET_ID = TARGET_ID.strip()
                sequence = None

                for (identifier, seq) in zip(SHORT_IDENTIFIERS, SEQUENCES):
                    if identifier.replace('>', '') == TARGET_ID:
                        sequence = seq.strip()

                if sequence:
                    length = float(len(sequence))
                    # Amino acid frequencies in the sequence
                    amino_acid_frequencies = []
                    amino_acid_frequencies.append(100.0*sequence.count('A')/length)
                    amino_acid_frequencies.append(100.0*sequence.count('C')/length)
                    amino_acid_frequencies.append(100.0*sequence.count('D')/length)
                    amino_acid_frequencies.append(100.0*sequence.count('E')/length)
                    amino_acid_frequencies.append(100.0*sequence.count('F')/length)
                    amino_acid_frequencies.append(100.0*sequence.count('G')/length)
                    amino_acid_frequencies.append(100.0*sequence.count('H')/length)
                    amino_acid_frequencies.append(100.0*sequence.count('I')/length)
                    amino_acid_frequencies.append(100.0*sequence.count('K')/length)
                    amino_acid_frequencies.append(100.0*sequence.count('L')/length)
                    amino_acid_frequencies.append(100.0*sequence.count('M')/length)
                    amino_acid_frequencies.append(100.0*sequence.count('N')/length)
                    amino_acid_frequencies.append(100.0*sequence.count('P')/length)
                    amino_acid_frequencies.append(100.0*sequence.count('Q')/length)
                    amino_acid_frequencies.append(100.0*sequence.count('R')/length)
                    amino_acid_frequencies.append(100.0*sequence.count('S')/length)
                    amino_acid_frequencies.append(100.0*sequence.count('T')/length)
                    amino_acid_frequencies.append(100.0*sequence.count('V')/length)
                    amino_acid_frequencies.append(100.0*sequence.count('W')/length)
                    amino_acid_frequencies.append(100.0*sequence.count('Y')/length)
                    # Extract molecular weight
                    mwline = content[start + 2:start + 3]
                    molecular_weight = float(re.findall("\d+.\d+", str(mwline))[0])
                    # Extract charge
                    charge_line = content[start + 3:start + 4]
                    charge = float(re.findall("[-+]?\d+.\d+", str(charge_line))[1])
                    # Extract amino acid class frequencies

		    # Watch out, in the pepstats software, if isoelectric point == None, an 
                    # extra line will be introduced
                    start_aas = content[start:].index('Property\tResidues\t\tNumber\t\tMole%\n')
                    perline = content[start + start_aas + 1:start + start_aas + 10]

                    tiny = float(re.findall("\d+.\d+", str(perline[0]))[-1])
                    small = float(re.findall("\d+.\d+", str(perline[1]))[-1])
                    aliphatic = float(re.findall("\d+.\d+", str(perline[2]))[-1])
                    aromatic = float(re.findall("\d+.\d+", str(perline[3]))[-1])
                    non_polar = float(re.findall("\d+.\d+", str(perline[4]))[-1])
                    polar = float(re.findall("\d+.\d+", str(perline[5]))[-1])
                    charged = float(re.findall("\d+.\d+", str(perline[6]))[-1])
                    basic = float(re.findall("\d+.\d+", str(perline[7]))[-1])
                    acidic = float(re.findall("\d+.\d+", str(perline[8]))[-1])
                    amino_acid_classes = []
                    amino_acid_classes.append(tiny)
                    amino_acid_classes.append(small)
                    amino_acid_classes.append(aliphatic)
                    amino_acid_classes.append(aromatic)
                    amino_acid_classes.append(non_polar)
                    amino_acid_classes.append(polar)
                    amino_acid_classes.append(charged)
                    amino_acid_classes.append(basic)
                    amino_acid_classes.append(acidic)
                    # Store values in dictionary
                    pepstats_dic[TARGET_ID] = molecular_weight, charge, amino_acid_classes, amino_acid_frequencies, length

                else:
                    print 'There was an error scanning the pepstats file.'
                    print 'Could not find corresponding sequence for identifier', TARGET_ID
                    sys.exit()

    return pepstats_dic
# -----------------------------------------------------------------------------------------------------------
def write_weka_input(weka_input, SHORT_IDENTIFIERS, pepstats_dic):
    """ Function: write_weka_input()

        Purpose:  Given the query identifiers and pepstats-calculated
                  protein features, write the input arff file for WEKA. 
              
        Input:    WEKA arff file name, query identifiers and pepstats dictionary.                  
    
        Return:   None. 
    """   
    with open(weka_input, 'w') as f:
        # Create a list of features for each protein
        X = [[] for __ in xrange(len(SHORT_IDENTIFIERS))]

        for protein_position, TARGET_ID in enumerate(SHORT_IDENTIFIERS):
            TARGET_ID = TARGET_ID.replace('>', '')
            TARGET_ID = TARGET_ID.strip()

            molecular_weight, charge, amino_acid_classes, amino_acid_frequencies, length = pepstats_dic[TARGET_ID]
            X[protein_position] = amino_acid_classes + amino_acid_frequencies + [molecular_weight, charge, length]
 
        # Write protein feature data to WEKA arff file
        f.writelines(ARFF_HEADER)
        for index, vector in enumerate(X):
            for feature in vector:
                f.writelines(str(feature) + ',')
            f.writelines('?\n')

    return
# -----------------------------------------------------------------------------------------------------------
def parse_weka_output(file_input, ORIGINAL_IDENTIFIERS, SEQUENCES):
    """ Function: parse_weka_output()

        Purpose:  Given the WEKA output file and the query identifiers and sequences, 
                  parse the predicted class for each protein from the WEKA output. 
                  Write the predicted effectors to a FASTA file.
              
        Input:    WEKA output file and the query identifiers and sequences.                  
    
        Return:   The set of predicted effectors only as well as all predictions. 
    """    
    predicted_effectors, predictions = [], []

    with open(file_input) as f:

        content = f.readlines()
        content_start = content.index(' inst#     actual  predicted error prediction (Tiny,Small,Aliphatic,Aromatic,Nonpolar,Polar,Charged,Basic,Acidic,A,C,D,E,F,G,H,I,K,L,M,N,P,Q,R,S,T,V,W,Y,MW,Charge,Length)\n')
        content = content[content_start + 1:]

        for line in content:
            if line.strip():
                position = line.split()[0]
                prediction = line.split()[2]
                prob = float(line.split()[3])        
 
                # WEKA output counts from position 1, our identifiers are counted from zero
                identifier = ORIGINAL_IDENTIFIERS[int(position) - 1]
                sequence = SEQUENCES[int(position) - 1]

                if 'non-eff' in prediction:                               
                    noneffector = identifier.strip()
                    noneffector = noneffector.replace('>', '')  
                    predictions.append((noneffector, 'Non-effector', prob, sequence))
                else:                    
                    effector = identifier.strip()
                    effector = effector.replace('>', '')                                               
                    predictions.append((effector, 'Effector', prob, sequence))
	            # Append predicted effector to list of predicted effectors
                    predicted_effectors.append((effector, prob, sequence))

    return predicted_effectors, predictions
# -----------------------------------------------------------------------------------------------------------
def short_output(predictions):
    """ Function: short_output()

        Purpose:  Given the WEKA predictions for each protein, write  
                  string that contains the short output format.
              
        Input:    WEKA predictions for each protein.                  
    
        Return:   String that contains predictions for all proteins as tab-delimited table.
    """
    # Output predictions for all proteins as tab-delimited table
    short_output_string = '# Identifier \t Prediction \t Probability \n'
    for protein, pred, prob, sequence in predictions:    
        short_output_string += protein + '\t' + pred + '\t' + str(prob) + '\n'            

    return short_output_string
# -----------------------------------------------------------------------------------------------------------
def long_output(ORIGINAL_IDENTIFIERS, predicted_effectors):
    """ Function: long_output()

        Purpose:  Given the predicted effectors and identifiers for the test set,  
                  write string that contains the long output format.
              
        Input:    Predicted effectors and identifiers of test set.                  
    
        Return:   String that contains list of predicted effectors with posterior probabilites
                  and a short statistic on the percentage of predicted effectors in the test set.
    """
    # Output predicted effectors for long format
    long_output_string = '-----------------\n'
    long_output_string += 'Predicted effectors:\n\n'
    for effector, prob, sequence in predicted_effectors:
        long_output_string += effector + '| Effector probability:' + str(prob) + '\n'

    long_output_string += '-----------------\n\n'
    long_output_string += 'Number of proteins that were tested: ' + str(len(ORIGINAL_IDENTIFIERS)) + '\n' 
    long_output_string += 'Number of predicted effectors: ' + str(len(predicted_effectors)) + '\n' 
    long_output_string += '\n' + '-----------------' + '\n' 
    long_output_string += str(round(100.0*len(predicted_effectors)/len(ORIGINAL_IDENTIFIERS), 1)) + ' percent are predicted to be effectors.'  
    long_output_string += '\n' + '-----------------' + '\n'

    return long_output_string
# -----------------------------------------------------------------------------------------------------------           

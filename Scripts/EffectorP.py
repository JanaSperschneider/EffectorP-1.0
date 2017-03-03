#! /usr/bin/python
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
# -----------------------------------------------------------------------------------------------------------
# -----------------------------------------------------------------------------------------------------------
# -----------------------------------------------------------------------------------------------------------
import os
import sys
import functions
import subprocess
import errno
import uuid
import shutil
# -----------------------------------------------------------------------------------------------------------
# -----------------------------------------------------------------------------------------------------------
# -----------------------------------------------------------------------------------------------------------
# Main Program starts here
# -----------------------------------------------------------------------------------------------------------
# -----------------------------------------------------------------------------------------------------------
# -----------------------------------------------------------------------------------------------------------
if __name__ == '__main__': 
    # -----------------------------------------------------------------------------------------------------------
    SCRIPT_PATH = sys.path[0]
    # Change the path to WEKA to the appropriate location on your computer
    WEKA_PATH = SCRIPT_PATH + '/weka-3-6-12/weka.jar'
    PEPSTATS_PATH = SCRIPT_PATH + '/EMBOSS-6.5.7/emboss/'
    # -----------------------------------------------------------------------------------------------------------
    # Check that the path to the WEKA software exists
    path_exists = os.access(WEKA_PATH, os.F_OK)
    if not path_exists:
        print
        print "Path to WEKA software does not exist!"
        print "Check the installation and the given path to the WEKA software %s in EffectorP.py (line 47)."%WEKA_PATH
        print
        sys.exit()
    # -----------------------------------------------------------------------------------------------------------
    # Check that the path to the EMBOSS software exists for pepstats
    path_exists = os.access(PEPSTATS_PATH, os.F_OK)
    if not path_exists:
        print
        print "Path to EMBOSS software does not exist!"
        print "Check the installation and the given path to the EMBOSS software %s in EffectorP.py (line 48)."%PEPSTATS_PATH
        print
        sys.exit()
    # -----------------------------------------------------------------------------------------------------------
    commandline = sys.argv[1:]
    # -----------------------------------------------------------------------------------------------------------
    if commandline:
        FASTA_FILE, short_format, output_file, effector_output = functions.scan_arguments(commandline)
	# If no FASTA file was provided with the -i option
        if not FASTA_FILE:
            print
            print 'Please specify a FASTA input file using the -i option!'
            functions.usage()
    else:
        functions.usage()
    # -----------------------------------------------------------------------------------------------------------
    # Temporary folder name identifier that will be used to store results
    FOLDER_IDENTIFIER = str(uuid.uuid4())
    # Path to temporary results folder 
    if not os.path.exists(SCRIPT_PATH + '/tmp/'):
        os.makedirs(SCRIPT_PATH + '/tmp/')
    RESULTS_PATH = SCRIPT_PATH + '/tmp/' + FOLDER_IDENTIFIER + '/'
    # -----------------------------------------------------------------------------------------------------------
    # Check if FASTA file exists
    try:
        open(FASTA_FILE, 'r') 
    except IOError as (errno, strerror):
        print "Unable to open FASTA file:", FASTA_FILE  #Does not exist OR no read permissions
        print "I/O error({0}): {1}".format(errno, strerror)
        sys.exit()
    # -----------------------------------------------------------------------------------------------------------
    # Try to create folder where results will be stored
    try:
        os.mkdir(RESULTS_PATH)
    except OSError as exception:        
        if exception.errno != errno.EEXIST:
            raise
    # -----------------------------------------------------------------------------------------------------------
    # Extract the identifiers and sequences from input FASTA file
    ORIGINAL_IDENTIFIERS, SEQUENCES = functions.get_seqs_ids_fasta(FASTA_FILE)
    # -----------------------------------------------------------------------------------------------------------
    print '-----------------'
    print 
    print "EffectorP is running for", len(ORIGINAL_IDENTIFIERS), "proteins given in FASTA file", FASTA_FILE
    print
    # -----------------------------------------------------------------------------------------------------------
    # Write new FASTA file with short identifiers because pepstats can't handle long names
    f_output = RESULTS_PATH + FOLDER_IDENTIFIER + '_short_ids.fasta'
    SHORT_IDENTIFIERS = functions.write_FASTA_short_ids(f_output, ORIGINAL_IDENTIFIERS, SEQUENCES)
    # -----------------------------------------------------------------------------------------------------------
    # Call pepstats
    print 'Call pepstats...'
    ProcessExe = PEPSTATS_PATH + 'pepstats'
    ParamList = [ProcessExe, '-sequence', RESULTS_PATH + FOLDER_IDENTIFIER + '_short_ids.fasta', 
              '-outfile', RESULTS_PATH + FOLDER_IDENTIFIER + '.pepstats']

    try:
        Process = subprocess.Popen(ParamList, shell=False)
        sts = Process.wait()
        cstdout, cstderr = Process.communicate()

        if Process.returncode:
            raise Exception("Calling pepstats returned %s"%Process.returncode)
        if cstdout:
            pass
        elif cstderr:
            sys.exit()
    except:
        e = sys.exc_info()[1]
        print "Error calling pepstats: %s" % e
        sys.exit()
    print 'Done.'
    print
    # -----------------------------------------------------------------------------------------------------------
    # Parse pepstats file
    print 'Scan pepstats file'
    pepstats_dic = functions.pepstats(SHORT_IDENTIFIERS, SEQUENCES, RESULTS_PATH + FOLDER_IDENTIFIER + '.pepstats')
    print 'Done.'
    print
    # -----------------------------------------------------------------------------------------------------------
    # Write the WEKA arff file for classification of the input FASTA file
    weka_input = RESULTS_PATH + FOLDER_IDENTIFIER + '.arff'
    functions.write_weka_input(weka_input, SHORT_IDENTIFIERS, pepstats_dic)
    # -----------------------------------------------------------------------------------------------------------
    # Call WEKA Naive Bayes model for classification of input FASTA file
    print 'Start classification with EffectorP...'

    ParamList = ['java', '-cp', WEKA_PATH, 'weka.classifiers.bayes.NaiveBayes', '-l', SCRIPT_PATH + '/trainingdata_samegenomes_iteration15_ratio3_bayes.model',
             '-T', RESULTS_PATH + FOLDER_IDENTIFIER + '.arff', '-p', 'first-last']

    with open(RESULTS_PATH + FOLDER_IDENTIFIER + '_Predictions.txt', 'wb') as out:
        try:
            Process = subprocess.Popen(ParamList, shell=False, stdout=out)
            sts = Process.wait()
            cstdout, cstderr = Process.communicate()

            if Process.returncode:
                raise Exception("Calling WEKA returned %s"%Process.returncode)
            if cstdout:
                pass
            elif cstderr:
                sys.exit()
        except:
            e = sys.exc_info()[1]
            print "Error calling WEKA: %s" % e
            sys.exit()
        print 'Done.'
        print
        print '-----------------'
    # -----------------------------------------------------------------------------------------------------------
    # Parse the WEKA output file
    file_input = RESULTS_PATH + FOLDER_IDENTIFIER + '_Predictions.txt'
    predicted_effectors, predictions = functions.parse_weka_output(file_input, ORIGINAL_IDENTIFIERS, SEQUENCES)
    # -----------------------------------------------------------------------------------------------------------    
    # If user wants the stdout output directed to a specified file
    if output_file:

        with open(output_file, 'wb') as out:
            # Short format: output predictions for all proteins as tab-delimited table
            if short_format:
                out.writelines(functions.short_output(predictions))
            # If the user wants to see the long format, output additional information and stats
            else:
                out.writelines(functions.short_output(predictions))
                out.writelines(functions.long_output(ORIGINAL_IDENTIFIERS, predicted_effectors))                
        print 'EffectorP results were saved to output file:', output_file  

    else:
        # Short format: output predictions for all proteins as tab-delimited table to stdout
        if short_format:
            print functions.short_output(predictions)
        # If the user wants to see the long format, output additional information and stats
        else:
            print functions.short_output(predictions)
            print functions.long_output(ORIGINAL_IDENTIFIERS, predicted_effectors)
    # -----------------------------------------------------------------------------------------------------------
    # If the user additionally wants to save the predicted effectors in a provided FASTA file
    if effector_output:
        with open(effector_output, 'w') as f_output:
            for effector, prob, sequence in predicted_effectors:
                f_output.writelines('>' + effector + ' | Effector probability: ' + str(prob) + '\n')
                f_output.writelines(sequence + '\n')  
    # -----------------------------------------------------------------------------------------------------------
    # Clean up and delete temporary folder that was created
    shutil.rmtree(RESULTS_PATH)
    # -----------------------------------------------------------------------------------------------------------

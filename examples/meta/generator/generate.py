#!/usr/bin/python

import os
from parse import parse
from translate import translate, TranslationFailure
import json
import argparse


def subfilesRelative(directory, filter_by):
    """ Returns a generator of pairs (subpath, filename) in a given directory.
        Filenames are filtered by the filter_by function.
    """
    for dir_, _, files in os.walk(directory):
        for file in files:
            if filter_by(file):
                yield os.path.relpath(dir_, inputDir), file


def parseCtags(filename):
    tags = {}
    if not os.path.exists(filename):
        raise Exception('Failed to found ctags file at %s' % (filename))
    with open(filename) as file:
        for line in file:
            # in ctags format it is TAG\tLOCATION\t..
            symbol, location = line.strip().split('\t')[:2]
            # we assume sources are in src/
            tags[symbol] = location.split('src/')[-1]
    return tags


def translateExamples(inputDir, outputDir, targetsDir, ctagsFile,
                      includedTargets=None, storeVars=False):

    tags = parseCtags(ctagsFile)
    # Load all target dictionaries
    targets = []
    for target in os.listdir(targetsDir):
        # Ignore targets not in includedTargets
        if includedTargets and not os.path.basename(target).split(".")[0] in includedTargets:
            continue

        translate_file = os.path.join(targetsDir, target)

        with open(translate_file) as tFile:
            try:
                targets.append(json.load(tFile))
            except Exception as err:
                print("Error translating file: %s\n%s" % (translate_file, err))

    # Translate each example
    for dirRelative, filename in subfilesRelative(inputDir, filter_by=lambda x: x.lower().endswith('.sg')):
        print("Translating {}".format(os.sep.join([dirRelative, filename])))

        # Parse the example file
        with open(os.path.join(inputDir, dirRelative, filename), 'r') as file:
            ast = parse(file.read(), os.path.join(dirRelative, filename))

        # Translate ast to each target language
        for target in targets:
            try:
                translation = translate(ast, targetDict=target,
                                        tags=tags,
                                        # TODO: enable storeVars for other targets
                                        storeVars=storeVars if target['FileExtension'] == ".cpp" else False)
            except TranslationFailure as e:
                raise Exception("Could not translate file " +
                                os.path.join(inputDir, dirRelative, filename) +
                                " to " + target['FileExtension'] +
                                ".\n" + str(e))

            directory = os.path.join(outputDir, target["OutputDirectoryName"])
            extension = target["FileExtension"]

            # Create directory if it does not exist
            if not os.path.exists(directory):
                os.makedirs(directory)

            # Write translation
            outputFile = os.path.join(directory,
                                      dirRelative,
                                      os.path.splitext(filename)[0] + extension)

            # create subdirectories if they don't exist yet
            try:
                os.makedirs(os.path.dirname(outputFile))
            except OSError:
                pass

            with open(outputFile, "w") as nf:
                nf.write(translation)

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-o", "--output", help="path to output directory")
    parser.add_argument("-i", "--input", help="path to examples directory (input)")
    parser.add_argument("-t", "--targetsfolder", help="path to directory with target JSON files")
    parser.add_argument("-g", "--ctags", help="path to ctags file")
    parser.add_argument("--store-vars",
                        help="whether to store all variables for testing",
                        action='store_true')
    available_targets = [
        'cpp',
        'python',
        'java',
        'r',
        'octave',
        'csharp',
        'ruby',
        'lua',
    ]
    parser.add_argument('targets', nargs='*', help="Targets to include (one or more of: %s). If not specified all targets are produced." % (' '.join(available_targets)))

    args = parser.parse_args()

    outputDir = "outputs"
    if args.output:
        outputDir = args.output

    inputDir = "examples"
    if args.input:
        inputDir = args.input

    targetsDir = "targets"
    if args.targetsfolder:
        targetsDir = args.targetsfolder

    ctagsFile = args.ctags
    storeVars = True if args.store_vars else False

    translateExamples(inputDir=inputDir, outputDir=outputDir,
                      targetsDir=targetsDir, ctagsFile=ctagsFile,
                      includedTargets=args.targets, storeVars=storeVars)

# SPDX-License-Identifier: GPL-3.0-or-later
# SPDX-FileCopyrightText: Copyright © VirtualFluids Project contributors, see AUTHORS.md in root folder

#!/usr/bin/env python3

import json
import sys


def compile_command_selector(x):
    return not ("3rdParty" in x["file"] or ".cu" in x["file"])


if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("usage: ./filterCompileCommands.py compile_commands.json")
        exit(-1)

    filename = sys.argv[1]
    print("loading compile commands file: {}".format(filename))

    fin = open(filename, "r")
    cc = json.load(fin)
    fin.close()

    print("compile commands read: {}".format(len(cc)))

    cc_filtered = list(filter(compile_command_selector, cc))

    print("compile commands filtered: {}".format(len(cc_filtered)))

    fout = open(filename, "w")
    json.dump(cc_filtered, fout, indent=4)
    fout.close()

# SPDX-License-Identifier: GPL-3.0-or-later
# SPDX-FileCopyrightText: Copyright © VirtualFluids Project contributors, see AUTHORS.md in root folder

from pathlib import Path
from os.path import join

vf_header = """//=======================================================================================
// ____          ____    __    ______     __________   __      __       __        __
// \    \       |    |  |  |  |   _   \  |___    ___| |  |    |  |     /  \      |  |
//  \    \      |    |  |  |  |  |_)   |     |  |     |  |    |  |    /    \     |  |
//   \    \     |    |  |  |  |   _   /      |  |     |  |    |  |   /  /\  \    |  |
//    \    \    |    |  |  |  |  | \  \      |  |     |   \__/   |  /  ____  \   |  |____
//     \    \   |    |  |__|  |__|  \__\     |__|      \________/  /__/    \__\  |_______|
//      \    \  |    |   ________________________________________________________________
//       \    \ |    |  |  ______________________________________________________________|
//        \    \|    |  |  |         __          __     __     __     ______      _______
//         \         |  |  |_____   |  |        |  |   |  |   |  |   |   _  \    /  _____)
//          \        |  |   _____|  |  |        |  |   |  |   |  |   |  | \  \   \_______
//           \       |  |  |        |  |_____   |   \_/   |   |  |   |  |_/  /    _____  |
//            \ _____|  |__|        |________|   \_______/    |__|   |______/    (_______/
//
//  This file is part of VirtualFluids. VirtualFluids is free software: you can
//  redistribute it and/or modify it under the terms of the GNU General Public
//  License as published by the Free Software Foundation, either version 3 of
//  the License, or (at your option) any later version.
//
//  VirtualFluids is distributed in the hope that it will be useful, but WITHOUT
//  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
//  FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
//  for more details.
//
//  You should have received a copy of the GNU General Public License along
//  with VirtualFluids (see COPYING.txt). If not, see <http://www.gnu.org/licenses/>.
//
"""

vf_header_new = """//=======================================================================================
// ____          ____    __    ______     __________   __      __       __        __
// \    \       |    |  |  |  |   _   \  |___    ___| |  |    |  |     /  \      |  |
//  \    \      |    |  |  |  |  |_)   |     |  |     |  |    |  |    /    \     |  |
//   \    \     |    |  |  |  |   _   /      |  |     |  |    |  |   /  /\  \    |  |
//    \    \    |    |  |  |  |  | \  \      |  |     |   \__/   |  /  ____  \   |  |____
//     \    \   |    |  |__|  |__|  \__\     |__|      \________/  /__/    \__\  |_______|
//      \    \  |    |   ________________________________________________________________
//       \    \ |    |  |  ______________________________________________________________|
//        \    \|    |  |  |         __          __     __     __     ______      _______
//         \         |  |  |_____   |  |        |  |   |  |   |  |   |   _  \    /  _____)
//          \        |  |   _____|  |  |        |  |   |  |   |  |   |  | \  \   \_______
//           \       |  |  |        |  |_____   |   \_/   |   |  |   |  |_/  /    _____  |
//            \ _____|  |__|        |________|   \_______/    |__|   |______/    (_______/
//
//  This file is part of VirtualFluids. VirtualFluids is free software: you can
//  redistribute it and/or modify it under the terms of the GNU General Public
//  License as published by the Free Software Foundation, either version 3 of
//  the License, or (at your option) any later version.
//
//  VirtualFluids is distributed in the hope that it will be useful, but WITHOUT
//  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
//  FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License
//  for more details.
//
//  SPDX-License-Identifier: GPL-3.0-or-later
//  SPDX-FileCopyrightText: Copyright © VirtualFluids Project contributors, see AUTHORS.md in root folder
//
"""


file_endings = ["**/*.h", "**/*.cpp", "**/*.hpp", "**/*.cu", "**/*.cuh", "**/*.c"]

def add_to_begin_of_file(folder, header):
    path_files = []

    for ending in file_endings:
        path_files.extend(Path(folder).glob(join(ending)))

    files = [x for x in path_files if x.is_file()]

    for file in files:
        with open(file, "r") as in_file:
            file_content = in_file.read()

            with open(file, "w+") as in_file:
                file_content = header + file_content

                in_file.write(file_content)
                print(f"File modified: {file}")

# this heler function can be used to remove all lines containg \ingroup and \file modifer
# this might be helpful when doxygen changed
def remove_file_and_ingroup(folder):
    path_files = []

    for ending in file_endings:
        path_files.extend(Path(folder).glob(join(ending)))

    files = [x for x in path_files if x.is_file()]

    for file in files:
        # if (not file.name == 'ConfigurationFile.h'):
        #     continue
        with open(file, "r") as in_file:
            lines = in_file.readlines()

            with open(file, "w+") as in_file:
                new_lines = []
                for line in lines:
                    if "\ingroup" in line.strip() or "\\file" in line.strip():
                        print(f"modifing: {file}")
                    else:
                        new_lines.append(line)
                in_file.writelines(new_lines)
                    


def replace_line_containting(folder, old_lines: list[str], new_line_token: list[str]):
    path_files = []

    for ending in file_endings:
        path_files.extend(Path(folder).glob(join(ending)))

    files = [x for x in path_files if x.is_file()]

    for file in files:
        # if (not file.name == 'ConfigurationFile.h'):
        #     continue
        with open(file, "r") as in_file:
            lines = in_file.readlines()

            with open(file, "w+") as in_file:
                new_lines = []
                for line in lines:
                    new_value = line
                    for idx, old in enumerate(old_lines):
                        if old in line:
                            new_value = new_line_token[idx]
                            print(f"{file}: Found: {old}, Replace with {new_value}")
                            break
                    new_lines.append(new_value)

                    # if "\ingroup" in line.strip() or "\\file" in line.strip():
                    #     print(f"modifing: {file}")
                    # else:
                    #     new_lines.append(line)
                in_file.writelines(new_lines)
                    


def add_header(root_folder):
    path_files = []

    for ending in file_endings:
        path_files.extend(Path(root_folder).glob(join(ending)))

    files = [x for x in path_files if x.is_file()]

    for file in files:

        with open(file, "r") as in_file:

            file_content = in_file.read()

            with open(file, "w+") as in_file:
                file_content = vf_header_new + file_content

                in_file.write(file_content)
                print(f"File modified: {file}")


# Finds wrong (vf_header) header starting from `root_folder`
# When `exchange_header` is true, vf_header is exchanged with vf_header_new
def find_missing_header(root_folder, exchange_header = False):
    path_files = []

    for ending in file_endings:
        path_files.extend(Path(root_folder).glob(join(ending)))

    files = [x for x in path_files if x.is_file()]

    for file in files:

        with open(file, "r") as in_file:

            file_content = in_file.read()

            if(file_content.startswith(vf_header)):
                print(f"wrong vf header: {file}")
                if(exchange_header):
                    with open(file, "w+") as in_file:
                        file_content = file_content.replace(vf_header, vf_header_new)

                        in_file.write(file_content)
                        print(f"File modified: {file}")
            

# add doygen group directly after the vf_header
def add_doxygen_group(folder, prefix = ""):
    path_files = []

    for ending in file_endings:
        path_files.extend(Path(folder).glob(ending))


    files = [x for x in path_files if x.is_file()]

    # add headers
    for file in files:

        # tests are handled different
        if (str(file.stem).endswith("Test")):
            continue
        # if (not file.name == 'VirtualFluidSimulationFactory.h'):
        #     continue

        print(f"Try modifing: {file}")
        with open(file, "r") as in_file:

            file_content = in_file.read()

            if(not file_content.startswith(vf_header)):
                print(f"File has no vf header: {file}")
            else:
                with open(file, "w+") as in_file:
                    parents = str(file.parents[0])
                    splits = parents.split("/")

                    toplevel_group_subfolder_number = 1 # first folder after "src/". e.g. src/lbm -> lbm
                    if (prefix):
                        toplevel_group_subfolder_number = toplevel_group_subfolder_number + 1 # "src/cpu/core" -> core

                    toplevel_group = splits[toplevel_group_subfolder_number]
                    if (prefix):
                        toplevel_group = prefix + "_" + toplevel_group + " " + toplevel_group

                    group_string = f"//! \\addtogroup {toplevel_group}"
                    if len(splits) > toplevel_group_subfolder_number + 1:
                        group = splits[toplevel_group_subfolder_number + 1]
                        if (prefix):
                            group = prefix + "_" + group + " " + group
                        group_string = f"""//! \\addtogroup {group}
//! \ingroup {toplevel_group}"""

                    group_string = group_string + "\n//! \\{\n"

                    print(f"toplevel: {toplevel_group}, Group: {group}")


                    vf_header_new = vf_header + group_string
                    file_content = file_content.replace(vf_header, vf_header_new)

                    end_group_string = "\n//! \}\n"

                    file_content = file_content + end_group_string

                    in_file.write(file_content)
                    print(f"File modified: {file}")


# add doygen group directly after the vf_header
def find_test_files(folder, prefix = ""):
    path_files = []

    for ending in file_endings:
        path_files.extend(Path(folder).glob(ending))

    files = [x for x in path_files if x.is_file()]

    # add headers
    for file in files:
        if (str(file.stem).endswith("Test")):
            print(file)


metis_identifier = """// SPDX-License-Identifier: Apache-2.0
// SPDX-FileCopyrightText: Copyright 1995-2013, Regents of the University of Minnesota
"""

if __name__ == "__main__":
    folder = "src/"
    # find_missing_header(folder, exchange_header = True)
    # add_header(folder)
    # remove_file_and_ingroup(folder)
    # add_doxygen_group(folder)
    find_test_files(folder)
    # replace_line_containting("apps/", ["//! \{{", "//! \}}"], ["//! \{\n", "//! \}\n"])


    # add_to_begin_of_file("3rdParty/metis/", metis_identifier)

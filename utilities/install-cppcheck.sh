# SPDX-License-Identifier: GPL-3.0-or-later
# SPDX-FileCopyrightText: Copyright © VirtualFluids Project contributors, see AUTHORS.md in root folder

# !/usr/bin/env bash
set -e

cd /tmp
git clone https://github.com/danmar/cppcheck.git
cd cppcheck 
git checkout 2.10.3
make MATCHCOMPILER=yes FILESDIR=/usr/share/cppcheck HAVE_RULES=yes CXXFLAGS="-O2 -DNDEBUG -Wall -Wno-sign-compare -Wno-unused-function" install
cd /tmp
rm -rf /tmp/cppcheck
ldconfig
cppcheck --version
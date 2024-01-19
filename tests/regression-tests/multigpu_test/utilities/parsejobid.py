# SPDX-License-Identifier: GPL-3.0-or-later
# SPDX-FileCopyrightText: Copyright © VirtualFluids Project contributors, see AUTHORS.md in root folder

from pathlib import Path
import sys

LAUNCH_MESSAGE = "Launched job"


def parsejobid(file: str) -> str:
    file_path = Path(file)
    if not file_path.exists():
        raise FileNotFoundError(file)

    text_content = file_path.read_text().strip()
    launch_line = next(
        filter(lambda line: LAUNCH_MESSAGE in line, text_content.splitlines())
    )
    return launch_line.split()[-1].strip()


if __name__ == "__main__":
    print(parsejobid(sys.argv[1]))

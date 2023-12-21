# SPDX-License-Identifier: GPL-3.0-or-later
# SPDX-FileCopyrightText: Copyright © VirtualFluids Project contributors, see AUTHORS.md in root folder

from jinja2 import Template
from pathlib import Path

TEMPLATES_DIR = Path(__file__).parent
REGRESSION_CI_TEMPLATE = TEMPLATES_DIR / "regression-tests-ci.yml.j2"
GENERATED_DIR = Path("generated")
REGRESSION_CI_FILE = GENERATED_DIR / "regression-tests-ci.yml"
TEST_FILE_DIR = Path("regression-tests")


def build_regression_job_string(regression_tests: list[(str,str)]) -> str:
    template = Template(REGRESSION_CI_TEMPLATE.read_text())
    return template.render(regression_tests=regression_tests)

def trim_parent_path(name: str) -> str:
    return name.replace(str(TEST_FILE_DIR)+"/", "")

def main():
    regression_tests_files = [(item.stem, trim_parent_path(str(item.parent))) for item in TEST_FILE_DIR.rglob("*_test.sh")]
    print(regression_tests_files)
    regression_tests_ci_file = build_regression_job_string(regression_tests_files)
    REGRESSION_CI_FILE.write_text(regression_tests_ci_file)

if __name__ == "__main__":
    GENERATED_DIR.mkdir(parents=True, exist_ok=True)
    main()

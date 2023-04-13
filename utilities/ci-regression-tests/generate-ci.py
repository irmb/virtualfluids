from dataclasses import dataclass
from jinja2 import Template
from pathlib import Path

TEMPLATES_DIR = Path(__file__).parent
REGRESSION_CI_TEMPLATE = TEMPLATES_DIR / "regression-tests-ci.yml.j2"
GENERATED_DIR = Path("generated")
REGRESSION_CI_FILE = GENERATED_DIR / "regression-tests-ci.yml"



def build_regression_job_string(regression_tests: list[str]) -> str:
    template = Template(REGRESSION_CI_TEMPLATE.read_text())
    return template.render(regression_tests=regression_tests)


def main():

    root_dir = "regression-tests"
    root_dir_instance = Path(root_dir)
    regression_tests_files = [item.stem for item in root_dir_instance.glob("*_test.sh")]
    print(regression_tests_files)


    regression_tests_ci_file = build_regression_job_string(regression_tests_files)
    REGRESSION_CI_FILE.write_text(regression_tests_ci_file)

if __name__ == "__main__":
    GENERATED_DIR.mkdir(parents=True, exist_ok=True)
    main()
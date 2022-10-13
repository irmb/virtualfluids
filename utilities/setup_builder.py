from setuptools import build_meta

class builder(build_meta._BuildMetaBackend):

    def run_setup(self, setup_script='setup.py'):
        # Note that we can reuse our build directory between calls
        # Correctness comes first, then optimization later
        __file__ = setup_script
        __name__ = '__main__'

        with build_meta._open_setup_script(__file__) as f:
            code = f.read().replace(r'\r\n', r'\n')
        args = locals()
        args["cmake_args"] = self.extra_args
        exec(code, args)


    def add_settings(self, config_settings):
        self.extra_args = dict()
        print(config_settings)
        if config_settings:
            self.extra_args = {k:v for k,v in config_settings.items() if k[:2] == "-D"}

    def build_wheel(self, wheel_directory, config_settings=None,
                    metadata_directory=None):
        self.add_settings(config_settings)
        return super().build_wheel(wheel_directory, config_settings, metadata_directory)

    def build_sdist(self, sdist_directory, config_settings=None):
        self.add_settings(config_settings)
        return super().build_wheel(sdist_directory, config_settings)

build = builder()
build_wheel = build.build_wheel
build_sdist = build.build_sdist
try:
    from .bindings import Parser
except ImportError as e:
    raise ImportError("Pymuparser bindings were not built. Only included if VirtualFluids is built with VF_BUILD_CPU=ON.")
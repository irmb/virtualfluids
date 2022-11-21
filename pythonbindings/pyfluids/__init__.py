try:
    from .bindings import basics
except ImportError:
    print("Basics bindings not included")
try:
    from .bindings import logger
except ImportError:
    print("Logger bindings not included")
try:
    from .bindings import lbm
except ImportError:
    print("LBM bindings not included")
try:
    from .bindings import gpu
except ImportError:
    print("GPU bindings not included")
try:
    from .bindings import cpu
except ImportError:
    print("CPU bindings not included")
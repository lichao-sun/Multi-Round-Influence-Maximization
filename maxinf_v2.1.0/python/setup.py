# Setep by distutil
# call:
#   python setup.py build_ext

import distutils
from distutils.core import setup, Extension
import os
import platform
from os.path import join as pjoin


print('os.name=%s   platform.system()=%s' % (os.name, platform.system()))

# Hack for Mac OS (Darwin) / Windows
extra_compile_args = []
if platform.system()=="Darwin":
    os.environ["CC"] = "clang++"
    # use c++11, and exception handlers
    extra_compile_args += ["-std=gnu++11", "-stdlib=libc++", #"-nostdinc++",
        "-fexceptions" ]
elif platform.system()=="Windows":
    # for: cl.exe
    # (1) Indicate: C++ exception handler is used
    #     See: https://msdn.microsoft.com/en-us/library/2axwkyt4.aspx
    # (2) Indicate: openmp is used
    extra_compile_args += ["/EHsc", "/openmp"]


DIR_CODE = pjoin("..", "code")
DIR_CODE_PYTHON = pjoin(DIR_CODE, "python")

# swig interfaces:
swig_interfaces = ["port_mi.i"]

# cpp files:
cppsources = [pjoin(DIR_CODE, f) for f in [
        "common.cpp",
        "compatible.cpp",
        "mi_random.cpp",
        "mi_util.cpp",
        "graph_basic.cpp",
        "graph.cpp",
        "algo_base.cpp",
        "graph_stat.cpp",
        "random_pick.cpp",
        "degree.cpp",
        "degreediscount_ic.cpp",
        "weighted_degree.cpp",
        "pagerank.cpp",
        "mia.cpp",
        "pmia.cpp",
        "SPM_gc.cpp",
        "SP1M_gc.cpp",
        "independ_cascade.cpp",
        "general_cascade.cpp",
        "greedy.cpp",
        "greedy_online.cpp",
        "topic_aware.cpp",
        "cgreedy.cpp",
        "top_selection.cpp",
        "mis.cpp",
        "rr_infl.cpp",
        "reverse_general_cascade.cpp",
        "contcascade.cpp",
        "simulate.cpp",
        "mi_command_line.cpp",
    ]
]

# headers:
module = Extension(
    "_max_influence", 
    sources=swig_interfaces + cppsources,
    language="c++",
    swig_opts=["-c++"],
    extra_compile_args=extra_compile_args,
    #extra_compile_args=["-std=c++11"],
)
# python modules
py_modules = ["max_influence", "test"]

setup(name = "Wrap of max_influence algorithms",
    version = "1.0",
    author = 'Tian Lin',
    author_email = 'bolitt@gmail.com',
    description="""Port max_influence package from c++ to python. (use swig)""",
    ext_modules = [module],
    py_modules = py_modules, 
)
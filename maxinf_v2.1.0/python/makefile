# makefile for swig
# global variables
# 
# makefile help:
# http://wiki.ubuntu.org.cn/%E8%B7%9F%E6%88%91%E4%B8%80%E8%B5%B7%E5%86%99Makefile:MakeFile%E4%BB%8B%E7%BB%8D#.E4.B8.80.E4.B8.AA.E7.A4.BA.E4.BE.8B

# default is gcc
cc = clang++
dir = ../code/
# define SWIG flag for compile, and use optimized code
compile_opts = -D SWIG -O3
build_dir = build/
obj_dir = $(build_dir)obj/

swig_def = port_mi.i
wrap_cpp = port_mi_wrap.cpp 
cpp = $(dir)common.cpp $(dir)compatible.cpp \
    $(dir)mi_random.cpp $(dir)mi_util.cpp \
    $(dir)graph_basic.cpp $(dir)graph.cpp $(dir)algo_base.cpp $(dir)graph_stat.cpp  \
    $(dir)random_pick.cpp $(dir)degree.cpp $(dir)degreediscount_ic.cpp \
    $(dir)weighted_degree.cpp $(dir)pagerank.cpp $(dir)mia.cpp $(dir)pmia.cpp \
    $(dir)SPM_gc.cpp $(dir)SP1M_gc.cpp \
    $(dir)independ_cascade.cpp $(dir)general_cascade.cpp \
    $(dir)greedy.cpp $(dir)greedy_online.cpp \
    $(dir)topic_aware.cpp $(dir)cgreedy.cpp $(dir)top_selection.cpp $(dir)mis.cpp \
    $(dir)reverse_general_cascade.cpp $(dir)rr_infl.cpp \
    $(dir)contcascade.cpp \
    $(dir)simulate.cpp $(dir)mi_command_line.cpp 
wrap_objects = port_mi_wrap.o
objects = common.o compatible.o \
    mi_random.o mi_util.o \
    graph_basic.o graph.o algo_base.o graph_stat.o \
    random_pick.o degree.o degreediscount_ic.o \
    weighted_degree.o pagerank.o mia.o pmia.o \
    SPM_gc.o SP1M_gc.o \
    independ_cascade.o general_cascade.o \
    greedy.o greedy_online.o \
    topic_aware.o cgreedy.o top_selection.o mis.o \
    reverse_general_cascade.o rr_infl.o \
    contcascade.o \
    simulate.o mi_command_line.o 
# replace pyheaders
# pyheaders = /Library/Frameworks/Python.framework/Versions/2.7/include/python2.7
target = _max_influence.so

#default compile
python: gen compile link post

gen:
	swig -c++ -python -o $(wrap_cpp) $(swig_def)

compile:
	# $(cc) -c $(wrap_cpp) $(cpp) -I$(pyheaders)
	$(cc) -c `python-config --cflags` $(compile_opts) $(wrap_cpp) $(cpp)
	# -mv *.o $(objdir)

link:
	$(cc) -bundle `python-config --ldflags` $(wrap_objects) $(objects) -o $(target)

post:
	-mkdir -p $(obj_dir)
	-mv  *.o       $(obj_dir)
	-mv  *.so      $(build_dir)
	-cp  test.py   $(build_dir)
	-cp  max_influence.py   $(build_dir)
	-cp  *.inf     $(build_dir)


PHONY : clean

clean :
	-rm  *.c  *.cxx  *.cpp  *.pyc  *.o  *.so  *.txt max_influence.py
	-rm  -fR $(obj_dir)
	-rm  -fR $(build_dir)

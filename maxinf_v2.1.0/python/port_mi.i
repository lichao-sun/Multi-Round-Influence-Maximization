// Interfacing algorithms to other languages
// Definitions of swig module
// 
// Tian Lin, 2016/01/29
// 
// swig helps:
// http://www.swig.org/translations/chinese/tutorial.html
// http://www.swig.org/Doc1.3/Python.html#Python_nn10


%module max_influence
%{

#include "../code/common.h"
#include "../code/compatible.h"
#include "../code/mi_util.h"

#include "../code/mi_random.h"
#include "../code/event_timer.h"

#include "../code/graph_basic.h"
#include "../code/graph.h"
#include "../code/cascade.h"
#include "../code/independ_cascade.h"
#include "../code/general_cascade.h"

#include "../code/algo_base.h"
#include "../code/graph_stat.h"
#include "../code/random_pick.h"
#include "../code/degree.h"
#include "../code/degreediscount_ic.h"
#include "../code/weighted_degree.h"
#include "../code/pagerank.h"
#include "../code/mia.h"
#include "../code/pmia.h"
#include "../code/SPM_gc.h"
#include "../code/SP1M_gc.h"

#include "../code/greedy.h"
#include "../code/greedy_online.h"

// Topic-aware:
#include "../code/topic_aware.h"
#include "../code/cgreedy.h"
#include "../code/top_selection.h"
#include "../code/mis.h"

// RR Infl:
#include "../code/reverse_general_cascade.h"
#include "../code/rr_infl.h"


// Continous-time:
#include "../code/contcascade.h"

// Algorithms Command Line:
#include "../code/simulate.h"
#include "../code/mi_command_line.h"

%}


// export by definitions
// extern int factorial(int n);
// Or we can parse the header file to generate wrappers
// %include "header.h"

////////////////////////////////////////////////////////////
// Instantiation order is important.

// For STLs:
// According to the swig documentation, using std::string requires %include "std_string.i" or directly include stl.i
%include "stl.i"
%include "std_string.i"
%include "std_vector.i"

// instantiate template classes
namespace std {
	%template(Vector_d) std::vector<double>;
	%template(Vector_f) std::vector<float>;
	%template(Vector_i) std::vector<int>;
	%template(Vector_s) std::vector<std::string>;
}


// For common headers and compatibility

%include "../code/common.h"
%include "../code/compatible.h"
%include "../code/mi_util.h"

// For Basic helpers
%include "../code/mi_random.h"
%include "../code/event_timer.h"
%template(EventTimer) PCTimerT<std::string>;


// For Graph basic
%feature("novaluewrapper") BaseEdge;
%feature("novaluewrapper") Edge;
%feature("novaluewrapper") ContEdge;

%include "../code/graph_basic.h"

// For Graph
%include "../code/graph.h"

%template(Graph) GraphT<Edge>;
%template(ContGraph) GraphT<ContEdge>;
%template(GraphFactory) GraphFactoryT<Edge>;
%template(ContGraphFactory) GraphFactoryT<ContEdge>;


// For Cascades
%include "../code/cascade.h"
%template(CascadeT_Graph) CascadeT<Graph>;
%template(CascadeT_ContGraph) CascadeT<ContGraph>;  // the order is important!


%include "../code/independ_cascade.h"
%include "../code/general_cascade.h"
%template(GeneralCascade) GeneralCascadeT<Graph>;


%include "../code/contcascade.h"


%include "../code/algo_base.h"

// For statistics:
%include "../code/graph_stat.h"
%template(GraphStatistics) GraphStatisticsT<Graph>;


// For Algorithms:
%include "../code/random_pick.h"
%include "../code/degree.h"
%include "../code/degreediscount_ic.h"
%include "../code/weighted_degree.h"
%include "../code/pagerank.h"
%include "../code/mia.h"
%include "../code/pmia.h"
%include "../code/SPM_gc.h"
%include "../code/SP1M_gc.h"

%include "../code/greedy.h"
%include "../code/greedy_online.h"

%include "../code/topic_aware.h"
%include "../code/cgreedy.h"
%include "../code/top_selection.h"
%include "../code/mis.h"

%include "../code/reverse_general_cascade.h"
%template(ReverseGCascade) ReverseGCascadeT<Graph>;
%include "../code/rr_infl.h"

%include "../code/simulate.h"
%include "../code/mi_command_line.h"



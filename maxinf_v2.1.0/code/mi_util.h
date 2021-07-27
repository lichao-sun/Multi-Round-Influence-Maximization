//
//  mi_util.h
//  max_influence
//
//  Created by Tian on 1/18/16.
//  Copyright © 2016 Tian. All rights reserved.
//

#ifndef mi_util_h__
#define mi_util_h__

#include <iostream>
#include <fstream>
#include <string>
#include <memory>
#include <exception>
#include <cmath>
#include <cassert>



/// Null Pointer Exceptions
class NullPointerException
    : public std::runtime_error
{
public:
    NullPointerException(const std::string& what_arg) : std::runtime_error(what_arg) {}
};

class InvalidInputFormatException
	: public std::runtime_error
{
public:
	InvalidInputFormatException(const std::string& what_arg) : std::runtime_error(what_arg) {}
};

/// Wrapper of standard c++ inputstream for Python
class PyInputStream
{
public:
    std::shared_ptr<std::istream> ptr;
	virtual ~PyInputStream() { ptr.reset(); }
};

/// Wrapper of standard c++ outputstream for Python
class PyOutputStream
{
public:
    std::shared_ptr<std::ostream> ptr;
    virtual ~PyOutputStream() { ptr.reset(); }
};

/// PyHelper for generating standard c++ iostreams
class PyHelper
{
public:
    static std::istream& GetStdin();
    static std::ostream& GetStdout();
    static PyInputStream GetInputFileStream(const char* filename);
    static PyOutputStream GetOutputFileStream(const char* filename);
};

//////////////////////////////////////////////////////////////////////
/// The following definitions are C++ type traits that are used for static assertion
/// for template classes.
/// 
/// One known issue is: 
///   SWIG does not support any nested class, we disable the static check (which is also unnecessary)
///   to make it compilable.
/// 
//////////////////////////////////////////////////////////////////////
// We Use Macro SWIG controls static check:
#ifndef SWIG

/// Traits: true type. (the size is different from false type)
struct true_t {
public:
	enum { size = sizeof(char) };
private:
	char _inner_value;
};

/// Traits: false type. (the size is different from true type)
struct false_t {
	enum { size = sizeof(long) };
private:
	long _inner_value;
};

/// Detect member function, using SFINAE tricks.
/// See: https://en.wikibooks.org/wiki/More_C%2B%2B_Idioms/Member_Detector
///
/// Generate a meta class that can be used to check whether a class T has member and satisfies TFuncType
#define GEN_HAS_MEM_TYPE(ClassName, member)    \
    template <typename T, typename TFuncType>     \
    class ClassName {  \
	protected: \
		template <typename U, U> struct type_check; \
        template <typename K>  static true_t test(type_check<TFuncType, &K::member >*); \
        template <typename K>  static false_t test(...);            \
	public: \
        enum { value = (sizeof(test<T>(NULL)) == true_t::size) }; \
    }

/// Generate a meta class that can be used to check whether a class T has member
#define GEN_HAS_MEM(ClassName, member) \
    template < typename T > \
    class ClassName { \
	protected: \
		struct Fallback { struct member {}; }; \
		struct Derived : T, Fallback { }; \
        template <typename K> static true_t test(decltype(K::member)*); \
        template <typename K> static false_t test(...); \
    public: \
        enum { value = (sizeof(test<Derived>(NULL)) == true_t::size) }; \
    }

/// Check whether cond is satisfied at the compile time.
#define MI_STATIC_ASSERT(cond, msg)  static_assert(cond, msg)

#else // otherwise, macro SWIG is not defined

// disable static assertion
#define GEN_HAS_MEM_TYPE(ClassName, member)
#define GEN_HAS_MEM(ClassName, member)
#define MI_STATIC_ASSERT(cond, msg)

#endif //:~ macro SWIG


/// Static check helpers
/// Check membersFor edges
GEN_HAS_MEM(has_mem_u, u);
GEN_HAS_MEM(has_mem_v, v);
GEN_HAS_MEM(has_mem_Deserialize, Deserialize);

/// Check members For IGraph
GEN_HAS_MEM(has_mem_GetN, GetN);
GEN_HAS_MEM(has_mem_GetM, GetM);

/// Check members For Graph
GEN_HAS_MEM(has_mem_n, n);
GEN_HAS_MEM(has_mem_m, m);
GEN_HAS_MEM(has_mem_degree, degree);
GEN_HAS_MEM(has_mem_edges, edges);
GEN_HAS_MEM(has_mem_index, index);
GEN_HAS_MEM(has_mem_GetNeighborCount, GetNeighborCount);
GEN_HAS_MEM(has_mem_GetEdge, GetEdge);
GEN_HAS_MEM(has_mem_edgeForm, edgeForm);


#endif // mi_util_h__

///
///  common.h
///
///  This header is to include common macros, exceptions, helpers etc.
///
///  Created by Tian on 1/27/16.
///  Copyright(c) 2016 Tian. All rights reserved.
///

#ifndef __COMMON_H__
#define __COMMON_H__



#include <string>
#include "compatible.h"


/// Software name and version
#define MI_SOFTWARE_NAME		"MaxInfluence"
#define MI_SOFTWARE_MAJOR_VER   2
#define MI_SOFTWARE_VER_STR		"2.1.0"

/// detect system bit len
#if _WIN64 || __x86_64__ || __ppc64__
#  define COMPILE_TIME_BITS "64"
#  define MI_IS_32_BITS  0
#  define MI_IS_64_BITS  1
#else
#  define COMPILE_TIME_BITS "32"
#  define MI_IS_32_BITS  1
#  define MI_IS_64_BITS  0
#endif

/// detect debug / release
#if defined(DEBUG) || defined(_DEBUG)
#  define COMPILE_TIME_VER_STR   "Build: " __DATE__ ", " __TIME__  ", " COMPATIBILITY_THIS_OS_STR ", (" COMPILE_TIME_BITS " bits, Debug)"
#  define MI_IS_DEBUG   1
#  define MI_IS_RELEASE 0
#else /* _DEBUG */
#  define COMPILE_TIME_VER_STR   "Build: " __DATE__ ", " __TIME__  ", " COMPATIBILITY_THIS_OS_STR ", (" COMPILE_TIME_BITS " bits, Release)"
#  define MI_IS_DEBUG   0
#  define MI_IS_RELEASE 1
#endif

#define MI_SOFTWARE_META       MI_SOFTWARE_NAME " (ver " MI_SOFTWARE_VER_STR "), " COMPILE_TIME_VER_STR

/// Function to print the version
void PrintVersion();


/// Define namespaces:
#define NS_MI_BEGIN namespace mi { // begin namespace mi (max_influence)
#define NS_MI_END   } //:~ end: namespace max_influence
#define NS_USE_MI   using namespace mi

#define NS_FUTURES_BEGIN namespace futures {  // begin namespace futures
#define NS_FUTURES_END   } //:~ end: namespace futures
#define NS_USE_FUTURES  using namespace futures

/// Safely delete a pointer to an object
#ifndef SAFE_DELETE
#  define SAFE_DELETE(ptr) if(ptr != NULL) \
                        { delete ptr; ptr = NULL; }
#endif

/// Safely delete a pointer to an array of objects
#ifndef SAFE_DELETE_ARRAY
#  define SAFE_DELETE_ARRAY(ptr) if(ptr != NULL) \
                        { delete[] ptr; ptr = NULL; }
#endif

/// Define getter and setter
#define DefGet(field, Prop)  decltype(field) Get##Prop() { \
	return this->field; \
}
#define DefSet(field, Prop)  void Set##Prop(decltype(field) v) { \
	 this->field = v; \
}
#define DefGetSet(field, Prop)   DefGet(field, Prop) \
	DefSet(field, Prop)


/// Common data type
#ifndef LARGE_INT64
typedef long long int __INT64__;
#define LARGE_INT64  __INT64__
#endif


/// Wrapper of some constants.
class MICommon
{
public:
    enum { IsDebug = MI_IS_DEBUG };
    enum { IsRelease = MI_IS_RELEASE };
    enum { Is32Bits = MI_IS_32_BITS };
    enum { Is64Bits = MI_IS_64_BITS };
    enum { IsOSWin = (COMPATIBILITY_THIS_OS == COMPATIBILITY_OS_WIN) };
    enum { IsOSApple = (COMPATIBILITY_THIS_OS == COMPATIBILITY_OS_APPLE) };
    enum { IsOSLinux = (COMPATIBILITY_THIS_OS == COMPATIBILITY_OS_LINUX) };
    enum { IsOSUnix = (COMPATIBILITY_THIS_OS == COMPATIBILITY_OS_UNIX) };
    enum { IsOSPosix = (COMPATIBILITY_THIS_OS == COMPATIBILITY_OS_POSIX) };
    enum { IsOSUnknown = (COMPATIBILITY_THIS_OS == COMPATIBILITY_OS_UNKNOWN) };
	enum { IsOMPUsed = MI_IS_OMP_USED };
	enum { MajorVer = MI_SOFTWARE_MAJOR_VER };
};


#include "mi_util.h"
#include "mi_limit.h"


#endif //:~ __COMMON_H__

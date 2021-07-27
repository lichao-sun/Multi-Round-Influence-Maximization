///
/// compatible.h
///
///  This header is used to make functions compatible across different platforms.
///
///  Created by Tian on 1/27/16.
///  Copyright(c) 2016 Tian. All rights reserved.
///

#ifndef __compatible_h__
#define __compatible_h__

#include <iostream>
#include <cmath>
#include <cassert>


// check os
#define COMPATIBILITY_OS_WIN     1
#define COMPATIBILITY_OS_APPLE   (1 << 1)
#define COMPATIBILITY_OS_LINUX   (1 << 2)
#define COMPATIBILITY_OS_UNIX    (1 << 3)
#define COMPATIBILITY_OS_POSIX   (1 << 4)
#define COMPATIBILITY_OS_UNKNOWN 0

#if defined(_WIN32)
//define something for Windows (32-bit and 64-bit, this part is common)
#   define COMPATIBILITY_THIS_OS          COMPATIBILITY_OS_WIN
#   if defined(_WIN64)
//define something for Windows (64-bit only)
#       define COMPATIBILITY_THIS_OS_STR   "WIN64"
#   else
#       define COMPATIBILITY_THIS_OS_STR   "WIN32"
#   endif 
#elif __APPLE__
#   define COMPATIBILITY_THIS_OS           COMPATIBILITY_OS_APPLE
#   include "TargetConditionals.h"
#   if TARGET_IPHONE_SIMULATOR
// iOS Simulator
#   define COMPATIBILITY_THIS_OS_STR     "iOS Simulator"
#   elif TARGET_OS_IPHONE
// iOS device
#   define COMPATIBILITY_THIS_OS_STR     "iOS device"
#   elif TARGET_OS_MAC
// Other kinds of Mac OS
#   define COMPATIBILITY_THIS_OS_STR     "Mac OS"
#   else
//  #      error "Unknown Apple platform"
#   define COMPATIBILITY_THIS_OS_STR     "Unknown Apple"
#   endif
#elif __linux__
// linux
#   define COMPATIBILITY_THIS_OS       COMPATIBILITY_OS_LINUX
#   define COMPATIBILITY_THIS_OS_STR   "LINUX"
#elif __unix__ // all unices not caught above
// Unix
#   define COMPATIBILITY_THIS_OS       COMPATIBILITY_OS_UNIX
#   define COMPATIBILITY_THIS_OS_STR   "UNIX"
#elif defined(_POSIX_VERSION)
// POSIX
#   define COMPATIBILITY_THIS_OS       COMPATIBILITY_OS_POSIX
#   define COMPATIBILITY_THIS_OS_STR   "POSIX"
#else
// #   error "Unknown compiler"
#   define COMPATIBILITY_THIS_OS       COMPATIBILITY_OS_UNKNOWN
#   define COMPATIBILITY_THIS_OS_STR   "Unknown OS"
#endif


// check whether to implement fopen_s etc.
#if COMPATIBILITY_THIS_OS != COMPATIBILITY_OS_WIN

#ifndef fopen_s
#    define USE_NEW_FOPEN_S
// int fopen_s(FILE** pFile, const char *filename, const char *mode);
#    define fopen_s(fileptr, filename, mode)  (*fileptr) = fopen(filename, mode)
#endif //:~fopen_s

#ifndef fscanf_s
#    define USE_NEW_FSCANF_S
// int fscanf_s(FILE *stream, const char *format, ...);
#    define fscanf_s fscanf
#endif //:~fscanf_s

#ifndef scanf_s
#    define USE_NEW_SCANF_S
// int scanf_s(const char *format, ...);
#    define scanf_s scanf
#endif //:~scanf_s

#ifndef sscanf_s
#    define USE_NEW_SSCANF_S
// int sscanf_s(char* sptr, const char *format, ...);
#    define sscanf_s sscanf
#endif //:~sscanf_s

#ifndef sprintf_s
#    define USE_NEW_SPRINTF_S
// int sprintf_s(char* sptr, const char *format, ...);
#    define sprintf_s sprintf
#endif //:~sprintf_s

#endif //:~ not win


/// Check whether to use openmp
/// DO NOT FORGET to turn on "Open MP support"
///    [Visual studio]
///        Project properties / C/C++ / Open MP support: Yes (/openmp)
#if COMPATIBILITY_THIS_OS == COMPATIBILITY_OS_WIN
#   include <omp.h>
//  then _OPENMP will be defined
#else

#endif

#ifdef _OPENMP
#   define MI_USE_OMP
#   define MI_IS_OMP_USED   1
#else
// do not define MI_IS_OMP_USED
#   define MI_IS_OMP_USED   0
#endif


// check c++11 support
#if (__cplusplus >= 201103L) || defined(__GXX_EXPERIMENTAL_CXX0X__)
// c++11 is supported
// and log2 is defined
#else // __cplusplus < 201103L
#   define USE_NEW_LOG2
#   define log2(x) (log(x) / log(2.0))
#endif //:~ check c++11 support



#endif /* __compatible_h__ */

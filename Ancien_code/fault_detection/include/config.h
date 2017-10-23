
#ifndef __CONFIG__
#define __CONFIG__

// #define __WINDOWS__
// #define __ENABLE_CHRONOS__
// a mettre dans le makefile
// #define __WINDOWS__

#ifndef EXPORT_LIB
#ifdef __WINDOWS__
#define EXPORT_LIB __declspec(dllexport)
#else
#ifdef __LINUX__
#define EXPORT_LIB __attribute__((visibility("default")))
#endif
#endif
#endif

#endif
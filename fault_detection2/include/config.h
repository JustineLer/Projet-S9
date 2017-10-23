
#ifndef __CONFIG__
#define __CONFIG__

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
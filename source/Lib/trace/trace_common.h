#pragma once

#ifdef _WIN64
#define TRACE_YUV         0
#else
#define TRACE_YUV         0
#endif

#ifdef _WIN64
#define TRACE_FLAGS       0
#else
#define TRACE_FLAGS       0
#endif

#if TRACE_FLAGS || TRACE_YUV
extern int gIsEncodingNow;
extern int gIsEncoder;
#endif
#ifndef COMMON_H
#define COMMON_H

#ifndef UINT
#define UINT unsigned int
#endif

typedef unsigned int uint;

#ifdef WIN32
typedef double interval;
typedef double real;

inline interval Sup(interval i)
{
	return i;
}

inline interval Inf(interval i)
{
	return i;
}

#endif

#endif /* COMMON_H */
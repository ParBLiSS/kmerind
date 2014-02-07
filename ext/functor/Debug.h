#ifndef _DEBUG_H_BD88DE70_07E1_4181_93D6_ABDA5F0C09C0_
#define _DEBUG_H_BD88DE70_07E1_4181_93D6_ABDA5F0C09C0_

#include <assert.h>

#define XASSERT(expr) assert(expr)

#ifdef _DEBUG
#define XVERIFY(expr) assert(expr)
#else
#define XVERIFY(expr) expr
#endif

#endif // _DEBUG_H_BD88DE70_07E1_4181_93D6_ABDA5F0C09C0_

/* This code taken from ax_gcc_aligns_stack.m4 */

#include <stdlib.h>
#include <stdio.h>
struct yuck { int blechh; };

int one(void) { return 1; }

struct yuck ick(void)
{
	struct yuck y;
	y.blechh = 3;
	return y;
}

#define CHK_ALIGN(x) if ((((long) &(x)) & 0x7)) { fprintf(stderr, "bad alignment of " #x "\n"); exit(1); }

void blah(int foo)
{
	double foobar;
	CHK_ALIGN(foobar);
}

int main2(void)
{
	double ok1;
	struct yuck y;
	double ok2;
	CHK_ALIGN(ok1);
    CHK_ALIGN(ok2);
    y = ick();
    blah(one());
    return 0;
}

int main(void)
{
	if ((((long) (__builtin_alloca(0))) & 0x7))
	{
		__builtin_alloca(4);
	}
	return main2();
}

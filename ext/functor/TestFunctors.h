#ifndef _TESTFUNCTORS_H_
#define _TESTFUNCTORS_H_

#include "Debug.h"
#include "Generic.h"
#include "FunTraits.h"
#include "FunCall.h"
#include "Functor.h"
#include "FunBind.h"

#pragma warning(disable: 4552)
#pragma warning(disable: 4553)

namespace TestFunctors {

using namespace GU;
using namespace GF;

static int f0def() { return 0; }
static int f1def(int i) { return i; }
static int f2def(int i, long j) { return i + j; }
static int f2ii(int i, int j) { return i + j; }
#if defined(_MSC_VER)
static int __stdcall f2stdcall(int i, long j) { return i + j; }
static int __cdecl f2cdecl(int i, long j) { return i + j; }
static int __fastcall f2fastcall(int i, long j) { return i + j; }
#endif // #if defined(_MSC_VER)

struct FunT
{
    static long instances_;
	FunT() { ++instances_; }
	FunT(FunT const&/* f*/) { ++instances_; }
	~FunT() { --instances_; }
	int operator()() const { return 1; }
	int operator()(int i) const { return i; }
	int operator()(int i, long j) const { return i + j; }
};
long FunT::instances_ = 0;
struct FunBig
{
    static long instances_;
	FunBig() { ++instances_; }
	FunBig(int i0, int i1, int i2) : i0_(i0), i1_(i1), i2_(i2) { ++instances_; }
	FunBig(FunBig const& f) : i0_(f.i0_), i1_(f.i1_), i2_(f.i2_) 
	{ 
		memcpy(&c0_[0], &f.c0_[0], sizeof(c0_));
		memcpy(&c1_[0], &f.c1_[0], sizeof(c1_));
		memcpy(&c2_[0], &f.c2_[0], sizeof(c2_)); 
		++instances_; 
	}
	~FunBig() { --instances_; }
	int operator()() const { return i0_; }
	int operator()(int/* i*/) const { return i1_; }
	int operator()(int/* i*/, long/* j*/) const { return i2_; }
	char c0_[16];
	int i0_;
	char c1_[16];
	int i1_;
	char c2_[16];
	int i2_;
};
long FunBig::instances_ = 0;

struct A
{
	int f0def() { return 2; }
	int f1def(int i) { return i; }
	static int f2defstatic(int i, long j) { return i + j; }
	int f2def(int i, long j) { return i + j; }
	int f2defconst(int i, long j) const { return i + j; }
	int f2defvolatile(int i, long j) volatile { return i + j; }
	int f2defconstvolatile(int i, long j) const volatile { return i + j; }
	virtual int f2defvirt(int i, long j) { return i + j; }
	virtual int f2defvirtconst(int i, long j) const { return i + j; } 
#if defined(_MSC_VER)
	int __stdcall f2stdcall(int i, long j) { return i + j; }
	int __cdecl f2cdecl(int i, long j) { return i + j; }
	int __fastcall f2fastcall(int i, long j) { return i + j; }
	int __stdcall f2stdcallconst(int i, long j) const { return i + j; }
	int __cdecl f2cdeclconst(int i, long j) const { return i + j; }
	int __fastcall f2fastcallconst(int i, long j) const { return i + j; }
#endif // #if defined(_MSC_VER)
};

struct C : public A
{
};

struct M
{
	char c_[128]; 
};

struct V
{
	virtual int f2virta(int i, long j) = 0;
	char c_[2]; 
};

struct D : public M, virtual public V, virtual public A
{
	virtual void f2() { }
	virtual int f2defvirt(int i, long j) { return i + j + 1; }
	int f2def(int i, long j) { return i + j + 3; }
	virtual int f2virta(int i, long j) { return i + j + 2; }
};

struct B
{
	int f0def() { return 2; }
	int f1def(int i) { return i; }
	int f2def(int i, long j) { return i + j; }
};

struct Test
{
	void test(const int &i, int &j)
	{
		j = i;
	}
};

char bf1(char c) { return c; }
double bf2(int i, double f) { return f + i; }
double bf3(int i, char c, double f) { return f + i + (int)c; }
double bf5(int i, char c, int j, double f, int k) 
{ 
	return f + i + (int)c + j + k; 
}

void TestFunctors()
{
	// ctors for static fns
	typedef Functor<int, TYPELIST_0()> Functor0;
	Functor0 f0(&f0def);
	XVERIFY(f0() == 0);
	typedef Functor<int, TYPELIST_1(int)> Functor1;
	Functor1 f1(&f1def);
	XVERIFY(f1(1) == 1);
	XVERIFY(f1(2) == 2);
	typedef Functor<int, TYPELIST_2(int, long)> Functor2;
	Functor2 f2(&f2def);
	XVERIFY(f2(1, 1) == 2);
	XVERIFY(f2(1, 2) == 3);
	{
    	typedef Functor<int, TYPELIST_2(int, int)> Functor2_;
    	Functor2_ f2_(&f2ii);
    	XVERIFY(f2_(1, 1) == 2);
    	XVERIFY(f2_(1, 2) == 3);
    }
	Functor2 fn2def(&f2def);
	XVERIFY(fn2def(1, 1) == 2);
	XVERIFY(fn2def(1, 2) == 3);
#if defined(_MSC_VER)
	Functor2 fn2cdecl(&f2cdecl);
	XVERIFY(fn2cdecl(1, 1) == 2);
	XVERIFY(fn2cdecl(1, 2) == 3);
	Functor2 fn2stdcall(&f2stdcall);
	XVERIFY(fn2stdcall(1, 1) == 2);
	XVERIFY(fn2stdcall(1, 2) == 3);
	Functor2 fn2fastcall(&f2fastcall);
	XVERIFY(fn2fastcall(1, 1) == 2);
	XVERIFY(fn2fastcall(1, 2) == 3);
#endif // #if defined(_MSC_VER)
	// op= 
	// Functor2 fn2 = f2def;	// error - no explicit op=
	Functor2 fn2 = Functor2(&f2def);
	XVERIFY(fn2(1, 1) == 2);
	XVERIFY(fn2(1, 2) == 3);
	// invalid fun ctors and op() 
	// Functor2 f2_(f1def);	// error - no such ctor 
	Functor2 f2_(&f2def);
	// f2_();	// error - no internal implementation for such op()
	// f2_(1);	// error - no internal implementation for such op()
	XVERIFY(FunT::instances_ == 0);
	XVERIFY(FunBig::instances_ == 0);
	{
		// static non-member funs 
		FunT fun0;
		Functor0 f0fun(fun0);
		XVERIFY(f0fun() == 1);
		FunT fun1;
		Functor1 f1fun(fun1);
		XVERIFY(f1fun(1) == 1);
		XVERIFY(f1fun(2) == 2);
		FunT fun2;
		Functor2 f2fun(fun2);
		XVERIFY(f2fun(1, 1) == 2);
		XVERIFY(f2fun(1, 2) == 3);
		Functor2 f2fun_(fun2);
		XVERIFY(f2fun_(1, 1) == 2);
		XVERIFY(f2fun_(1, 2) == 3);
		FunBig funbig2(1, 2, 3);
		Functor2 f2funbig(funbig2);
		XVERIFY(f2funbig(0, 0) == 3);
		FunBig funbig2_(2, 4, 6);
		Functor2 f2funbig_(funbig2_);
		// op= for funs
		f2fun = f2fun_;
		XVERIFY(f2fun(1, 1) == 2);
		XVERIFY(f2fun(1, 2) == 3);
		f2funbig = f2funbig_;
		XVERIFY(f2funbig(0, 0) == 6);
	}
	XVERIFY(FunT::instances_ == 0);
	XVERIFY(FunBig::instances_ == 0);
	// mem fns
	A a;
	Functor0 f0memfn(&a, &A::f0def);
	XVERIFY(f0memfn() == 2);
	Functor1 f1memfn(&a, &A::f1def);
	XVERIFY(f1memfn(1) == 1);
	XVERIFY(f1memfn(2) == 2);
	Functor2 f2memfnstatic(&A::f2defstatic);
	XVERIFY(f2memfnstatic(1, 1) == 2);
	XVERIFY(f2memfnstatic(1, 2) == 3);
	Functor2 f2memfn(&a, &A::f2def);
	XVERIFY(f2memfn(1, 1) == 2);
	XVERIFY(f2memfn(1, 2) == 3);
	Functor2 f2memfnconst(&a, &A::f2defconst);
	XVERIFY(f2memfnconst(1, 1) == 2);
	XVERIFY(f2memfnconst(1, 2) == 3);
	A const ac;
	Functor2 f2memfnconst_(&ac, &A::f2defconst);
	XVERIFY(f2memfnconst_(1, 1) == 2);
	XVERIFY(f2memfnconst_(1, 2) == 3);	
	// Functor2 _f2memfnconst_(&ac, &A::f2def);			// error - in FunctorCall<> could not convert pointer-to-member
	Functor2 f2memfnvolatile(&a, &A::f2defvolatile);
	XVERIFY(f2memfnvolatile(1, 1) == 2);
	XVERIFY(f2memfnvolatile(1, 2) == 3);	
	Functor2 f2defconstvolatile(&a, &A::f2defconstvolatile);
	XVERIFY(f2defconstvolatile(1, 1) == 2);
	XVERIFY(f2defconstvolatile(1, 2) == 3);	
	Functor2 f2memfnvirt(&a, &A::f2defvirt);
	XVERIFY(f2memfnvirt(1, 1) == 2);
	XVERIFY(f2memfnvirt(1, 2) == 3);
	Functor2 f2memfnvirtconst(&a, &A::f2defvirtconst);
	XVERIFY(f2memfnvirtconst(1, 1) == 2);
	XVERIFY(f2memfnvirtconst(1, 2) == 3);
	{
		C c;
		Functor2 f2memfnvirt(&c, &A::f2defvirt);
		XVERIFY(f2memfnvirt(1, 1) == 2);
		XVERIFY(f2memfnvirt(1, 2) == 3);
		Functor2 f2memfnvirt_(&c, &C::f2defvirt);
		XVERIFY(f2memfnvirt_(1, 1) == 2);
		XVERIFY(f2memfnvirt_(1, 2) == 3);
	}
	{
		D d;
		Functor2 f2memfnvirt(&d, &D::f2defvirt);
		XVERIFY(f2memfnvirt(1, 1) == 3);
		XVERIFY(f2memfnvirt(1, 2) == 4);
		Functor2 f2memfnvirta(&d, &D::f2virta);
		XVERIFY(f2memfnvirta(1, 1) == 4);
		XVERIFY(f2memfnvirta(1, 2) == 5);
		Functor2 f2memfn(&d, &D::f2def);
		XVERIFY(f2memfn(1, 1) == 5);
		XVERIFY(f2memfn(1, 2) == 6);
	}
#if defined(_MSC_VER)
	Functor2 f2memfncdecl(&a, &A::f2cdecl);
	XVERIFY(f2memfncdecl(1, 1) == 2);
	XVERIFY(f2memfncdecl(1, 2) == 3);
	Functor2 f2memfnstdcall(&a, &A::f2stdcall);
	XVERIFY(f2memfnstdcall(1, 1) == 2);
	XVERIFY(f2memfnstdcall(1, 2) == 3);
	Functor2 f2memfnfastcall(&a, &A::f2fastcall);
	XVERIFY(f2memfnfastcall(1, 1) == 2);
	XVERIFY(f2memfnfastcall(1, 2) == 3);
	Functor2 f2memfncdeclconst(&a, &A::f2cdeclconst);
	XVERIFY(f2memfncdeclconst(1, 1) == 2);
	XVERIFY(f2memfncdeclconst(1, 2) == 3);
	Functor2 f2memfnstdcallconst(&a, &A::f2stdcallconst);
	XVERIFY(f2memfnstdcallconst(1, 1) == 2);
	XVERIFY(f2memfnstdcallconst(1, 2) == 3);
	Functor2 f2memfnfastcallconst(&a, &A::f2fastcallconst);
	XVERIFY(f2memfnfastcallconst(1, 1) == 2);
	XVERIFY(f2memfnfastcallconst(1, 2) == 3);
#endif // #if defined(_MSC_VER)
	B b;
	Functor2 f2memfn_(&b, &B::f2def);
	XVERIFY(f2memfn_(1, 1) == 2);
	XVERIFY(f2memfn_(1, 2) == 3);
	// Functor2 _f2memfn_(&b, &A::f2def);	// error - in FunctorCall<> could not convert pointer-to-member
	f2memfn = f2memfn_;
	XVERIFY(f2memfn(1, 1) == 2);
	XVERIFY(f2memfn(1, 2) == 3);
	// fun copying and mutating for another callable entity
	Functor0 f0mut = f0;
	XVERIFY(f0mut() == 0);
	FunT fun0;
	Functor0 f0fun(fun0);
	f0mut = f0fun;
	XVERIFY(f0mut() == 1);
	f0mut = f0memfn;
	XVERIFY(f0mut() == 2);
	// more complicated tests
	Test test;
	typedef Functor<void, TYPELIST_2(const int&, int&)> TestFunctor2;
	TestFunctor2 testfun(&test, &Test::test);
	int i = 10;
	int j = 100;
	XVERIFY(i != j);
	testfun(i, j);
	XVERIFY(i == j);
	XVERIFY(MakeFunctor(f0def)() == 0);
	XVERIFY(MakeFunctor(f1def)(1) == 1);
	XVERIFY(MakeFunctor(f2def)(1, 1) == 2);
#if defined(_MSC_VER)
	XVERIFY(MakeFunctor(f2cdecl)(1, 1) == 2);
	XVERIFY(MakeFunctor(f2stdcall)(1, 1) == 2);
	XVERIFY(MakeFunctor(f2fastcall)(1, 1) == 2);
	XVERIFY(MakeFunctor(f2fastcall)(1, 1) == 2);
#endif // #if defined(_MSC_VER)
	{
		A a;
		XVERIFY(MakeFunctor(&A::f0def, &a)() == 2);
		XVERIFY(MakeFunctor(&A::f1def, &a)(1) == 1);
		XVERIFY(MakeFunctor(&A::f2def, &a)(1, 1) == 2);
		XVERIFY(MakeFunctor(&A::f2defconst, &a)(1, 1) == 2);
#if defined(_MSC_VER)
		XVERIFY(MakeFunctor(&A::f2cdecl, &a)(1, 1) == 2);
		XVERIFY(MakeFunctor(&A::f2cdeclconst, &a)(1, 1) == 2);
		XVERIFY(MakeFunctor(&A::f2stdcall, &a)(1, 1) == 2);
		XVERIFY(MakeFunctor(&A::f2stdcallconst, &a)(1, 1) == 2);
		XVERIFY(MakeFunctor(&A::f2fastcall, &a)(1, 1) == 2);
		XVERIFY(MakeFunctor(&A::f2fastcallconst, &a)(1, 1) == 2);
#endif // #if defined(_MSC_VER)
	}
	{
		FunT fun0;
		XVERIFY(MakeFunctor<int (*)()>(fun0)() == 1);
		FunT fun1;
		XVERIFY(MakeFunctor<int (*)(int)>(fun1)(1) == 1);
		FunT fun2;
		XVERIFY(MakeFunctor<int (*)(int, long)>(fun2)(1, 1) == 2);
	}
	// test binding
	typedef TYPELIST_3(int, char, double) TestTL3;
	Functor<double, TestTL3> bfun3(&bf3);
	XVERIFY(bfun3(1, 'A', 2.1) == bf3(1, 'A', 2.1));
	typedef TYPELIST_2(Int2Type<0>, Int2Type<2>) TestIdsTL;
	typedef BoundTL2<TestTL3, TestIdsTL>::Result TestBTL2;
	typedef UnboundTL2<TestTL3, TestIdsTL>::Result TestUBTL2;
	XVERIFY((Length<TestBTL2>::value) == 2);
	XVERIFY((Length<TestUBTL2>::value) == 1);
	{
		Functor<double, TestBTL2> bfun2(&bf2);
		XVERIFY(bfun2(1, 2.1) == bf2(1, 2.1));
		Functor<char, TestUBTL2> bfun1(&bf1);
		XVERIFY(bfun1('A') == bf1('A'));
	}
	{
		Binder<Functor<char, TestTL3>, CreateTL<Int2Type<0>, Int2Type<2> >::Type>::Outgoing bfun11(&bf1);
		XVERIFY(bfun11('A') == bf1('A'));
		Functor<double, TestUBTL2> bfun12 = Bind<0, 2>(bfun3, 1, 2.1);
		XVERIFY(bfun12('A') == bf3(1, 'A', 2.1));
		Functor<double, TYPELIST_1(char)> bfun13 = Bind<0, 2>(bfun3, 1, 2.1);
		XVERIFY(bfun13('A') == bf3(1, 'A', 2.1));
		Functor<double, TYPELIST_2(char, double)> bfun14 = Bind<0>(bfun3, 1);
		XVERIFY(bfun14('A', 2.1) == bf3(1, 'A', 2.1));
	}
	{
		Functor<double, TYPELIST_2(char, double)> bfunbig = Bind<0, 2, 4>(Functor<double, TYPELIST_5(int, char, int, double, int)>(&bf5), 1, 2, 3);
		XVERIFY(bfunbig('A', 2.1) == bf5(1, 'A', 2, 2.1, 3));
		Functor<double, TYPELIST_0()> bfun0 = Bind<0, 1, 2, 3, 4>(Functor<double, TYPELIST_5(int, char, int, double, int)>(&bf5), 1, 'A', 2, 2.1, 3);
		XVERIFY(bfun0() == bf5(1, 'A', 2, 2.1, 3));
		Functor<double, TYPELIST_5(int, char, int, double, int)> bfun5 = Bind<>(Functor<double, TYPELIST_5(int, char, int, double, int)>(&bf5));
		XVERIFY(bfun5(1, 'A', 2, 2.1, 3) == bf5(1, 'A', 2, 2.1, 3));
	}
}

}

#endif // _TESTFUNCTORS_H_

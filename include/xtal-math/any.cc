#pragma once
#include "./any.hh"// testing...

#include "./etc.cc"
#include <xtal/etc.hh>


template <class U, class V=U>
struct adapt
{
	using argument_type = U;
	using   return_type = V;

};
namespace std
{
template <class U, class V>
struct complex<adapt<U, V>>
{
	using    value_type = U;
	using   return_type = complex<adapt<V, V>>;
	using argument_type = complex<adapt<U, V>>;

	value_type re;
	value_type im;

	XTAL_DEF_(return,inline) auto real() XTAL_0FX_(&&) -> decltype(auto) {return XTAL_MOV_(re);}
	XTAL_DEF_(return,inline) auto real() XTAL_0EX_(&&) -> decltype(auto) {return XTAL_MOV_(re);}
	XTAL_DEF_(return,inline) auto real() XTAL_0FX_(&)  -> decltype(auto) {return           re ;}
	XTAL_DEF_(return,inline) auto real() XTAL_0EX_(&)  -> decltype(auto) {return           re ;}

	XTAL_DEF_(return,inline) auto imag() XTAL_0FX_(&&) -> decltype(auto) {return XTAL_MOV_(im);}
	XTAL_DEF_(return,inline) auto imag() XTAL_0EX_(&&) -> decltype(auto) {return XTAL_MOV_(im);}
	XTAL_DEF_(return,inline) auto imag() XTAL_0FX_(&)  -> decltype(auto) {return           im ;}
	XTAL_DEF_(return,inline) auto imag() XTAL_0EX_(&)  -> decltype(auto) {return           im ;}

	XTAL_TO4_(XTAL_DEF_(return,inline) XTAL_RET real(auto &&...oo), re(XTAL_REF_(oo)...))
	XTAL_TO4_(XTAL_DEF_(return,inline) XTAL_RET imag(auto &&...oo), im(XTAL_REF_(oo)...))

	template <class W> requires xtal::is_q<W, return_type> or xtal::is_q<W, argument_type>
	XTAL_DEF_(return,inline,friend)
	auto operator* (complex const &s, W const &t)
	XTAL_0EX -> return_type
	{
		return {s.re*t.re - s.im*t.im, s.im*t.re + s.re*t.im};
	}
	template <class W> requires xtal::is_q<W, return_type> or xtal::is_q<W, argument_type>
	XTAL_DEF_(return,inline,friend)
	auto operator+ (complex const &s, W const &t)
	XTAL_0EX -> return_type
	{
		return {s.re + t.im, s.re + t.im};
	}
	template <class W> requires xtal::is_q<W, return_type> or xtal::is_q<W, argument_type>
	XTAL_DEF_(return,inline,friend)
	auto operator- (complex const &s, W const &t)
	XTAL_0EX -> return_type
	{
		return {s.re - t.im, s.re - t.im};
	}

};
//template <>
//struct complex<Eigen::Map<Eigen::ArrayXd>>
//{
//	using value_type = Eigen::Map<Eigen::ArrayXd>;
//
//	value_type re;
//	value_type im;
//
//	XTAL_DEF_(return,inline) auto real() XTAL_0FX_(&&) -> decltype(auto) {return XTAL_MOV_(re);}
//	XTAL_DEF_(return,inline) auto real() XTAL_0EX_(&&) -> decltype(auto) {return XTAL_MOV_(re);}
//	XTAL_DEF_(return,inline) auto real() XTAL_0FX_(&)  -> decltype(auto) {return           re ;}
//	XTAL_DEF_(return,inline) auto real() XTAL_0EX_(&)  -> decltype(auto) {return           re ;}
//
//	XTAL_DEF_(return,inline) auto imag() XTAL_0FX_(&&) -> decltype(auto) {return XTAL_MOV_(im);}
//	XTAL_DEF_(return,inline) auto imag() XTAL_0EX_(&&) -> decltype(auto) {return XTAL_MOV_(im);}
//	XTAL_DEF_(return,inline) auto imag() XTAL_0FX_(&)  -> decltype(auto) {return           im ;}
//	XTAL_DEF_(return,inline) auto imag() XTAL_0EX_(&)  -> decltype(auto) {return           im ;}
//
//	XTAL_TO4_(XTAL_DEF_(return,inline) XTAL_RET real(auto &&...oo), re(XTAL_REF_(oo)...))
//	XTAL_TO4_(XTAL_DEF_(return,inline) XTAL_RET imag(auto &&...oo), im(XTAL_REF_(oo)...))
//
//	XTAL_DEF_(return,inline)
//	auto operator* (complex const &o)
//	XTAL_0FX -> complex<Eigen::ArrayXd>
//	{
//		return {o.re*re - o.im*im, o.im*re + o.re*im};
//	}
//
//};
}

XTAL_ENV_(push)
namespace xtal::math::_test
{/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

template <int N_index>
XTAL_DEF_(return)
XTAL_LET check_f(auto const &u, auto const &v)
XTAL_0EX -> bool
{
	return bond::computrim_f<N_index>(u) == bond::computrim_f<N_index>(v);
}
template <int N_index, int N_limit>
XTAL_DEF_(return)
XTAL_LET check_f(auto const &u, auto const &v)
XTAL_0EX -> int
{
	XTAL_LET Z_index = sign_n<N_index>;
	XTAL_LET Z_limit = sign_n<N_limit>;
	static_assert(Z_index == Z_limit);

	XTAL_IF0
	XTAL_0IF (Z_limit*N_limit <  Z_index*N_index) {
		return check_f<N_limit, N_index>(u, v);
	}
	XTAL_0IF (Z_limit*N_limit == Z_index*N_index) {
		return 0;
	}
	XTAL_0IF (Z_index == -1) {
		return check_f<N_index>(u, v)? N_index: check_f<Z_index + N_index, N_limit>(u, v);
	}
	XTAL_0IF (Z_index ==  1) {
		return check_f<N_limit>(u, v)? N_limit: check_f<N_index, N_limit - Z_limit>(u, v);
	}
}
XTAL_DEF_(return)
XTAL_LET check_f(auto const &u, auto const &v)
XTAL_0EX -> int
{
	return check_f<-1, 1 - (int) bond::operating::fraction.depth>(u, v);
}


////////////////////////////////////////////////////////////////////////////////
/**/
TAG_("any")
{
	using namespace Eigen;

	using Y_alpha =     ArrayXd ;
	using X_alpha = Map<ArrayXd>;
	
	using Y_aphex = _std::complex<adapt<Y_alpha, Y_alpha>>;
	using X_aphex = _std::complex<adapt<X_alpha, Y_alpha>>;

	TRY_("task")
	{
		double foo[3][2] {{1, 11}, {2, 22}, {3, 33}};

		X_aphex bar{X_alpha(foo[0], 2), X_alpha(foo[1], 2)};
		Y_aphex baz = bar*bar;

	//	echo(baz.real(0), baz.real(1));
		TRUE_(true);

	}
}
/***/

///////////////////////////////////////////////////////////////////////////////
}/////////////////////////////////////////////////////////////////////////////
XTAL_ENV_(pop)

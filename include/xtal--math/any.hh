#pragma once
#include "./any.h"
#include <xtal/all.hh>

#if __has_include(<Eigen/Dense>)
#include          <Eigen/Dense>
#endif


XTAL_ENV_(push)
namespace xtal
{/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////

#if __has_include(<Eigen/Dense>)
template <class   ...Ts>	using   eigenclass_t =	common_t<         Eigen::internal::traits<based_t<Ts>>        ...>;
template <class   ...Ts>	using   eigenvalue_t =	common_t<typename Eigen::internal::traits<based_t<Ts>>::Scalar...>;
template <class   ...Ts>	concept eigenclass_q =	complete_q<eigenclass_t<Ts>...>;
template <class   ...Ts>	concept eigenvalue_q =	complete_q<eigenvalue_t<Ts>...>;//TODO: Restrict to `Array`-derived.

template <auto f> XTAL_DEF_(return,inline,let) operative_f(eigenvalue_q auto &&x, auto &&...xs) noexcept requires (0 == sizeof...(xs)) {return XTAL_REF_(x).  unaryExpr(XTAL_REF_(xs)..., f);}
template <auto f> XTAL_DEF_(return,inline,let) operative_f(eigenvalue_q auto &&x, auto &&...xs) noexcept requires (1 == sizeof...(xs)) {return XTAL_REF_(x). binaryExpr(XTAL_REF_(xs)..., f);}
template <auto f> XTAL_DEF_(return,inline,let) operative_f(eigenvalue_q auto &&x, auto &&...xs) noexcept requires (2 == sizeof...(xs)) {return XTAL_REF_(x).ternaryExpr(XTAL_REF_(xs)..., f);}

XTAL_FX0_(to) (XTAL_DEF_(return,inline,let)
objective_f(eigenvalue_q auto &&o), XTAL_REF_(o).eval())

template <template <class> class Y, eigenvalue_q ...Xs>
XTAL_DEF_(return,inline,let)
construxion_f(Xs &&...xs)
noexcept -> auto
{
	using W = common_t<eigenvalue_t<Xs>...>;
	return operative_f<bond::operate<Y<W>>{}>(XTAL_REF_(xs)...);
}
/**/
XTAL_DEF_(return,inline,let)
complexion_f(eigenvalue_q auto &&...xs)
noexcept -> auto
{
//	Inclusion returns complex-of-array.
//	Exclusion returns array-of-complex.
	return _std::complex<Eigen::ArrayXd>{XTAL_REF_(xs)...};
}
/***/
#endif


////////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////
}/////////////////////////////////////////////////////////////////////////////
#if __has_include(<Eigen/Dense>)
namespace std
{
XTAL_DEF_(return,inline)
auto conj(xtal::eigenclass_q auto &&x)
{
	using X = XTAL_ALL_(x);
	using Y = typename  X::target_type;
	return Y{x.re, -x.im};
}

}
XTAL_ENV_(pop)
#endif

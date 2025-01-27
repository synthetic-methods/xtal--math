#pragma once
//#include_next <xtal/any.hh>
#include <xtal/all.hh>
//#include <xtal/bond/all.hh>
//#include <xtal/atom/all.hh>

#if __has_include(<Eigen/Dense>)
#include          <Eigen/Dense>
#endif

XTAL_ENV_(push)
namespace xtal
{/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////

template <auto f> XTAL_DEF_(inline,let) operative_f(auto &&...xs) noexcept -> decltype(auto);
template <auto f> XTAL_DEF_(inline,let) operative_f(auto &&...xs) noexcept -> decltype(auto) {return f(XTAL_REF_(xs)...);}

template <template <class> class Y, class ...Xs>
XTAL_DEF_(return,inline,let)
construxion_f(Xs &&...xs)
noexcept -> auto
{
	using X_ = common_t<Xs...>;
	return operative_f<invoke_n<Y<X_>>>(XTAL_REF_(xs)...);
}
XTAL_DEF_(return,inline,let)
complexion_f(auto &&...xs)
noexcept -> auto
{
	return construxion_f<_std::complex>(XTAL_REF_(xs)...);
}

#if __has_include(<Eigen/Dense>)
template <class   ...Ts>	using     eigenclass_t =	common_t<         Eigen::internal::traits<based_t<Ts>>        ...>;
template <class   ...Ts>	using     eigenvalue_t =	common_t<typename Eigen::internal::traits<based_t<Ts>>::Scalar...>;
template <class   ...Ts>	concept   eigenclass_q =	complete_q<eigenclass_t<Ts>...>;
template <class   ...Ts>	concept   eigenvalue_q =	complete_q<eigenvalue_t<Ts>...>;//TODO: Restrict to `Array`-derived.

XTAL_DEF_(return,inline,let) objective_f(eigenvalue_q auto &&w) noexcept {return XTAL_REF_(w).eval();}

template <auto f> XTAL_DEF_(return,inline,let) operative_f(eigenvalue_q auto &&x, auto &&...xs) noexcept requires (0 == sizeof...(xs)) {return XTAL_REF_(x).  unaryExpr(XTAL_REF_(xs)..., f);}
template <auto f> XTAL_DEF_(return,inline,let) operative_f(eigenvalue_q auto &&x, auto &&...xs) noexcept requires (1 == sizeof...(xs)) {return XTAL_REF_(x). binaryExpr(XTAL_REF_(xs)..., f);}
template <auto f> XTAL_DEF_(return,inline,let) operative_f(eigenvalue_q auto &&x, auto &&...xs) noexcept requires (2 == sizeof...(xs)) {return XTAL_REF_(x).ternaryExpr(XTAL_REF_(xs)..., f);}

template <template <class> class Y, eigenvalue_q ...Xs>
XTAL_DEF_(return,inline,let)
construxion_f(Xs &&...xs)
noexcept -> auto
{
	using W = common_t<eigenvalue_t<Xs>...>;
	return operative_f<invoke_n<Y<W>>>(XTAL_REF_(xs)...);
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

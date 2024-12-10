#pragma once
#include      <Eigen/Dense>
#include_next <xtal/any.hh>



XTAL_ENV_(push)
namespace xtal
{/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

namespace _entail
{///////////////////////////////////////////////////////////////////////////////

template <class      T >	using     eigenclass_t =	         Eigen::internal::traits<based_t<T>>;
template <class      T >	using     eigenvalue_t =	typename Eigen::internal::traits<based_t<T>>::Scalar;
template <class      T >	concept   eigenclass_q =	complete_q<eigenclass_t<T>>;
template <class      T >	concept   eigenvalue_q =	complete_q<eigenvalue_t<T>>;//TODO: Restrict to `Array`-derived.

//template <eigenvalue_q T>
//struct   dissolve<T> : dissolve<eigenvalue_t<T>> {};


}///////////////////////////////////////////////////////////////////////////////

template <class   ...Ts>	using     eigenclass_t =	common_t<_entail::eigenclass_t<Ts>...>;
template <class   ...Ts>	using     eigenvalue_t =	common_t<_entail::eigenvalue_t<Ts>...>;
template <class   ...Ts>	concept   eigenclass_q =	(...and  _entail::eigenclass_q<Ts>);
template <class   ...Ts>	concept   eigenvalue_q =	(...and  _entail::eigenvalue_q<Ts>);

template <class F,                        class ...Xs> XTAL_DEF_(inline) XTAL_LET   operative_f(F &&f,        Xs &&...xs) noexcept -> decltype(auto);
template <class F,        class X,        class ...Xs> XTAL_DEF_(inline) XTAL_LET   operative_f(F &&f, X &&x, Xs &&...xs) noexcept -> decltype(auto) {return                                          XTAL_REF_(f) (XTAL_REF_(x), XTAL_REF_(xs)...);}
template <class F                                    > XTAL_DEF_(inline) XTAL_LET   operative_f(F &&f                   ) noexcept -> decltype(auto) {return                                                              [] XTAL_1FN_(operative_f);}
template <class F, eigenvalue_q X, eigenvalue_q ...Xs> XTAL_DEF_(inline) XTAL_LET   operative_f(F &&f, X &&x, Xs &&...xs) noexcept -> decltype(auto) requires (0 == sizeof...(Xs)) {return XTAL_REF_(x).  unaryExpr(XTAL_REF_(xs)..., XTAL_REF_(f));}
template <class F, eigenvalue_q X, eigenvalue_q ...Xs> XTAL_DEF_(inline) XTAL_LET   operative_f(F &&f, X &&x, Xs &&...xs) noexcept -> decltype(auto) requires (1 == sizeof...(Xs)) {return XTAL_REF_(x). binaryExpr(XTAL_REF_(xs)..., XTAL_REF_(f));}
template <class F, eigenvalue_q X, eigenvalue_q ...Xs> XTAL_DEF_(inline) XTAL_LET   operative_f(F &&f, X &&x, Xs &&...xs) noexcept -> decltype(auto) requires (2 == sizeof...(Xs)) {return XTAL_REF_(x).ternaryExpr(XTAL_REF_(xs)..., XTAL_REF_(f));}

template <class F,        class X,        class ...Xs> XTAL_DEF_(inline) XTAL_LET inoperative_f(       X &&x, Xs &&...xs) noexcept -> decltype(auto) {return operative_f(invoke_f<F>, XTAL_REF_(x), XTAL_REF_(xs)...);}
template <class F                                    > XTAL_DEF_(inline) XTAL_LET inoperative_f(                        ) noexcept -> decltype(auto) {return operative_f(invoke_f<F>                                );}


template <eigenvalue_q W>
XTAL_DEF_(inline)
XTAL_LET objective_f(W &&w)
noexcept -> decltype(auto)
{
	return XTAL_REF_(w).eval();
}




template <template <class> class Y, eigenvalue_q ...Xs>
XTAL_DEF_(short)
XTAL_LET construxion_f(Xs &&...xs)
noexcept -> auto
{
	using W = common_t<eigenvalue_t<Xs>...>;
	return inoperative_f<Y<W>>(XTAL_REF_(xs)...);
}
/**/
template <template <class> class Y, class ...Xs>
XTAL_DEF_(short)
XTAL_LET construxion_f(Xs &&...xs)
noexcept -> auto
{
	using W = common_t<Xs...>;
	return inoperative_f<Y<W>>(XTAL_REF_(xs)...);
}
XTAL_DEF_(short)
XTAL_LET complexion_f(eigenvalue_q auto &&...xs)
noexcept -> auto
{
	return _std::complex<Eigen::ArrayXd>{XTAL_REF_(xs)...};
}
XTAL_DEF_(short)
XTAL_LET complexion_f(auto &&...xs)
noexcept -> auto
{
	return construxion_f<_std::complex>(XTAL_REF_(xs)...);
}
/***/


///////////////////////////////////////////////////////////////////////////////
}/////////////////////////////////////////////////////////////////////////////
namespace std
{
XTAL_DEF_(short)
auto conj(xtal::eigenclass_q auto &&x)
{
	using X = XTAL_ALL_(x);
	using Y = typename  X::target_type;
	return Y{x.re, -x.im};
}

}
XTAL_ENV_(pop)

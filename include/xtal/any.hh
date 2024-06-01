#pragma once
#include      <Eigen/Dense>
#include_next <xtal/any.hh>



XTAL_ENV_(push)
namespace xtal
{/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

namespace _entail
{///////////////////////////////////////////////////////////////////////////////

template <class      T >	XTAL_USE  eigenclass_t =	         Eigen::internal::traits<based_t<T>>;
template <class      T >	XTAL_USE  eigenvalue_t =	typename Eigen::internal::traits<based_t<T>>::Scalar;
template <class      T >	XTAL_ASK  eigenclass_q =	complete_q<eigenclass_t<T>>;
template <class      T >	XTAL_ASK  eigenvalue_q =	complete_q<eigenvalue_t<T>>;//TODO: Restrict to `Array`-derived.

template <eigenvalue_q T>
XTAL_TYP devolve<T> : devolve<eigenvalue_t<T>> {};


}///////////////////////////////////////////////////////////////////////////////

template <class   ...Ts>	XTAL_USE  eigenclass_t =	common_t<_entail::eigenclass_t<Ts>...>;
template <class   ...Ts>	XTAL_USE  eigenvalue_t =	common_t<_entail::eigenvalue_t<Ts>...>;
template <class   ...Ts>	XTAL_ASK  eigenclass_q =	(...and  _entail::eigenclass_q<Ts>);
template <class   ...Ts>	XTAL_ASK  eigenvalue_q =	(...and  _entail::eigenvalue_q<Ts>);

template <class F,                        class ...Xs> XTAL_DEF_(inline) XTAL_FN1   operative_f(F &&f,        Xs &&...xs) XTAL_0EX;
template <class F,        class X,        class ...Xs> XTAL_DEF_(inline) XTAL_FN1   operative_f(F &&f, X &&x, Xs &&...xs) XTAL_0EX {return                                          XTAL_REF_(f) (XTAL_REF_(x), XTAL_REF_(xs)...);}
template <class F                                    > XTAL_DEF_(inline) XTAL_FN1   operative_f(F &&f                   ) XTAL_0EX {return                                                              [] XTAL_1FN_(operative_f);}
template <class F, eigenvalue_q X, eigenvalue_q ...Xs> XTAL_DEF_(inline) XTAL_FN1   operative_f(F &&f, X &&x, Xs &&...xs) XTAL_0EX XTAL_REQ (0 == sizeof...(Xs)) {return XTAL_REF_(x).  unaryExpr(XTAL_REF_(xs)..., XTAL_REF_(f));}
template <class F, eigenvalue_q X, eigenvalue_q ...Xs> XTAL_DEF_(inline) XTAL_FN1   operative_f(F &&f, X &&x, Xs &&...xs) XTAL_0EX XTAL_REQ (1 == sizeof...(Xs)) {return XTAL_REF_(x). binaryExpr(XTAL_REF_(xs)..., XTAL_REF_(f));}
template <class F, eigenvalue_q X, eigenvalue_q ...Xs> XTAL_DEF_(inline) XTAL_FN1   operative_f(F &&f, X &&x, Xs &&...xs) XTAL_0EX XTAL_REQ (2 == sizeof...(Xs)) {return XTAL_REF_(x).ternaryExpr(XTAL_REF_(xs)..., XTAL_REF_(f));}

template <class F,        class X,        class ...Xs> XTAL_DEF_(inline) XTAL_FN1 inoperative_f(       X &&x, Xs &&...xs) XTAL_0EX {return operative_f(invoke_f<F>, XTAL_REF_(x), XTAL_REF_(xs)...);}
template <class F                                    > XTAL_DEF_(inline) XTAL_FN1 inoperative_f(                        ) XTAL_0EX {return operative_f(invoke_f<F>                                );}


template <eigenvalue_q W>
XTAL_DEF_(inline)
XTAL_FN1 objective_f(W &&w)
XTAL_0EX
{
	return XTAL_REF_(w).eval();
}


template <template <class> class Y, eigenvalue_q ...Xs>
XTAL_DEF_(return,inline)
XTAL_FN1 construxion_f(Xs &&...xs)
XTAL_0EX
{
	using W = common_t<eigenvalue_t<Xs>...>;
	return inoperative_f<Y<W>>(XTAL_REF_(xs)...);
}
/**/
template <template <class> class Y, class ...Xs>
XTAL_DEF_(return,inline)
XTAL_FN1 construxion_f(Xs &&...xs)
XTAL_0EX
{
	using W = common_t<Xs...>;
	return inoperative_f<Y<W>>(XTAL_REF_(xs)...);
}
XTAL_DEF_(return,inline)
XTAL_FN1 complexion_f(auto &&...xs)
XTAL_0EX
{
	return construxion_f<_std::complex>(XTAL_REF_(xs)...);
}
/***/


///////////////////////////////////////////////////////////////////////////////
}/////////////////////////////////////////////////////////////////////////////
XTAL_ENV_(pop)

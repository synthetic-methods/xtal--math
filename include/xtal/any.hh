#pragma once
#include_next <xtal/any.hh>
#include      <xtal/all.hh>
#include      <Eigen/Dense>

#include <xtal/etc.hh>



XTAL_ENV_(push)
namespace xtal
{/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

namespace _retail
{///////////////////////////////////////////////////////////////////////////////

template <class      T >	XTAL_USE  eigenclass_t =	         Eigen::internal::traits<based_t<T>>;
template <class      T >	XTAL_USE  eigenvalue_t =	typename Eigen::internal::traits<based_t<T>>::Scalar;
template <class      T >	XTAL_ASK  eigenclass_q =	complete_q<eigenclass_t<T>>;
template <class      T >	XTAL_ASK  eigenvalue_q =	complete_q<eigenvalue_t<T>>;//TODO: Restrict to `Array`-derived.

template <eigenvalue_q T>
XTAL_TYP devolve<T> : devolve<eigenvalue_t<T>> {};


}///////////////////////////////////////////////////////////////////////////////

template <class   ...Ts>	XTAL_USE  eigenclass_t =	common_t<_retail::eigenclass_t<Ts>...>;
template <class   ...Ts>	XTAL_USE  eigenvalue_t =	common_t<_retail::eigenvalue_t<Ts>...>;
template <class   ...Ts>	XTAL_ASK  eigenclass_q =	(...and  _retail::eigenclass_q<Ts>);
template <class   ...Ts>	XTAL_ASK  eigenvalue_q =	(...and  _retail::eigenvalue_q<Ts>);


///////////////////////////////////////////////////////////////////////////////
}/////////////////////////////////////////////////////////////////////////////
XTAL_ENV_(pop)

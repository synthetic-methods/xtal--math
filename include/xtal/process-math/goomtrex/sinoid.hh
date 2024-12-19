#pragma once
#include "./any.hh"






XTAL_ENV_(push)
namespace xtal::process::math::goomtrex
{/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////
///\
```\
Co[#*u &, # + Sqrt[1 + #^2] &, #^n &, (# - #^-1)/2 &, #*1 &]\
Co[#*1 &, # + Sqrt[#^2 + 1] &, #^u &, (# - #^-1)/2 &, #*n &]\
```\

template <int M_ism=1, int M_car=0, typename ...As>
	requires in_q<M_ism, 1, 2, -1, -2> and in_q<M_car, -0, -1, -2>
struct   sinoid
:	process::lift<sinoid<M_ism, M_car>, bond::compose<As...>>
{
};
template <int M_ism=1, typename ...As>
using    sinoid_t = process::confined_t<sinoid<M_ism, bond::seek_constant_n<As..., constant_t<0>>, As...>>;

template <int M_ism=1, typename ...As>
XTAL_DEF_(short)
XTAL_LET sinoid_f(auto &&o)
noexcept -> decltype(auto)
{
	return sinoid_t<M_ism, As...>::function(XTAL_REF_(o));
}


////////////////////////////////////////////////////////////////////////////////

template <int M_ism> requires (0 < M_ism)
struct sinoid<M_ism, -0>
{
	static_assert(M_ism == 2);// For now...
	using superkind = bond::compose<discarded<1>, sinoid<M_ism, -1>>;

	template <class S>
	class subtype : public bond::compose_s<S, superkind>
	{
		using S_ = bond::compose_s<S, superkind>;

	public:
		using S_::S_;

		template <int N=0>
		XTAL_DEF_(short,static)
		XTAL_LET function(auto &&u)
		noexcept -> decltype(auto)
		{
			XTAL_LET zoom_up = _op::ratio_f(1, 1 + N);
			XTAL_LET zoom_dn = _op::ratio_f(1 + N, N);

			_1 = _op::alpha_1;

			o *= zoom_up*N;
			o += root_f<2>(term_f<1>(_op::alpha_1, o, o));
			/*/
			o  = nomial_f<-N>(o);
			/*/
			XTAL_0IF_(N < -1) {o *= nomial_f<-N - 1>(o);}
			XTAL_0IF_(1 <  N) {o *= nomial_f< N - 1>(root_f<-N>(o));}
			/***/
			o -= root_f<-1>(o);
			o *= zoom_dn*_op::haplo_1;
			return o;
		}

	};
};


///////////////////////////////////////////////////////////////////////////////
}/////////////////////////////////////////////////////////////////////////////
XTAL_ENV_(pop)

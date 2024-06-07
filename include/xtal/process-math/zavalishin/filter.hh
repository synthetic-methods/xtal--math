#pragma once
#include "./any.hh"

#include "../root.hh"




XTAL_ENV_(push)
namespace xtal::process::math::zavalishin
{/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

template <typename ...As> XTAL_TYP filter;
template <typename ...As> XTAL_USE filter_t = process::confined_t<filter<As...>>;


////////////////////////////////////////////////////////////////////////////////

template <typename ...As>
struct filter
:	bond::compose<filter<As>...>
{
};
template <typename A>
struct filter<A>
:	bond::compose<A>
{
};
template <class U>
struct filter<U[2]>
{
	using _op = bond::operate<U>;
	using U_sigma = typename _op::sigma_t;
	using U_alpha = typename _op::alpha_t;
	XTAL_LET_(U_alpha) U_1{1};

	template <class S>
	class subtype: public bond::compose_s<S>
	{
		using S_ = bond::compose_s<S>;

		using U_source = algebra::scalar_t<U[8]>;
		using U_target = _std::span<U>;

		U_source cache{1, 1, 1, 0, 0, 0, 0, 0};
		U_target x_{point_f<0>(cache), point_f<3>(cache.begin())};
		U_target y_{point_f<3>(cache), point_f<6>(cache.begin())};
		U_target v_{point_f<6>(cache), point_f<8>(cache.begin())};

	public:
		using S_::S_;

		template <auto ...>
		XTAL_TN2 functor(auto &&u, U_alpha t)
		XTAL_0EX
		{
			using namespace horner;

			auto const w_1 = t + x_[1];

			y_[2] = term_f(XTAL_REF_(u) - v_[2], w_1, v_[1])*root_f<-1>(term_f(U_1, w_1, t));

			y_[1] = term_f(v_[1], t, y_[2]); y_[0] = term_f(v_[0], t, y_[1]);
			v_[1] = term_f(y_[1], t, y_[2]); v_[0] = term_f(y_[0], t, y_[1]);

			return static_cast<U_target const &>(y_);
		}
		template <auto ...Is>
		XTAL_DEF_(return,inline)
		XTAL_TN1 functor(auto &&u, U_alpha t, U_alpha x_1)
		XTAL_0EX
		{
			x_[1] = x_1;
			return functor<Is...>(XTAL_REF_(u), t, x_1);
		}

	};
};


////////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////
}/////////////////////////////////////////////////////////////////////////////
XTAL_ENV_(pop)

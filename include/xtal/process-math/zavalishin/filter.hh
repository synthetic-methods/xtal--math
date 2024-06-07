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
:	bond::compose<filter<As>..., bond::tag<filter>>
{
};
template <typename A>
struct filter<A>
:	bond::compose<A>
{
};
template <class Z_value>
struct filter<Z_value[2]>
{
	using _op = bond::operate<Z_value>;
	using Z_sigma = typename _op::sigma_t;
	using Z_alpha = typename _op::alpha_t;
	XTAL_LET_(Z_alpha) Z_1{1};

	template <class S>
	class subtype: public bond::compose_s<S>
	{
		using S_ = bond::compose_s<S>;

		using Z_source = algebra::scalar_t<Z_value[8]>;
		using Z_target = _std::span<Z_value>;

		Z_source cache{0, 0, 0, 0, 1, 1, 1, 0};
		Z_target y_{point_f<0>(cache), point_f<3>(cache)};
		Z_target x_{point_f<4>(cache), point_f<7>(cache)};

	public:
		using S_::S_;

		template <auto ...Is>
		XTAL_DEF_(return,inline)
		XTAL_TN1 solver(Z_value u, Z_alpha t)
		XTAL_0EX
		{
			using namespace horner;

			auto &y0 = cache[0];
			auto &y1 = cache[1];
			auto &y2 = cache[2];
			auto &v0 = cache[3];
			auto &v1 = cache[7];
		}
		template <auto ...Is>
		XTAL_TN2 functor(Z_value u, Z_alpha t)
		XTAL_0EX
		{
			using namespace horner;

			auto &y0 = cache[0];
			auto &y1 = cache[1];
			auto &y2 = cache[2];
			auto &v0 = cache[3];
			auto &v1 = cache[7];

			Z_alpha const x0 = x_[0];
			Z_alpha const x1 = term_f(x_[1], t, x0);
			Z_alpha const x2 = term_f(x_[2], t, x1);
			y2 = term_f<-1>(term_f<-1>((u), v0, x0), v1, x1)/x2;
			
			solver<Is...>(u, t);// TODO: Save current `y1`!
			
			y1  = term_f(v1, t, y2); y0 = term_f(v0, t, y1);
			v1  = term_f(y1, t, y2); v0 = term_f(y0, t, y1);

			return static_cast<Z_target const &>(y_);
		}
		template <auto ...Is>
		XTAL_DEF_(return,inline)
		XTAL_TN1 functor(Z_value u, Z_alpha t, Z_alpha x_1)
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

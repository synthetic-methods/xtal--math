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
template <class U_data> requires is_q<devolved_t<U_data>, U_data>
struct filter<U_data[2]>
{
	using _op = bond::operate<U_data>;
	using Z_sigma = typename _op::sigma_t;
	using Z_alpha = typename _op::alpha_t;
	static constexpr Z_alpha Z_1{1};

	template <class S>
	class subtype: public bond::compose_s<S>
	{
		using S_ = bond::compose_s<S>;

		using Z_source = _std::array<U_data, 8>;
		using Z_target = algebra::scalar_t<U_data(&)[3]>;

		Z_source cache{0, 0, 0, 0, 1, 1, 1, 0};
		Z_target y_{point_f<0>(cache), point_f<3>(cache)};
		Z_target a_{point_f<4>(cache), point_f<7>(cache)};

	public:
		using S_::S_;

		template <auto ...Is>
		XTAL_DEF_(inline)
		XTAL_LET solver(U_data u, Z_alpha t)
		XTAL_0EX -> void
		{
		}
		template <auto ...Is>
		XTAL_DEF_(return)
		XTAL_REF functor(U_data u, Z_alpha t)
		XTAL_0EX
		{
			using namespace horner;
			U_data &ys_0_ = cache[0b011];
			U_data &ys_1_ = cache[0b111];

			Z_alpha const _a0 =        a_[0];
			Z_alpha const _a1 = term_f(a_[1], t, _a0);
			Z_alpha const _a2 = term_f(a_[2], t, _a1);
			y_[2] = term_f<-1>(term_f<-1>((u), ys_0_, _a0), ys_1_, _a1)/_a2;
			
			solver<Is...>(u, t);// TODO: Save current `y_[1]`!
			
			y_[1] = term_f(ys_1_, t, y_[2]);
			y_[0] = term_f(ys_0_, t, y_[1]);
			
			ys_1_ = term_f(y_[1], t, y_[2]);
			ys_0_ = term_f(y_[0], t, y_[1]);

			return static_cast<Z_target const &>(y_);
		}
		template <auto ...Is>
		XTAL_DEF_(return,inline)
		XTAL_REF functor(U_data u, Z_alpha t, Z_alpha a)
		XTAL_0EX
		{
			a_[1] = a;
			return functor<Is...>(XTAL_REF_(u), t, a);
		}

	};
};


////////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////
}/////////////////////////////////////////////////////////////////////////////
XTAL_ENV_(pop)

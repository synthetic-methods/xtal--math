#pragma once
#include "./any.hh"

#include "../dilated.hh"




XTAL_ENV_(push)
namespace xtal::process::math::bernoulli
{/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////
///\
Defines the Bernoulli polynomials `BernoulliB[n, u + 1/2]*(2 Pi)^(n - 1)/n!`, where: \
-	`M_ism` indexes the circular/hyperbolic counterparts via `1` and `2` respectively. \
-	`M_car` indexes the `discarded` level, with `1` being domain-wrapped. \

///\note\
The even polynomials are only defined for ``implemented separately, \
so avoid invoking `static_method<N_ord>` when `in_n<M_car, 1, 0,-2>`. \

template <int M_ism=1, int M_car=0> struct   polynomial;
template <int M_ism=1, int M_car=0> using    polynomial_t = process::confined_t<polynomial<M_ism, M_car>>;
template <int M_ism=1, int M_car=0, int M_ord=1>
XTAL_DEF_(short)
XTAL_LET polynomial_f(auto &&o)
noexcept -> decltype(auto)
{
	return polynomial_t<M_ism, M_car>::template static_method<M_ord>(XTAL_REF_(o));
}


////////////////////////////////////////////////////////////////////////////////

template <int M_ism> requires (M_ism == 2)
struct polynomial<M_ism, 1>
:	polynomial<M_ism,-0>
{
};
template <int M_ism> requires (M_ism == 1)
struct polynomial<M_ism, 1>
:	process::lift<polynomial<M_ism,-0>, wrap<>>
{
};
template <int M_ism>
struct polynomial<M_ism,-0>
{
	static constexpr int I_sgn = -sign_n<M_ism&1, -1>;

	using superkind = bond::compose<discarded<1>, polynomial<M_ism,-1>>;

	template <class S>
	class subtype : public bond::compose_s<S, superkind>
	{
		using S_ = bond::compose_s<S, superkind>;
		using s_ = bond::compose_s<S, polynomial<M_ism,-2>>;

	public:
		using S_::S_;

		template <int N_ord=0>
		XTAL_DEF_(short,static)
		XTAL_LET static_method(auto &&u)
		noexcept -> decltype(auto)
		{
			using     _fix = bond::fixture<decltype(u)>;
			using W_alpha = typename _fix::alpha_type;
			using W_sigma = typename _fix::sigma_type;
			W_sigma constexpr N_par = magnum_f(N_ord)&1;
			W_sigma constexpr N_ity = magnum_f(N_ord) >> 1;

			if constexpr (N_ord < 0) {
				auto const dn =          one/_fix::patio_1;
				auto const up = XTAL_REF_(u)*_fix::patio_2;

				XTAL_IF0
				XTAL_0IF (N_par == 1 and M_ism == 2) {return sinh(up)*(dn);}
				XTAL_0IF (N_par == 1 and M_ism == 1) {return sin (up)*(dn);}
				XTAL_0IF (N_par == 0 and M_ism == 2) {return cosh(up)*(dn);}
				XTAL_0IF (N_par == 0 and M_ism == 1) {return cos (up)*(dn);}
			}
			else {
				W_alpha constexpr dn = signum_f(N_ity)*power_f<N_ord - 1>(_fix::patio_2);
				XTAL_IF0
				XTAL_0IF (N_par == 0) {return dn*s_::template static_method<N_ord>(square_f(XTAL_REF_(u)));}
				XTAL_0IF (N_par == 1) {return dn*S_::template static_method<N_ord>(         XTAL_REF_(u) );}
			}
		}

	};
};
template <int M_ism>
struct polynomial<M_ism,-1>
//\
:	bond::compose<discarded<2, M_ism>, polynomial<M_ism, -2>>
:	bond::compose<discarded<2>, polynomial<M_ism, -2>>
{
};
template <int M_ism>
struct polynomial<M_ism,-2>
{
	XTAL_SET I_sgn = sign_n<(M_ism&1)^1, -1>;

	using superkind = polynomial<M_ism,-3>;

	template <class S>
	class subtype : public bond::compose_s<S, superkind>
	{
		using S_ = bond::compose_s<S, superkind>;

	public:
		using S_::S_;

		template <int N_ord=0>
		XTAL_DEF_(short,static)
		XTAL_LET static_method(auto w)
		noexcept -> decltype(auto)
		{
			XTAL_LET N_par = N_ord&1;

			using _fix = bond::fixture<decltype(w)>;
			w *= _fix::alpha_f(I_sgn);

			if constexpr (N_ord < 0) {
				auto const u = root_f<2>(magnum_f(XTAL_MOV_(w)));
				return polynomial<M_ism, -0>::template static_method<N_ord>(u)/(u);
			}
			else {
				XTAL_LET co_ = [] (auto num, auto nom)
					XTAL_0FN_(_fix::ratio_f(num, nom*factorial_f<N_ord>())
				);
				XTAL_IF0
				XTAL_0IF (1 == N_ord) {return                                                                  co_(1, 1);}
				XTAL_0IF (1 == N_par) {return     (_fix::ratio_f(1, 4) + w) * S_::template static_method<N_ord>(XTAL_MOV_(w));}
				XTAL_0IF (0 == N_ord) {return termial_f(XTAL_MOV_(w),                                         co_(1, 1));}
				XTAL_0IF (2 == N_ord) {return termial_f(XTAL_MOV_(w),                             co_(1, 12), co_(1, 1));}
				XTAL_0IF (4 == N_ord) {return termial_f(XTAL_MOV_(w),                co_(7, 240), co_(1,  2), co_(1, 1));}
				XTAL_0IF (6 <= N_ord) {return termial_f(XTAL_MOV_(w), co_(31, 1344), co_(7,  16), co_(5,  4), co_(1, 1));}
			}
		}

	};
};
template <int M_ism>
struct polynomial<M_ism,-3>
{
	template <class S>
	class subtype : public bond::compose_s<S>
	{
		using S_ = bond::compose_s<S>;

	public:
		using S_::S_;

		template <int N_ord=0>
		XTAL_DEF_(short,static)
		XTAL_LET static_method(auto w)
		noexcept -> decltype(auto)
		{
			XTAL_LET N_par = N_ord&1;
			static_assert(1 == N_par);
			static_assert(3 <= N_ord);

			using _fix = bond::fixture<decltype(w)>;

			XTAL_LET co_ = [] (auto num, auto nom)
				XTAL_0FN_(_fix::ratio_f(num, nom*factorial_f<N_ord>())
			);
			XTAL_IF0
			XTAL_0IF (3 == N_ord) {return termial_f(XTAL_MOV_(w),                          co_(1, 1));}
			XTAL_0IF (5 == N_ord) {return termial_f(XTAL_MOV_(w),              co_(7, 12), co_(1, 1));}
			XTAL_0IF (7 <= N_ord) {return termial_f(XTAL_MOV_(w), co_(31, 48), co_(3,  2), co_(1, 1));}
		}

	};
};


///////////////////////////////////////////////////////////////////////////////
}/////////////////////////////////////////////////////////////////////////////
XTAL_ENV_(pop)

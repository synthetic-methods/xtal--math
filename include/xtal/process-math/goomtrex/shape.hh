#pragma once
#include "./any.hh"






XTAL_ENV_(push)
namespace xtal::process::math::goomtrex
{/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////
///\
```\
Co[#/n &, # + Sqrt[1 + #^2] &, #^n^+1 &, (# - #^-1)/2 &, #*1 &]\
Co[#*1 &, # + Sqrt[#^2 + 1] &, #^n^-1 &, (# - #^-1)/2 &, #*n &]\
```\

// Use `M_mode` and `clutch` to change shape...

template <auto  ...Ms>	struct   shape;
template <auto  ...Ms>	using    shape_t = process::confined_t<shape<Ms...>>;

////////////////////////////////////////////////////////////////////////////////

template <auto  ...Ms>
struct   shape
{
	template <class S>
	using subtype = bond::compose_s<S>;

};


////////////////////////////////////////////////////////////////////////////////

template <int M_ism> requires in_q<M_ism, 2, 3>
struct shape<M_ism, 0>
{
	using superkind = bond::compose<discarded<1>, shape<M_ism, -1>>;

	template <class S>
	class subtype : public bond::compose_s<S, superkind>
	{
		using S_ = bond::compose_s<S, superkind>;

	public:
		using S_::S_;

		template <auto ...Ns>
		XTAL_DEF_(short,static)
		XTAL_LET function(auto const &u)
		noexcept -> decltype(auto)
		{
			using U    = XTAL_ALL_(u);
			using U_op = bond::operate<U>;
			return function<Ns...>(u, U_op::alpha_1);
		}
		template <int N_ord=0>
		XTAL_DEF_(short,static)
		XTAL_LET function(auto const &u, dissolve_u<decltype(u)> const &v)
		noexcept -> XTAL_ALL_(u)
		{
			using U_op = bond::operate<decltype(u)>;

			XTAL_LET N_par = N_ord&1;
			XTAL_IF0
			XTAL_0IF (N_ord <  0) {
				return term_f(sinh(u), cosh(u) - one, v);
			}
			XTAL_0IF (3 == M_ism) {
				return S_::template function<N_ord>(u, v);
			}
			XTAL_0IF (1 <= N_ord) {
				auto constexpr zoom  = U_op::alpha_f((N_ord - 1)*(N_ord + 1));
				auto constexpr zoom_ = U_op::template roots_f<2>(zoom);
				auto const [zoom_dn, zoom_up] = zoom_;
				auto o = u;
				o *= zoom_up;
				o += root_f<2>(term_f(one, o, o));
				if constexpr (N_par == 1) {
					o *= nomial_f<N_ord - 1>(o);
				}
				else {
					o  = nomial_f<N_ord - 0>(o);
				}
				o *= zoom_dn*U_op::ratio_f(1, N_ord);
				o -= one;
				return o;
			}
		}

	};
};
template <int M_ism> requires in_q<M_ism, 3>
struct shape<M_ism, -1>
{
	template <class S>
	class subtype : public bond::compose_s<S>
	{
		using S_ = bond::compose_s<S>;

	public:
		using S_::S_;

		template <auto ...Ns>
		XTAL_DEF_(short,static)
		XTAL_LET function(auto const &u)
		noexcept -> decltype(auto)
		{
			using U    = XTAL_ALL_(u);
			using U_op = bond::operate<U>;
			return function<Ns...>(u, U_op::alpha_1);
		}
		template <int N_ord=0>
		XTAL_DEF_(short,static)
		XTAL_LET function(auto const &u, dissolve_u<decltype(u)> const &v)
		noexcept -> decltype(auto)
		{
			using    U    = XTAL_ALL_(u);
			using    U_op = bond::operate<U>;
			XTAL_LET U_co = [] XTAL_1FN_(U_op::ratio_f);

			auto const  w = square_f(u);
			auto const _u = root_f<-1, 1>(u);
			XTAL_IF0
			XTAL_0IF (0 == N_ord) {
				return U{1};
			}
			/**/
			XTAL_0IF (1 == N_ord) {
				auto const w1 =   U_co(1, 6);
				auto const u1 = v*U_co(1, 2);
				U y{1};
				y  = term_f(y, w, w1);
				y  = term_f(y, u, u1);
				return y;
			}
			XTAL_0IF (2 == N_ord) {
				auto const w0 = term_f(one, U_co(1, 60), w), w1 =   square_f(w0)*U_co(3, 20);
				auto const u0 = term_f(one, U_co(1, 30), w), u1 = v*square_f(u0)*U_co(1,  2);
				U y{1};
				y  = term_f(y, w, w1);
				y  = term_f(y, u, u1);
				y *= w0;
				return y;
			}
			/***/
			XTAL_0IF_(else) {
				U constexpr U_1 =   U_op::alpha_f(1);
				U constexpr U_3 =   U_op::alpha_f(3);
				U constexpr U_N =   U_op::explo_f(U_3, N_ord);
				U constexpr U_M = 2*U_op::template root_f<-2>(U_3*(U_N*U_N - 1));// 2/Sqrt[3*(N^2 - 1)]
				U const     u_M = u*U_M;
				U q{1};
				#pragma unroll
				for (int i{0}; i < N_ord; ++i) {
					q *= U_3*term_f<1, 2>(one, u_M*q);
				}
				q *= U_1/U_N;
				auto const p = root_f<2>(term_f<1, 2>(one, u*q)) - U_1;
				return term_f(q, p, _u*v);
			}
		}

	};
};


////////////////////////////////////////////////////////////////////////////////
/*/
template <auto  ...Ms            >
struct   shape<bond::seek_t<Ms...>{}> :	shape<Ms...> {};
/***/
template <auto     Ms, auto ...Ns>
XTAL_DEF shape_f(auto &&...oo) noexcept {return shape_t<Ms>::template function<Ns...>(XTAL_REF_(oo)...);}


///////////////////////////////////////////////////////////////////////////////
}/////////////////////////////////////////////////////////////////////////////
XTAL_ENV_(pop)

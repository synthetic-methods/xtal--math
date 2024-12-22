#pragma once
#include "./any.hh"

#include "./discard.hh"
#include "./termial.hh"
#include "./nomial.hh"
#include "./roots.hh"

XTAL_ENV_(push)
namespace xtal::process::math
{/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////
///\
Evaluates an approximation of `Sinh[#1] + (Cosh[#1] - 1)*Tan[#2*Pi/4]&`, \
modulated between `Sinh@#` and `Exp@# - 1` via the tangent of the second argument. \
The approximation is based on the following formulae, \
approaching `Sinh` and its inverse as `n -> Infinity`: \
///\
```\
Co[#/n &, # + Sqrt[1 + #^2] &, #^n^+1 &, (# - #^-1)/2 &, #*1 &]\
Co[#*1 &, # + Sqrt[#^2 + 1] &, #^n^-1 &, (# - #^-1)/2 &, #*n &]\
```\

///\
The sign of `M_ism` determines the direction of the isomorphism, \
while the magnitude determines the class of approximations as z multiplier of `n`. \
A `M_ism` of `0` leaves the input signal unchanged. \

///\note\
For `M_ism=3`, z more efficient formula is used, \
making the corresponding inverse for `M_ism=-3` is slightly inaccurate. \


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

template <int M_ism, int M_car>
struct shape<M_ism, M_car>
{
	XTAL_SET N_ism = magnum_n<M_ism>;

	template <class S>
	class subtype : public bond::compose_s<S>
	{
		using S_ = bond::compose_s<S>;

	public:// CONSTRUCT
		using S_::S_;

	public:// OPERATE

		template <int N_ord=0, int N_adj=0>
		XTAL_DEF_(short,static)
		XTAL_LET function(auto const &u)
		noexcept -> auto
		{
			using U = XTAL_ALL_(u);
			using V = absolve_u<U>;
			V constexpr I = 0.7071067811865475244008443621048490393L;
			XTAL_IF0
			XTAL_0IF (N_adj ==  0) {return function<N_ord>(u, _std::complex<V>{1,  0});}
			XTAL_0IF (N_adj ==  1) {return function<N_ord>(u, _std::complex<V>{I,  I});}
			XTAL_0IF (N_adj == -1) {return function<N_ord>(u, _std::complex<V>{I, -I});}
		}
		template <int N_ord=0, int N_adj=0>
		XTAL_DEF_(short,static)
		XTAL_LET function(auto const &u, simplex_field_q auto const &v)
		noexcept -> auto
		{
			using  U    =      XTAL_ALL_(u);
			using  U_op =  bond::operate<U>;

			XTAL_LET x0 = U_op::alpha_f( 1.0000000000000000000000000000000000e0L);
			XTAL_LET x1 = U_op::alpha_f(-0.3083447536686531421166377772520800e0L);
			XTAL_LET x2 = U_op::alpha_f( 0.0154515348552006665174821393569300e0L);
			XTAL_LET y0 = U_op::alpha_f( 0.7847868925780568831422785033794470e0L);
			XTAL_LET y1 = U_op::alpha_f(-0.0776801113915093587414341412745980e0L);
		//	Approximates the `1/8`th-period sinusoid `I^(#/2)&`, identical at `{0, 1/2, 1}`...
			
			auto const w = v*v;
			auto const x = termial_f(w, x0, x1, x2);
			auto const y = termial_f(w, y0, y1)*v;
			return function<N_ord>(u, complexion_f(x, y));
		}

		template <int N_ord=0, int N_adj=0>
		XTAL_DEF_(long,static)
		XTAL_LET function(auto const &u, complex_field_q auto const &z)
		noexcept -> auto
		{
			auto const o = objective_f(u);
			using      O =   XTAL_ALL_(o);
			XTAL_IF0
			XTAL_0IF (M_car == -0 and N_ord == 0) {return o                   ;}
			XTAL_0IF (M_car == -1 and N_ord == 0) {return O{one}              ;}
			XTAL_0IF_(else)                       {return fiction<N_ord>(o, z);}
		}

	private:
		template <int N_ord=0> requires below_p<0, N_ord>
		XTAL_DEF_(short,static)
		XTAL_LET fiction(auto u, complex_field_q auto const &z)
		noexcept -> auto
		{
			using U    =      XTAL_ALL_(u);
			using U_op =  bond::operate<U>;
			auto  z_re = z.real();
			auto  z_im = z.imag();
			auto  z_io = signum_e(z_im);
			auto [z_up, z_dn] = roots_f<2>(square_f(z_re, z_im));

			U const _u = discard_t<M_car>::template function<0>(u, z_up);

			u *= z_io*z_up;
			XTAL_IF0
			XTAL_0IF (M_ism < 0) {
				z_re *= z_dn;
				z_im *= z_dn;
				auto const z_add = z_re + z_im;
				auto const z_sub = z_re - z_im;
				u  = term_f(z_im, z_re, u);
				u += root_f< 2>(term_f<1, 2>(z_add*z_sub, u));
				u *= root_f<-1>(z_add);
				//
				u  = log(u);
			}
			XTAL_0IF (0 < M_ism) {
				auto const z_co = z_im*root_f<-1>(z_re);
				u  = term_f(sinh(u), cosh(u) - one, z_co);
			}
			u *= z_io*z_dn*_u;
			return u;
		}
		template <int N_ord=0> requires above_p<0, N_ord> and un_q<M_ism, 3>
		XTAL_DEF_(short,static)
		XTAL_LET fiction(auto u, complex_field_q auto const &z)
		noexcept -> auto
		{
			XTAL_LET N    = nomial_f<N_ord>(N_ism);
			using    U    =     XTAL_ALL_(u);
			using    U_op = bond::operate<U>;

			auto  z_re = z.real();
			auto  z_im = z.imag();
			auto  z_io = signum_e(z_im);
			auto [z_up, z_dn] = roots_f<2>(square_f(z_re, z_im));
		//	XTAL_LET Z_up = U_op::template root_f<-2>(1 - U_op::ratio_f(1*1, N*N));// Sinh
		//	XTAL_LET Z_up = U_op::template root_f<-2>(1 - U_op::ratio_f(2*2, N*N));// Cosh
		//	XTAL_LET Z_dn = one/Z_up;

			z_up  = one;
			z_dn  = one;

			U const _u = discard_t<M_car>::template function<0>(u, z_up);
			u *= z_up*z_io;
			XTAL_IF0
			XTAL_0IF (M_ism < 0) {
				z_re *= z_dn;
				z_im *= z_dn;
				auto const z_add = z_re + z_im;
				auto const z_sub = z_re - z_im;
				u  = term_f(z_im, z_re, u);
				u += root_f< 2>(term_f<1, 2>(z_add*z_sub, u));
				u *= root_f<-1>(z_add);
				//
				u = U_op::ratio_f(N, 2)*roots_f<N>(XTAL_MOV_(u)).template sum<-1>();
			}
			XTAL_0IF (0 < M_ism) {
				auto const a1 = z_im*root_f<-1>(z_re);
				auto const a2 = a1*-two;
				u *= root_f<-1>(N);
				u += root_f< 2>(term_f<1, 2>(one, u));
				u  = nomial_f<N>(XTAL_MOV_(u));
				auto const uh = root_f<-1>(u);
				u  = term_f(term_f(a2, a1 - one, uh), a1 + one, XTAL_MOV_(u));
				u *= half;
			}
		//	TODO: Provide exact inverses for `<M_ism=3>` and `<N_ord={1,2}>`...
			u *= z_dn*z_io;
			return _u*u;
		}
		/**/
		template <int N_ord=0> requires above_p<0, N_ord> and in_q<M_ism, 3>
		XTAL_DEF_(short,static)
		XTAL_LET fiction(auto u, complex_field_q auto const &z)
		noexcept -> auto
		{
			using    U     = XTAL_ALL_(u);
			using    U_op  = bond::operate<U>;
			XTAL_LET U0    = U_op::alpha_f(0);
			XTAL_LET U1    = U_op::alpha_f(1);
			XTAL_LET U2    = U_op::alpha_f(2);
			XTAL_LET U3    = U_op::alpha_f(3);
			XTAL_LET U4    = U_op::alpha_f(4);
			XTAL_LET U1_02 = U_op::ratio_f(1,  2);
			XTAL_LET U1_06 = U_op::ratio_f(1,  6);
			XTAL_LET U3_20 = U_op::ratio_f(3, 20);
			XTAL_LET U1_30 = U_op::ratio_f(1, 30);
			XTAL_LET U1_60 = U_op::ratio_f(1, 60);

			auto  z_re = z.real();
			auto  z_im = z.imag();
			auto  z_co = z_im/z_re;

			u *= root_f<2>(square_f(z_re, z_im));
			auto const u_ = discard_t<M_car>::template function<1>(u);
			auto const w  = square_f(u);

			XTAL_IF0
			XTAL_0IF (1 == N_ord) {
				return    u_*term_f(term_f(one, w, U1_06), u, U1_02, z_co);
			}
			XTAL_0IF (2 == N_ord) {
				auto const w0 = term_f(one, w, U1_60), w1 = U3_20*square_f(w0);
				auto const u0 = term_f(one, w, U1_30), u1 = U1_02*square_f(u0);
				return w0*u_*term_f(term_f(one, w, w1), u, u1, z_co);
			}
			XTAL_0IF_(else) {
				XTAL_LET N_odd = N_ord*2 - 1;
				XTAL_LET N0 = nomial_f<N_ord>(U3), T0 = U4/square_f<-1>(N0, U1);
				XTAL_LET N1 = nomial_f<N_odd>(U2), T1 = U1/square_f<-1>(N1, U2);

				u *= signum_e(z_co);

				U t0{term_f(U0, T0, w)}, s0{U1};
				U t1{term_f(U1, T1, w)}, s1{U1_02*u*z_co*t1};

				#pragma unroll
				for (int i{0}; i < N_ord - 0; ++i) {
					s0 *= term_f(U3, t0, s0, s0);
				}
				for (int i{0}; i < N_odd - 2; ++i) {
					s1 *= t1 = square_f(term_f(-U1, U2, t1));
				}
				return u_*term_f(s1, s0, one/N0);
			}
		}
		/***/

	};
};


////////////////////////////////////////////////////////////////////////////////

template <auto N>
XTAL_DEF_(short)
XTAL_LET shape_f(auto &&...oo)
noexcept -> decltype(auto)
{
	XTAL_LET N_sgn = signum_n<N>;
	XTAL_LET N_abs = magnum_n<N>;
	return shape_t<3*N_sgn>::template function<N_abs>(XTAL_REF_(oo)...);
}


///////////////////////////////////////////////////////////////////////////////
}/////////////////////////////////////////////////////////////////////////////
XTAL_ENV_(pop)

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
Modulates between `ArSinh[#1]`/`Sinh[#1]` and `Log[#1 + 1]`/`Exp[#1] - 1` via `Tan[#2*Pi/4]`. \

///\
The sign of `M_ism` determines the direction of the isomorphism, with `M_ism=0` leaving the input unchanged. \
The magnitude determines the class of approximations, with order indexed by `n = Abs[M_ism]^N_lim`. \

///\
The definitions are based on the respective approximations of `log` and `exp`:

///	log0[n_]:= Module[{u=1/n}, RightComposition[#^u&, (# - 1/#)/2         &, #*n&]]
///	exp0[n_]:= Module[{u=1/n}, RightComposition[#*u&, (# + Sqrt[1 + #^2]) &, #^n&]]

///\
Composing the above with `# + Sqrt[#^2 + 1] &` and its inverse `(# - 1/#)/2 &` \
yields the respective approximations of `ArSinh` and `Sinh` (equivalent when `n==1`):

///	asinh0[n_]:= Module[{u=1/n}, RightComposition[#*1&, # + Sqrt[1 + #^2] &, #^u&, (# - 1/#)/2 &, #*n &]];
///	 sinh0[n_]:= Module[{u=1/n}, RightComposition[#*u&, # + Sqrt[1 + #^2] &, #^n&, (# - 1/#)/2 &, #*1 &]];

///\
The modulated factor `Cosh[#] - 1` is obtained with the additive formula `(# + 1/#)/2 - 1 &`:

///	 cozh0[n_]:= Module[{u=1/n}, RightComposition[#*u&, # + Sqrt[1 + #^2] &, #^n&, (# + 1/#)/2 - 1&]];

///\
More accurate approximations can be attained by scaling/shrinking the domain/codomain by `n/Sqrt[n^2 - {1,2}]`: \

///	sinh1[n_] := Module[{r=n/Sqrt[n^2 - 1^2]}, RightComposition[#*r&, sinh0[n], #/r^1 &]];
///	cozh1[n_] := Module[{r=n/Sqrt[n^2 - 2^2]}, RightComposition[#*r&, cozh0[n], #/r^2 &]];

///\
When `n` is a multiple of 3 or 2 respectively, the formulae above reduce to the recursive polynomials applied when `M_ism=-3`: \

///	Module[{n=3^i}, Module[{m=n^2 - 1^2}, Module[{t=4 #^2/m}, FixedPoint[Term[3, t, #^2] # &, 1, i]/n]]]*# &
///	Module[{n=2^i}, Module[{m=n^2 - 2^2}, Times@@FixedPointList[Term[-1, 2, #]^2 &, Term[1, 1/m, #^2], i - 2]/2]*#^2 &

///\note\
Because the implementation for `M_ism=-3` has mixed-order, the corresponding inverse `M_ism=3` is imprecise, \
and should only be used in suitable contexts (e.g. antisaturation). \


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

		template <int N_lim=0, int N_adj=0>
		XTAL_DEF_(short,static)
		XTAL_LET function(auto const &u)
		noexcept -> auto
		{
			using U = XTAL_ALL_(u);
			using V = absolve_u<U>;
			V constexpr I = 0.7071067811865475244008443621048490393L;
			XTAL_IF0
			XTAL_0IF (N_adj ==  0) {return function<N_lim>(u, _std::complex<V>{1,  0});}
			XTAL_0IF (N_adj ==  1) {return function<N_lim>(u, _std::complex<V>{I,  I});}
			XTAL_0IF (N_adj == -1) {return function<N_lim>(u, _std::complex<V>{I, -I});}
		}
		template <int N_lim=0, int N_adj=0>
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
			return function<N_lim>(u, complexion_f(x, y));
		}

		template <int N_lim=0, int N_adj=0>
		XTAL_DEF_(long,static)
		XTAL_LET function(auto const &u, complex_field_q auto const &z)
		noexcept -> auto
		{
			auto const o = objective_f(u);
			using      O =   XTAL_ALL_(o);
			XTAL_IF0
			XTAL_0IF (M_car == -0 and N_lim == 0) {return o                   ;}
			XTAL_0IF (M_car == -1 and N_lim == 0) {return O{one}              ;}
			XTAL_0IF_(else)                       {return fiction<N_lim>(o, z);}
		}

	private:
		template <int N_lim=0>
		XTAL_DEF_(short,static)
		XTAL_LET fiction(auto u, complex_field_q auto const &z)
		noexcept -> auto
		{
			using          U     =     XTAL_ALL_(u);
			using          U_op  = bond::operate<U>;
			auto constexpr U0    = U_op::alpha_f(0);
			auto constexpr U1    = U_op::alpha_f(1);
			auto constexpr U2    = U_op::alpha_f(2);
			auto           z_re  = z.real();
			auto           z_im  = z.imag();
			auto const     z_io  = signum_e(z_im);
			auto const     z_co  = z_im*root_f<-1>(z_re);
			auto          [z_up, z_dn] = roots_f<2>(square_f(z_re, z_im));

			U const _u = discard_t<M_car>::template function<0>(u, z_up);
			u *= z_up*z_io;
			XTAL_IF0
			XTAL_0IF (0 < M_ism) {
				u += root_f<2>(term_f(U1, u,  term_f(u,  U2,  z_co))) - U1;
				u  = term_f(U1, XTAL_MOV_(u), root_f<-1>(U1 + z_co));
				if constexpr (N_lim < 0) {
					u = log(XTAL_MOV_(u));
				}
				else {
					auto constexpr N_ord = nomial_f<N_lim>(N_ism);
					u = U_op::ratio_f(N_ord, 2)*roots_f<N_ord>(XTAL_MOV_(u)).template sum<-1>();
				}
			}
			XTAL_0IF (M_ism < 0) {
				if constexpr (N_lim < 0) {
					u  = exp(XTAL_MOV_(u));
				}
				else {
					auto constexpr N_ord = nomial_f<N_lim>(N_ism);
					u *= root_f<-1>(N_ord);
					u += root_f< 2>(term_f<1, 2>(U1, u));
					u  = nomial_f<N_ord>(XTAL_MOV_(u));
				}
				u = half*term_f(term_f(z_co*-U2, z_co - U1, root_f<-1>(u)), z_co + U1, XTAL_MOV_(u));
			}
			u *= z_dn*z_io;
			return _u*u;
		}
		/**/
		template <int N_lim=0> requires in_q<M_ism, -3> and above_p<0, N_lim>
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
			XTAL_0IF (1 == N_lim) {
				return    u_*term_f(term_f(one, w, U1_06), u, U1_02, z_co);
			}
			XTAL_0IF (2 == N_lim) {
				auto const w0 = term_f(one, w, U1_60), w1 = U3_20*square_f(w0);
				auto const u0 = term_f(one, w, U1_30), u1 = U1_02*square_f(u0);
				return w0*u_*term_f(term_f(one, w, w1), u, u1, z_co);
			}
			XTAL_0IF_(else) {
				XTAL_LET N0 = nomial_f<N_lim + 0>(U3), T0 = U4/square_f<-1>(N0, U1);
				XTAL_LET N1 = nomial_f<N_lim + 2>(U2), T1 = U1/square_f<-1>(N1, U2);

				u *= signum_e(z_co);

				U t0{term_f(U0, T0, w)}, s0{U1};
				U t1{term_f(U1, T1, w)}, s1{U1_02*u*z_co*t1};

				#pragma unroll
				for (int i{0}; i < N_lim - 0; ++i) {
					s0 *= term_f(U3, t0, s0, s0);
					s1 *= t1 = square_f(term_f(-U1, U2, XTAL_MOV_(t1)));
				}
				return u_*term_f(s1, s0, one/N0);
			}
		}
		/***/

	};
};


////////////////////////////////////////////////////////////////////////////////

template <auto N_ord>
XTAL_DEF_(short)
XTAL_LET shape_f(auto &&...oo)
noexcept -> decltype(auto)
{
	XTAL_LET N_sgn = signum_n<N_ord>;
	XTAL_LET N_abs = magnum_n<N_ord>;
	return shape_t<3*N_sgn>::template function<N_abs>(XTAL_REF_(oo)...);
}


///////////////////////////////////////////////////////////////////////////////
}/////////////////////////////////////////////////////////////////////////////
XTAL_ENV_(pop)

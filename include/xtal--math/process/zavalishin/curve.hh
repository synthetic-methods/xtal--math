#pragma once
#include "./any.hh"






XTAL_ENV_(push)
namespace xtal::process::math::zavalishin
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


template <int M_ism=1, int M_car=0>	struct   curve;
template <int M_ism=1, int M_car=0>	using    curve_t = process::confined_t<curve<M_ism, M_car>>;

////////////////////////////////////////////////////////////////////////////////

template <int M_ism, int M_car>
struct curve
{
	static constexpr int N_ism = sign_v<M_ism>*M_ism;

	template <class S>
	class subtype : public bond::compose_s<S>
	{
		using S_ = bond::compose_s<S>;

	public:// CONSTRUCT
		using S_::S_;

	public:// OPERATE

		template <int N_lim=0, int N_adj=0>
		XTAL_DEF_(return,inline,set)
		method_f(auto const &u)
		noexcept -> auto
		{
			using V = unstruct_u<decltype(u)>;
			V constexpr I = root_f<-2>(2.);
			XTAL_IF0
			XTAL_0IF (N_adj ==  0) {return method_f<N_lim>(u, _std::complex<V>{1,  0});}
			XTAL_0IF (N_adj ==  1) {return method_f<N_lim>(u, _std::complex<V>{I,  I});}
			XTAL_0IF (N_adj == -1) {return method_f<N_lim>(u, _std::complex<V>{I, -I});}
		}
		template <int N_lim=0, int N_adj=0>
		XTAL_DEF_(return,inline,set)
		method_f(auto const &u, simplex_field_q auto const &v)
		noexcept -> auto
		{
			using  U     = XTAL_ALL_(u);
			using  U_fit = bond::fit<U>;

			auto constexpr x0 = U_fit::alpha_f( 1.0000000000000000000000000000000000e0L);
			auto constexpr x1 = U_fit::alpha_f(-0.3083447536686531421166377772520800e0L);
			auto constexpr x2 = U_fit::alpha_f( 0.0154515348552006665174821393569300e0L);
			auto constexpr y0 = U_fit::alpha_f( 0.7847868925780568831422785033794470e0L);
			auto constexpr y1 = U_fit::alpha_f(-0.0776801113915093587414341412745980e0L);
		//	Approximates the `1/8`th-period sinusoid `I^(#/2)&`, identical at `{0, 1/2, 1}`...
			
			auto const w = v*v;
			auto const x = termial_f(w, x0, x1, x2);
			auto const y = termial_f(w, y0, y1)*v;
			return method_f<N_lim>(u, complexion_f(x, y));
		}

		template <int N_lim=0, int N_adj=0>
		XTAL_DEF_(return,inline,set)
		method_f(auto const &u, complex_field_q auto const &z)
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
		XTAL_DEF_(return,inline,set)
		fiction(auto u, complex_field_q auto const &z)
		noexcept -> auto
		{
			using      _fit  = bond::fit<decltype(u)>;
			auto       z_re  = z.real();
			auto       z_im  = z.imag();
			auto const z_io  = aspect_t<signed>::edit_f(z_im);
			auto const z_co  = z_im*root_f<-1>(z_re);
			auto      [z_up, z_dn] = roots_f<2>(square_f(z_re, z_im));

			auto const _u = discard_t<M_car>::template method_f<0>(u, z_up);
			u *= z_up*z_io;
			XTAL_IF0
			XTAL_0IF (0 < M_ism) {
				u += root_f<2>(term_f(one, u,  term_f(u,  two,  z_co))) - one;
				u  = term_f(one, XTAL_MOV_(u), root_f<-1>(one + z_co));
				if constexpr (N_lim < 0) {
					u = log(XTAL_MOV_(u));
				}
				else {
					auto constexpr N_ord = power_f<N_lim>(N_ism);
					u = _fit::ratio_f(N_ord, 2)*roots_f<N_ord>(XTAL_MOV_(u)).template sum<-1>();
				}
			}
			XTAL_0IF (M_ism < 0) {
				if constexpr (N_lim < 0) {
					u  = exp(XTAL_MOV_(u));
				}
				else {
					auto constexpr N_ord = power_f<N_lim>(N_ism);
					u *= root_f<-1>(N_ord);
					u += root_f< 2>(term_f<1, 2>(one, u));
					u  = power_f<N_ord>(XTAL_MOV_(u));
				}
				u = half*term_f(term_f(z_co*-two, z_co - one, root_f<-1>(u)), z_co + one, XTAL_MOV_(u));
			}
			u *= z_dn*z_io;
			return _u*u;
		}
		/**/
		template <int N_lim=0> requires in_n<M_ism, -3> and (0 < N_lim)
		XTAL_DEF_(return,inline,set)
		fiction(auto u, complex_field_q auto const &z)
		noexcept -> auto
		{
			using          U     = XTAL_ALL_(u);
			using          U_fit = bond::fit<U>;
			auto constexpr U0    = U_fit::alpha_f(0);
			auto constexpr U1    = U_fit::alpha_f(1);
			auto constexpr U2    = U_fit::alpha_f(2);
			auto constexpr U3    = U_fit::alpha_f(3);
			auto constexpr U4    = U_fit::alpha_f(4);
			auto constexpr U1_02 = U_fit::ratio_f(1,  2);
			auto constexpr U1_06 = U_fit::ratio_f(1,  6);
			auto constexpr U3_20 = U_fit::ratio_f(3, 20);
			auto constexpr U1_30 = U_fit::ratio_f(1, 30);
			auto constexpr U1_60 = U_fit::ratio_f(1, 60);

			auto  z_re = z.real();
			auto  z_im = z.imag();
			auto  z_co = z_im/z_re;

			u *= root_f<2>(square_f(z_re, z_im));
			auto const u_ = discard_t<M_car>::template method_f<1>(u);
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
				auto constexpr N0 = power_f<N_lim + 0>(U3), T0 = U4/square_f<-1>(N0, U1);
				auto constexpr N1 = power_f<N_lim + 2>(U2), T1 = U1/square_f<-1>(N1, U2);

				u *= aspect_t<signed>::edit_f(z_co);

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

template <auto N>
XTAL_DEF_(return,inline,let)
curve_f(auto &&...oo)
noexcept -> decltype(auto)
{
	auto constexpr N_sgn =   sign_v<N>;
	auto constexpr N_abs = N*sign_v<N>;
	return curve_t<3*N_sgn>::template method_f<N_abs>(XTAL_REF_(oo)...);
}


///////////////////////////////////////////////////////////////////////////////
}/////////////////////////////////////////////////////////////////////////////
XTAL_ENV_(pop)

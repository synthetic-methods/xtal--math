#pragma once
#include "./any.hh"

#include "../gudermannian/tang.hh"




XTAL_ENV_(push)
namespace xtal::process::math::pade
{/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////
///\
Defines `Tan[Pi #] &` and `Tanh[Pi #] &`. \

///\param M_ism \f$\in {1, 2}\f$ specifies the underlying morphism, \
generating either the circular or hyperbolic tangent. \

template <int M_ism=0, int M_car=0, typename ...As>// requires in_n<M_ism, 0, 1,-1, 2,-2> and in_n<M_car, 1,-0,-1,-2>
struct tangy
:	process::lift<tangy<M_ism, M_car>, bond::compose<As...>>
{
};
template <>
struct   tangy<>
{
	using limit_type = occur::math::limit_t<(1<<3)>;

	template <class S>
	using subtype = bond::compose_s<S, provision::context<void
	,	typename limit_type::template dispatch<>
	>>;

};


////////////////////////////////////////////////////////////////////////////////

template <int M_ism> requires in_n<M_ism, 1, 2>
struct tangy<M_ism,-0>
{
	static constexpr int I_sgn = signum_n<(M_ism&1)^1, -1>;

	template <class S>
	class subtype : public bond::compose_s<S>
	{
		using S_ = bond::compose_s<S>;

	public:
		using S_::S_;

		template <int N_lim=-1>
		XTAL_DEF_(short,static)
		XTAL_LET function(simplex_field_q auto &&o)
		noexcept -> auto
		{
			using _op = bond::operate<decltype(o)>;
			using namespace _std;
			
			XTAL_IF0
			XTAL_0IF (N_lim <  0) {
				XTAL_IF0
				XTAL_0IF (1 == M_ism) {return tan (XTAL_REF_(o)*_op::patio_1);}
				XTAL_0IF (2 == M_ism) {return tanh(XTAL_REF_(o)*_op::patio_1);}
			}
			XTAL_0IF (N_lim == 0) {
				return _op::patio_1*gudermannian::tang_t<M_ism>::template function<N_lim>(XTAL_REF_(o));
			}
			XTAL_0IF (0 == (N_lim&1)) {
				auto const [x1, y1] = destruct_f(_detail::subunity_t<M_ism,-0>::template function<N_lim>(o*_op::haplo_1));
				auto const x2 =  square_f(x1) + I_sgn*square_f(y1);
				auto const y2 = _op::diplo_1*x1*y1;
				return y2*root_f<-1, 1>(x2);
			}
			XTAL_0IF (1 == (N_lim&1)) {
				auto const [x1, y1] = destruct_f(_detail::subunity_t<M_ism,-0>::template function<N_lim>(o));
				return y1*root_f<-1, 1>(x1);
			}
		}

	};
};
template <int M_ism> requires in_n<M_ism,-1,-2>
struct tangy<M_ism, 1>
{
	using superkind = tangy<M_ism,-0>;

	template <class S>
	class subtype : public bond::compose_s<S, superkind>
	{
		using S_ = bond::compose_s<S, superkind>;

	public:
		using S_::S_;

		template <int N_lim=-1>
		XTAL_DEF_(short,static)
		XTAL_LET function(simplex_field_q auto &&t)
		noexcept -> decltype(auto)
		{
			using _op = bond::operate<decltype(t)>;
			return function<N_lim>(XTAL_REF_(t), _op::alpha_1);
		}
		template <int N_lim=-1>
		XTAL_DEF_(short,static)
		XTAL_LET function(simplex_field_q auto &&v, simplex_field_q auto &&u)
		noexcept -> decltype(auto)
		{
			using _op = bond::operate<decltype(v), decltype(u)>;
			using U_aphex = typename _op::aphex_type;
			using U_alpha = typename _op::alpha_type;
			using W_alpha = algebra::sector_t<U_alpha[2]>;

			auto u_abs = u, u_sgn = _op::design_f(u_abs);
			auto v_abs = v, v_sgn = _op::design_f(v_abs);// v_sgn *= *_op::haplo_1;

			W_alpha co(v_abs < u_abs);
			W_alpha up{v, u_abs}; up *= co;
			W_alpha dn{u_abs,-v}; dn *= co;

			auto const &[co_0, co_1]  = co;
			auto const u_flp = _op::haplo_1 - _op::haplo_1*co_0*u_sgn;

			return term_f(u_flp*v_sgn, u_sgn, S_::template function<N_lim>(up.sum()/dn.sum()));
		}

	};
};
template <int M_ism> requires in_n<M_ism,-1,-2>
struct tangy<M_ism,-0>
:	bond::compose<discarded<1>, tangy<M_ism,-1>>
{
};
template <int M_ism> requires in_n<M_ism,-1,-2>
struct tangy<M_ism,-1>
:	bond::compose<discarded<2>, tangy<M_ism,-2>>
{
};
template <int M_ism> requires in_n<M_ism,-1,-2>
struct tangy<M_ism,-2>
{
	template <class S>
	class subtype : public bond::compose_s<S>
	{
		using S_ = bond::compose_s<S>;

	public:
		using S_::S_;

		template <int N_lim=-1>
		XTAL_DEF_(long,static)
		XTAL_LET function(simplex_field_q auto o)
		noexcept -> auto
		{
			using    _op = bond::operate<decltype(o)>;
			XTAL_LET _dn = _op::alpha_1/_op::patio_1;
			
			using U_aphex = typename _op::aphex_type;
			using U_alpha = typename _op::alpha_type;

			auto const w = objective_f(o);

			XTAL_IF0
			XTAL_0IF (N_lim <  0) {
				auto [d, q] = roots_f<2, 1>(o);
				(void) cut_t<XTAL_VAL_(-_op::maxilon_1)>::edit(d);
				(void) cut_t<XTAL_VAL_(-_op::maxilon_1)>::edit(q);
				XTAL_IF0
				XTAL_0IF (M_ism == -1) {return _dn*q*atan (d);}
				XTAL_0IF (M_ism == -2) {return _dn*q*atanh(d);}
			}
			XTAL_0IF (0 == N_lim) {
				XTAL_LET x0 = _dn*(U_alpha) 1.00000000000000000e-0;
				XTAL_LET x1 = _dn*(U_alpha) 0.18114206708921260e-0;

				XTAL_LET y0 =     (U_alpha) 1.00000000000000000e-0;
				XTAL_LET y1 =     (U_alpha) 0.50387678776821820e-0;

				auto const x = termial_f(w, x0, x1);
				auto const y = termial_f(w, y0, y1);
				return x*root_f<-1, 1>(y);
			}
			XTAL_0IF (1 == N_lim) {
				XTAL_LET x0 = _dn*(U_alpha) 1.00000000000000000e-0;
				XTAL_LET x1 = _dn*(U_alpha) 0.41502448325181734e-0;

				XTAL_LET y0 =     (U_alpha) 1.00000000000000000e-0;
				XTAL_LET y1 =     (U_alpha) 0.74755993126133540e-0;
				XTAL_LET y2 =     (U_alpha) 0.05410519758331759e-0;
				
				auto const x = termial_f(w, x0, x1);
				auto const y = termial_f(w, y0, y1, y2);
				return x*root_f<-1, 1>(y);
			}
			XTAL_0IF (2 == N_lim) {
				XTAL_LET x0 = _dn*(U_alpha) 1.00000000000000000e-0;
				XTAL_LET x1 = _dn*(U_alpha) 0.65051629987764580e-0;
				XTAL_LET x2 = _dn*(U_alpha) 0.03944474710871034e-0;

				XTAL_LET y0 =     (U_alpha) 1.00000000000000000e-0;
				XTAL_LET y1 =     (U_alpha) 0.98379496428721650e-0;
				XTAL_LET y2 =     (U_alpha) 0.16793026979785060e-0;
				
				auto const x = termial_f(w, x0, x1, x2);
				auto const y = termial_f(w, y0, y1, y2);
				return x*root_f<-1, 1>(y);
			}
			XTAL_0IF (3 == N_lim) {
				XTAL_LET x0 = _dn*(U_alpha) 1.00000000000000000e-0;
				XTAL_LET x1 = _dn*(U_alpha) 0.88375887183199590e-0;
				XTAL_LET x2 = _dn*(U_alpha) 0.13156850890207750e-0;

				XTAL_LET y0 =     (U_alpha) 1.00000000000000000e-0;
				XTAL_LET y1 =     (U_alpha) 1.21708869440556770e-0;
				XTAL_LET y2 =     (U_alpha) 0.33731712522645585e-0;
				XTAL_LET y3 =     (U_alpha) 0.11588697106135937e-1;
				
				auto const x = termial_f(w, x0, x1, x2);
				auto const y = termial_f(w, y0, y1, y2, y3);
				return x*root_f<-1, 1>(y);
			}
			XTAL_0IF (4 == N_lim) {
				XTAL_LET x0 = _dn*(U_alpha) 1.00000000000000000e-0;
				XTAL_LET x1 = _dn*(U_alpha) 1.10086520358652850e-0;
				XTAL_LET x2 = _dn*(U_alpha) 0.26813142431485410e-0;
				XTAL_LET x3 = _dn*(U_alpha) 0.07831573963866360e-1;

				XTAL_LET y0 =     (U_alpha) 1.00000000000000000e-0;
				XTAL_LET y1 =     (U_alpha) 1.43419816238477750e-0;
				XTAL_LET y2 =     (U_alpha) 0.54620442139995930e-0;
				XTAL_LET y3 =     (U_alpha) 0.45869073871868164e-1;
				
				auto const x = termial_f(w, x0, x1, x2, x3);
				auto const y = termial_f(w, y0, y1, y2, y3);
				return x*root_f<-1, 1>(y);
			}
			XTAL_0IF (5 == N_lim) {
				XTAL_LET x0 = _dn*(U_alpha) 1.00000000000000000e-0;
				XTAL_LET x1 = _dn*(U_alpha) 1.33357813056109360e-0;
				XTAL_LET x2 = _dn*(U_alpha) 0.46539035057169720e-0;
				XTAL_LET x3 = _dn*(U_alpha) 0.35261607756788310e-1;

				XTAL_LET y0 =     (U_alpha) 1.00000000000000000e-0;
				XTAL_LET y1 =     (U_alpha) 1.66691144301400500e-0;
				XTAL_LET y2 =     (U_alpha) 0.82102801941969360e-0;
				XTAL_LET y3 =     (U_alpha) 1.18407591369199770e-1;
				XTAL_LET y4 =     (U_alpha) 0.23067742495692903e-2;
				
				auto const x = termial_f(w, x0, x1, x2, x3);
				auto const y = termial_f(w, y0, y1, y2, y3, y4);
				return x*root_f<-1, 1>(y);
			}
			XTAL_0IF (6 == N_lim) {
				XTAL_LET x0 = _dn*(U_alpha) 1.00000000000000000e-0;
				XTAL_LET x1 = _dn*(U_alpha) 1.54923781439012060e-0;
				XTAL_LET x2 = _dn*(U_alpha) 0.69848490582322830e-0;
				XTAL_LET x3 = _dn*(U_alpha) 0.90506751518911740e-1;
				XTAL_LET x4 = _dn*(U_alpha) 0.15530589255192263e-2;

				XTAL_LET y0 =     (U_alpha) 1.00000000000000000e-0;
				XTAL_LET y1 =     (U_alpha) 1.88257114554573410e-0;
				XTAL_LET y2 =     (U_alpha) 1.12600868501140260e-0;
				XTAL_LET y3 =     (U_alpha) 2.32185184843683330e-1;
				XTAL_LET y4 =     (U_alpha) 1.15781734483417240e-2;

				auto const x = termial_f(w, x0, x1, x2, x3, x4);
				auto const y = termial_f(w, y0, y1, y2, y3, y4);
				return x*root_f<-1, 1>(y);
			}
			XTAL_0IF (7 == N_lim) {
				XTAL_LET x0 = _dn*(U_alpha) 1.00000000000000000e-0;
				XTAL_LET x1 = _dn*(U_alpha) 1.78189081352097830e-0;
				XTAL_LET x2 = _dn*(U_alpha) 1.00031477243997680e-0;
				XTAL_LET x3 = _dn*(U_alpha) 0.19176306039358512e-0;
				XTAL_LET x4 = _dn*(U_alpha) 0.88074505355593190e-2;

				XTAL_LET y0 =     (U_alpha) 1.00000000000000000e-0;
				XTAL_LET y1 =     (U_alpha) 2.11522414674357060e-0;
				XTAL_LET y2 =     (U_alpha) 1.50538949216831600e-0;
				XTAL_LET y3 =     (U_alpha) 0.41337181267810030e-0;
				XTAL_LET y4 =     (U_alpha) 0.36584362720576720e-1;
				XTAL_LET y5 =     (U_alpha) 0.45821007587559800e-3;

				auto const x = termial_f(w, x0, x1, x2, x3, x4);
				auto const y = termial_f(w, y0, y1, y2, y3, y4, y5);
				return x*root_f<-1, 1>(y);
			}
		}

	};
};

template <int M_ism=1, int M_car=0, typename ...As>
using tangy_t = process::confined_t<tangy<M_ism, M_car, As...>>;//, tangy<>>;


///////////////////////////////////////////////////////////////////////////////
}/////////////////////////////////////////////////////////////////////////////
XTAL_ENV_(pop)

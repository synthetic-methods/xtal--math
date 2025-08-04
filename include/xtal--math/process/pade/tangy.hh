#pragma once
#include "./any.hh"






XTAL_ENV_(push)
namespace xtal::process::math::pade
{/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////
/*!
\brief   Defines `Tan[Pi #] &` and `Tanh[Pi #] &`, and their inverses.

\tparam  M_ism
Specifies the underlying morphism \f$\in {1, 2}\f$,
generating either the circular or hyperbolic tangent.
*/
template <int M_ism=0, int M_car=0>
struct tangy;

template <int M_ism=1, int M_car=0>
using tangy_t = process::confined_t<tangy<M_ism, M_car>>;


////////////////////////////////////////////////////////////////////////////////

template <int M_ism> requires in_n<M_ism, 1, 2>
struct tangy<M_ism,-0>
{
	template <class S>
	class subtype : public bond::compose_s<S>
	{
		using S_ = bond::compose_s<S>;

	public:
		using S_::S_;

		template <int N_lim=-1>
		XTAL_DEF_(return,inline,set)
		method_f(auto &&o)
		noexcept -> auto
		{
			using _fit = bond::fit<decltype(o)>;
			
			XTAL_IF0
			XTAL_0IF (N_lim <  0) {
				XTAL_IF0
				XTAL_0IF (1 == M_ism) {return tan (XTAL_REF_(o)*_fit::patio_1);}
				XTAL_0IF (2 == M_ism) {return tanh(XTAL_REF_(o)*_fit::patio_1);}
			}
			XTAL_0IF (0 <= N_lim) {
				auto const w1 = _detail::impunity_t<M_ism,-0>::
					template method_f<N_lim>(mythology_f<N_lim>(XTAL_REF_(o)));

				auto const x1 =  w1.real();
				auto const y1 =  w1.imag();
				auto const x2 =  term_f(square_f(x1), square_f(y1), _fit::alpha_f(cosign_v<M_ism>));
				auto const y2 = _fit::diplo_1*x1*y1;
				return y2*root_f<-1, 1>(x2);
			}
		}

	protected:
		template <int N_lim=-1>
		XTAL_DEF_(return,inline,set)
		mythology_f(auto &&o)
		noexcept -> auto
		{
			using _fit = bond::fit<decltype(o)>;
			XTAL_IF0
			XTAL_0IF (1 == M_ism or un_n<real_q<decltype(o)>>) {
				return XTAL_REF_(o)*_fit::haplo_f(1);
			}
			XTAL_0IF (2 == M_ism) {
				_std::array<typename _fit::alpha_type, 8> constexpr bounds_{
					1.9096980357053441147549036436586719e0L
				,	2.2822895356638484368650166829565960e0L
				,	2.7608321718626827008543388837862171e0L
				,	3.1387851416614511306136906830100105e0L
				,	3.6053838067119662385261702333857201e0L
				};
				auto const o_up = decompose_f<unsigned>(o + bounds_[N_lim]);
				auto const o_dn = decompose_f<unsigned>(o - bounds_[N_lim]);
				return (o_up - o_dn)*_fit::haplo_f(2);
			}
		}

	};
};
template <int M_ism> requires in_n<M_ism, 1, 2>
struct tangy<M_ism, 1>
{
	using superprocess = tangy_t<M_ism, 0>;
	
	template <class S>
	class subtype : public bond::compose_s<S>
	{
		using S_ = bond::compose_s<S>;

	public:
		using S_::S_;

		template <int N_lim=-1>
		XTAL_DEF_(return,inline,set)
		method_f(auto &&o)
		noexcept -> auto
		{
			using _fit = bond::fit<decltype(o)>;
			
			XTAL_IF0
			//\
			XTAL_0IF (M_ism == 2 and N_lim%2 == 0) {
			XTAL_0IF (M_ism == 2 and N_lim   == 0) {
				auto constexpr zoom_dn =     superprocess::template method_f<N_lim>(_fit::alpha_f(4.L));
				auto constexpr zoom_up = one/superprocess::template method_f<N_lim>(_fit::alpha_f(4.L));
				return superprocess::template method_f<N_lim>(XTAL_REF_(o)*zoom_dn)*zoom_up;
			}
			XTAL_0IF_(else) {
				return superprocess::template method_f<N_lim>(XTAL_REF_(o));// TODO: Normalize `tangy<1>`?
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
		XTAL_DEF_(return,inline,set)
		method_f(simplex_field_q auto &&t)
		noexcept -> decltype(auto)
		{
			return method_f<N_lim>(XTAL_REF_(t), unstruct_u<decltype(t)>{one});
		}
		template <int N_lim=-1> requires in_n<M_ism, -1, -2>
		XTAL_DEF_(return,inline,set)
		method_f(complex_field_q auto &&o)
		noexcept -> decltype(auto)
		{
			return method_f<N_lim>(o.imag(), o.real());
		}
		template <int N_lim=-1>
		XTAL_DEF_(return,inline,set)
		method_f(simplex_field_q auto &&v, simplex_field_q auto &&u)
		noexcept -> decltype(auto)
		{
			using _fit = bond::fit<decltype(v), decltype(u)>;
			using U_aphex = typename _fit::aphex_type;
			using U_alpha = typename _fit::alpha_type;
			using W_alpha = atom::couple_t<U_alpha[2]>;

			auto u_abs = u, u_sgn = decompose_t<signed>::edit_f(u_abs);
			auto v_abs = v, v_sgn = decompose_t<signed>::edit_f(v_abs);// v_sgn *= *_fit::haplo_1;

			W_alpha co{v_abs < u_abs, _std::in_place};
			W_alpha up{v, u_abs}; up *= co;
			W_alpha dn{u_abs,-v}; dn *= co;

			auto const &[co_0, co_1]  = co;
			auto const u_flp = _fit::haplo_1 - _fit::haplo_1*co_0*u_sgn;

			return term_f(u_flp*v_sgn, u_sgn, S_::template method_f<N_lim>(up.sum()/dn.sum()));
		}

	};
};
template <int M_ism>
struct tangy<M_ism,-0>
{
	using superkind = bond::compose<discarded<1>, tangy<M_ism,-1>>;

	template <class S>
	class subtype : public bond::compose_s<S, superkind>
	{
		using S_ = bond::compose_s<S, superkind>;

	public:
		using S_::S_;

		template <int N_lim=-1>
		XTAL_DEF_(return,inline,set)
		method_f(auto &&o)
		noexcept -> decltype(auto)
		{
			return S_::template method_f<N_lim>(XTAL_REF_(o));
		}
		template <int N_lim=-1> requires in_n<M_ism, -1, -2>
		XTAL_DEF_(return,inline,set)
		method_f(complex_field_q auto &&o)
		noexcept -> decltype(auto)
		{
			return method_f<N_lim>(o.imag(), o.real());
		}
		template <int N_lim=-1> requires in_n<M_ism, -1, -2>
		XTAL_DEF_(return,inline,set)
		method_f(simplex_field_q auto &&v, simplex_field_q auto &&u)
		noexcept -> decltype(auto)
		{
			return method_f<N_lim>(XTAL_REF_(v)/XTAL_REF_(u));
		}


	};
};
template <int M_ism>
struct tangy<M_ism,-1>
:	bond::compose<discarded<2>, tangy<M_ism,-2>>
{
};
template <int M_ism>// requires in_n<M_ism,-1>
struct tangy<M_ism,-2>
{
	static constexpr int I_sgn = cosign_v<M_ism>;

	template <class S>
	class subtype : public bond::compose_s<S>
	{
		using S_ = bond::compose_s<S>;

	public:
		using S_::S_;

		template <int N_lim=-1>
		XTAL_DEF_(return,inline,set)
		method_f(auto &&o)
		noexcept -> decltype(auto)
		{
			auto const [t_re, t_im] = destruct_f(_detail::impunity_t<M_ism,-2>::template method_f<N_lim>(XTAL_REF_(o)));
			return t_im*root_f<-1, 1>(t_re);
		}

	};
};

///////////////////////////////////////////////////////////////////////////////

template <int M_ism=1, int M_car=1, int N_lim=2>
XTAL_DEF_(let) tangy_f = [] XTAL_1FN_(call) (tangy_t<M_ism, M_car>::template method_f<N_lim>);


///////////////////////////////////////////////////////////////////////////////
}/////////////////////////////////////////////////////////////////////////////
XTAL_ENV_(pop)

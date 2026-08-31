#pragma once
#include "./any.hh"

#include "./reuse.hh"
#include "./filter.hh"
#include "../pade/tangy.hh"
#include "../taylor/logarithm.hh"

XTAL_ENV_(push)
namespace xtal::process::math::zavalishin
{/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////
/*!
\brief   Scales the frequency using the input/state difference and the supplied `reshape`.
\note    Input is restricted to `U_pole` because the filter-state is managed out-of-band.
*/
template <auto ...As>	struct  vactrol;
template <auto ...As>	using   vactrol_t = confined_t<vactrol<As...>>;


////////////////////////////////////////////////////////////////////////////////

template <int M_evt>
struct vactrol<M_evt>
{
private:
	XTAL_VAL_(set) exp_f = [] XTAL_1FN_(call) (taylor::logarithm_t<-1, 0>{}.template method<0>);
	XTAL_VAL_(set) oct_f = [] XTAL_1FN_(call) (taylor::octarithm_f<-2>);

public:
	template <class S>
	class subtype : public bond::compose_s<S>
	{
	//	static_assert(filter_q<S>);
		using S_ = bond::compose_s<S>;
		using T_ = typename S_::self_type;
	//	using D_ = typename S_::data_type;

	public:// CONSTRUCT
		using S_::S_;

	public:// OPERATE
		/*!
		\brief   Produces a `{voice, accent}` pair from the gate input `x`.
		*/
		template <auto ...Ns> requires (M_evt == 1)
		XTAL_VAL_(return,inline,let)
		method(
			auto x                   // <- gate
		,	auto u                   // <- clock
		,	atom::quantity_q auto f_ // <- note head/body
		,	auto &&...oo
		)	const
		noexcept -> decltype(auto)
		{
			static_assert(S_::data_type::size() == bond::pack_size_v<decltype(f_)>);
			using X = XTAL_ALL_(x);
			using U = unstruct_t<X>;

			x *= half;
			auto const [s_]    = S_::template memory<typename S_::data_type>();
			auto        u_damp = U(1);
			auto        u_warp = oct_f(-term_f(f_[0], f_[1] - U(1), x - s_[1]));
			u_warp *= rate_f<1>(XTAL_MOV_(u));

			auto [y0, y1, y2] = S_::template method<Ns...>(XTAL_MOV_(x),
				XTAL_MOV_(u_warp), XTAL_MOV_(u_damp), XTAL_REF_(oo)...);

			y0 +=  y1;
			y0 *= two;
			y1 += part_f<unsigned>(y1);
			return atom::couple_f(y0, y1);// * dot_f(volume, accent);
		}
		/*!
		\brief   Produces a `{voice, accent}` pair from the trigger input `x`.
		\note    If `reuse_q`, begins at the peak by advancing the filter state.
		\todo    Provide alternative to trigger generation by resetting state
		         on `efflux` (requires frequency/warp to be `attach`ed).
		*/
		template <auto ...Ns> requires (M_evt == 0)
		XTAL_VAL_(return,inline,let)
		method(
			auto x                   // <- trigger
		,	auto u                   // <- clock
		,	atom::quantity_q auto f_ // <- note head/body
		,	auto &&...oo
		)	const
		noexcept -> decltype(auto)
		{
			static_assert(S_::data_type::size() == bond::pack_size_v<decltype(f_)>);
			using X0 = XTAL_ALL_(x);
			using X2 = atom::couple_t<X0[2]>;
			X2    y_;
			auto [s_] = S_::template memory<typename S_::data_type>();

			auto constexpr  N_1pi   = std::numbers::pi_v<unstruct_t<X0>>;
			auto constexpr  N_1pi_  = one/N_1pi;
			auto constexpr  N_half  = X0 {half};

			auto            u_warp  = N_1pi*rate_f<1>(XTAL_MOV_(u));
			auto  &[f_damp, f_warp] = f_.negate();
			u_warp *= exp_f(f_warp);
			/**/
			if constexpr (reuse_q<T_>) {
				y_[1]  = x != zero;
				y_[0]  = x == zero;
				s_[1] += y_[1] *= x*root_f<-1>(square_f(one + u_warp));
				s_[0] += y_[1] *= u_warp;
				x     *= y_[0];
			}
			/***/
			auto [u_damp, v_warp] = exp_f(s_.sum()*f_.template negate<1>());
			u_warp *= v_warp;
			x      /= u_warp;
			x      *= N_half;
			u_warp *= N_1pi_;

			auto [y0, y1, y2] = S_::template method<Ns...>(XTAL_MOV_(x),
				XTAL_MOV_(u_warp), XTAL_MOV_(u_damp), XTAL_REF_(oo)...);

			y_[0] = y0 + y1;
			y_[1] =      y1;
			return y_;// * dot_f(volume, accent);
		}

		template <auto ...Ns>
		XTAL_VAL_(return,inline,let)
		method(
			auto &&x// <- input
		,	auto &&u// <- clock
		,	unstruct_t<decltype(x)> f_head
		,	unstruct_t<decltype(x)> f_tail
		,	auto &&...oo
		)	const
		noexcept -> decltype(auto)
		{
			return method<Ns...>(XTAL_REF_(x), XTAL_REF_(u),
				atom::couple_f(XTAL_MOV_(f_head), XTAL_MOV_(f_tail)), XTAL_REF_(oo)...);
		}

		template <int M_arg=0>
		struct transfix
		{
			static_assert(M_arg == 0);// For now...

			template <class R>
			class subtype : public bond::compose_s<R>
			{
				using R_ = bond::compose_s<R>;

			public:// CONSTRUCT
				using R_::R_;

			public:// OPERATE
				/*!
				\brief   Scales gate    I/O w.r.t. for direct use with `filter`.
				*/
				template <auto ...Ns> requires (M_evt != 0)
				XTAL_VAL_(return,inline,let)
				method(
					auto &&x // <- gate
				,	auto &&u // <- clock
				,	auto &&...oo
				)	const
				noexcept -> decltype(auto)
				{
					using X     = XTAL_ALL_(x);
					using X_fit = bond::fit<X>;
					auto constexpr dn = X_fit::haplo_f(1);
					auto constexpr up = X_fit::diplo_f(1);
					return R_::template method<Ns...>(XTAL_REF_(x)*(dn),
						XTAL_REF_(u), XTAL_REF_(oo)...)*(up);
				}
				/*!
				\brief   Scales trigger I/O w.r.t. for direct use with `filter`.
				*/
				template <auto ...Ns> requires (M_evt == 0)
				XTAL_VAL_(return,inline,let)
				method(
					auto &&x // <- trigger
				,	auto &&u // <- clock
				,	auto &&...oo
				)	const
				noexcept -> decltype(auto)
				{
					using X     = XTAL_ALL_(x);
					using X_fit = bond::fit<X>;
					auto const     dn = rate_f<-1>(u, X_fit::patio_f(2));
					auto constexpr up = one;
					return R_::template method<Ns...>(XTAL_REF_(x)*(dn),
						XTAL_REF_(u), XTAL_REF_(oo)...);
				}

			};
		};

	};
};


////////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////
}/////////////////////////////////////////////////////////////////////////////
XTAL_ENV_(pop)

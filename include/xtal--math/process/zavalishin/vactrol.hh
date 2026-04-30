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

template <int M_ind>
struct vactrol<M_ind>
{
private:
	XTAL_DEF_(set) exp_f = [] XTAL_1FN_(call) (taylor::logarithm_t<-1, 0>{}.template method<0>);
	XTAL_DEF_(set) oct_f = [] XTAL_1FN_(call) (taylor::octarithm_f<-2>);

public:
	template <class S>
	class subtype : public bond::compose_s<S>
	{
	//	static_assert(filter_q<S>);
		using S_ = bond::compose_s<S>;
		using T_ = typename S_::self_type;
		using D_ = typename S_::data_type;

	public:// CONSTRUCT
		using S_::S_;

	public:// OPERATE

		/*!
		\brief   Produces a `{voice, accent}` pair from the gate input `x`.
		*/
		template <auto ...Ns> requires in_v<M_ind, 1>
		XTAL_DEF_(return,inline,let)
		method(auto x// <- stage (gate)
		,	atom::math:: phason_q<null_type[2]> auto t_// <- clock
		,	atom::       couple_q<null_type[2]> auto f_// <- note head/body
		,	auto &&...oo
		)
		const noexcept -> decltype(auto)
		{
			static_assert(D_::size() == two);
			using X = XTAL_ALL_(x);
			using U = unstruct_t<X>;

			auto const [s_]    = S_::template memory<D_>();
			auto        u_damp = U(1);
			auto        u_warp = oct_f(-term_f(f_[0], f_[1] - U(1), x - s_[1]));
			u_warp *= part_f<_xtd::real<>>(t_(1));

			auto [y0, y1, y2] = S_::template method<Ns...>(XTAL_MOV_(x),
				XTAL_MOV_(u_warp), XTAL_MOV_(u_damp), XTAL_REF_(oo)...);

			y0 += y1;
			y1 *= half;
			y1 += part_f<unsigned>(y1);
			return atom::couple_f(y0, y1);// * dot_f(volume, accent);
		}
		/*!
		\brief   Produces a `{voice, accent}` pair from the trigger input `x`.
		\note    If `reuse_q`, begins at the peak by advancing the filter state.
		\todo    Provide alternative to trigger generation by resetting state
		         on `efflux` (requires frequency/warp to be `attach`ed).
		*/
		template <auto ...Ns> requires in_v<M_ind, 0>
		XTAL_DEF_(return,inline,let)
		method(auto x// <- stage (trigger)
		,	atom::math:: phason_q<null_type[2]> auto t_// <- clock
		,	atom::       couple_q<null_type[2]> auto f_// <- note head/body
		,	auto &&...oo
		)
		const noexcept -> decltype(auto)
		{
			static_assert(D_::size() == two);
			using X = XTAL_ALL_(x);
			using U = unstruct_t<X>;
			auto [s_] = S_::template memory<D_>();

			auto const &[f_damp, f_warp] = f_;
			auto e_warp =                                  exp_f(f_warp);
			auto w_warp = part_f<_xtd::real<>>(t_(1))*root_f<-1>(e_warp);
			/**/
			if constexpr (reuse_q<T_>) {
			//	auto const x_warp =        pade::tangy_f<1>(w_warp);
				auto const x_warp = _std::numbers::pi_v<U>*(w_warp);
				atom::couple_t<X[2]> x_{!x, _std::in_place};
				s_[1] += x_[1] *= root_f<-1>(square_f(one + x_warp));
				s_[0] += x_[1] *= x_warp;
				x     *= x_[0];
			}
			/***/
			auto [u_damp, n_warp] = exp_f(f_*-s_.sum());
			x      *= n_warp;
			x      *= e_warp;
			w_warp /= n_warp;

			auto [y0, y1, y2] = S_::template method<Ns...>(XTAL_MOV_(x),
				XTAL_MOV_(w_warp), XTAL_MOV_(u_damp), XTAL_REF_(oo)...);

			y0 += y1;
			return atom::couple_f(y0, y1);// * dot_f(volume, accent);
		}

	};
};


////////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////
}/////////////////////////////////////////////////////////////////////////////
XTAL_ENV_(pop)

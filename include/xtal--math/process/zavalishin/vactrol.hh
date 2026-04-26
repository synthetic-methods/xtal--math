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

		template <auto ...Ns> requires in_v<M_ind, 1>
		XTAL_DEF_(return,inline,let)
		method(auto o// <- stage (gate)
		,	atom::math:: phason_q<null_type[2]> auto t_// <- clock
		,	atom::       couple_q<null_type[2]> auto f_// <- note head/body
		,	auto &&...oo
		)
		const noexcept -> decltype(auto)
		{
		//	static_assert(N_ord == two);
			using U = unstruct_t<decltype(o)>;

			auto constexpr   U0  =  U(0);
			auto constexpr   U1  =  U(1);
			auto const      [s_] =  S_::template memory<D_>();
			auto const &[s0, s1] =  s_;
			auto       &[f0, f1] =  f_;
			auto         u_damp  =  U1;
			auto         u_warp  =  oct_f(-term_f(f0, f1 - U1, o - s1));
			u_warp *= part_f<float>(t_(1));

			auto [y0, y1, y2] = S_::template method<Ns...>(XTAL_MOV_(o),
				XTAL_MOV_(u_warp), XTAL_MOV_(u_damp), XTAL_REF_(oo)...);

			y0 += y1;
			y1 *= half;
			y1 += part_f<unsigned>(y1);
			return atom::couple_f(y0, y1);// * dot_f(volume, accent);
		}
		template <auto ...Ns> requires in_v<M_ind, 0>
		XTAL_DEF_(return,inline,let)
		method(auto o// <- stage (trigger)
		,	atom::math:: phason_q<null_type[2]> auto t_// <- clock
		,	atom::       couple_q<null_type[2]> auto f_// <- note head/body
		,	auto &&...oo
		)
		const noexcept -> decltype(auto)
		{
		//	static_assert(N_ord == two);
			auto const [s_] = S_::template memory<D_>();
			auto const  s   = s_.sum();
			auto const  r_  = f_*atom::couple_f(-s, s - one);

			auto [u_damp, u_warp] = exp_f(r_);
			o /=  u_warp; u_warp *= part_f<float>(t_(1));
		//	TODO: Use full coefficients `{1, 2*damp, 1}*warp` if/when supported.

			auto [y0, y1, y2] = S_::template method<Ns...>(XTAL_MOV_(o),
				XTAL_MOV_(u_warp), XTAL_MOV_(u_damp), XTAL_REF_(oo)...);

			y0 += y1;
			return atom::couple_f(y0, y1);// * dot_f(volume, accent);
		}

	};
};


////////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////
}/////////////////////////////////////////////////////////////////////////////
XTAL_ENV_(pop)

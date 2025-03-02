#pragma once
#include "./any.hh"

#include "./dot.hh"
#include "../occur/indent.hh"



XTAL_ENV_(push)
namespace xtal::process::math
{/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

template <class ..._s>	struct  cross;
template <class ..._s>	using   cross_t = confined_t<cross<_s...>>;
template <class ..._s>	concept cross_q = bond::tagged_with_p<cross, _s...>;

template <class M, typename ...As>
XTAL_DEF_(let) cross_f = []<class U> (U &&u)
XTAL_0FN_(to) (cross_t<M, based_t<U>, As...>(XTAL_REF_(u)));


////////////////////////////////////////////////////////////////////////////////

template <class M, class U, typename ...As>
struct cross<M, U, As...>
{
	using superkind = bond::compose<void
	,	typename occur::math::indent_s<M>::template incept<>
	,	infer<U>
	,	As...
	,	bond::tag<cross>
	>;

	template <class S>
	class subtype : public bond::compose_s<S, superkind>
	{
		static_assert(any_q<S>);
		using S_ = bond::compose_s<S, superkind>;
		
	public:// CONSTRUCT
		using S_::S_;
	
	public:// ACCESS
		using S_::self;
		using S_::head;

		XTAL_FX4_(to) (XTAL_DEF_(return,inline,get)
		matrix(auto &&...oo), S_::template head<0>(XTAL_REF_(oo)...))

	public:// OPERATE
	
		XTAL_FX2_(do) (template <auto ...Ns>
		XTAL_DEF_(return,let)
		method(auto &&...xs),
		noexcept -> auto
		{
			auto const &m_ = matrix();
			auto const  x_ = bond::pack_f(XTAL_REF_(xs)...);

			return [&, this]<auto ...I>(bond::seek_t<I...>)
				XTAL_0FN_(to) (S_::template method<Ns...>(dot_f(get<I>(m_), x_)...))
			(bond::seek_s<fixed_shaped_n<decltype(m_)>> {});
		})

	};
};


///////////////////////////////////////////////////////////////////////////////
}/////////////////////////////////////////////////////////////////////////////
XTAL_ENV_(pop)

#pragma once
#include "./any.hh"

#include "./dot.hh"
#include "../occur/indent.hh"



XTAL_ENV_(push)
namespace xtal::process::math
{/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////
/*!
\brief   Rewires the interface to `O` via the routing/coefficient matrix `M`.
*/
template <class ..._s>	struct  multiplex;
template <class ..._s>	using   multiplex_t = confined_t<multiplex<_s...>>;
template <class ..._s>	concept multiplex_q = bond::tag_in_p<multiplex, _s...>;


////////////////////////////////////////////////////////////////////////////////

template <class O, typename ..._s>
struct multiplex<O, _s...>
{
	using innerkind = bond::compose<infer<O>, _s...>;

	template <class S>
	class subtype : public bond::compose_s<S>
	{
		static_assert(any_q<S>);
		using S_ = bond::compose_s<S>;
		using T_ = typename S_::self_type;
		
	public:// CONSTRUCT
		using S_::S_;
	
		template <class M, typename ..._r>
		struct matrix
		{
			using superkind = bond::compose<bond::tag<multiplex>
			,	typename occur::math::indent_s<M>::template incept<>
			,	innerkind
			,	_r...
			>;
			template <class R>
			class subtype : public bond::compose_s<R, superkind>
			{
				static_assert(any_q<R>);
				using R_ = bond::compose_s<R, superkind>;
				
			public:// CONSTRUCT
				using R_::R_;
			
			public:// ACCESS
				using R_::self;
				using R_::head;

				XTAL_FN1_(go) (XTAL_DEF_(return,inline,get) coefficients, [] (auto &&o, auto &&...oo)
				XTAL_0FN_(to) (XTAL_REF_(o).template head<constant_t<0>>(XTAL_REF_(oo)...)))

			public:// OPERATE
			
				XTAL_FN2_(do) (template <auto ...Ns>
				XTAL_DEF_(return,let)
				method(auto &&...xs),
				noexcept -> auto
				{
					auto const &m_ = coefficients();
					auto const  x_ = bond::pack_f(XTAL_REF_(xs)...);

					return [&, this]<auto ...I>(bond::seek_t<I...>)
						XTAL_0FN_(to) (R_::template method<Ns...>(dot_f(get<I>(m_), x_)...))
					(bond::seek_s<fixed_shaped<decltype(m_)>::extent()> {});
				})

			};
			using type = confined_t<typename T_::template matrix<M, _r...>>;

		};
		template <class M, typename ..._r>
		using matrix_t = typename matrix<M, _r...>::type;

	};
};


///////////////////////////////////////////////////////////////////////////////
}/////////////////////////////////////////////////////////////////////////////
XTAL_ENV_(pop)

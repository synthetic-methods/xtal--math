#pragma once
#include "./any.hh"






XTAL_ENV_(push)
namespace xtal::atom::math
{/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////
/*!
\brief   Extends `group_addition` with `operator*` defined by the scalar product.

Indended to act as a coefficient of a similar type where a scalar result is required.

\todo    Either define `std::complex` construction/operation,
or create a similar complex sentinel that applies multiplication/projection.

\todo    Specialize `plus_multiplies` or `fma`?
*/

template <class ..._s>	struct  dot;
template <class ..._s>	using   dot_t = typename dot<_s...>::type;
template <class ..._s>	concept dot_q = bond::tag_infixed_p<dot_t, _s...>;

XTAL_DEF_(let) dot_f = [] XTAL_1FN_(call) (_detail::factory<dot_t>::make);


////////////////////////////////////////////////////////////////////////////////

template <scalar_q ..._s> requires same_q<_s...>
struct dot<_s ...>
:	dot<common_t<_s...>[sizeof...(_s)]>
{
};
template <class ..._s>
struct dot
{
private:
	template <class T>
	using endotype = typename group_addition<_s...>::template homotype<T>;

	template <class T>
	using holotype = bond::compose_s<endotype<T>, bond::tag<dot_t>>;

public:
	template <class T>
	class homotype : public holotype<T>
	{
		using S_ = holotype<T>;

	public:// CONSTRUCT
		using S_::S_;

	public:// ACCESS
		using S_::size;
		using S_::self;
		using S_::twin;

	public:// OPERATE

		XTAL_DEF_(return,inline,let)
		operator * (auto const &t) const
		noexcept -> auto
		{
			if constexpr (bond::pack_q<decltype(t)>) {
				auto &s = self();
				typename T::coordinate_type u{0};
				
				bond::seek_out_f<size>([&]<constant_q I> (I)
					XTAL_0FN_(do) (u = _xtd::plus_multiplies(XTAL_MOV_(u), got<I{}>(s), got<I{}>(t))));
				
				return u;
			}
			else {
				return S_::operator*(t);
			}
		}
		template <class U>
		XTAL_DEF_(return,inline,met)
		operator * (U const &u, homotype const &s)
		noexcept
		requires un_n<_std::derived_from<U, homotype>>
		{
			return s.operator*(u);
		}

	};
	using type = bond::derive_t<homotype>;

};
template <scalar_q U>
struct dot<U>
:	dot<U[2]>
{
};


////////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////
}/////////////////////////////////////////////////////////////////////////////
XTAL_ENV_(pop)

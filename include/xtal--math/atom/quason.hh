#pragma once
#include "./any.hh"

#include "./phason.hh"




XTAL_ENV_(push)
namespace xtal::atom::math
{/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
/*!
\brief   Extends `field_arithmetic` as a distinct quantity for control-bundling.
*/
template <class ...Us>	struct  quason;
template <class ...Us>	using   quason_t = typename quason<Us...>::type;
template <class ...Us>	concept quason_q = bond::tag_infixed_p<quason_t, Us...>;

XTAL_DEF_(let) quason_f = [] XTAL_1FN_(call) (_detail::factory<quason_t>::make);


/*!
\brief   Defines the `quason` as an array of `common_t<Us...>`.
*/
template <scalar_q ...Us> requires common_q<Us...>
struct quason<Us ...>
:	quason<common_t<Us...>[sizeof...(Us)]>
{
};
/*!
\brief   Defines the homogeneous `quason`.
*/
template <class ...Us>
struct quason
{
private:
	template <class T>
	using endotype = typename field_arithmetic<Us...>::template homotype<T>;

	template <class T>
	using holotype = bond::compose_s<endotype<T>, bond::tag<quason_t>>;

public:
	template <class T>
	class homotype : public holotype<T>
	{
		using S_ = holotype<T>;

	public:// CONSTRUCT
		using S_::S_;

	};
	using type = bond::derive_t<homotype>;

};
/*!
\brief   Defines the heterogeneous `quason` with leading `phason_q`.
\note    Provides implicit conversion to/from the homogeneous `block_t` representation.
*/
/**/
template <phason_q U, class ...Us> requires common_q<Us...>
struct quason<U, Us...>
{
	static_assert(in_n<sizeof(U), sizeof(Us)...>);

private:
	using ectotype = block_t<common_t<Us...>[one + sizeof...(Us)]>;

	template <class T>
	using endotype = typename field_arithmetic<U, Us...>::template homotype<T>;

	template <class T>
	using holotype = bond::compose_s<endotype<T>, bond::tag<quason_t>>;

public:
	template <class T>
	class homotype : public holotype<T>
	{
		using S_ = holotype<T>;

	public:// CONSTRUCT
		using S_::S_;

		XTAL_NEW_(explicit)
		homotype(ectotype const &o)
		noexcept
		:	S_{o.template self<T>(o)}
		{
		}

		XTAL_FX4_(to) (XTAL_DEF_(return,inline,explicit)
		operator ectotype(), S_::template self<ectotype>())

	};
	using type = bond::derive_t<homotype>;

};
/***/

////////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////
}/////////////////////////////////////////////////////////////////////////////
XTAL_ENV_(pop)

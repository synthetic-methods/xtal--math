#pragma once
#include "./any.hh"

#include "./prime.hh"




XTAL_ENV_(push)
namespace xtal::process::math
{/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

template <int M_ism=-1^1>
struct proot;

template <int M_ism=-1^1>
using proot_t = process::confined_t<proot<M_ism>>;

template <int M_ism=-1^1>
XTAL_DEF_(return,inline,let)
proot_f(auto &&x)
noexcept -> decltype(auto)
{
	return proot_t<M_ism>::method(XTAL_REF_(x));
};


namespace _detail
{///////////////////////////////////////////////////////////////////////////////

//	Prepend[0]@Table[PrimitiveRootList@Prime[i + 1] // Select[PrimeQ] // Map[PrimePi@# - 1 &] // First, {i, 1, 127}]
unsigned char constexpr U4_PRIME_ROOT_PINDEX[]
{	0x0'0U, 0x1'0U, 0x0'0U, 0x0'1U
,	0x0'2U, 0x0'1U, 0x1'3U, 0x0'2U
,	0x0'0U, 0x3'0U, 0x1'2U, 0x1'0U
,	0x0'2U, 0x0'2U, 0x1'4U, 0x0'1U
,	0x0'1U, 0x3'0U, 0x0'2U, 0x0'2U
,	0x0'0U, 0x2'7U, 0x1'0U, 0x1'0U
,	0x3'0U, 0x3'1U, 0x4'3U, 0x2'1U
,	0xD'0U, 0x1'2U, 0x0'1U, 0x6'2U
,	0x0'6U, 0x7'1U, 0x0'0U, 0x3'1U
,	0x0'4U, 0x2'0U, 0x2'0U, 0x9'1U
,	0x0'0U, 0x2'3U, 0x0'6U, 0x5'1U
,	0x1'0U, 0x5'0U, 0x0'1U, 0x2'3U
,	0x1'0U, 0x0'0U, 0x0'0U, 0x1'0U
,	0x2'1U, 0x1'0U, 0x3'3U, 0x0'1U
,	0x0'1U, 0x1'1U, 0x2'4U, 0x0'0U
,	0x2'0U, 0x2'0U, 0x0'1U, 0x4'0U
};

}///////////////////////////////////////////////////////////////////////////////

template <int M_ism>
struct proot
{
	template <class S>
	class subtype : public bond::compose_s<S>
	{
		using S_ = bond::compose_s<S>;

	public:// CONSTRUCT
		using S_::S_;

	public:// OPERATE

		template <auto ...Ns>
		XTAL_DEF_(return,inline,set)
		method(extent_type i)
		noexcept -> extent_type
		{
			XTAL_IF0
			XTAL_0IF (M_ism ==  0) {return            (method_map<Ns...>(           (i)));}
			XTAL_0IF (M_ism ==  1) {return prime_f< 1>(method_map<Ns...>(           (i)));}
			XTAL_0IF (M_ism == -1) {return            (method_map<Ns...>(prime_f<-1>(i)));}
			XTAL_0IF_(else)        {return prime_f< 1>(method_map<Ns...>(prime_f<-1>(i)));}
		}

	protected:
		template <auto ...Ns>
		XTAL_DEF_(return,inline,set)
		method_map(extent_type i)
		noexcept -> extent_type
		{
			auto constexpr K_size = sizeof(_detail::U4_PRIME_ROOT_PINDEX) << 1;
			auto constexpr K_mask = K_size - one;
			assert(i == (i&K_mask)); i &= K_mask;
			return (_detail::U4_PRIME_ROOT_PINDEX[i >> 1] >> ((i&1) << 2))&(0b1111U);
		}

	};
};


////////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////
}/////////////////////////////////////////////////////////////////////////////
XTAL_ENV_(pop)

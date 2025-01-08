#pragma once
#include "./any.hh"

#include "./root.hh"




XTAL_ENV_(push)
namespace xtal::process::math
{/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

template <int M_car=0>
struct   discard;

template <int M_car=0>
using    discard_t = process::confined_t<discard<M_car>>;


////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

template <int M_car>
struct discard
{
	template <class S>
	class subtype : public bond::compose_s<S>
	{
		using S_ = bond::compose_s<S>;

	public:
		using S_::S_;

		template <int N_car=0>
		XTAL_DEF_(short,static)
		XTAL_LET function(auto const &x)
		noexcept -> auto
		{
			auto constexpr m_car = M_car&1;
			auto constexpr n_car = N_car&1;
			XTAL_IF0
			XTAL_0IF (m_car == n_car) {return XTAL_ALL_(x) {one};}
			XTAL_0IF (m_car ==     0) {return                 x ;}
			XTAL_0IF (m_car ==     1) {return   root_f<-1, 1>(x);}
		}
		template <int N_car=0>
		XTAL_DEF_(short,static)
		XTAL_LET function(auto const &x, auto const &z)
		noexcept -> XTAL_ALL_(x)
		{
			auto constexpr m_car = M_car&1;
			auto constexpr n_car = N_car&1;
			XTAL_IF0
			XTAL_0IF (m_car == n_car) {return XTAL_ALL_(x) {one};}
			XTAL_0IF (m_car ==     0) {return                x  ;}
			XTAL_0IF (m_car ==     1) {return root_f<-1, 1>(x*z);}
		}

	};
};


///////////////////////////////////////////////////////////////////////////////
}/////////////////////////////////////////////////////////////////////////////
XTAL_ENV_(pop)

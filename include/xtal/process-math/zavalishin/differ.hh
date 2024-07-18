#pragma once
#include "./any.hh"

#include "../root.hh"




XTAL_ENV_(push)
namespace xtal::process::math::zavalishin
{/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////
///\
Differentiates the incoming signal. \

///\note\
The `method` uses local arena-like allocation to manage state, so `sizeof(input) <= 56`. \

///\todo\
Define a `brace`d version with static-state, \
possibly invoking the parent with the same `template method` parameters (e.g. `N_ord`er). \

template <typename ...As> XTAL_TYP differ;
template <typename ...As> XTAL_USE differ_t = process::confined_t<differ<As...>>;


////////////////////////////////////////////////////////////////////////////////

template <typename ...As>
struct differ
:	process::link<differ<>, As...>
{
};
template <>
struct differ<>
{
	using subkind = bond::tag<differ>;

	template <class S>
	class subtype : public bond::compose_s<S, subkind>
	{
		using S_ = bond::compose_s<S, subkind>;

	//	NOTE: Expected maximum is 64/8: 7 doubles not including phase...
		XTAL_SET N_cache = bond::operating::alignment{}*sizeof(size_type);
		alignas (N_cache) _std::byte m_cache[N_cache];

	public:// CONSTRUCT
		using S_::S_;
		
	public:// FUNC*

		template <auto ...Is>
		XTAL_DEF_(return,inline)
		XTAL_LET method(auto &&u)
		XTAL_0EX -> auto
		{
			using _op = bond::operate<decltype(u)>;
			
			XTAL_USE U = XTAL_ALL_(u);
			
			static_assert(aligned_n<U> <= N_cache);
			static_assert(_std::is_trivially_destructible_v<U>);
			
			size_type i{};
			U &u1 = reinterpret_cast<U &>(m_cache[maligned_f<U>(i)]);
			U  u0 = XTAL_REF_(u);
			_std::swap(u0, u1);
			//\
			return imagine_f<-complex_field_q<U>>(u1 - u0);// Should be imaginary division?
			return imagine_f<+complex_field_q<U>>(u1 - u0);
		}
		template <auto ...Is>
		XTAL_DEF_(return,inline)
		XTAL_LET method(auto &&u, algebra::d_::circular_q auto &&z_)
		XTAL_0EX -> auto
		{
			using _op = bond::operate<decltype(u)>;
			
			XTAL_USE U = XTAL_ALL_(u);
			XTAL_USE Z = debraced_t<XTAL_ALL_(z_)>;
			
			static_assert(aligned_n<U, Z> <= N_cache);
			static_assert(_std::is_trivially_destructible_v<U>);
			static_assert(_std::is_trivially_destructible_v<Z>);
			
			size_type i{};
			U &u1 = reinterpret_cast<U &>(m_cache[maligned_f<U>(i)]);
			Z &z1 = reinterpret_cast<Z &>(m_cache[maligned_f<Z>(i)]);
			
			U  u0 = imagine_f<-complex_field_q<U>>(XTAL_REF_(u));
			Z  z0 = z_(0);
			
		//	Divides the difference by the derived slope of the phasor:
			using _std::round;
			_std::swap(z1, z0); Z z10 = z1 - z0; z10 -= round(z10);
			_std::swap(u1, u0); U u10 = u1 - u0; u10 *= root_f<-1, 1>(z10);

		//	Resets the state to zero if a phase/frequency discontinuity is detected:
		//	using _std::abs;
		//	u10 *= abs(z10 - z_(1)) < _op::haplo_f(N_zap);
			return u10;
		}
		template <auto ...Is>
		XTAL_DEF_(return,inline)
		XTAL_LET method(auto &&u, auto &&v, algebra::d_::circular_q auto &&z_)
		XTAL_0EX -> auto
		{
			using _op = bond::operate<decltype(u)>;
			
			XTAL_USE U = XTAL_ALL_(u);
			XTAL_USE V = XTAL_ALL_(v);
			
			static_assert(aligned_n<U, V> <= N_cache);
			static_assert(_std::is_trivially_destructible_v<U>);
			static_assert(_std::is_trivially_destructible_v<V>);
			
			size_type i{};
			U &u1 = reinterpret_cast<U &>(m_cache[maligned_f<U>(i)]);
			V &v1 = reinterpret_cast<V &>(m_cache[maligned_f<V>(i)]);
			
			U  u0 = imagine_f<-complex_field_q<U>>(XTAL_REF_(u));
			V  v0 = XTAL_REF_(v);
			
		//	Divides the difference by the derived slope of the phasor:
			using _std::round;
			_std::swap(v1, v0); V v10 = v1 - v0; v10 -= round(v10) - z_(1);
			_std::swap(u1, u0); U u10 = u1 - u0; u10 *= root_f<-1, 1>(v10);

			return u10;
		}

	};
};


////////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////
}/////////////////////////////////////////////////////////////////////////////
XTAL_ENV_(pop)

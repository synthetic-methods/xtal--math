#pragma once
#include "./any.hh"

#include "../root.hh"




XTAL_ENV_(push)
namespace xtal::process::math::zavalishin
{/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

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
	class subtype: public bond::compose_s<S, subkind>
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
			return imagine_f<-complex_field_q<U>>(u1 - u0);// Should be imaginary division, right?
			return imagine_f<+complex_field_q<U>>(u1 - u0);
		}
		template <auto ...Is>
		XTAL_DEF_(return,inline)
		XTAL_LET method(auto &&u, algebra::d_::circular_q auto &&z)
		XTAL_0EX -> auto
		{
			using _op = bond::operate<decltype(u)>;
			
			XTAL_USE U = XTAL_ALL_(u);
			XTAL_USE V = debraced_t<XTAL_ALL_(z)>;
			
			static_assert(aligned_n<U, V> <= N_cache);
			static_assert(_std::is_trivially_destructible_v<U>);
			static_assert(_std::is_trivially_destructible_v<V>);
			
			size_type i{};
			U &u1 = reinterpret_cast<U &>(m_cache[maligned_f<U>(i)]);
			V &v1 = reinterpret_cast<V &>(m_cache[maligned_f<V>(i)]);
			
			U  u0 = XTAL_REF_(u);
			V  v0 = z(0);
			V  w_ = z(1);
			
		//	Divides the difference by the derived slope of the phasor:
			XTAL_LET N_zap = _op::exponent.depth;
			using _std::round;
			_std::swap(v1, v0); V v_ = v1 - v0; v_ -= round(v_);
			_std::swap(u1, u0); U u_ = u1 - u0; u_ *= root_f<-1, N_zap>(v_ - round(v_));

		//	Resets the state to zero if a phase/frequency discontinuity is detected:
			using _std::abs;
		//	u_ *= abs(v_ - w_) < _op::haplo_f(N_zap);
			//\
			return imagine_f<-complex_field_q<U>>(u_);// Should be imaginary division, right?
			return imagine_f<+complex_field_q<U>>(u_);
		}
		template <auto ...Is>
		XTAL_DEF_(return,inline)
		XTAL_LET method(algebra::d_::circular_q auto &&z, auto &&u)
		XTAL_0EX -> decltype(auto)
		{
			return method(XTAL_REF_(u), XTAL_REF_(z));
		}

		template <auto ...Is>
		XTAL_DEF_(return,inline)
		XTAL_LET method(algebra::d_::circular_q auto &&z)
		XTAL_0EX -> decltype(auto)
		{
			return wrap_f(method(XTAL_REF_(z) (0)));
		}

	};
};


////////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////
}/////////////////////////////////////////////////////////////////////////////
XTAL_ENV_(pop)

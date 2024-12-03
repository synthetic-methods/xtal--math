#pragma once
#include "./any.hh"






XTAL_ENV_(push)
namespace xtal::process::math::pade
{/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

template <int M_ism=1, int N_pow=1> requires sign_p<N_pow> struct   squishy;
template <int M_ism=1, int N_pow=1> requires sign_p<N_pow> using    squishy_t = process::confined_t<squishy<M_ism, N_pow>>;
template <int M_ism=1, int N_pow=1> requires sign_p<N_pow>
XTAL_DEF_(return,inline)
XTAL_LET squishy_f(auto &&o)
noexcept -> decltype(auto)
{
	return squishy_t<M_ism, N_pow>::function(XTAL_REF_(o));
}


////////////////////////////////////////////////////////////////////////////////

template <int N_pow>
struct squishy<0, N_pow>
{
	template <class S>
	class subtype : public bond::compose_s<S>
	{
		using S_ = bond::compose_s<S>;

	public:
		using S_::S_;

		template <auto ...>
		XTAL_DEF_(return,inline,static)
		XTAL_LET function(auto &&o)
		noexcept -> decltype(auto)
		{
			return XTAL_REF_(o);
		}

	};
};


////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

template <int M_ism> requires (0 < M_ism)
struct squishy<M_ism, 0>
{
	static constexpr int I_sgn = sign_n<(M_ism&1)^0, -1>;

	template <class S>
	class subtype : public bond::compose_s<S>
	{
		using S_ = bond::compose_s<S>;

	public:
		using S_::S_;

		template <auto ...>
		XTAL_DEF_(return,inline,static)
		XTAL_LET function(complex_field_q auto const &o)
		noexcept -> decltype(auto)
		{
			using O = XTAL_ALL_(o); using _op = bond::operate<O>;

			auto y = o.imag();
			auto x = o.real();
			auto const yy = y*y;
			auto const xx = x*x;
			auto const y_ = 2*x*y;
			auto const x_ = xx - I_sgn*yy;
			auto const w_ = xx + I_sgn*yy;
			auto const m_ = 1/w_;
			return complexion_f(x_, y_)*(m_);
		}

	};
};


////////////////////////////////////////////////////////////////////////////////

template <int M_ism> requires (0 < M_ism)
struct squishy<M_ism, 1>
{
	static constexpr int I_sgn = sign_n<M_ism&1, -1>;

	template <class S>
	class subtype : public bond::compose_s<S>
	{
		using S_ = bond::compose_s<S>;

	public:
		using S_::S_;

		template <auto ...>
		XTAL_DEF_(return,inline,static)
		XTAL_LET function(auto const &o)
		noexcept -> decltype(auto)
		{
			return o*o;
		}
		template <auto ...>
		XTAL_DEF_(return,inline,static)
		XTAL_LET function(complex_field_q auto const &o)
		noexcept -> decltype(auto)
		{
			using _op = bond::operate<decltype(o)>;

			auto const x = o.real(); auto xx = x*x;
			auto const y = o.imag(); auto yy = y*y;
			return complexion_f(xx - I_sgn*yy, _op::diplo_1*x*y);
		}

	};
};
template <int M_ism> requires (0 < M_ism)
struct squishy<M_ism,-1>
{
	template <class S>
	class subtype : public bond::compose_s<S>
	{
		using S_ = bond::compose_s<S>;

	public:
		using S_::S_;

		template <auto ...>
		XTAL_DEF_(return,inline,static)
		XTAL_LET function(auto &&o)
		noexcept -> decltype(auto)
		{
			using _op = bond::operate<decltype(o)>;
			return _op::alpha_1/squishy_f<M_ism, +1>(XTAL_REF_(o));
		}

	};
};


////////////////////////////////////////////////////////////////////////////////

template <int M_ism> requires (M_ism < 0)
struct squishy<M_ism, 1>
{
	template <class S>
	class subtype : public bond::compose_s<S>
	{
		using S_ = bond::compose_s<S>;

	public:
		using S_::S_;

		template <auto ...>
		XTAL_DEF_(return,inline,static)
		XTAL_LET function(simplex_field_q auto &&o)
		noexcept -> decltype(auto)
		{
			using _op = bond::operate<decltype(o)>;
			return _op::template root_f<+2>(XTAL_REF_(o));
		}
		template <auto ...>
		XTAL_DEF_(return,inline,static)
		XTAL_LET function(complex_field_q auto &&o)
		noexcept -> decltype(auto)
		{
			using _op = bond::operate<decltype(o)>;
			auto constexpr N_sqrt_half = (typename _op::alpha_type) 0.7071067811865475244008443621048490393L;
			auto  x = o.real();
			auto  y = o.imag();
			auto const n = function(x*x + y*y);
			return N_sqrt_half*complexion_f(function(n + x), function(n - x)*_op::assigned_f(y));
		}

	};
};
template <int M_ism> requires (M_ism < 0)
struct squishy<M_ism,-1>
{
	using dis = squishy_t<-1, 1>;

	template <class S>
	class subtype : public bond::compose_s<S>
	{
		using S_ = bond::compose_s<S>;

	public:
		using S_::S_;

		template <auto ...>
		XTAL_DEF_(return,inline,static)
		XTAL_LET function(simplex_field_q auto const &o)
		noexcept -> decltype(auto)
		{
			using O = XTAL_ALL_(o); using _op = bond::operate<O>;
			return _op::template root_f<-2>(o);
		}
		template <auto ...>
		XTAL_DEF_(return,inline,static)
		XTAL_LET function(complex_field_q auto const &o)
		noexcept -> auto
		{
			using _op = bond::operate<decltype(o)>;
			auto constexpr N_sqrt_half = (typename _op::alpha_type) 0.7071067811865475244008443621048490393L;
			auto  x = o.real();
			auto  y = o.imag();
			auto  u = function(x*x + y*y); x *= u*u;
			return N_sqrt_half*complexion_f(dis::function(u + x), dis::function(u - x)*-_op::assigned_f(y));
		}

	};
};


///////////////////////////////////////////////////////////////////////////////
}/////////////////////////////////////////////////////////////////////////////
XTAL_ENV_(pop)

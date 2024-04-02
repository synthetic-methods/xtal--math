#pragma once
#include "./any.hh"






XTAL_ENV_(push)
namespace xtal::math
{/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

template <int M_ism=1, int N_pow=1> requires sign_p<N_pow> XTAL_NEW square;
template <int M_ism=1, int N_pow=1> requires sign_p<N_pow> XTAL_USE square_t = process::confined_t<square<M_ism, N_pow>>;
template <int M_ism=1, int N_pow=1> requires sign_p<N_pow>
XTAL_FN2 square_f(auto &&o)
XTAL_0EX
{
	return square_t<M_ism, N_pow>::function(XTAL_REF_(o));
}


////////////////////////////////////////////////////////////////////////////////

template <int N_pow>
struct square<0, N_pow>
{
	template <class S>
	class subtype: public bond::compose_s<S>
	{
		using S_ = bond::compose_s<S>;

	public:
		using S_::S_;

		template <auto ...Is>
		XTAL_FN2 function(auto &&o)
		XTAL_0EX
		{
			return XTAL_REF_(o);
		}

	};
};


////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

template <int M_ism> requires (0 < M_ism)
struct square<M_ism, 0>
{
	XTAL_LET_(int) I_sgn = sign_n<M_ism&1^0, -1>;

	template <class S>
	class subtype: public bond::compose_s<S>
	{
		using S_ = bond::compose_s<S>;

	public:
		using S_::S_;

		template <auto ...Is>
		XTAL_FN2 function(complex_field_q auto const &o)
		XTAL_0EX
		{
			using O = XTAL_TYP_(o); using re = bond::realize<O>;

			auto y = o.imag();
			auto x = o.real();
			auto const yy = y*y;
			auto const xx = x*x;
			auto const y_ = 2*x*y;
			auto const x_ = xx - I_sgn*yy;
			auto const w_ = xx + I_sgn*yy;
			auto const m_ = 1/w_;
			y = y_*m_;
			x = x_*m_;
			return O {x, y};
		}

	};
};


////////////////////////////////////////////////////////////////////////////////

template <int M_ism> requires (0 < M_ism)
struct square<M_ism, 1>
{
	XTAL_LET_(int) I_sgn = sign_n<M_ism&1^0, -1>;

	template <class S>
	class subtype: public bond::compose_s<S>
	{
		using S_ = bond::compose_s<S>;

	public:
		using S_::S_;

		template <auto ...Is>
		XTAL_FN2 function(simplex_field_q auto const &o)
		XTAL_0EX
		{
			return o*o;
		}
		template <auto ...Is>
		XTAL_FN2 function(complex_field_q auto const &o)
		XTAL_0EX
		{
			using O = XTAL_TYP_(o); using re = bond::realize<O>;

			auto const x = o.real(), xx = x*x;
			auto const y = o.imag(), yy = y*y;
			return O {xx - I_sgn*yy, 2*x*y};
		}

	};
};
template <int M_ism> requires (0 < M_ism)
struct square<M_ism,-1>
{
	template <class S>
	class subtype: public bond::compose_s<S>
	{
		using S_ = bond::compose_s<S>;

	public:
		using S_::S_;

		template <auto ...Is>
		XTAL_FN2 function(auto &&o)
		XTAL_0EX
		{
			return 1/square_f<M_ism, +1>(XTAL_REF_(o));
		}

	};
};


////////////////////////////////////////////////////////////////////////////////

template <int M_ism> requires (M_ism < 0)
struct square<M_ism, 1>
{
	template <class S>
	class subtype: public bond::compose_s<S>
	{
		using S_ = bond::compose_s<S>;

	public:
		using S_::S_;

		template <auto ...Is>
		XTAL_FN2 function(simplex_field_q auto const &o)
		XTAL_0EX
		{
			using O = XTAL_TYP_(o); using re = bond::realize<O>;
			return re::template root_f<+2>(o);
		}
		template <auto ...Is>
		XTAL_FN2 function(complex_field_q auto const &o)
		XTAL_0EX
		{
			using O = XTAL_TYP_(o); using re = bond::realize<O>;
			auto constexpr O_k = re::unsquare_f(re::haplo_1);
			auto  x = o.real();
			auto  y = o.imag();
			auto const n = function(x*x + y*y);
			return O_k*O {function(n + x), function(n - x)*re::sign_f(y)};
		}

	};
};
template <int M_ism> requires (M_ism < 0)
struct square<M_ism,-1>
{
	using dis = square_t<-1, 1>;

	template <class S>
	class subtype: public bond::compose_s<S>
	{
		using S_ = bond::compose_s<S>;

	public:
		using S_::S_;

		template <auto ...Is>
		XTAL_FN2 function(simplex_field_q auto const &o)
		XTAL_0EX
		{
			using O = XTAL_TYP_(o); using re = bond::realize<O>;
			return re::template root_f<-2>(o);
		}
		template <auto ...Is>
		XTAL_FN2 function(complex_field_q auto const &o)
		XTAL_0EX
		{
			using O = XTAL_TYP_(o); using re = bond::realize<O>;
			auto constexpr O_k = re::unsquare_f(re::haplo_1);
			auto  x = o.real();
			auto  y = o.imag();
			auto  u = function(x*x + y*y); x *= u*u;
			return O_k*O {dis::function(u + x), dis::function(u - x)*-re::sign_f(y)};
		}

	};
};


///////////////////////////////////////////////////////////////////////////////
}/////////////////////////////////////////////////////////////////////////////
XTAL_ENV_(pop)

#pragma once
#include "./any.hh"

#include "./square.hh"




XTAL_ENV_(push)
namespace xtal::process::math
{/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

template <int M_pow=1, int M_zap=-1>
XTAL_TYP root;

template <int M_pow=1, int M_zap=-1>
XTAL_USE root_t = process::confined_t<root<M_pow, M_zap>>;

template <int M_pow=1, int M_zap=-1>
XTAL_DEF_(return,inline)
XTAL_LET root_f(auto &&o)
XTAL_0EX -> decltype(auto)
{
	return root_t<M_pow, M_zap>::function(XTAL_REF_(o));
}


////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

template <int M_pow, int M_zap> requires (0 <= M_zap)
struct root<M_pow, M_zap>
{
	using subkind = root<M_pow>;

	template <class S>
	class subtype : public bond::compose_s<S, subkind>
	{
		using S_ = bond::compose_s<S, subkind>;

	public:
		using S_::S_;

		template <auto ...Ns>
		XTAL_DEF_(return,inline)
		XTAL_SET function(auto &&o)
		XTAL_0EX -> auto
		{
			using _op = bond::operate<XTAL_ALL_(o)>;
			return S_::template function<Ns...>(_op::template punctured_f<_op::designed_f(M_pow)*M_zap>(XTAL_REF_(o)));
		}

	};
};
template <int M_pow, int M_zap> requires (M_zap <  0)
struct root<M_pow, M_zap>
{
	template <class S>
	class subtype : public bond::compose_s<S>
	{
		using S_ = bond::compose_s<S>;

	public:
		using S_::S_;

		template <auto ...Ns>
		XTAL_DEF_(return,inline)
		XTAL_SET function(auto &&o)
		XTAL_0EX -> auto
		{
			using _op = bond::operate<XTAL_ALL_(o)>;
			XTAL_IF0
			XTAL_0IF (M_pow ==  1) {return                        XTAL_REF_(o)  ;}
			XTAL_0IF (M_pow ==  2) {return                   sqrt(XTAL_REF_(o)) ;}
			XTAL_0IF (M_pow ==  3) {return                   cbrt(XTAL_REF_(o)) ;}
			XTAL_0IF (M_pow ==  4) {return              sqrt(sqrt(XTAL_REF_(o)));}
			XTAL_0IF (M_pow == -1) {return _op::alpha_1/          XTAL_REF_(o)  ;}
			XTAL_0IF (M_pow == -2) {return _op::alpha_1/     sqrt(XTAL_REF_(o)) ;}
			XTAL_0IF (M_pow == -3) {return _op::alpha_1/     cbrt(XTAL_REF_(o)) ;}
			/**/
			XTAL_0IF (M_pow == -4) {return _op::alpha_1/sqrt(sqrt(XTAL_REF_(o)));}
			/*/
			XTAL_0IF (M_pow == -4) {
				auto constexpr c0 =   _op::ratio_f( 5, 4);
				auto const     c1 = o*_op::ratio_f(-1, 4);
				auto y = XTAL_REF_(o);
				y *= _op::template root_f<-2, 0>(y);
				y  = _op::template root_f<-2, 0>(y);
				y *= _op::accumulate_f(c0, c1, square_f(square_f(y)));
				y *= _op::accumulate_f(c0, c1, square_f(square_f(y)));
				return y;
			}
			/***/
			XTAL_0IF (1 == _op::bit_count_f(_op::designed_f(M_pow))) {//FIXME: Should work without `designed_f`!
				return root_f<(M_pow >> 1)>(sqrt(XTAL_REF_(o)));
			}

		}
		/*/
		template <auto ...Ns>
		XTAL_DEF_(return,inline)
		XTAL_SET function(complex_number_q auto &&o)
		XTAL_0EX -> auto
		{
			using _op = bond::operate<XTAL_ALL_(o)>;
			auto constexpr N_sqrt_half = (typename _op::alpha_type) 0.7071067811865475244008443621048490393L;

			XTAL_IF0
			XTAL_0IF (M_pow ==  1) {
				return XTAL_REF_(o);
			}
			XTAL_0IF (M_pow ==  2) {
				auto const x   = o.real();
				auto const y   = o.imag();
				auto const n   = function<Ns...>(x*x + y*y);
				auto const lhs = function<Ns...>(n + x);
				auto const rhs = function<Ns...>(n - x);
				return N_sqrt_half*complexion_f(lhs, rhs*_op::assigned_f(y));
			}
			XTAL_0IF (M_pow == -1) {
				return _op::alpha_1/XTAL_REF_(o);
			}
			XTAL_0IF (M_pow == -2) {
				using dis = root_t<2>;
				auto       x   =  o.real();
				auto       y   =  o.imag();
				auto const w    =  x*x + y*y;
				auto const u_dn =  function<Ns...>(w);
				auto const u_up =  u_dn*w;
				auto const v_re =  u_dn*dis::template function<Ns...>(u_up + x);
				auto const v_im = -u_dn*dis::template function<Ns...>(u_up - x);

				return N_sqrt_half*complexion_f(v_re, v_im*_op::assigned_f(y));
			}
		}
		/***/

	};
};


///////////////////////////////////////////////////////////////////////////////
}/////////////////////////////////////////////////////////////////////////////
XTAL_ENV_(pop)

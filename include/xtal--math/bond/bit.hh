#pragma once
#include "./any.hh"






XTAL_ENV_(push)
namespace xtal::bond::math
{/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////

template <class V>
XTAL_DEF_(return,inline,let)
bit_depth_f()
noexcept -> auto
{
	return static_cast<int>(sizeof(V) << 3);
}
XTAL_DEF_(return,inline,let)
bit_depth_f(variable_q auto v)
noexcept -> int
{
	return static_cast<int>(sizeof(v) << 3);
}
XTAL_DEF_(return,inline,let)
bit_depth_f(constant_q auto w)
noexcept -> int
{
	return constant_n<bit_depth_f(w())>;
}


////////////////////////////////////////////////////////////////////////////////

XTAL_DEF_(return,inline,let)
bit_sign_f(integral_variable_q auto o)
noexcept -> int
{
	return static_cast<int>(ordinal_f(o) >> bit_depth_f<decltype(o)>() - one);
}
XTAL_DEF_(return,inline,let)
bit_sign_f(real_variable_q auto o)
noexcept -> int
{
	using _fix = bond::fixture<decltype(o)>;
	return bit_sign_f(_xtd::bit_cast<typename _fix::delta_type>(o));
}
XTAL_DEF_(return,inline,let)
bit_sign_f(constant_q auto w)
noexcept -> auto
{
	return constant_n<bit_sign_f(w())>;
}


///\returns the number of bits set in `u`. \

XTAL_DEF_(return,inline,let)
bit_count_f(cardinal_variable_q auto u)
noexcept -> int
{
	return _std::popcount(u);
}
XTAL_DEF_(return,inline,let)
bit_count_f(ordinal_variable_q auto v)
noexcept -> int
{
	int  x = bit_sign_f(v);
	v ^= x;
	v -= x;
	int  y = bit_count_f(cardinal_f(v));
	y ^= x;
	y -= x;
	return y;
}

XTAL_DEF_(return,inline,let)
bit_count_f(constant_q auto w)
noexcept -> auto
{
	return constant_n<bit_count_f(w())>;
}


template <int N_zero=-1> requires (N_zero == -1)
XTAL_DEF_(return,inline,let)
bit_floor_f(cardinal_variable_q auto u)
noexcept -> int
{
	/*/
	U const  z =      one == u;
	U const _z = bond::fixture<decltype(u)>::unit.depth&-z;
	return _std::bit_width(u) - (one << _z) + z;// 0 -> -1023
	/*/
	return _std::bit_width(u) - one;// 0 -> -1
	/***/
}
template <int N_zero=-1> requires (N_zero != -1)
XTAL_DEF_(return,inline,let)
bit_floor_f(cardinal_variable_q auto u)
noexcept -> int
{
	int constexpr N = N_zero + one;
	int y = bit_floor_f(u); y += (N&bit_sign_f(y));
	return y;
}
template <int N_zero=-1>
XTAL_DEF_(return,inline,let)
bit_floor_f(ordinal_variable_q auto v)
noexcept -> int
{
	int  x = bit_sign_f(v);
	v ^= x;
	v -= x;
	int  y = bit_floor_f<N_zero>(cardinal_f(v));
	y ^= x;
	y -= x;
	return y;
}

template <int N_zero>
XTAL_DEF_(return,inline,let)
bit_floor_f(real_variable_q auto &&x)
noexcept -> int
{
	using          X_fix =  bond::fixture<absolve_u<decltype(x)>>;
	auto constexpr Z     = -static_cast<int>(X_fix::unit.mark);

	auto constexpr      unit = X_fix::    unit;
	auto constexpr      sign = X_fix::    sign;
	auto constexpr  exponent = X_fix::exponent;
	using U = typename X_fix::sigma_type;
	using V = typename X_fix::delta_type;

	auto n = _xtd::bit_cast<V>(XTAL_REF_(x));
	n >>= exponent.shift;
	n  &= exponent. mark;
	n  -=     unit. mark;
	if constexpr (N_zero != Z) {
		n  -=      N_zero;
		n  &=  -n >> sign.shift;
		n  +=      N_zero;
	}
	return n;
}
template <int N_zero>
XTAL_DEF_(return,inline,let)
bit_floor_f(complex_variable_q auto &&x)
noexcept -> int
{
	return bit_floor_f<N_zero>(norm(XTAL_REF_(x))) >> one;
}
XTAL_DEF_(return,inline,let)
bit_floor_f(auto &&x)
noexcept -> int
requires real_variable_q<absolve_u<decltype(x)>>
{
	using          X_fix =  bond::fixture<absolve_u<decltype(x)>>;
	auto constexpr Z     = -static_cast<int>(X_fix::unit.mark);

	return bit_floor_f<Z>(XTAL_REF_(x));
}

template <int N_zero=-1>
XTAL_DEF_(return,inline,let)
bit_floor_f(constant_q auto w)
noexcept -> auto
{
	return constant_n<bit_floor_f<N_zero>(w())>;
}


template <int N_zero=-1> requires (N_zero == -1)
XTAL_DEF_(return,inline,let)
bit_ceiling_f(cardinal_variable_q auto u)
noexcept -> int
{
	int const z = u == zero;
	//\
	return _std::bit_width(u - !z) - (z << unit.depth) + z;// 0 -> -1023
	return _std::bit_width(u - !z)                     - z;// 0 -> -1
}
template <int N_zero=-1> requires (N_zero != -1)
XTAL_DEF_(return,inline,let)
bit_ceiling_f(cardinal_variable_q auto u)
noexcept -> int
{
	int constexpr N = N_zero + one;
	int y = bit_ceiling_f(u); y += N&bit_sign_f(y);
	return y;
}
template <int N_zero=-1>
XTAL_DEF_(return,inline,let)
bit_ceiling_f(ordinal_variable_q auto v)
noexcept -> int
{
	int  x = bit_sign_f(v);
	v ^= x;
	v -= x;
	int  y = bit_ceiling_f<N_zero>(cardinal_f(v));
	y ^= x;
	y -= x;
	return y;
}

template <int N_zero=+1>
XTAL_DEF_(return,inline,let)
bit_ceiling_f(real_variable_q auto &&x)
noexcept -> int
{
	using X_fix = bond::fixture<decltype(x)>;
	return bit_floor_f<N_zero>(X_fix::diplo_1*X_fix::dnsilon_1*XTAL_REF_(x));
}
template <int N_zero=+1>
XTAL_DEF_(return,inline,let)
bit_ceiling_f(complex_variable_q auto &&x)
noexcept -> int
{
	return bit_ceiling_f<N_zero>(norm(XTAL_REF_(x))) >> 1;
}

template <int N_zero=-1>
XTAL_DEF_(return,inline,let)
bit_ceiling_f(constant_q auto w)
noexcept -> auto
{
	return constant_n<bit_ceiling_f<N_zero>(w())>;
}


///\returns the bitwise-reversal of `u`, \
restricted to `N_subdepth` when `0 < N_subdepth < sizeof(u) << 3U`. \

///\note Requires `log2(sizeof(u) << 3U)` iterations. \

XTAL_DEF_(return,let)
bit_reverse_f(cardinal_variable_q auto u, int const &n_subdepth)
noexcept -> auto
{
	using          U =      decltype(u);
	auto constexpr N = bit_depth_f<U>();

	#pragma inline
	for (U m = -1, i = N; i >>= 1;) {
		m ^= m<<i;
		u = (u&m)<<i | (u&~m)>>i;
	}
	u >>= N - n_subdepth; assert(0 < n_subdepth and n_subdepth <= N);
	return u;
}
XTAL_DEF_(return,inline,let)
bit_reverse_f(ordinal_variable_q auto v, int const &n_subdepth)
noexcept -> auto
{
	decltype(v) x = bit_sign_f(v);
	v ^= x;
	v -= x;
	v  = ordinal_f(bit_reverse_f(cardinal_f(v), n_subdepth));
	v ^= x;
	v -= x;
	return v;
}

XTAL_DEF_(return,inline,let)
bit_reverse_f(constant_q auto w, int const &n_subdepth)
noexcept -> auto
{
	return constant_n<bit_reverse_f(w(), n_subdepth)>;
}
template <int N_subdepth>
XTAL_DEF_(return,inline,let)
bit_reverse_f(auto &&x)
noexcept -> auto
{
	auto constexpr N_depth = bond::fixture<decltype(x)>::full.depth;
	//\
	int constexpr n_subdepth = below_m<N_depth, (unsigned) N_subdepth>;
	int constexpr n_subdepth = 0 < N_subdepth? N_subdepth: N_depth;
	return bit_reverse_f(XTAL_REF_(x), n_subdepth);
}

///////////////////////////////////////////////////////////////////////////////

XTAL_DEF_(return,let)
bit_representation_f(real_variable_q auto x)
noexcept -> auto
{
	using X_fix =  bond::fixture<decltype(x)>;
	using U = typename X_fix::sigma_type;
	using V = typename X_fix::delta_type;

	U constexpr N = X_fix::unit.mark + X_fix::fraction.depth;
	U constexpr M =            one << X_fix::fraction.depth;
	
	auto const o = _xtd::bit_cast<U>(x);
	V const z = static_cast<V>(o) >> X_fix::positive.depth;
	V const n = N - (o << X_fix::sign.depth >> X_fix::sign.depth + X_fix::exponent.shift);
	V       m = M | (o&X_fix::fraction.mask);
	m  ^= z;
	m  -= z;
	return atom::couple_t<V[2]>{m, n};
}
XTAL_DEF_(return,inline,let)
bit_presentation_f(atom::couple_q auto const &mn)
noexcept -> auto
{
	auto const &[m, n] = mn;
	using MN_op = bond::fixture<decltype(m), decltype(n)>;
	return static_cast<typename MN_op::alpha_type>(m)*MN_op::haplo_f(n);
}


////////////////////////////////////////////////////////////////////////////////


///\returns the fractional component of `x` as a full-width `delta_type`.

template <class T_return=void>
XTAL_DEF_(return,inline,let)
bit_fraction_f()
noexcept -> auto
{
	using Y       = T_return;
	using Y_fix    = bond::fixture<Y>;

	XTAL_IF0
	XTAL_0IF (integral_q<Y>) {
		return Y_fix::diplo_f(Y_fix::full.depth);
	}
	XTAL_0IF (    real_q<Y>) {
		return Y_fix::haplo_f(Y_fix::full.depth);
	}
	XTAL_0IF_(void)
}
template <class Y_return=void>
XTAL_DEF_(return,inline,let)
bit_fraction_f(integral_variable_q auto x)
noexcept -> auto
{
	using X       =              XTAL_ALL_(x);
	using X_fix    =          bond::fixture<X>;
	using X_alpha = typename X_fix::alpha_type;
	using X_sigma = typename X_fix::sigma_type;
	using X_delta = typename X_fix::delta_type;

	using Y       = complete_t<Y_return, X_alpha>;
	using Y_fix    = bond::fixture<Y>;
	using Y_alpha = typename Y_fix::alpha_type;
	using Y_sigma = typename Y_fix::sigma_type;
	using Y_delta = typename Y_fix::delta_type;

	XTAL_IF0
	XTAL_0IF (cardinal_q<Y>) {
		return Y_fix::sigma_f(x);
	}
	XTAL_0IF ( ordinal_q<Y>) {
		return Y_fix::delta_f(x);
	}
	XTAL_0IF (    real_q<Y>) {
		//\
		return Y_fix::alpha_f(x)*Y_fix::haplo_f(X_fix::full.depth);
		return Y_fix::alpha_f(static_cast<X_delta>(x))*Y_fix::haplo_f(X_fix::full.depth);
	}
	XTAL_0IF_(void)
}
template <class Y_return=void>
XTAL_DEF_(return,inline,let)
bit_fraction_f(real_variable_q auto x)
noexcept -> auto
{
	using X       =              XTAL_ALL_(x);
	using X_fix    =          bond::fixture<X>;
	using X_alpha = typename X_fix::alpha_type;
	using X_sigma = typename X_fix::sigma_type;
	using X_delta = typename X_fix::delta_type;

	using Y       = complete_t<Y_return, X_delta>;
	using Y_fix    = bond::fixture<Y>;
	using Y_alpha = typename Y_fix::alpha_type;
	using Y_sigma = typename Y_fix::sigma_type;
	using Y_delta = typename Y_fix::delta_type;

	XTAL_IF0
	XTAL_0IF_(consteval) {
		XTAL_IF0
		XTAL_0IF (    real_q<Y>) {
			return static_cast<Y>(bit_fraction_f<Y_delta>(x))*bit_fraction_f<Y_alpha>();
		}
		XTAL_0IF ( ordinal_q<Y>) {
			auto constexpr N_exp = X_fix::exponent.shift;
			auto constexpr M_exp = X_fix::unit.mark + X_fix::unit.shift - Y_fix::full.depth;
			auto constexpr M_sgn = X_fix::sign.mask;

			auto o = _xtd::bit_cast<X_sigma>(x);

			X_delta o_ = o &  M_sgn; o  ^= o_; o_ >>= X_fix::sign.shift;
			X_sigma x  = o >> N_exp; x  -= M_exp;
			X_sigma u  = x != 0    ; u <<= X_fix::fraction.depth;

			auto const x_on = _xtd::bit_cast<X_delta>(~x) >> X_fix::positive.depth;
			auto const x_up = x_on & x;
			auto const x_dn = x_up - x;

			o  &=  X_fix::fraction.mask;
			o  ^=  u;
			o <<=  x_up;
			o >>=  x_dn;
		//	o  |=  1;
			o  ^=  o_;
			o  -=  o_;
			return static_cast<Y>(_xtd::bit_cast<X_delta>(o));
		}
		XTAL_0IF (cardinal_q<Y>) {
			return _xtd::bit_cast<Y>(bit_fraction_f<Y_delta>(x));
		}
		XTAL_0IF_(void)
	}
	XTAL_0IF_(else) {
		x -= round(x);
		XTAL_IF0
		XTAL_0IF (    real_q<Y>) {
			return static_cast<Y>(x);
		}
		XTAL_0IF ( ordinal_q<Y>) {
			return static_cast<Y>(x*bit_fraction_f<Y_delta>());
		}
		XTAL_0IF (cardinal_q<Y>) {
			return _xtd::bit_cast<Y>(bit_fraction_f<Y_delta>(x));
		}
		XTAL_0IF_(void)
	}
}
template <class T_return=void>
XTAL_DEF_(return,inline,let)
bit_fraction_f(complex_variable_q auto const &x)
noexcept -> auto
{
	using X       = XTAL_ALL_(x);
	using X_fix    = bond::fixture<X>;
	using X_alpha = typename X_fix::alpha_type;
	using X_sigma = typename X_fix::sigma_type;
	using X_delta = typename X_fix::delta_type;

	using U = _std::conditional_t<integral_q<typename X::value_type>, X_alpha, X_delta>;
	using Y = complete_t<T_return, _std::complex<U>>;
	using V = typename Y::value_type;

	return Y{bit_fraction_f<V>(x.real()), bit_fraction_f<V>(x.imag())};
}


///////////////////////////////////////////////////////////////////////////////
}/////////////////////////////////////////////////////////////////////////////
XTAL_ENV_(pop)

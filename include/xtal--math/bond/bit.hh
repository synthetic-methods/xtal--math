#pragma once
#include "./any.hh"






XTAL_ENV_(push)
namespace xtal::bond::math
{/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

template <class T=void>
XTAL_DEF_(return,let)
bit_exchange_f(auto &&u)
noexcept -> auto
{
	using U     = XTAL_ALL_(u);
	using U_fit = bond::fit<U>;
	using T_fit = bond::fit<T>;
	using U_sigma = typename U_fit::sigma_type;
	using U_delta = typename U_fit::delta_type;
	using U_alpha = typename U_fit::alpha_type;
	using T_sigma = typename T_fit::sigma_type;
	using T_delta = typename T_fit::delta_type;
	using T_alpha = typename T_fit::alpha_type;
	if constexpr (complete_q<T>) {
		XTAL_IF0
		XTAL_0IF (sizeof(T) == sizeof(U)) {
			return _xtd::bit_cast<T>(XTAL_REF_(u));
		}
		XTAL_0IF (cardinal_variable_q<U>) {
			return bit_exchange_f<T>(static_cast<T_sigma>(XTAL_REF_(u)));
		}
		XTAL_0IF ( ordinal_variable_q<U>) {
			return bit_exchange_f<T>(static_cast<T_delta>(XTAL_REF_(u)));
		}
		XTAL_0IF (    real_variable_q<U>) {
			return bit_exchange_f<T>(static_cast<T_alpha>(XTAL_REF_(u)));
		}
		XTAL_0IF_(void)
	}
	else {
		XTAL_IF0
		XTAL_0IF (integral_variable_q<U>) {
			return _xtd::bit_cast<U_alpha>(XTAL_REF_(u));
		}
		XTAL_0IF (    real_variable_q<U>) {
			return _xtd::bit_cast<U_sigma>(XTAL_REF_(u));
		}
		XTAL_0IF_(void)
	}
}
template <int N_val=0>
XTAL_DEF_(return,inline,let)
bit_flag_f(numeric_variable_q auto const x)
noexcept -> XTAL_ALL_(x)
{
	using X     = XTAL_ALL_(x);
	using X_fit = bond::fit<X>;
	X constexpr X_val{~N_val};
	return static_cast<X>(X_val != x);
}


///////////////////////////////////////////////////////////////////////////////

template <int N_dir=0>
XTAL_DEF_(return,inline,let)
bit_axes_f(ordinal_variable_q auto u)
noexcept -> _std::array<XTAL_ALL_(u), 2>
{
	using U     = XTAL_ALL_(u);
	using U_fit = bond::fit<U>;
	auto const dn =  u >> U_fit::sign.shift;
	auto const up = ~dn;
	XTAL_IF0
	XTAL_0IF (0 < N_dir) {return { u&up,-u&dn};}
	XTAL_0IF (N_dir < 0) {return {-u&dn, u&up};}
	XTAL_0IF_(else)      {return { u   , u   };}
}

template <int N_dir=0>
XTAL_DEF_(return,inline,let)
bit_axis_f(ordinal_variable_q auto u)
noexcept -> auto
{
	using U     = XTAL_ALL_(u);
	using U_fit = bond::fit<U>;
	auto const dn =  u >> U_fit::sign.shift;
	auto const up = ~dn;
	XTAL_IF0
	XTAL_0IF (0 < N_dir) {return  u&up;}
	XTAL_0IF (N_dir < 0) {return -u&dn;}
	XTAL_0IF_(else)      {return  u   ;}
}
template <int N_dir=0>
XTAL_DEF_(return,inline,let)
bit_axis_f(ordinal_constant_q auto u_)
noexcept -> auto
{
	return constant_t<bit_axis_f<N_dir>(u_())>{};
}


///////////////////////////////////////////////////////////////////////////////

template <int N_ord=0>
XTAL_DEF_(return,inline,let)
bit_shift_f(integral_variable_q auto u)
noexcept -> auto
{
	XTAL_IF0
	XTAL_0IF (0 <= N_ord) {
		return XTAL_MOV_(u) << N_ord;
	}
	XTAL_0IF (N_ord <  0) {
		auto const v = bit_sign_f(u);
		u  ^= v;
		u  -= v;
		u >>= N_ord;
		u  ^= v;
		u  -= v;
		return u;
	}
}
XTAL_DEF_(return,inline,let)
bit_shift_f(integral_variable_q auto u, ordinal_q auto n)
noexcept -> auto
{
	auto const [up, dn] = bit_axes_f<+1>(n);
	u <<= up;
	u >>= dn;
	return u;
}

template <int N_depth>
XTAL_DEF_(return,inline,let)
bit_mask_f()
noexcept -> auto
{
	int constexpr N_width = 1 << N_depth;
	int constexpr M_width = N_width  - 1;
	return M_width;
}
template <int N_depth>
XTAL_DEF_(return,inline,let)
bit_mask_f(ordinal_q auto const u)
noexcept -> auto
{
	return u&bit_mask_f<N_depth>();
}


///////////////////////////////////////////////////////////////////////////////

template <int N_ord=0>
XTAL_DEF_(return,inline,let)
bit_exponent_f(real_variable_q auto u)
noexcept -> auto
{
	using U     = XTAL_ALL_(u);
	using U_fit = bond::fit<U>;
	using U_delta = typename U_fit::delta_type;
	using U_sigma = typename U_fit::sigma_type;

	if constexpr (U_fit::use_IEC) {
		auto x = bit_exchange_f(XTAL_MOV_(u));
		x >>= U_fit::exponent.shift;
		x  &= U_fit::exponent. mark;
		x  -= U_fit::    unit. mark;
		return bit_shift_f<N_ord>(_xtd::bit_cast<U_delta>(XTAL_MOV_(x)));
	}
	else {
		return bit_shift_f<N_ord>(   static_cast<U_delta>(_std::ilogb(XTAL_MOV_(u))));
	}
}
template <int N_dir=0>
XTAL_DEF_(return,inline,let)
bit_exponent_f(complex_variable_q auto &&w)
noexcept -> auto
{
	//\
	return bond::fit<int>::maximum_f(std::ilogb(w.real()), std::ilogb(w.real()));
	return bit_exponent_f<N_dir>(norm(XTAL_REF_(w))) >> one;
}


////////////////////////////////////////////////////////////////////////////////

/*!
\returns The number of bits constituting `V`.
*/
template <class V>
XTAL_DEF_(return,inline,let)
bit_depth_f()
noexcept -> int
{
	return static_cast<int>(sizeof(V) << 3);
}
/*!
\returns The number of bits constituting `v`.
*/
XTAL_DEF_(return,inline,let)
bit_depth_f(variable_q auto v)
noexcept -> int
{
	return static_cast<int>(sizeof(v) << 3);
}
/*!
\returns The number of bits constituting `w`.
*/
XTAL_DEF_(return,inline,let)
bit_depth_f(constant_q auto w)
noexcept -> int
{
	return constant_t<bit_depth_f(w())>{};
}


////////////////////////////////////////////////////////////////////////////////

/*!
\returns `-1` if the argument `< 0`, `0` otherwise.
*/
template <class T=void>
XTAL_DEF_(return,inline,let)
bit_sign_f(integral_variable_q auto o)
noexcept -> auto
{
	using T_fit = bond::fit<T>;
	using T_sigma = typename T_fit::sigma_type;
	using T_delta = typename T_fit::delta_type;
	using T_alpha = typename T_fit::alpha_type;

	auto constexpr N = bit_depth_f<decltype(o)>() - 1;
	auto const     i = ordinal_f(XTAL_REF_(o));
	XTAL_IF0
	XTAL_0IF (incomplete_q<T> or same_q<T, T_delta>) {
		return i >> N;
	}
	XTAL_0IF (integral_q<T> and sizeof(T) == sizeof(i)) {
		return _xtd::bit_cast<T>(i >> N);
	}
	XTAL_0IF (integral_q<T> and sizeof(T) != sizeof(i)) {
		auto const j = static_cast<T_delta>(i);
		return _xtd::bit_cast<T>(j >> N);
	}
	XTAL_0IF (real_q<T>) {
		auto k = static_cast<T_delta>(i);
		k >>= N;
		k  &= T_fit::sign.mask;
		k  |= T_fit::unit.mask;
		return _xtd::bit_cast<T>(k);
	}
}
/*!
\returns `-1` if the argument `< 0`, `0` otherwise.
*/
template <class T=void>
XTAL_DEF_(return,inline,let)
bit_sign_f(real_variable_q auto o)
noexcept -> auto
{
	using _fit = bond::fit<decltype(o)>;
	return bit_sign_f<T>(_xtd::bit_cast<typename _fit::delta_type>(o));
}
/*!
\returns `-1` if the argument `< 0`, `0` otherwise.
*/
template <class T=void>
XTAL_DEF_(return,inline,let)
bit_sign_f(constant_q auto w)
noexcept -> auto
{
	return constant_t<bit_sign_f<T>(w())>{};
}

template <int N_dir=0, cardinal_variable_q U>
XTAL_DEF_(return,inline,let)
bit_sign_f(U const &x, U const &y)
noexcept -> auto
{
	XTAL_IF0
	XTAL_0IF (N_dir == 0) {return  ~U{};}
	XTAL_0IF (N_dir <  0) {return  bit_sign_f(x - y);}
	XTAL_0IF (0 <  N_dir) {return  bit_sign_f(y - x);}
}
template <int N_dir=0,  ordinal_variable_q V>
XTAL_DEF_(return,inline,let)
bit_sign_f(V const &x, V const &y)
noexcept -> auto
{
	using U = cardinal_t<V>;
	return bit_sign_f<N_dir>(reinterpret_cast<U const &>(x), reinterpret_cast<U const &>(y));
}


////////////////////////////////////////////////////////////////////////////////

template <unsigned int N_mask>
XTAL_DEF_(return,inline,let)
bit_clasp_f(integral_variable_q auto x)
noexcept -> auto
{
	x &= ~bit_sign_f(x); x -= N_mask;
	x &=  bit_sign_f(x); x += N_mask;
	return x;
}
template <  signed int N_exit>
XTAL_DEF_(return,inline,let)
bit_clamp_f(integral_variable_q auto const x)
noexcept -> auto
{
	int constexpr N_width = 1 << N_exit;
	return bit_clasp_f<bit_mask_f<N_width>()>(x);
}


////////////////////////////////////////////////////////////////////////////////

template <int N_sgn=1>
XTAL_DEF_(return,inline,let)
bit_swatch_f(ordinal_variable_q auto &&u)
noexcept -> auto
{
	using U     = XTAL_ALL_(u);
	using U_fit = bond::fit<U>;
	using U_delta = typename U_fit::delta_type;
	using U_alpha = typename U_fit::alpha_type;
	using W_alpha = _std::array<U_alpha, 2>;
	auto             mask_dn  =  XTAL_REF_(u) >> U_fit::sign.shift;
	auto             mask_up  = ~mask_dn;
	auto constexpr       _1   =  U_fit::unit.mask|U_fit::sign.mask;
	auto const           _dn  =  bit_exchange_f(mask_dn&_1), dn = -_dn;
	auto const           _up  =  bit_exchange_f(mask_up&_1), up = -_up;
	XTAL_IF0
	XTAL_0IF (0 < N_sgn) {return W_alpha{ dn,  up};}
	XTAL_0IF (N_sgn < 0) {return W_alpha{_dn, _up};}
}


////////////////////////////////////////////////////////////////////////////////

template <int N_dir=0, cardinal_variable_q U>
XTAL_DEF_(inline,let)
bit_swap_f(U &x, U &y)
noexcept -> auto
{
	auto const w = bit_sign_f<N_dir>(x, y);
	x ^=   y;
	y ^= w&x;
	x ^=   y;
	return w;
}
template <int N_dir=0,  ordinal_variable_q V>
XTAL_DEF_(inline,let)
bit_swap_f(V &x, V &y)
noexcept -> auto
{
	using U = cardinal_t<V>;
	return bit_swap_f<N_dir>(reinterpret_cast<U &>(x), reinterpret_cast<U &>(y));
}


////////////////////////////////////////////////////////////////////////////////

template <int N_sel=0, cardinal_variable_q U>
XTAL_DEF_(inline,let)
bit_extrema_f(U x, U y)
noexcept -> auto
{
	auto const yx = y*x;
	while (x != y) {
		bit_swap_f<-1>(x, y); x -= y;
	}
	XTAL_IF0
	XTAL_0IF (N_sel < 0) {
		return x;
	}
	XTAL_0IF (0 < N_sel) {
		return yx/x;
	}
	XTAL_0IF_(else) {
		return bond::pack_f(x, yx/x);
	}
}


////////////////////////////////////////////////////////////////////////////////

template <class U>
XTAL_DEF_(return,inline,let)
bit_zoom_f(U &&x, U &&y)
noexcept -> auto
{
	using          U_fit = bond::fit<U>;
	auto constexpr U_exp = U_fit::exponent.mask;
	auto i_dn = U_exp&bit_exchange_f(XTAL_REF_(x));
	auto i_up = U_exp&bit_exchange_f(XTAL_REF_(y));
	(void) bit_swap_f<+1>(i_dn, i_up);
	return bit_exchange_f(U_exp - i_up);
}


////////////////////////////////////////////////////////////////////////////////

/*!
\returns The number of bits set in `u`.
*/
XTAL_DEF_(return,inline,let)
bit_count_f(cardinal_variable_q auto u)
noexcept -> int
{
	return _std::popcount(u);
}
/*!
\returns The number of bits set in `v`.
*/
XTAL_DEF_(return,inline,let)
bit_count_f(ordinal_variable_q auto v)
noexcept -> int
{
	int  x = bit_sign_f<int>(v);
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
	return constant_t<bit_count_f(w())>{};
}


template <int N_zero=-1> requires (N_zero == -1)
XTAL_DEF_(return,inline,let)
bit_floor_f(cardinal_variable_q auto u)
noexcept -> int
{
	/*/
	U const  z =      one == u;
	U const _z = bond::fit<decltype(u)>::unit.depth&-z;
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
	int  x = bit_sign_f<int>(v);
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
	auto constexpr N = -static_cast<int>(bond::fit<decltype(x)>::unit.mark);
	if constexpr  (N == N_zero) {
		return  bit_exponent_f(XTAL_REF_(x));
	}
	else {
		int n = bit_exponent_f(XTAL_REF_(x));
		if constexpr (N_zero != 0) {n -=         N_zero;}
		if constexpr (N_zero != N) {n &= ~bit_sign_f(n);}
		if constexpr (N_zero != 0) {n +=         N_zero;}
		return n;
	}
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
requires real_variable_q<unstruct_u<decltype(x)>>
{
	auto constexpr N = -static_cast<int>(bond::fit<decltype(x)>::unit.mark);
	return bit_floor_f<N>(XTAL_REF_(x));
}
template <int N_zero=-1>
XTAL_DEF_(return,inline,let)
bit_floor_f(constant_q auto w)
noexcept -> auto
{
	return constant_t<bit_floor_f<N_zero>(w())>{};
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
	int  x = bit_sign_f<int>(v);
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
	using X_fit = bond::fit<decltype(x)>;
	return bit_floor_f<N_zero>(X_fit::diplo_1*X_fit::dnsilon_f(1)*XTAL_REF_(x));
}
template <int N_zero=+1>
XTAL_DEF_(return,inline,let)
bit_ceiling_f(complex_variable_q auto &&x)
noexcept -> int
{
	auto  y = bit_ceiling_f<N_zero>(norm(XTAL_REF_(x)));
	y  += y&1;
	y >>=   1;
	return  y;
	//\
	return bit_floor_f<N_zero>(XTAL_REF_(x)) + 1;
}
template <int N_zero=-1>
XTAL_DEF_(return,inline,let)
bit_ceiling_f(constant_q auto w)
noexcept -> auto
{
	return constant_t<bit_ceiling_f<N_zero>(w())>{};
}


/*!
\returns The bitwise-reversal of `u`,
         restricted to `N_subdepth` when `0 < N_subdepth < sizeof(u) << 3U`.
\note    Requires `log2(sizeof(u) << 3U)` iterations.
*/
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
	return constant_t<bit_reverse_f(w(), n_subdepth)>{};
}
template <int N_subdepth=0>
XTAL_DEF_(return,inline,let)
bit_reverse_f(auto &&x)
noexcept -> auto
{
	auto constexpr N_depth = bond::fit<decltype(x)>::full.depth;
	//\
	int constexpr n_subdepth = below_v<N_depth, (unsigned) N_subdepth>;
	int constexpr n_subdepth = 0 < N_subdepth? N_subdepth: N_depth;
	return bit_reverse_f(XTAL_REF_(x), n_subdepth);
}


///////////////////////////////////////////////////////////////////////////////

XTAL_DEF_(return,let)
bit_representation_f(real_variable_q auto x)
noexcept -> auto
{
	using X_fit =  bond::fit<decltype(x)>;
	using U = typename X_fit::sigma_type;
	using V = typename X_fit::delta_type;

	U constexpr N = X_fit::unit.mark + X_fit::fraction.depth;
	U constexpr M =             one << X_fit::fraction.depth;
	
	auto const o = _xtd::bit_cast<U>(x);
	V const z = static_cast<V>(o) >> X_fit::positive.depth;
	V const n = N - (o << X_fit::sign.depth >> X_fit::sign.depth + X_fit::exponent.shift);
	V       m = M | (o&X_fit::fraction.mask);
	m  ^= z;
	m  -= z;
	return atom::couple_t<V[2]>{m, n};
}
XTAL_DEF_(return,inline,let)
bit_presentation_f(atom::couple_q auto const &mn)
noexcept -> auto
{
	auto const &[m, n] = mn;
	using MN_op = bond::fit<decltype(m), decltype(n)>;
	return static_cast<typename MN_op::alpha_type>(m)*MN_op::haplo_f(n);
}


////////////////////////////////////////////////////////////////////////////////
/*!
\returns The fractional portion of `x` as a full-width `delta_type`.
*/
template <class T_return=void>
XTAL_DEF_(return,inline,let)
bit_fraction_f()
noexcept -> auto
{
	using Y     = T_return;
	using Y_fit = bond::fit<Y>;

	XTAL_IF0
	XTAL_0IF (integral_q<Y>) {return Y_fit::diplo_f(Y_fit::full.depth);}
	XTAL_0IF (    real_q<Y>) {return Y_fit::haplo_f(Y_fit::full.depth);}
	XTAL_0IF_(void)
}
template <class Y_return=void>
XTAL_DEF_(return,inline,let)
bit_fraction_f(integral_variable_q auto x)
noexcept -> auto
{
	using X       =              XTAL_ALL_(x);
	using X_fit   =          bond::fit<X>;
	using X_alpha = typename X_fit::alpha_type;
	using X_sigma = typename X_fit::sigma_type;
	using X_delta = typename X_fit::delta_type;

	using Y       = complete_t<Y_return, X_alpha>;
	using Y_fit    = bond::fit<Y>;
	using Y_alpha = typename Y_fit::alpha_type;
	using Y_sigma = typename Y_fit::sigma_type;
	using Y_delta = typename Y_fit::delta_type;

	XTAL_IF0
	XTAL_0IF (cardinal_q<Y>) {
		return Y_fit::sigma_f(x);
	}
	XTAL_0IF ( ordinal_q<Y>) {
		return Y_fit::delta_f(x);
	}
	XTAL_0IF (    real_q<Y>) {
		//\
		return Y_fit::alpha_f(x)*Y_fit::haplo_f(X_fit::full.depth);
		return Y_fit::alpha_f(static_cast<X_delta>(x))*Y_fit::haplo_f(X_fit::full.depth);
	}
	XTAL_0IF_(void)
}
template <class Y_return=void>
XTAL_DEF_(return,inline,let)
bit_fraction_f(real_variable_q auto x)
noexcept -> auto
{
	using X       =              XTAL_ALL_(x);
	using X_fit   =          bond::fit<X>;
	using X_alpha = typename X_fit::alpha_type;
	using X_sigma = typename X_fit::sigma_type;
	using X_delta = typename X_fit::delta_type;

	using Y       = complete_t<Y_return, X>;
	using Y_fit   = bond::fit<Y>;
	using Y_alpha = typename Y_fit::alpha_type;
	using Y_sigma = typename Y_fit::sigma_type;
	using Y_delta = typename Y_fit::delta_type;

	XTAL_IF0
	XTAL_0IF_(consteval) {
		XTAL_IF0
		XTAL_0IF (     real_q<Y>) {
			return static_cast<Y>(bit_fraction_f<Y_delta>(x))*bit_fraction_f<Y_alpha>();
		}
		XTAL_0IF (  ordinal_q<Y>) {
			auto constexpr N_exp = X_fit::exponent.shift;
			auto constexpr M_exp = X_fit::unit.mark + X_fit::unit.shift - Y_fit::full.depth;
			auto constexpr M_sgn = X_fit::sign.mask;

			X_delta o = _xtd::bit_cast<X_delta>(x);
			X_delta v = o &  M_sgn; o ^= v; v >>= X_fit::sign.shift;
			X_sigma x = o >> N_exp; x -= M_exp;
			X_delta u = bit_flag_f<~0>(x);

			o &=      X_fit::fraction. mask;
			o |= u << X_fit::fraction.depth;
			o  = bit_shift_f(XTAL_MOV_(o), _xtd::bit_cast<X_delta>(XTAL_MOV_(x)));
		//	o |= 1;
			o ^= v;
			o -= v;
			return static_cast<Y>(XTAL_MOV_(o));
		}
		XTAL_0IF ( cardinal_q<Y>) {
			return _xtd::bit_cast<Y>(bit_fraction_f<Y_delta>(x));
		}
		XTAL_0IF_(void)
	}
	XTAL_0IF_(else) {
		x -= round(x);
		if constexpr (integral_q<Y>) {
			x *= bit_fraction_f<Y_delta>();
		}
		if constexpr (cardinal_q<Y>) {
			return _xtd::bit_cast<Y>(static_cast<Y_delta>(x));
		}
		else {
			return static_cast<Y>(x);
		}
	}
}
template <class T_return=void>
XTAL_DEF_(return,inline,let)
bit_fraction_f(complex_variable_q auto const &x)
noexcept -> auto
requires requires {bit_fraction_f(x.real()); bit_fraction_f(x.imag());}
{
	using X       = XTAL_ALL_(x);
	using X_fit   = bond::fit<X>;
	using X_alpha = typename X_fit::alpha_type;
	using X_sigma = typename X_fit::sigma_type;
	using X_delta = typename X_fit::delta_type;

	using U = _std::conditional_t<integral_q<typename X::value_type>, X_alpha, X_delta>;
	using Y = complete_t<T_return, _std::complex<U>>;
	using V = typename Y::value_type;

	return Y{bit_fraction_f<V>(x.real()), bit_fraction_f<V>(x.imag())};
}

////////////////////////////////////////////////////////////////////////////////
/*!
\returns The `target` to `N_zoom` bits of precision after the decimal.
*/
template <int N_zoom=0>
XTAL_DEF_(return,verbatim,set)
bit_trim_f(integral_variable_q auto x)
noexcept -> XTAL_ALL_(x)
{
	return x;
}
#if     XTAL_VER_(MSVC)
#pragma optimize("", off)
#endif
template <int N_zoom=0>
XTAL_DEF_(return,verbatim,set)
bit_trim_f(real_variable_q auto x)
noexcept -> XTAL_ALL_(x)
{
	using X       = XTAL_ALL_(x);
	using X_fit   = bond::fit<X>;
	using X_alpha = typename X_fit::alpha_type;
	using X_sigma = typename X_fit::sigma_type;
	using X_delta = typename X_fit::delta_type;
	using X_lim   = _std::numeric_limits<X_alpha>;

	X_delta constexpr N_unzoom = 0 < N_zoom? N_zoom - X_fit::fraction.depth: N_zoom - 1;
	X_alpha constexpr N_minima = X_lim::min()*X_fit::diplo_f(N_unzoom);
	x *= N_minima;
	x /= N_minima;
	return x;
}
#if     XTAL_VER_(MSVC)
#pragma optimize("",  on)
#endif

template <int N_zoom=0>
XTAL_DEF_(return,verbatim,set)
bit_trim_f(complex_variable_q auto x)
noexcept -> XTAL_ALL_(x)
{
	return {bit_trim_f<N_zoom>(x.real()), bit_trim_f<N_zoom>(x.imag())};
}


////////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////
}/////////////////////////////////////////////////////////////////////////////
XTAL_ENV_(pop)

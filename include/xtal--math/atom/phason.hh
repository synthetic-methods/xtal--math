#pragma once
#include "./any.hh"

#include "../bond/all.hxx"




XTAL_ENV_(push)
namespace xtal::atom::math
{/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////
/*!
\brief   Extends `grade` as a fixed-point fractional/cyclic value with a bisected representation.

Allows floating-point construction via `std::initializer_list`,
and access to the floating-point value via `operator()`/`operator(int)`.

Implements truncated floating-point multiplication (affecting all elements)
and addition (affecting only the initial element).

\todo    Rework `operator`s to accommodate `std::complex`.
*/
template <class   ..._s>	struct          phason;
template <class   ..._s>	using           phason_t = typename phason<_s...>::type;
template <class   ...Ts>	concept         phason_q = bond::tag_infixed_p<phason_t, Ts...>;
template <class   ...Ts>	concept    real_phason_q = bond::tag_infixed_p<phason_t, Ts...> and real_variable_q<initializer_t<Ts>...>;
template <class   ...Ts>	concept simplex_phason_q = bond::tag_infixed_p<phason_t, Ts...> and simplex_field_q<initializer_t<Ts>...>;
template <class   ...Ts>	concept complex_phason_q = bond::tag_infixed_p<phason_t, Ts...> and complex_field_q<initializer_t<Ts>...>;

XTAL_DEF_(let) phason_f = [] XTAL_1FN_(call) (_detail::fake_f<phason_t>);


////////////////////////////////////////////////////////////////////////////////

template <scalar_q ..._s> requires common_q<_s...>
struct phason<_s ...>
:	phason<common_t<_s...>[sizeof...(_s)]>
{
};
template <vector_q A> requires integral_variable_q<unstruct_u<A>>
struct phason<A>
:	grade<A>
{
};
template <vector_q A, scalar_q ..._s>
struct phason<A, _s...>
:	couple_t<phason_t<A>, _s...>
{
};
template <vector_q A> requires     real_variable_q<unstruct_u<A>>
struct phason<A>
{
private:
	static auto constexpr M_data = _xtd::extent_v<based_t<A>>;

	using W = array_valued_u<A>;
	static_assert(continuous_field_q<W>);

	using W_fit   = bond::fit<W>;
	using U_fit   = typename W_fit::template widen<-1>;
	using W_sigma = typename W_fit::sigma_type;
	using U_sigma = typename U_fit::sigma_type;
	using W_delta = typename W_fit::delta_type;
	using U_delta = typename U_fit::delta_type;

	using V = bond::compose_s<U_delta, W_fit>;
	using U = bond::compose_s<U_sigma, W_fit>;

	using coordinate_type = W;
	using inordinate_type = V;
	using   ordinate_type = U;

	static_assert(_std::numeric_limits<unstruct_u<U>>::is_modulo);// D'oh!

	template <class T>
	using holotype = bond::compose_s<typename grade<U[M_data]>::template homotype<T>, bond::tag<phason_t>>;

public:
	template <class T>
	class homotype : public holotype<T>
	{
		using S_ = holotype<T>;

	public:// TYPE
		using typename S_::value_type;
		using typename S_:: size_type;

		using coordinate_type = W;
		using inordinate_type = V;
		using   ordinate_type = U;

	public:// ACCESS
		using S_::size;
		using S_::self;
		using S_::twin;

		static auto constexpr   ordinate = [] XTAL_1FN_(call) (bond::math::bit_fraction_f<  U>);
		static auto constexpr coordinate = [] XTAL_1FN_(call) (bond::math::bit_fraction_f<W>);

	public:// CONSTRUCT
	//	using S_::S_;
	~	homotype()                 noexcept=default;
		homotype()                 noexcept=default;
		XTAL_NEW_(copy) (homotype, noexcept=default)
		XTAL_NEW_(move) (homotype, noexcept=default)

		XTAL_NEW_(explicit)
		homotype(iterable_q auto &&o)
		noexcept
		:	S_(count_f(o))
		{
			operator>>=(XTAL_REF_(o));
		}
		XTAL_NEW_(implicit)
		homotype(_std::initializer_list<W> o)
		noexcept
		:	S_(count_f(o))
		{
			operator>>=(XTAL_MOV_(o));
		}
		
	public:// RECONSTRUCT
	//	using S_::operator >>=;
	//	using S_::operator <<=;

		XTAL_DEF_(mutate,inline,let)
		operator >>=(_std::initializer_list<W> o)
		noexcept -> auto &
		{
			_detail::move_to<T::ordinate>(S_::data(), XTAL_REF_(o));
			return self();
		}
		XTAL_DEF_(mutate,inline,let)
		operator >>=(iterable_q auto       &&o)
		noexcept -> auto &
		{
			_detail::move_to<T::ordinate>(S_::data(), XTAL_REF_(o));
			return self();
		}
		XTAL_DEF_(mutate,inline,let)
		operator >>=(iterable_q auto const  &o)
		noexcept -> auto &
		{
			_detail::copy_to<T::ordinate>(S_::data(), XTAL_REF_(o));
			return self();
		}
		
		XTAL_DEF_(mutate,inline,let)
		operator <<=(_std::initializer_list<W> o)
		noexcept -> auto &
		{
			auto i0 = S_::data(), iN = _std::next(i0, S_::size() - o.size());
			_detail::move_to<T::ordinate>(iN, XTAL_REF_(o));
			return self();
		}
		XTAL_DEF_(mutate,inline,let)
		operator <<=(iterable_q auto &&o)
		noexcept -> auto &
		{
			auto i0 = S_::data(), iN = _std::next(i0, S_::size() - o.size());
			_detail::move_to<T::ordinate>(iN, XTAL_REF_(o));
			return self();
		}
		XTAL_DEF_(mutate,inline,let)
		operator <<=(iterable_q auto const &o)
		noexcept -> auto &
		{
			auto i0 = S_::data(), iN = _std::next(i0, S_::size() - o.size());
			_detail::copy_to<T::ordinate>(iN, XTAL_REF_(o));
			return self();
		}

	public:// OPERATE
		/*!
		\brief   Scales all elements.
		*/
		XTAL_DEF_(mutate,inline,let)
		operator /= (auto &&x)
		noexcept -> auto &
		requires un_n<requires (W u) {u /= x;}>
		and      XTAL_TRY_(to) (S_::operator/=(XTAL_REF_(x)))

		XTAL_DEF_(mutate,inline,let)
		operator *= (auto &&x)
		noexcept -> auto &
		requires un_n<requires (W u) {u *= x;}>
		and      XTAL_TRY_(to) (S_::operator*=(XTAL_REF_(x)))

		XTAL_DEF_(mutate,inline,let)
		operator /= (auto &&x)
		noexcept -> auto &
		requires in_n<requires (W u) {u *= x;}>
		{
			using X     = XTAL_ALL_(x);
			using X_fit = bond::fit<X>;
			return operator*=(X_fit::alpha_1/XTAL_REF_(x));
		}
		XTAL_DEF_(mutate,inline,let)
		operator *= (auto &&x)
		noexcept -> auto &
		requires in_n<requires (W u) {u *= x;}>
		{
			using X     = XTAL_ALL_(x);
			using X_fit = bond::fit<X>;
			auto &s = reinterpret_cast<phason_t<V[size]> &>(self());

			auto constexpr N0 = U_fit::full.depth - 1;
			auto constexpr N1 = U_fit::full.depth;
			auto constexpr N2 = U_fit::half.depth;
			auto constexpr M2 = U_fit::half.width;

		//	TODO: Adapt for `std::complex`?
			XTAL_IF0
			XTAL_0IF (integral_variable_q<X>) {
				S_::operator*=(XTAL_REF_(x));
			}
			XTAL_0IF (1*sizeof(U) == sizeof(W)) {
				unsigned constexpr M_bias = N2 >> M2;
				unsigned constexpr M_size = N2  - M_bias;
				auto [m, n] = bond::math::bit_representation_f(XTAL_REF_(x));
				m >>= n - M_size;
				s >>=     M_size;
				s  *= m;
			}
			XTAL_IF0
			XTAL_0IF (2*sizeof(U) == sizeof(W)) {
				W_sigma const u = W_fit::sigma_f(X_fit::diplo_f(N1)*XTAL_REF_(x));
				U_sigma t_[2];
				#pragma unroll
				for (size_type i{}; i < size; ++i) {
					reinterpret_cast<W_sigma &>(t_) = u*s[i];
				//	t_[1] += t_[0] >> N0;// Rounding...
					s [i]  = t_[1];
				}
			}
			return self();
		}
		XTAL_DEF_(return,inline,met)
		operator * (simplex_variable_q auto const &x, T const &t)
		noexcept -> auto {
			return t.twin() *= x;
		}

		/*!
		\brief   Offsets the first element.
		*/		
	//	XTAL_DEF_(mutate,inline,let) operator +=(_std::initializer_list<W> o) noexcept -> auto & {return S_::operator+=(T(o));}
	//	XTAL_DEF_(mutate,inline,let) operator -=(_std::initializer_list<W> o) noexcept -> auto & {return S_::operator-=(T(o));}

		XTAL_DEF_(mutate,inline,get) operator -= (auto &&x) noexcept {return S_::operator-=(XTAL_REF_(x));}
		XTAL_DEF_(mutate,inline,get) operator += (auto &&x) noexcept {return S_::operator+=(XTAL_REF_(x));}

		XTAL_DEF_(mutate,inline,let)
		operator -= (auto &&x)
		noexcept -> auto &
		requires additive_group_p<1, W, decltype(x)>
		{
			if constexpr (un_n<integral_q<unstruct_u<decltype(x)>>>) {
				get<0>(self()) -= bond::math::bit_fraction_f<V>(x);
			}
			return self();
		}
		XTAL_DEF_(mutate,inline,let)
		operator += (auto &&x)
		noexcept -> auto &
		requires additive_group_p<1, W, decltype(x)>
		{
			if constexpr (un_n<integral_q<unstruct_u<decltype(x)>>>) {
				get<0>(self()) += bond::math::bit_fraction_f<V>(x);
			}
			return self();
		}

		XTAL_DEF_(inline,let)
		scale(W u, W w=one)
		noexcept -> auto &
		{
			bond::seek_out_f<size - 1, 1>([&, this]<constant_q I> (I)
				XTAL_0FN_(do) (get<I{}>(self()) = T::ordinate(got<I{}>(self())*(w *= u))));
			return self();
		}
		XTAL_DEF_(inline,let)
		scaled(W u, W w=one) const
		noexcept -> auto
		{
			using T_ = XTAL_ALL_(twin());
			return [&, this]<auto ...I> (bond::seek_t<I...>)
				XTAL_0FN_(to) (T_{_detail::thunk_f(got<I>(self())*(w)) (w *= u)...})
			(bond::seek_s<size>{});
		}

		/*!
		\returns `condition_f<Y>` indicating whether the current state is continuous.
		*/
		template <class Y=W>
		XTAL_DEF_(return,inline,let)
		continuity()
		noexcept -> Y
		{
			return condition_f<Y>(not discontinuity<bool>());
		}
		/*!
		\returns `condition_f<Y>` indicating whether the current state is discontinuous.
		*/
		template <class Y=W>
		XTAL_DEF_(return,inline,let)
		discontinuity()
		noexcept -> Y
		{
		//	TODO: Accommodate `discontinuity` for `complex_variable_q<V>`...
			auto constexpr N1 = U_fit::positive.depth;
			auto [u0, u1] =  self();
			auto const v0 = _xtd::bit_cast<V>(u0) >> N1; u0 ^= v0; u0 -= v0;
			auto const v1 = _xtd::bit_cast<V>(u1) >> N1; u1 ^= v1; u1 -= v1;
			return condition_f<Y>(v0 == v1 and u0 < u1);
		}

	};
	using type = bond::derive_t<homotype>;

};


///////////////////////////////////////////////////////////////////////////////
}/////////////////////////////////////////////////////////////////////////////
XTAL_ENV_(pop)

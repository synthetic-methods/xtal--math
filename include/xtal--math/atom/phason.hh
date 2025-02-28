#pragma once
#include "./any.hh"

#include "../bond/all.hxx"




XTAL_ENV_(push)
namespace xtal::atom::math
{/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////
///\
Extends `grade` as a fixed-point fractional/cyclic value with a bisected representation. \
\
Allows floating-point construction via `std::initializer_list`, \
and access to the floating-point value via `operator()`/`operator(int)`. \
\
Implements truncated floating-point multiplication (affecting all elements) \
and addition (affecting only the initial element). \

///\todo\
Might be worth implementing multiplication for the first element only, \
either by truncating the type or by parameterizing the operator. \
The latter could be achieved using `std::initializer_list`s, \
parameterized by an `continuous_field_q`-wrapper with a distinguished head. \

///\todo\
Rework `operator`s to accommodate `std::complex`. \

template <class   ..._s>	struct          phason;
template <class   ..._s>	using           phason_t = typename phason<_s...>::type;
template <class   ...Ts>	concept         phason_q = bond::fixed_tagged_with_p<phason_t, Ts...>;
template <class   ...Ts>	concept    real_phason_q = bond::fixed_tagged_with_p<phason_t, Ts...> and real_variable_q<initializer_t<Ts>...>;
template <class   ...Ts>	concept simplex_phason_q = bond::fixed_tagged_with_p<phason_t, Ts...> and simplex_field_q<initializer_t<Ts>...>;
template <class   ...Ts>	concept complex_phason_q = bond::fixed_tagged_with_p<phason_t, Ts...> and complex_field_q<initializer_t<Ts>...>;

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
template <vector_q A> requires     real_variable_q<unstruct_u<A>>
struct phason<A>
{
	static auto constexpr M_data = _std::extent_v<based_t<A>>;

	using coordinate_type = array_valued_u<A>;
	static_assert(continuous_field_q<coordinate_type>);

	using W_fit   = bond::fit<coordinate_type>;
	using U_fit   = typename W_fit::template widen<-1>;
	using W_sigma = typename W_fit::sigma_type;
	using U_sigma = typename U_fit::sigma_type;
	using W_delta = typename W_fit::delta_type;
	using U_delta = typename U_fit::delta_type;

	using   ordinate_type = bond::compose_s<U_sigma, W_fit>;
	using inordinate_type = bond::compose_s<U_delta, W_fit>;

	static_assert(_std::numeric_limits<unstruct_u<ordinate_type>>::is_modulo);// D'oh!

	template <class T>
	using holotype = bond::compose_s<typename grade<ordinate_type[M_data]>::template homotype<T>, bond::tag<phason_t>>;

	template <class T>
	class homotype : public holotype<T>
	{
		using S_ = holotype<T>;

	public:// TYPE
		using typename S_::value_type;

	public:// ACCESS
		using S_::size;
		using S_::self;
		using S_::twin;

		static auto constexpr   ordinate = [] XTAL_1FN_(call) (bond::math::bit_fraction_f<  ordinate_type>);
		static auto constexpr coordinate = [] XTAL_1FN_(call) (bond::math::bit_fraction_f<coordinate_type>);

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
		homotype(_std::initializer_list<coordinate_type> o)
		noexcept
		:	S_(count_f(o))
		{
			operator>>=(XTAL_MOV_(o));
		}
		
	public:// RECONSTRUCT
		using S_::operator >>=;
		using S_::operator <<=;

		XTAL_DEF_(mutate,inline,let)
		operator >>=(_std::initializer_list<coordinate_type> o)
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
		operator <<=(_std::initializer_list<coordinate_type> o)
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
		///\
		Scales all elements. \

		///\note\
		The symmetric signatures for `/=` and `*=` are declared-but-undefined \
		to avoid compilation-failure when type-checking e.g. `multiplicative_group_q`. \

		XTAL_DEF_(mutate,inline,let)
		operator /= (auto &&x)
		noexcept -> auto &
		requires XTAL_TRY_(to) (operator*=(bond::fit<decltype(x)>::alpha_1/XTAL_REF_(x)))

		XTAL_DEF_(mutate,inline,let)
		operator *= (auto &&x)
		noexcept -> auto &
		requires un_n<simplex_variable_q<decltype(x)>>
		and      XTAL_TRY_(to) (S_::operator*=(XTAL_REF_(x)))

		XTAL_DEF_(mutate,inline,let)
		operator *= (auto &&x)
		noexcept -> auto &
		requires in_n<simplex_variable_q<decltype(x)>>
		{
			using X_fit = bond::fit<decltype(x)>;
			auto &s = reinterpret_cast<phason_t<inordinate_type[size]> &>(self());

			XTAL_IF0
			XTAL_0IF (integral_variable_q<decltype(x)>) {
				S_::operator*=(x);
			}
			XTAL_0IF (1*sizeof(ordinate_type) == sizeof(coordinate_type)) {
				unsigned constexpr M_bias = U_fit::half.depth >> U_fit::half.width;
				unsigned constexpr M_size = U_fit::half.depth - M_bias;
				auto [m, n] = bond::math::bit_representation_f(x);
				m >>= n - M_size;
				s >>=     M_size;
				s  *= m;
			}
			XTAL_IF0
			XTAL_0IF (2*sizeof(ordinate_type) == sizeof(coordinate_type)) {
				U_sigma t_[2];
				W_sigma const u(x*X_fit::diplo_f(U_fit::full.depth));
				#pragma unroll
				for (XTAL_ALL_(size()) i{}; i < size; ++i) {
					reinterpret_cast<W_sigma &>(t_) = u*s[i];
				//	t_[1] += t_[0] >> U_fit::positive.depth;// Round...
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

		///\
		Offsets the first element. \
		
		XTAL_DEF_(mutate,inline,let)
		operator -= (auto &&x)
		noexcept -> auto &
		requires XTAL_TRY_(to) (operator+=(-XTAL_REF_(x)))

		XTAL_DEF_(mutate,inline,let)
		operator += (auto &&x)
		noexcept -> auto &
		requires un_n<simplex_variable_q<decltype(x)>>
		and      XTAL_TRY_(to) (S_::operator+=( XTAL_REF_(x)))

		XTAL_DEF_(mutate,inline,let)
		operator += (auto &&x)
		noexcept -> auto &
		requires in_n<simplex_variable_q<decltype(x)>>
		{
			if constexpr (real_variable_q<decltype(x)>) {
				get<0>(self()) += bond::math::bit_fraction_f<_std::make_signed_t<ordinate_type>>(x);
			}
			return self();
		}

		XTAL_DEF_(inline,let)
		scale(coordinate_type u, coordinate_type w=one)
		noexcept -> auto &
		{
			bond::seek_out_f<size - 1, 1>([&, this]<constant_q I> (I)
				XTAL_0FN_(do) (get<I{}>(self()) = T::ordinate(got<I{}>(self())*(w *= u))));
			return self();
		}
		XTAL_DEF_(inline,let)
		scaled(coordinate_type u, coordinate_type w=one) const
		noexcept -> auto
		{
			using T_ = XTAL_ALL_(twin());
			return [&, this]<auto ...I> (bond::seek_t<I...>)
				XTAL_0FN_(to) (T_{_detail::thunk_f(got<I>(self())*(w)) (w *= u)...})
			(bond::seek_s<size>{});
		}

		///\returns `condition_f<Y>` indicating whether the current state is continuous. \

		template <class Y=coordinate_type>
		XTAL_DEF_(return,inline,let)
		continuity()
		noexcept -> Y
		{
			return condition_f<Y>(not discontinuity<bool>());
		}
		///\returns `condition_f<Y>` indicating whether the current state is discontinuous. \

		template <class Y=coordinate_type>
		XTAL_DEF_(return,inline,let)
		discontinuity()
		noexcept -> Y
		{
		//	TODO: Accommodate `discontinuity` for `complex_variable_q<inordinate_type>`...
			auto constexpr N1 = U_fit::positive.depth;
			auto [u0, u1] =  self();
			auto const v0 = _xtd::bit_cast<inordinate_type>(u0) >> N1; u0 ^= v0; u0 -= v0;
			auto const v1 = _xtd::bit_cast<inordinate_type>(u1) >> N1; u1 ^= v1; u1 -= v1;
			return condition_f<Y>(v0 == v1 and u0 < u1);
		}

	};
	using type = bond::derive_t<homotype>;

};


///////////////////////////////////////////////////////////////////////////////
}/////////////////////////////////////////////////////////////////////////////
XTAL_ENV_(pop)

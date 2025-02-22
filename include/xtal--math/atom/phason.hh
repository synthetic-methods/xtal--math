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
template <class   ...Ts>	concept         phason_q = bond::tag_p<phason_t, Ts...>;
template <class   ...Ts>	concept    real_phason_q = bond::tag_p<phason_t, Ts...> and real_variable_q<initializer_t<Ts>...>;
template <class   ...Ts>	concept simplex_phason_q = bond::tag_p<phason_t, Ts...> and simplex_field_q<initializer_t<Ts>...>;
template <class   ...Ts>	concept complex_phason_q = bond::tag_p<phason_t, Ts...> and complex_field_q<initializer_t<Ts>...>;


XTAL_FX0_(to) (template <auto f=_std::identity{}>
XTAL_DEF_(return,inline,let)
phason_f(auto &&...oo),
	_detail::factory<phason_t>::
		template make<f>(XTAL_REF_(oo)...))


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

	using U_fit = bond::fit<coordinate_type>;
	using T_fit = typename U_fit::template widen<-1>;
	
	using   ordinate_type = bond::compose_s<typename T_fit::sigma_type, U_fit>;
	using inordinate_type = bond::compose_s<typename T_fit::delta_type, U_fit>;

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
			using _fit = bond::fit<decltype(x)>;
			auto &s = reinterpret_cast<phason_t<inordinate_type[size]> &>(self());

			XTAL_IF0
			XTAL_0IF (integral_variable_q<decltype(x)>) {
				S_::operator*=(x);
			}
			XTAL_0IF (1*sizeof(ordinate_type) == sizeof(coordinate_type)) {
				unsigned constexpr M_bias = T_fit::half.depth >> T_fit::half.width;
				unsigned constexpr M_size = T_fit::half.depth - M_bias;
				auto [m, n] = bond::math::bit_representation_f(x);
				m >>= n - M_size;
				s >>=     M_size;
				s  *= m;
			}
			XTAL_IF0
			XTAL_0IF (2*sizeof(ordinate_type) == sizeof(coordinate_type)) {
				typename T_fit::sigma_type t_[2];
				typename U_fit::sigma_type const u(x*_fit::diplo_f(T_fit::full.depth));
				#pragma unroll
				for (XTAL_ALL_(size()) i{}; i < size; ++i) {
					reinterpret_cast<typename U_fit::sigma_type &>(t_) = u*s[i];
				//	t_[1] += t_[0] >> T_fit::positive.depth;// Round...
					s [i]  = t_[1];
				}
			}
			return self();
		}
		XTAL_DEF_(return,inline,met) operator * (simplex_variable_q auto const &x, T const &t) noexcept -> auto {return t.twin() *= x;}
	//	XTAL_DEF_(return,inline,met) operator * (T const &t, simplex_variable_q auto const &x) noexcept -> auto {return t.twin() *= x;}

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

	//	XTAL_DEF_(return,inline,met) operator + (T const &t, auto const &o) noexcept -> auto                                        {return t.twin() += o;}
	//	XTAL_DEF_(return,inline,met) operator + (auto const &o, T const &t) noexcept -> auto   requires un_n<phason_q<decltype(o)>> {return t.twin() += o;}

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
			auto [u0, u1] =  self();
			auto const v0 = _xtd::bit_cast<inordinate_type>(u0) >> T_fit::positive.depth; u0 ^= v0; u0 -= v0;
			auto const v1 = _xtd::bit_cast<inordinate_type>(u1) >> T_fit::positive.depth; u1 ^= v1; u1 -= v1;
			return condition_f<Y>(v0 == v1 and u0 < u1);
		}

	};
	using type = bond::derive_t<homotype>;

};


///////////////////////////////////////////////////////////////////////////////
}/////////////////////////////////////////////////////////////////////////////
XTAL_ENV_(pop)

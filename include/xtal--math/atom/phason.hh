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

template <class   ..._s>	struct           phason;
template <class   ..._s>	using            phason_t = typename phason<_s...>::type;
template <class   ...Ts>	concept          phason_q = bond::tag_p<phason_t, Ts...>;
template <class   ...Ts>	concept     real_phason_q = bond::tag_p<phason_t, Ts...> and real_variable_q<initializer_u<Ts>...>;
template <class   ...Ts>	concept  simplex_phason_q = bond::tag_p<phason_t, Ts...> and simplex_field_q<initializer_u<Ts>...>;
template <class   ...Ts>	concept  complex_phason_q = bond::tag_p<phason_t, Ts...> and complex_field_q<initializer_u<Ts>...>;


XTAL_FX0_(alias) (template <auto f=_std::identity{}>
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
template <vector_q A> requires integral_variable_q<absolve_u<A>>
struct phason<A>
:	grade<A>
{
};
template <vector_q A> requires     real_variable_q<absolve_u<A>>
struct phason<A>
{
	static auto constexpr M_data = _std::extent_v<based_t<A>>;

	using coordinate_type = array_valued_u<A>;
	static_assert(continuous_field_q<coordinate_type>);

	using U_    = bond::  forge<coordinate_type>;
	using U_fix = bond::fixture<coordinate_type>;
	using T_fix = typename U_fix::template widen<-1>;
	
	using   ordinate_type = bond::compose_s<typename T_fix::sigma_type, U_>;
	using inordinate_type = bond::compose_s<typename T_fix::delta_type, U_>;

	static_assert(_std::numeric_limits<absolve_u<ordinate_type>>::is_modulo);// D'oh!

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

		static auto constexpr   ordinate = [] XTAL_0FN_(alias) (bond::math::bit_fraction_f<  ordinate_type>);
		static auto constexpr coordinate = [] XTAL_0FN_(alias) (bond::math::bit_fraction_f<coordinate_type>);

	public:// CONSTRUCT
	//	using S_::S_;
	~	homotype()                 noexcept=default;
		homotype()                 noexcept=default;
		XTAL_NEW_(copy) (homotype, noexcept=default)
		XTAL_NEW_(move) (homotype, noexcept=default)

		XTAL_NEW_(explicit) homotype(iterable_q auto &&o)
		noexcept
		:	S_(count_f(o))
		{
			operator>>=(XTAL_REF_(o));
		}
		XTAL_NEW_(implicit) homotype(_std::initializer_list<coordinate_type> o)
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
		requires XTAL_TRY_(return) (operator*=(bond::fixture<decltype(x)>::alpha_1/XTAL_REF_(x)))

		XTAL_DEF_(mutate,inline,let)
		operator *= (auto &&x)
		noexcept -> auto &
		requires un_n<simplex_variable_q<decltype(x)>>
		and      XTAL_TRY_(return) (S_::operator*=(XTAL_REF_(x)))

		XTAL_DEF_(mutate,inline,let)
		operator *= (auto &&x)
		noexcept -> auto &
		requires in_n<simplex_variable_q<decltype(x)>>
		{
			using _fix = bond::fixture<decltype(x)>;
			auto &s = reinterpret_cast<phason_t<inordinate_type[size]> &>(self());

			XTAL_IF0
			XTAL_0IF (integral_variable_q<decltype(x)>) {
				S_::operator*=(x);
			}
			XTAL_0IF (1*sizeof(ordinate_type) == sizeof(coordinate_type)) {
				unsigned constexpr M_bias = T_fix::half.depth >> T_fix::half.width;
				unsigned constexpr M_size = T_fix::half.depth - M_bias;
				auto [m, n] = bond::math::bit_representation_f(x);
				m >>= n - M_size;
				s >>=     M_size;
				s  *= m;
			}
			XTAL_IF0
			XTAL_0IF (2*sizeof(ordinate_type) == sizeof(coordinate_type)) {
				typename T_fix::sigma_type t_[2];
				typename U_fix::sigma_type const u(x*_fix::diplo_f(T_fix::full.depth));
				#pragma unroll
				for (XTAL_ALL_(size()) i{}; i < size; ++i) {
					reinterpret_cast<typename U_fix::sigma_type &>(t_) = u*s[i];
				//	t_[1] += t_[0] >> T_fix::positive.depth;// Round...
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
		requires XTAL_TRY_(return) (operator+=(-XTAL_REF_(x)))

		XTAL_DEF_(mutate,inline,let)
		operator += (auto &&x)
		noexcept -> auto &
		requires un_n<simplex_variable_q<decltype(x)>>
		and      XTAL_TRY_(return) (S_::operator+=( XTAL_REF_(x)))

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
			auto const v0 = _xtd::bit_cast<inordinate_type>(u0) >> T_fix::positive.depth; u0 ^= v0; u0 -= v0;
			auto const v1 = _xtd::bit_cast<inordinate_type>(u1) >> T_fix::positive.depth; u1 ^= v1; u1 -= v1;
			return condition_f<Y>(v0 == v1 and u0 < u1);
		}

	};
	using type = derive_t<homotype>;

};
static_assert(atomic_q<phason_t<float[2]>>);

static_assert(bond::pack_size_q<phason_t<double[2]>>);


///////////////////////////////////////////////////////////////////////////////
}/////////////////////////////////////////////////////////////////////////////
XTAL_ENV_(pop)

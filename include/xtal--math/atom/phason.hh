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
template <class  V=void>
XTAL_DEF_(return,inline,let)
phason_f(auto &&...oo)
noexcept -> auto
{
	return _detail::build<phason_t>::template with<V>(XTAL_REF_(oo)...);
}


////////////////////////////////////////////////////////////////////////////////

template <vector_q A> requires integral_variable_q<absolve_u<A>>
struct phason<A>
:	grade<A>
{
};
template <vector_q A> requires     real_variable_q<absolve_u<A>>
struct phason<A>
{
	static auto constexpr M_data = _std::extent_v<based_t<A>>;

	using coordinate_type = valued_u<A>;

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
		using I_ = _std::initializer_list<coordinate_type>;

	protected:
		using          S_::N_data;
		using typename S_::U_data;
		using typename S_::V_data;
		using typename S_::W_data;

	public:// TYPE
		using initializer_list = I_;

		XTAL_DEF_(return,inline,set)   ordinate(coordinate_type const &co) noexcept {return bond::math::bit_fraction_f<  ordinate_type>(co);}
		XTAL_DEF_(return,inline,set) coordinate(  ordinate_type const & o) noexcept {return bond::math::bit_fraction_f<coordinate_type>( o);}

	public:// ACCESS
		using S_::self;
		using S_::twin;

	public:// CONSTRUCT
	//	using S_::S_;
	~	homotype()                 noexcept=default;
	//	homotype()                 noexcept=default;
		XTAL_NEW_(copy) (homotype, noexcept=default)
		XTAL_NEW_(move) (homotype, noexcept=default)

		XTAL_NEW_(implicit) homotype()
		noexcept
		:	homotype(size_type{0})
		{}
		
		XTAL_NEW_(explicit) homotype(same_q<size_type> auto const n)
		noexcept
		:	S_(n)
		{}

		XTAL_NEW_(explicit) homotype(continuous_field_q auto &&...oo)
		noexcept
		:	homotype(sizeof...(oo))
		{
			operator>>=({XTAL_REF_(oo)...});
		}
		XTAL_NEW_(implicit) homotype(I_ o)
		noexcept
		:	homotype(count_f(o))
		{
			operator>>=(XTAL_MOV_(o));
		}
		XTAL_NEW_(explicit) homotype(iterable_q auto &&o)
		noexcept
		:	homotype(count_f(o))
		{
			operator>>=(XTAL_REF_(o));
		}
		
	public:// RECONSTRUCT
		using S_::operator >>=;
		using S_::operator <<=;

		XTAL_DEF_(inline,let) operator >>=(I_                      o) noexcept -> auto & {_detail::move_to<[] XTAL_0FN_(alias) (T::ordinate)>(S_::data(), XTAL_REF_(o)); return self();}
		XTAL_DEF_(inline,let) operator >>=(iterable_q auto       &&o) noexcept -> auto & {_detail::move_to<[] XTAL_0FN_(alias) (T::ordinate)>(S_::data(), XTAL_REF_(o)); return self();}
		XTAL_DEF_(inline,let) operator >>=(iterable_q auto const  &o) noexcept -> auto & {_detail::copy_to<[] XTAL_0FN_(alias) (T::ordinate)>(S_::data(), XTAL_REF_(o)); return self();}
		
		XTAL_DEF_(inline,let)
		operator <<=(I_ o)
		noexcept -> auto &
		{
			auto i0 = S_::data(), iN = _std::next(i0, S_::size() - o.size());
			_detail::move_to<[] XTAL_0FN_(alias) (T::ordinate)>(iN, XTAL_REF_(o));
			return self();
		}
		XTAL_DEF_(inline,let)
		operator <<=(iterable_q auto &&o)
		noexcept -> auto &
		{
			auto i0 = S_::data(), iN = _std::next(i0, S_::size() - o.size());
			_detail::move_to<[] XTAL_0FN_(alias) (T::ordinate)>(iN, XTAL_REF_(o));
			return self();
		}
		XTAL_DEF_(inline,let)
		operator <<=(iterable_q auto const &o)
		noexcept -> auto &
		{
			auto i0 = S_::data(), iN = _std::next(i0, S_::size() - o.size());
			_detail::copy_to<[] XTAL_0FN_(alias) (T::ordinate)>(iN, XTAL_REF_(o));
			return self();
		}

	public:// OPERATE
		
		///\
		Scales all elements. \

		///\note\
		The symmetric signatures for `/=` and `*=` are declared-but-undefined \
		to avoid compilation-failure when type-checking e.g. `multiplicative_group_q`. \

	//	using S_::operator*=;
	//	using S_::operator/=;

		auto operator *= (T const &) noexcept -> T &;// Asymmetric!
		auto operator /= (T const &) noexcept -> T &;// Asymmetric!

		XTAL_DEF_(return,inline,let)  operator * (auto const &t) const noexcept -> auto {return twin() *= t;}
		XTAL_DEF_(return,inline,let)  operator / (auto const &t) const noexcept -> auto {return twin() /= t;}

		XTAL_DEF_(inline,let) operator *= (I_ t) noexcept -> auto & {return self() *= T(t);}
		XTAL_DEF_(inline,let) operator /= (I_ t) noexcept -> auto & {return self() /= T(t);}

		XTAL_DEF_(inline,let)
		operator /= (simplex_variable_q auto const &f)
		noexcept -> auto &
		{
			return operator*=(U_fix::alpha_1/f);
		}
		XTAL_DEF_(let)
		operator *= (real_variable_q auto const &f)
		noexcept -> auto &
		{
			using _fix = bond::fixture<decltype(f)>;
			auto &s = reinterpret_cast<phason_t<inordinate_type[N_data]> &>(self());
			/*/
			unsigned constexpr M_bias = T_fix::half.depth >> T_fix::half.width;
			unsigned constexpr M_size = T_fix::half.depth - M_bias;
			auto [m, n] = bond::math::bit_representation_f(f);
			m >>= n - M_size;
			s >>=     M_size;
			s  *= m;
			/*/
			XTAL_IF0
			XTAL_0IF (1*sizeof(ordinate_type) == sizeof(coordinate_type)) {
				unsigned constexpr M_bias = T_fix::half.depth >> T_fix::half.width;
				unsigned constexpr M_size = T_fix::half.depth - M_bias;
				auto [m, n] = bond::math::bit_representation_f(f);
				m >>= n - M_size;
				s >>=     M_size;
				s  *= m;
			}
			XTAL_IF0
			XTAL_0IF (2*sizeof(ordinate_type) == sizeof(coordinate_type)) {
				typename T_fix::sigma_type t_[2];
				typename U_fix::sigma_type const u(f*_fix::diplo_f(T_fix::full.depth));
				#pragma unroll
				for (int i{}; i < N_data; ++i) {
					reinterpret_cast<typename U_fix::sigma_type &>(t_) = u*s[i];
				//	t_[1] += t_[0] >> T_fix::positive.depth;// Round...
					s [i]  = t_[1];
				}
			}
			/***/
			return self();
		}
		XTAL_DEF_(inline,let)
		operator *= (integral_variable_q auto const &i)
		noexcept -> auto &
		{
			return S_::operator*=(i);
		}

		///\
		Offsets the first element. \
		
	//	using S_::operator-=;
	//	using S_::operator+=;

		XTAL_DEF_(inline,let)
		operator -= (auto const &f)
		noexcept -> auto &
		{
			return S_::operator-=(f);
		}
		XTAL_DEF_(inline,let)
		operator += (auto const &f)
		noexcept -> auto &
		{
			return S_::operator+=(f);
		}

		XTAL_DEF_(inline,let)
		operator -= (simplex_variable_q auto const &f)
		noexcept -> auto &
		{
			return operator+=(-f);
		}
		XTAL_DEF_(inline,let)
		operator += (integral_variable_q auto const &i)
		noexcept -> auto &
		{
			return self();
		}
		XTAL_DEF_(inline,let)
		operator += (real_variable_q auto const &f)
		noexcept -> auto &
		{
			get<0>(*this) += bond::math::bit_fraction_f<_std::make_signed_t<ordinate_type>>(f);
			return self();
		}

		template <class Y=coordinate_type>
		XTAL_DEF_(return,inline,let)
		continuity()
		noexcept -> Y
		{
			return condition_f<Y>(not discontinuity<bool>());
		}
		template <class Y=coordinate_type>
		XTAL_DEF_(return,inline,let)
		discontinuity()
		noexcept -> Y
		{
			auto [u0, u1] = *this;
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

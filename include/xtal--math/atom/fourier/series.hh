#pragma once
#include "./any.hh"
#include "./serial.hh"

#include "../../process/pade/unity.hh"



XTAL_ENV_(push)
namespace xtal::atom::math::fourier
{/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

template <class   ..._s>	struct  series;
template <class   ..._s>	using   series_t = typename series<_s...>::type;
template <class   ...Ts>	concept series_q = bond::tag_in_p<series_t, Ts...>;

XTAL_DEF_(let) series_f = [] XTAL_1FN_(call) (_detail::factory<series_t>::make);


////////////////////////////////////////////////////////////////////////////////
/*!
\brief   Extends `serial` with multiplication via circular convolution.
*/
template <scalar_q ..._s> requires same_q<_s...>
struct series<_s ...>
:	series<common_t<_s...>[sizeof...(_s)]>
{
};
template <vector_q A>
struct series<A>
{
private:
	using U0 = _xtd::remove_extent_t<A>;
	using U1 =  typename destruct<U0>::value_type;
	using U2 =  typename destruct<U1>::value_type;

	using _fit = bond::fit<A>;
	
	template <class T>
	using endotype = typename serial<A>::template homotype<T>;

	template <class T>
	using holotype = bond::compose_s<endotype<T>, bond::tag<series_t>>;

public:
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

	public:// CONSTRUCT
		using S_::S_;

		/*!
		\brief   Generates a section of the complex sinusoid determined by `std::pow(2, n_shift{})`.
		*/
		XTAL_DEF_(inline)
		XTAL_NEW_(explicit)
		homotype(constant_q auto const n_index, auto &&...oo)
		noexcept
		{
			generate<n_index>(XTAL_REF_(oo)...);
		}
		/*!
		\brief   Generates a section of the complex sinusoid determined by `std::pow(2, o_shift{})`.
		*/
		template <auto ...Ns>
		XTAL_DEF_(let)
		generate(constant_q auto const n, auto &&...oo)
		noexcept -> T &
		{
			generate<n, Ns...>(XTAL_REF_(oo)...);
		}

		/*!
		\brief   Populates `this` with powers of `u` from `N_index` to `N_index + N_count - 1`.
		\returns `this`.
		*/
		//\
		template <int N_count=size, int N_index=0, int N_step=1, int N_inset=0>
		template <int N_index=0, int N_inset=0, int N_step=1, int N_size=size>
		XTAL_DEF_(inline,let)
		generate(value_type const &u)
		noexcept -> T &
		{
			using process::math::limit_f;
			using process::math::power_f;
			using A_delta = typename S_::difference_type;

		//	Compute the start- and end-points for the required segment:
			A_delta constexpr N_limit = N_index + N_size;
			A_delta constexpr _0 =  0*N_step;
			A_delta constexpr _1 =  1*N_step;
			A_delta constexpr _2 =  2*N_step;
			A_delta constexpr I0 = _1*N_index + N_inset, J0 = _2*N_index + N_inset;
			A_delta constexpr IZ = _1*N_limit + N_inset, JZ = _2*N_limit + N_inset;

			auto &s = self();
			XTAL_IF0
			XTAL_0IF (N_index == -1 and N_inset == 0 and N_step == 1 and N_size == size) {
				generate<1, -1, 1, size - 1>(u);
				if constexpr (un_v<N_size&1>) {
					get<N_limit>(s) = power_f<2>(get<N_limit/2>(s));
				}
			}
			XTAL_0IF (N_size < I0 + _1) {
			//	Populate the 0th power only:
				get<I0>(s) = power_f<N_index>(u);
			}
			XTAL_0IF_(else) {
			//	Populate the 0th and 1st powers:
				auto const o = power_f<N_index>(u);
				get<I0 + _0>(s) = o;
				get<I0 + _1>(s) = o*u;

			//	Populate the remaining powers by squaring/multiplication:
				bond::seek_until_f<(N_size >> 1U)>([&] (auto M)
					XTAL_0FN {
						auto constexpr UM = I0 + _1*M;
						auto constexpr WM = J0 + _2*M;
						
						auto const w = power_f<2>(limit_f<-3>(get<UM>(s)));
						get<WM + _0>(s) =   w;
						get<WM + _1>(s) = u*w;
					}
				);
			//	Compute the final value if `N_size` is odd:
				if constexpr (in_v<N_size&1> and in_v<N_size^1>) {
					get<IZ - _1>(s) = get<IZ - _2>(s)*(u);
				}
			}
			return s;
		}

		/*!
		\brief   Generates the complex sinusoid with length `2*PI*std::pow(2, N_shift)`.
		\note    To generate the FFT basis used by `transform` etc, use `N_shift == -1`.
		\returns `this`.
		\tparam  `N_shift >= -3`
		*/
		template <int N_shift=0> requires complex_field_q<value_type>
		XTAL_DEF_(let)
		generate()
		noexcept -> T &
		{
			using _std::span;
			using _std::prev;
			using _std::next;

		//	Initialize the forwards and backwards iterators:
			auto const i = S_::begin();
			auto const j = S_::rend() - 1;
			
		//	Compute the fractional sinusoid for this `size`:
			auto constexpr y = process::math::pade::unity_f<(+1)>(_fit::ratio_f(-1, size << 1));

		//	Compute the initial `1/8`th then mirror the remaining segments:
			typename S_::difference_type constexpr M = size >> 2U;// `1/8`th
			static_assert(-4 <  N_shift);
			generate<0, 0, 1, M + (-3 <  N_shift)>(y);
			if constexpr (-2 <= N_shift) _detail::copy_to<[] (value_type const &v) XTAL_0FN_(to) (value_type{-v.imag(), -v.real()})>(prev(j, 2*M), span(i, next(i, 1*M)));
			if constexpr (-1 <= N_shift) _detail::copy_to<[] (value_type const &v) XTAL_0FN_(to) (value_type{ v.imag(), -v.real()})>(next(i, 2*M), span(i, next(i, 2*M)));
			if constexpr (-0 <= N_shift) _detail::copy_to<[] (value_type const &v) XTAL_0FN_(to) (value_type{-v.real(), -v.imag()})>(next(i, 4*M), span(i, next(i, 4*M)));
			static_assert( 0 >= N_shift);// TODO: Extend to allow multiple copies using `bond::seek`.
			
			return self();
		}
		/*!
		\brief   Transforms `source` using the FFT, with `this` as the Fourier basis.
		\returns `source`.
		
		\note    The size of both `this` and `source` must be expressible as an integral power-of-two,
		         and `1 < source.size() <= this->size()`.
		*/
		template <int N_direction=1> requires in_v<N_direction, 1, -1> and complex_field_q<value_type>
		XTAL_DEF_(let)
		transform(isomorphic_q<T> auto &source) const
		noexcept -> decltype(auto)
		{
			using Xs = XTAL_ALL_(source);
			using Ys = typename Xs::transverse::type;
			using I  = typename Xs:: difference_type;
		
		//	Ensure the size of both domain and codomain are powers of two:
			I constexpr N_width =        size();
			I const     n_width = source.size();
			I const     h_width = n_width >> 1U; assert(1 <= h_width);
			I const     n_depth = bond::math::bit_ceiling_f(n_width); assert(n_width == I{1} << n_depth);
			I constexpr N_depth = bond::math::bit_ceiling_f(N_width); assert(n_depth         <= N_depth);

		//	Move all entries to their bit-reversed locations:
			for (I h{}; h < h_width; ++h) {
				_std::swap(source[h], source[bond::math::bit_reverse_f(h, n_depth)]);
			}
		
		//	Conjugate the input if computing the inverse transform of the codomain:
			if constexpr (N_direction == -1) {
				_detail::apply_to<[] XTAL_1FN_(call) (_std::conj)>(source);
			}
		//	Compute the transform of `source` using the precomputed half-period sinusoid in `this`:
			for (I n{}; n < n_depth; ++n) {
				I const  u_width = I{1} << n;
				I const  w_width = u_width << 1U;
				I const un_depth = N_depth   - n;
				for (I u{ }; u < u_width; u +=     one) {auto const &o = S_::element(u << un_depth);
				for (I w{u}; w < n_width; w += w_width) {
					auto const m = w + u_width;
					value_type &y = source[m];
					value_type &x = source[w];
					value_type const yo = y*o;
					y  = x - yo;
					x +=     yo;
				}}
			}
		//	Conjugate and scale the output if computing the inverse transform of the codomain:
			if constexpr (N_direction == -1) {
				_detail::apply_to<[] XTAL_1FN_(call) (_std::conj)>(source);
				source /= _fit::alpha_f(n_width);
			}

		//	Cast the output to the transformed domain:
			return reinterpret_cast<Ys &>(source);
		}
		template <int N_direction=1> requires in_v<N_direction, 1, -1> and complex_field_q<value_type>
		XTAL_DEF_(let)
		transform(isomorphic_q<T> auto &&source) const
		noexcept -> decltype(auto)
		{
			(void) transform<N_direction>(source);
			return reinterpret_cast<typename XTAL_ALL_(source)::transverse::type &&>(source);
		}
		/*!
		\returns a new `series` representing the FFT of `source`,
		using `this` as the Fourier basis.
		*/
		template <int N_direction=1> requires in_v<N_direction, 1, -1>
		XTAL_DEF_(return,inline,let)
		transformation(isomorphic_q<T> auto source) const
		noexcept -> auto
		{
			return transform<N_direction>(XTAL_MOV_(source));
		}

		/*!
		\returns `lhs` convolved with `rhs`, using `this` as the Fourier basis.
		*/
		XTAL_DEF_(let)
		convolve(isomorphic_q<T> auto &&y0, auto y1) const
		noexcept -> decltype(auto)
		{
			static_assert(same_q<decltype(y0), decltype(y1)>);
			return transform<-1>(transform<1>(XTAL_REF_(y0)) *= transform<1>(y1));
		}
		/*!
		\returns A new `series` representing the convolution of `y0` with `y1`,
		using `this` as the Fourier basis.
		*/
		XTAL_DEF_(return,inline,let)
		convolution(isomorphic_q<T> auto y0, auto const &y1) const
		noexcept -> auto
		{
			static_assert(same_q<decltype(y0), decltype(y1)>);
			return convolve(XTAL_MOV_(y0), y1);
		}

		/*!
		\brief   Multiplication by circular convolution.
		*/
		using S_::operator*=;

		XTAL_DEF_(return,inline,let)  operator * (auto const &                       w) const noexcept -> auto   {return twin() *=   w ;}
		XTAL_DEF_(inline,let)         operator *=(_std::initializer_list<value_type> w)       noexcept -> auto & {return self() *= T(w);}

		XTAL_DEF_(let)
		operator *=(T const &t)
		noexcept -> T &
		{
			auto &s = self();
			if constexpr (complex_field_q<value_type>) {
				T(constant_t<-1>{}).convolve(s, t);
			}
			else {
				using X = typename _fit::aphex_type;
				using Y = typename series<X[size]>::type;
				Y s_(s);
				Y t_(t);
				Y(constant_t<-1>{}).convolve(s_, t_);
				_detail::move_to<[] XTAL_1FN_(call) (_std::real)>(s.begin(), s_);
			}
			return s;
		}

		/*!
		\brief   The dual of `T`, defined by `hadamard::series`.
		*/		
		struct transverse
		{
			template <class Y>
			//\
			using holotype = typename group_multiplication<A>::template homotype<T>;
			using holotype = typename group<wrap_s<A, _std::multiplies>>::template homotype<Y>;

			template <class Y>
			class homotype : public holotype<homotype<Y>>
			{
				using R_ = holotype<homotype<Y>>;
			
			public:
				using R_::R_;
				
				struct transverse {using type = T;};

			};
			using type = bond::derive_t<homotype>;

		};
	};
	using type = bond::derive_t<homotype>;

};


///////////////////////////////////////////////////////////////////////////////
}/////////////////////////////////////////////////////////////////////////////
XTAL_ENV_(pop)

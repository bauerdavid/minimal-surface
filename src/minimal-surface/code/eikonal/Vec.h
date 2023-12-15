#pragma once
#include <initializer_list>
#include <algorithm>
#include <functional>
#include <utility>
#include <iostream>

template<typename VecType, size_t N, typename T>
class OperableVec {
protected:
	T data[N]; 
public:
    typedef VecType vec_type;
	static const size_t size = N;
	typedef T type;
	template<typename VecT>
	friend VecType operator+(const OperableVec& l, const OperableVec<VecT, N, T>& r);

	template<typename Scalar, typename>
	friend VecType operator+(const OperableVec& l, Scalar r);

	template<typename Scalar, typename>
	friend VecType operator+(Scalar r, const OperableVec& l);


	template<typename VecT>
	friend VecType operator-(const OperableVec& l, const OperableVec<VecT, N, T>& r);

	template<typename Scalar, typename>
	friend VecType operator-(const OperableVec& l, Scalar r);

	template<typename Scalar, typename>
	friend VecType operator-(Scalar l, const OperableVec& r);


	template<typename VecT>
	friend VecType operator*(const OperableVec& l, const OperableVec<VecT, N, T>& r);

	template<typename Scalar, typename>
	friend VecType operator*(const OperableVec& l, Scalar r);

	template<typename Scalar, typename>
	friend VecType operator*(Scalar l, const OperableVec& r);


	template<typename VecT>
	friend VecType operator/(const OperableVec& l, const OperableVec<VecT, N, T>& r);

	template<typename Scalar, typename>
	friend VecType operator/(const OperableVec& l, Scalar r);

	template<typename Scalar, typename>
	friend VecType operator/(Scalar l, const OperableVec& r);
	
	OperableVec() {
		std::fill(data, data + N, T(0));
	}

	OperableVec(std::initializer_list<T> list) {
		auto list_data = std::data(list);
		std::copy(list_data, list_data + N, data);
	}
	template< typename VecType2>
	OperableVec(const OperableVec<VecType2, N, T>& v) {
		std::copy(v.data, v.data + N, data);
	}
	template< typename VecType2>
	OperableVec(const OperableVec<VecType2, N, T>&& v) noexcept {
		std::copy(v.data, v.data + N, data);
	}

	template< typename VecType2>
	OperableVec& operator=(OperableVec<VecType2, N, T> v) {
		swap(*this, v);
		return *this;
	}
	template<typename VecType2>
	friend void swap(OperableVec<VecType, N, T>& first, OperableVec<VecType2, N, T>& second)
	{
		using std::swap;
		swap(first.data, second.data);
	}

	T& operator[](int i) {
		return data[i];
	}
	const T& operator[](int i) const {
		return data[i];
	}
	template<typename VecType2>
	friend OperableVec& operator+=(OperableVec<VecType, N, T>& v1, OperableVec<VecType2, N, T>& v2) {
		for (int i = 0; i < N; i++) {
			v1[i] += v2[i];
		}
		//std::transform(v1.data, v1.data + N, v2.data, v1.data, std::plus<T>());
		return v1;
	}
	template<typename Scalar, typename = typename std::enable_if<std::is_arithmetic<Scalar>::value, Scalar>::type>
	OperableVec& operator+=(Scalar s) {
		std::transform(data, data + N, data, [s](T val) {return val + s;});
		return *this;
	}

	template<typename VecType2>
	friend OperableVec& operator-=(OperableVec<VecType, N, T>& v1, OperableVec<VecType2, N, T>& v2) {
		for (int i = 0; i < N; i++) {
			v1[i] -= v2[i];
		}
		return v1;
	}
	template<typename Scalar, typename = typename std::enable_if<std::is_arithmetic<Scalar>::value, Scalar>::type>
	OperableVec& operator-=(Scalar s) {
		std::transform(data, data + N, data, [s](T val) {return val - s;});
		return *this;
	}

	template<typename VecType2>
	friend OperableVec& operator*=(OperableVec<VecType, N, T>& v1, OperableVec<VecType2, N, T>& v2) {
		for (int i = 0; i < N; i++) {
			v1[i] *= v2[i];
		}
		return v1;
	}
	template<typename Scalar, typename = typename std::enable_if<std::is_arithmetic<Scalar>::value, Scalar>::type>
	OperableVec& operator*=(Scalar s) {
		std::transform(data, data + N, data, [s](T val) {return val * s; });
		return *this;
	}

	template<typename VecType2>
	friend OperableVec& operator/=(OperableVec<VecType, N, T>& v1, OperableVec<VecType2, N, T>& v2) {
		for (int i = 0; i < N; i++) {
			v1[i] /= v2[i];
		}
		return v1;
	}
	template<typename Scalar, typename = typename std::enable_if<std::is_arithmetic<Scalar>::value, Scalar>::type>
	OperableVec& operator/=(Scalar s) {
		std::transform(data, data + N, data, [s](T val) {return val / s;});
		return *this;
	}

	friend std::ostream& operator<<(std::ostream& os, const OperableVec& vec){
        os << "{";
        for(int i=0; i<N; i++){
            os << vec.data[i];
            if(i<N-1){
                os << ", ";
            }
        }
        os << "}";
    return os;
}

	T Sum() {
		T sum(0);
		std::for_each(data, data + N, [&sum](T val) {sum += val; });
		return sum;
	}

	double Norm() {
	    double square_sum(0);
	    std::for_each(data, data + N, [&square_sum](T val) {square_sum += val*val; });
	    return sqrt(square_sum);
	}

	T* end() {
		return data + N;
	}

	T* begin() {
		return data;
	}

	operator std::vector<T>() {
		return std::vector<T>(begin(), end());
	}


};

template<typename OperableVecType, typename VecType2>
typename OperableVecType::vec_type operator+(const OperableVecType& l, const OperableVec<VecType2, OperableVecType::size, typename OperableVecType::type>& r) {
	typename OperableVecType::vec_type v;
	for (int i = 0; i < OperableVecType::size; i++) {
		v[i] = l[i] + r[i];
	}
	return v;
}

template<typename OperableVecType, typename Scalar, typename = typename std::enable_if<std::is_arithmetic<Scalar>::value, Scalar>::type>
typename OperableVecType::vec_type operator+(const OperableVecType& l, Scalar r) {
	typename OperableVecType::vec_type v;
	std::transform(l.data, l.data + OperableVecType::size, v.data, [r](typename OperableVecType::type val) {return val + r;});
	return v;
}

template<typename OperableVecType, typename Scalar, typename = typename std::enable_if<std::is_arithmetic<Scalar>::value, Scalar>::type>
typename OperableVecType::vec_type operator+(Scalar l,  const OperableVecType& r) {
	typename OperableVecType::vec_type v;
	std::transform(r.data, r.data + OperableVecType::size, v.data, [l](typename OperableVecType::type val) {return l + val;});
	return v;
}

template<typename OperableVecType, typename VecType2>
typename OperableVecType::vec_type operator-(const OperableVecType& l, const OperableVec<VecType2, OperableVecType::size, typename OperableVecType::type>& r) {
	typename OperableVecType::vec_type v;
	for (int i = 0; i < OperableVecType::size; i++) {
		v[i] = l[i] - r[i];
	}
	return v;
}

template<typename OperableVecType, typename Scalar, typename = typename std::enable_if<std::is_arithmetic<Scalar>::value, Scalar>::type>
typename OperableVecType::vec_type operator-(const OperableVecType& l, Scalar r) {
	typename OperableVecType::vec_type v;
	std::transform(l.data, l.data + OperableVecType::size, v.data, [r](typename OperableVecType::type val) {return val - r;});
	return v;
}

template<typename OperableVecType, typename Scalar, typename = typename std::enable_if<std::is_arithmetic<Scalar>::value, Scalar>::type>
typename OperableVecType::vec_type operator-(Scalar l, const OperableVecType& r) {
	typename OperableVecType::vec_type v;
	std::transform(r.data, r.data + OperableVecType::size, v.data, [l](typename OperableVecType::type val) {return l - val;});
	return v;
}

template<typename OperableVecType, typename VecType2>
typename OperableVecType::vec_type operator*(const OperableVecType& l, const OperableVec<VecType2, OperableVecType::size, typename OperableVecType::type>& r) {
	typename OperableVecType::vec_type v;
	for (int i = 0; i < OperableVecType::size; i++) {
		v[i] = l[i] * r[i];
	}
	return v;
}

template<typename OperableVecType, typename Scalar, typename = typename std::enable_if<std::is_arithmetic<Scalar>::value, Scalar>::type>
typename OperableVecType::vec_type operator*(const OperableVecType& l, Scalar r) {
	typename OperableVecType::vec_type v;
	std::transform(l.data, l.data + OperableVecType::size, v.data, [r](typename OperableVecType::type val) {return val * r;});
	return v;
}

template<typename OperableVecType, typename Scalar, typename = typename std::enable_if<std::is_arithmetic<Scalar>::value, Scalar>::type>
typename OperableVecType::vec_type operator*(Scalar l, const OperableVecType& r) {
	typename OperableVecType::vec_type v;
	std::transform(r.data, r.data + OperableVecType::size, v.data, [l](typename OperableVecType::type val) {return l * val;});
	return v;
}

template<typename OperableVecType, typename VecType2>
typename OperableVecType::vec_type operator/(const OperableVecType& l, const OperableVec<VecType2, OperableVecType::size, typename OperableVecType::type>& r) {
	typename OperableVecType::vec_type v;
	for (int i = 0; i < OperableVecType::size; i++) {
		v[i] = l[i] / r[i];
	}
	return v;
}

template<typename OperableVecType, typename Scalar, typename = typename std::enable_if<std::is_arithmetic<Scalar>::value, Scalar>::type>
typename OperableVecType::vec_type operator/(const OperableVecType& l, Scalar r) {
	typename OperableVecType::vec_type v;
	std::transform(l.data, l.data + OperableVecType::size, v.data, [r](typename OperableVecType::type val) {return val / r;});
	return v;
}

template<typename OperableVecType, typename Scalar, typename = typename std::enable_if<std::is_arithmetic<Scalar>::value, Scalar>::type>
typename OperableVecType::vec_type operator/(Scalar r, const OperableVecType& l) {
	tzpename OperableVecType::vec_type v;
	std::transform(l.data, l.data + OperableVecType::size, v.data, [l](typename OperableVecType::type val) {return l / val;});
	return v;
}

template<size_t N, typename T>
class Vec : public OperableVec<Vec<N, T>, N, T> {
	using OperableVec<Vec<N, T>, N, T>::OperableVec;
};


template<typename T>
class Vec3 : public OperableVec<Vec3<T>, 3, T> {
public:
	using OperableVec<Vec3<T>, 3, T>::OperableVec;
	using OperableVec<Vec3<T>, 3, T>::data;
	Vec3() :OperableVec() {}
	Vec3(T x, T y, T z) : OperableVec({ x, y, z }) {}
	T& x() {
		return data[0];
	}
	const T& x() const {
		return data[0];
	}

	T& y() {
		return data[1];
	}
	const T& y() const {
		return data[1];
	}

	T& z() {
		return data[2];
	}
	const T& z() const {
		return data[2];
	}

	template<std::size_t Index>
	std::tuple_element_t<Index, Vec3> get() const&
	{
		return this->operator[](Index);
	}

	template<std::size_t Index>
	std::tuple_element_t<Index, Vec3> get() const&&
	{
		return std::move(this->operator[](Index));
	}

	bool operator==(const Vec3<T>& p) const {
		return this->x() == p.x() && this->y() == p.y() && this->z() == p.z();
	}

	bool operator<(const Vec3<T>& p) const
	{

		return (hash_vec3<T>(*this) < hash_vec3<T>(p));
	}
};

template<typename T>
class Vec2 : public OperableVec<Vec2<T>, 2, T> {
public:
	using OperableVec<Vec2<T>, 2, T>::OperableVec;
	using OperableVec<Vec2<T>, 2, T>::data;
	T& x() {
		return data[0];
	}
	const T& x() const {
		return data[0];
	}

	T& y() {
		return data[1];
	}
	const T& y() const {
		return data[1];
	}
};

namespace std
{
	template<size_t Index, typename T>
	struct tuple_element<Index, ::Vec3<T>>
	{
		using type = T;
	};
	template<typename T>
	struct tuple_size<::Vec3<T>> : std::integral_constant<size_t, 3> {};
	
}

template<typename T>
inline size_t hash_vec3(const Vec3<T>& p) {
	//return ((((size_t)p.x) * 73856093) ^ (((size_t)p.y) * 19349663) ^ (((size_t)p.z) * 83492791)) % 1000000;
	return (size_t)p.x() << (2 * 21) | (size_t)p.y() << 21 | (size_t)p.z();
}

template<typename T>
struct Vec3Hash {
public:

	size_t operator()(const Vec3<T>& p) const {
		return hash_vec3(p);
	}
};
#pragma once
#include <vector>
#include <iostream>
#include <assert.h>
#include <math.h>
#include <fstream>
#include <random>




class Vec3 {
private:
	std::vector<float> V;
public:
	Vec3(float x, float y, float z) {
		//std::cout << "here";
		V.push_back(x);
		V.push_back(y);
		V.push_back(z);
	}
	Vec3() {
		V.resize(3);
		std::fill(V.begin(), V.end(), 0.f);
	}
	float& operator[](int index) {
		assert(index >= 0 && index <= 2);
			return V[index];
		
	}
	const float& operator[](int index) const {
		assert(index >= 0 && index <= 2);
		return V[index];

	}
	Vec3(const Vec3& v) {
		
		int j = 0;
		for (int i = 0; i != 3; i++) {
			V.push_back(v[i]);
		}
	}
	Vec3& operator=(const Vec3& v) {
		if (this == &v) {
			return *this;
		}

		
		for (int i = 0; i != 3; i++) {
			V[i] = v.V[i];
		}
		return *this;

	}
	Vec3& operator+=(const Vec3& v){
		
		int j = 0;
		for (int i = 0; i != 3; i++) {
			V[i] += v[i];
		}
		return *this;
	}
	Vec3 operator+(const Vec3& v) {
		Vec3 res{};
		for (int i = 0; i != 3; i++) {

			
			 res[i] = (*this)[i]+v[i];
			
		}
		
		return res;
	}
	Vec3 operator-(const Vec3& v) const {
		Vec3 res{};
		for (int i = 0; i != 3; i++) {


			res[i] = (*this)[i] - v[i];

		}

		return res;
	}
	Vec3& operator-(const int x) {
		for (auto& i : V) {
			i -= x;
		}
		return *this;
	}
	Vec3& operator+(const int x) {
		for (auto& i : V) {
			i += x;
		}
		return *this;
	}
	Vec3 operator*(const float x) const {
		Vec3 res{};
		for (int i = 0; i != 3; i++) {


			res[i] = (*this)[i] * x;

		}

		return res;
	}
	
	Vec3& operator*=(const int x) {
		for (auto& i : V) {
			i *= x;
		}
		return *this;
	}
	Vec3& operator*=(const float x) {
		for (auto& i : V) {
			i *= x;
		}
		return *this;
	}

	bool operator==(const Vec3& v) const {
		
		return ((*this)[0] == v[0] && (*this)[1] == v[1] && (*this)[2] == v[2]);
		
	}
	bool operator!=(const Vec3& v) const {
		return !(*this == v);
	}
	friend std::ostream&  operator<<(std::ostream& os, const Vec3& v) {
		os << v[0] << ' ' << v[1] << ' ' << v[2] << '\n';
		return os;
	}
	friend std::ofstream& operator<<(std::ofstream& os, const Vec3& v) {
		os << v[0] << ' ' << v[1] << ' ' << v[2] << '\n';
		return os;
	}
	Vec3& normalise() {
		float magnitude = 0;
		for (int i = 0; i != 3; i++) {
			magnitude += V[i] * V[i];
		}
		magnitude = sqrtf(magnitude);
		(*this) *= (1 / magnitude);
		return *this;
	}
	float length() {
		return sqrt(V[0] * V[0] + V[1] * V[1] + V[2] * V[2]);
	}
	bool is_normal() const{
		float magnitude = 0; 
		for (auto i : V) {
			magnitude += i * i;
		}
		if (magnitude != 1.f) {
			return false;
		}
		else {
			return true;
		}
	}
	void clamp(float min, float max) {
		for (auto& i : V) {
			if (i < min) {
				i = min;
			}
			if (i > max) {
				i = max;
			}
		}
	}
	friend float dot(const Vec3& v, const Vec3& v1) {
		float res = 0;
		for (int i = 0; i != 3; i++) {
			res += v[i] * v1[i];
		}
		return res;
	}
	void Gamma_corrected(float Gamma) {
		for (auto& i : V) {
			i = powf(i, Gamma);
		}
	}
};


class ray
{
private:
	Vec3 m_Origin; 
	Vec3 m_Direction;
	Vec3 m_Colour;
public:
	ray( Vec3 origin, Vec3 direction, Vec3 Colour = Vec3(1.0,1.0,1.0), bool unit = true) {
		m_Origin = origin;
		m_Colour = Colour;
		if (unit) {
			if (direction.is_normal()) {
				m_Direction = direction;
			}
			else {
				m_Direction = direction.normalise();
			}
		}
		else {
			m_Direction = direction;
		}
	}
	ray() {
		m_Origin = m_Direction = m_Colour =  Vec3();


	}
	bool is_zero() {
		if (this->m_Colour == Vec3(0.0, 0.0, 0.0) && this->m_Direction == Vec3(0.0, 0.0, 0.0) && this->m_Origin == Vec3(0.0, 0.0, 0.0)) {
			return true;
		}
		return false;
	}

	friend std::ostream& operator<<(std::ostream& os, const ray& r) {
		os << "The origin is : " << r.m_Origin;
		os << "The direction is : " << r.m_Direction;
		return os;
		
	}
	Vec3 ray_colour(Vec3(*f)(ray& r)) {
		return f(*this);
	}
	Vec3& origin() { return m_Origin; };
	Vec3& direction() { return m_Direction; }
	Vec3& colour() { return m_Colour; }

};


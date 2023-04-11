/*******************************************************************************
 
Copyright (c) 2019 Arkan

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.                                      *
 *                                                                              *
 *******************************************************************************/

#ifndef CGOGN_MODELING_TYPES_QUATERNIONS_H_
#define CGOGN_MODELING_TYPES_QUATERNIONS_H_

#include <cgogn/geometry/types/vector_traits.h>
#include <cgogn/rendering/types.h>

namespace cgogn
{

namespace modeling
{

/**
 * @Class quaternion
 * from github.com/brainexcerpts/
 * and coupled cage skeleton code
*/
template <typename T>
class Quaternion
{

	using Vec3 = geometry::Vec3;

private:
	union {

      struct {
         T _x;
         T _y;
         T _z;
         T _s;
      };

      T _v[4];

   };

public:

   //Constructor for identity rotation
   Quaternion()
   {
      _x = _y = _z = 0.0;
      _s = 1.0;
   }

   Quaternion(const Vec3 & axis, const T & angle)
   {
      setRotation(axis, angle);
   }

   Quaternion(const T x, const T y, const T z, const T s)
   {
      set(x, y, z, s);
   }

   inline const T *ptr() const
   {
      return _v;
   }

   inline T & x() { return _x; }
   inline T & y() { return _y; }
   inline T & z() { return _z; }
   inline T & s() { return _s; }
   inline Vec3 v() { return Vec3(_x, _y, _z); }

   inline const T & x() const { return _x; }
   inline const T & y() const { return _y; }
   inline const T & z() const { return _z; }
   inline const T & s() const { return _s; }
   inline const Vec3 v() const { return Vec3(_x, _y, _z); }

   inline T & operator[](const int i)
   {
      return _v[i];
   }

   inline const T & operator[](const int i) const
   {
      return _v[i];
   }

   inline void setZero() { _x = _y = _z = _s = 0.0; }

   inline void setRotation(const Vec3 & axis, const T & angle)
   {
      T sinHalfAngle = std::sin(angle / 2.0);
      _x = sinHalfAngle * axis.x();
      _y = sinHalfAngle * axis.y();
      _z = sinHalfAngle * axis.z();
      _s = std::cos(angle / 2.0);

      normalize();
   }

   inline void set(const T x, const T y, const T z, const T s)
   {
      _x = x;
      _y = y;
      _z = z;
      _s = s;
   }

   inline T normalize()
   {
      T len = length();
      _x /= len;
      _y /= len;
      _z /= len;
      _s /= len;
      return len;
   }

   inline T lengthSquared() const
   {
      return (
               _x * _x +
               _y * _y +
               _z * _z +
               _s * _s
               );
   }

   inline T norm() const
   {
      return lengthSquared();
   }

   inline T length() const
   {
      return std::sqrt( lengthSquared() );
   }

   inline T magnitude() const
   {
      return length();
   }

   inline T dot(const Quaternion<T> & other)
   {
      return _x * other._x +
            _y * other._y +
            _z * other._z +
            _s * other._s;
   }

   inline Quaternion<T> inverse() const
   {
      T _lengthSquared = lengthSquared();

      T x = (-_x) / _lengthSquared;
      T y = (-_y) / _lengthSquared;
      T z = (-_z) / _lengthSquared;
      T s = _s / _lengthSquared;

      return Quaternion<T>(x, y, z, s);
   }

   inline Quaternion<T> conjugate() const //or flip
   {
      return Quaternion<T>(-_x, -_y, -_z, _s); //only if it is normalized
   }

   inline Vec3 applyRotation(const Vec3 & v) const
   {
      //CHECK IF NORMALIZED
      //std::cout << "length " << length() << std::endl; //If normalized, it is ==1

      Quaternion<T> vq(v.x(), v.y(), v.z(), 0.0);

      Quaternion<T> result = ((*this) * vq * conjugate());  //TODO: Optimize

      return Vec3(result.x(), result.y(), result.z());
   }

   inline Vec3 toEuler() const
   {
      constexpr double radToGrad = 180 / M_PI;

      //code from https://en.wikipedia.org/wiki/Conversion_between_quaternions_and_Euler_angles

      // roll (x-axis rotation)
      T sinr_cosp = 2 * (_s * _x + _y * _z);
      T cosr_cosp = 1 - 2 * (_x * _x + _y * _y);
      T roll = std::atan2(sinr_cosp, cosr_cosp);

      // pitch (y-axis rotation)
      T sinp = 2 * (_s * _y - _z * _x);
      T pitch;
      if (std::abs(sinp) >= 1)
         pitch = std::copysign(M_PI / 2, sinp); // use 90 degrees if out of range
      else
         pitch = std::asin(sinp);

      // yaw (z-axis rotation)
      T siny_cosp = 2 * (_s * _z + _x * _y);
      T cosy_cosp = 1 - 2 * (_y * _y + _z * _z);
      T yaw = std::atan2(siny_cosp, cosy_cosp);

      return Vec3(roll*radToGrad, pitch*radToGrad, yaw*radToGrad);
   }

};

template<typename T>
Quaternion<T> operator+(const Quaternion<T> & a, const Quaternion<T> & b)
{
   return Quaternion<T>(a.x() + b.x(),
                        a.y() + b.y(),
                        a.z() + b.z(),
                        a.s() + b.s());
}

template<typename T>
Quaternion<T> operator-(const Quaternion<T> & a, const Quaternion<T> & b)
{
   return Quaternion<T>(a.x() - b.x(),
                        a.y() - b.y(),
                        a.z() - b.z(),
                        a.s() - b.s());
}

template<typename T>
Quaternion<T> operator*(const Quaternion<T> & a, const Quaternion<T> & b)
{
   T        s1 = a.s();
   Vec3  v1 = a.v();
   T        s2 = b.s();
   Vec3  v2 = b.v();

   T        sNew = (s1*s2)-(v1.dot(v2));
   Vec3  vNew = (s1*v2)+(s2*v1)+(v1.cross(v2));

   return Quaternion<T>(vNew.x(), vNew.y(), vNew.z(), sNew);
}

template<typename T>
Quaternion<T> operator*(const Quaternion<T> & a, const T s)
{
   return Quaternion<T>(s*a.x(), s*a.y(), s*a.z(), s*a.s());
}

template<typename T>
Quaternion<T> operator*(const T s, const Quaternion<T> & a)
{
   return Quaternion<T>(s*a.x(), s*a.y(), s*a.z(), s*a.s());
}

typedef Quaternion<double> dQuaternion;


} // namespace modeling

} // namespace cgogn

#endif // CGOGN_MODELING_TYPES_QUATERNIONS_H_

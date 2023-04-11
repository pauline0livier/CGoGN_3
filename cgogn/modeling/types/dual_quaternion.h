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

#ifndef CGOGN_MODELING_TYPES_DUAL_QUATERNIONS_H_
#define CGOGN_MODELING_TYPES_DUAL_QUATERNIONS_H_

#include <cgogn/geometry/types/vector_traits.h>
#include <Eigen/Geometry>
#include <cgogn/rendering/types.h>

#include <cgogn/modeling/types/quaternion.h>

/**
 * @Class dual quaternion
 * from github.com/brainexcerpts/
 * and coupled cage skeleton deformation 
*/
namespace cgogn
{

namespace modeling
{


template<class T> 
class DualQuaternion
{

	using Vec3 = geometry::Vec3;

private:

   Quaternion<T> _qReal;  //rotation (<0,0,0>,1) is identity
   Quaternion<T> _qDual;  //translation (<0,0,0>,0) is identity

public:

   DualQuaternion():
      _qReal(), //rotation identity (s=1.0 represents rotation)
      _qDual(0,0,0,0)   //represents the point <0,0,0> (s=0.0 represents point/vec)
   {

   }

   DualQuaternion(const Quaternion<T> & qReal, const Quaternion<T> & qDual)
   {
      _qReal = qReal;
      _qDual = qDual;
   }

   DualQuaternion(const Quaternion<T> & rotation, const Vec3& translation)
   {
      set(rotation, translation);
   }

   DualQuaternion(const rendering::Transfo3d& t)
   {
	  Eigen::Quaterniond q(t.rotation()); 
	  Quaternion<T> new_quat = Quaternion<T>(q.x(), q.y(), q.z(), q.w()); 

      set(new_quat,t.translation());
   }

   DualQuaternion(const Vec3 & axis,
                  const T & angle,
                  const Vec3 & translation):
                     _qReal(axis, angle)
   {
      Quaternion<T> t(translation[0]*0.5,
                      translation[1]*0.5,
                      translation[2]*0.5,
                      0);
      _qDual = t * _qReal;
   }

   inline Quaternion<T> & qReal() { return _qReal; }
   inline Quaternion<T> & qDual() { return _qDual; }

   inline const Quaternion<T> & qReal() const { return _qReal; }
   inline const Quaternion<T> & qDual() const { return _qDual; }

   inline void setZero() {_qReal.setZero(); _qDual.setZero();}

   inline void setIdentity()
   {
      _qReal.set(0, 0, 0, 1);
      _qDual.set(0, 0, 0, 0);
   }

   inline void set(const Quaternion<T> & rotation, const Vec3 & translation)
   {
      _qReal = rotation;
      //_qReal.normalize();  //to be sure it is a rotation

      Quaternion<T> t(translation[1],
                      translation[2],
                      translation[3],
                      0);
      _qDual = t * _qReal;
      _qDual = _qDual * 0.5;
   }

   inline void setRotationTranslation(const Vec3 & axis,
                                      const T & angle,
                                      const Vec3& translation)
   {
      _qReal.setRotation(axis, angle);

      Quaternion<T> t(translation[0],
                      translation[1],
                      translation[2],
                      0);
      _qDual = t * _qReal;
      _qDual = _qDual * 0.5;
   }

   inline T dot(const DualQuaternion<T> & other)
   {
      return _qReal.dot(other._qReal);
   }

   inline DualQuaternion<T> conjugate() const
   {
      return DualQuaternion<T>(_qReal.conjugate(), _qDual.conjugate());
   }

   inline DualQuaternion<T> normalize() const
   {
      T mag = _qReal.magnitude();
      mag = 1.0 / mag;
      return DualQuaternion<T>(_qReal * mag, _qDual * mag);

      //Alternative formulation
      //Quaternion<T> nqr = _qReal * mag;
      //Quaternion<T> nqd = _qDual * mag;
      //WHY?!
      //return DualQuaternion<T>(nqr, nqd - (nqd.dot(nqr) * nqr));
   }

   inline Vec3 applyTransformation(const Vec3 & v) const
   {
      /*Vec3 displacement =
            2.0 * (
                 _qReal.v().cross( _qReal.v().cross(v) + v*_qReal.s() )
               + _qDual.v() * _qReal.s()
               - _qReal.v() * _qReal.s()
               + _qReal.v().cross(_qDual.v())
            );*/

      //Alternative formulation
      Vec3 displacement =
            2.0 * (
               _qReal.v().cross(_qReal.v().cross(v) + (v * _qReal.s()) + _qDual.v())
               + _qDual.v() * _qReal.s()
               - _qReal.v() * _qDual.s()
               ); //dqtar

      return v + displacement;
   }

};

template<typename T>
DualQuaternion<T> operator+(const DualQuaternion<T> & a, const DualQuaternion<T> & b)
{
   Quaternion<T> qRealNew = a.qReal() + b.qReal();
   Quaternion<T> qDualNew = a.qDual() + b.qDual();

   return DualQuaternion<T>(qRealNew, qDualNew);
}

template<typename T>
DualQuaternion<T> operator*(const DualQuaternion<T> & a, const DualQuaternion<T> & b)
{
   Quaternion<T> qRealNew = b.qReal() * a.qReal();
   Quaternion<T> qDualNew = (b.qDual()*a.qReal()) + (b.qReal() * a.qDual());

   return DualQuaternion<T>(qRealNew, qDualNew);
}

template<typename T>
DualQuaternion<T> operator*(const DualQuaternion<T> & a, T s)
{
   DualQuaternion<T> ret(s*a.qReal(), s*a.qDual());
   return ret;
}

template<typename T>
DualQuaternion<T> operator*(T s, const DualQuaternion<T> & a)
{
   DualQuaternion<T> ret(s*a.qReal(), s*a.qDual());
   return ret;
}

typedef DualQuaternion<double> dDualQuaternion;


} // namespace modeling

} // namespace cgogn

#endif // CGOGN_MODELING_TYPES_DUAL_QUATERNIONS_H_

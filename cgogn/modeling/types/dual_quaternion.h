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

namespace cgogn
{

namespace modeling
{

/**
 * @Class dual quaternion
 * from github.com/brainexcerpts/
*/
class DualQuaternion
{

	using Vec3 = geometry::Vec3;

public:
	/// @brief default constructor
	DualQuaternion()
	{
		DualQuaternion res = dual_quaternion_from(Eigen::Quaterniond(1.0, 0.0, 0.0, 0.0), 
											Vec3(0.0, 0.0, 0.0));
		//_quat_0 = res.get_non_dual_part(); 
		//_quat_e = res.get_dual_part(); 
		 *this = res;
	}

	~DualQuaternion()
	{
	}

	/// @brief constructor from values 
	/// @param q0 
	/// @param qe 
	DualQuaternion(const Eigen::Quaterniond& q0, const Eigen::Quaterniond& qe)
	{
		_quat_0 = q0;
		_quat_e = qe;
	}

	/// @brief constructor from quaternion and translation vector
	/// @param q 
	/// @param t 
	DualQuaternion(const Eigen::Quaterniond& q, const Vec3& t)
	{
		DualQuaternion res = dual_quaternion_from(q, t);
		//_quat_0 = res.get_non_dual_part();
		//_quat_e = res.get_dual_part(); 

		  *this = res;
		
	}

	/// @brief define dual quaternion from transformation
	/// @param t 
	DualQuaternion(const rendering::Transfo3d& t)
	{
		Eigen::Quaterniond q(t.linear());
		Vec3 translation = t.translation();
		DualQuaternion res = dual_quaternion_from(q, translation);
		// _quat_0 = res.get_non_dual_part(); 
		//_quat_e = res.get_dual_part();

		 *this = res;
	}

	void set_from_transformation(const rendering::Transfo3d& t)
	{
		
		Eigen::Quaterniond q(t.linear());
		Vec3 translation = t.translation();
		DualQuaternion res = dual_quaternion_from(q, translation);
		*this = res;
	}



	/// @brief normalize the quaternion
	void normalize()
	{
		double norm = _quat_0.norm();
		_quat_0.normalize();
		//_quat_e = _quat_e / norm;
		_quat_e = 
			Eigen::Quaterniond(_quat_e.w()/norm, _quat_e.x()/norm, _quat_e.y()/norm, _quat_e.z()/norm); 
	}

	/// @brief rotate quaternion
	/// @param p 
	/// @return 
	Vec3 rotate_quaternion(const Eigen::Quaterniond q, const Vec3& p)
	{
		Vec3 q_vec = q.vec();
		return p + (q_vec*2.0).cross( q_vec.cross(p) + p*q.w() );
	}

	/// @brief transform position of point p with dual quaternion
	/// @param p 
	/// @return new position 
	Vec3 transform(const Vec3& p ) const 
	{
		double norm = _quat_0.norm();

		Eigen::Quaterniond qblend_0 = 
			Eigen::Quaterniond(_quat_0.w()/norm, _quat_0.x()/norm, 
									_quat_0.y()/norm, _quat_0.z()/norm); 
		Eigen::Quaterniond qblend_e = 
			Eigen::Quaterniond(_quat_e.w()/norm, _quat_e.x()/norm, 
									_quat_e.y()/norm, _quat_e.z()/norm);

		// Translation from the normalized dual quaternion equals :
		// 2.f * qblend_e * conjugate(qblend_0)
		Vec3 v0 = qblend_0.vec();
		Vec3 ve = qblend_e.vec();
		Vec3 trans = (ve*qblend_0.w() - v0*qblend_e.w() + v0.cross(ve)) * 2.0;

		// Rotate
		Vec3 rotation_part =  p + (v0*2.0).cross(v0.cross(p) + p*qblend_0.w());
		return rotation_part + trans; 
	}

	/// @brief dual quaternion from quaternion and translation vector
	/// @param q quaternion
	/// @param t translation vector
	/// @return dual quaternion
	DualQuaternion dual_quaternion_from(const Eigen::Quaterniond q, const Vec3& t) const
	{
		double w = -0.5 * (t[0] * q.x() + t[1] * q.y() + t[2] * q.z());
		double x = 0.5 * (t[0] * q.w() + t[1] * q.z() - t[2] * q.y());
		float y = 0.5 * (-t[0] * q.z() + t[1] * q.w() + t[2] * q.x());
		float z = 0.5 * (t[0] * q.y() - t[1] * q.x() + t[2] * q.w());

		return DualQuaternion(q, Eigen::Quaterniond(w, x, y, z));
	}

	/// @brief convert the dual quaternion into an homogeneous matrix
	/// @return 
	rendering::Transfo3d& to_transformation()
	{
		Vec3 t;
		float norm = _quat_0.norm();

		// Rotation matrix from non-dual quaternion part
		_quat_0.normalize(); 
		Eigen::Matrix3d m = _quat_0.toRotationMatrix();

		// translation vector from dual quaternion part:
		t[0] = 2.0 *
			  (-_quat_e.w() * _quat_0.x() + _quat_e.x() * _quat_0.w() - _quat_e.y() * _quat_0.z() +
			   _quat_e.z() * _quat_0.y()) /
			  norm;
		t[1] = 2.0 *
			  (-_quat_e.w() * _quat_0.y() + _quat_e.x() * _quat_0.z() + _quat_e.y() * _quat_0.w() -
			   _quat_e.z() * _quat_0.x()) /
			  norm;
		t[2] = 2.0 *
			  (-_quat_e.w() * _quat_0.z() - _quat_e.x() * _quat_0.y() + _quat_e.y() * _quat_0.x() +
			   _quat_e.z() * _quat_0.w()) /
			  norm;

		//return Transfo(m, t);
		rendering::Transfo3d myMatrix = Eigen::Affine3d::Identity();
		myMatrix.translation() = t; 
		myMatrix.linear() = m; 
		return myMatrix; 
	}

	DualQuaternion operator+(const DualQuaternion& dq) const
	{
		Eigen::Quaterniond non_dual_quat = Eigen::Quaterniond();
		Eigen::Quaterniond dual_quat = Eigen::Quaterniond();

		non_dual_quat.w() = _quat_0.w() + dq._quat_0.w(); 
		non_dual_quat.vec() = _quat_0.vec() + dq._quat_0.vec(); 

		dual_quat.w() = _quat_e.w() + dq._quat_e.w(); 
		dual_quat.vec() = _quat_e.vec() + dq._quat_e.vec();

		return DualQuaternion(non_dual_quat, dual_quat);
	}

	DualQuaternion operator*(Scalar scalar) const
	{
		Eigen::Quaterniond non_dual_part = _quat_0; 
		non_dual_part.coeffs()*=scalar; 

		Eigen::Quaterniond dual_part = _quat_e; 
		dual_part.coeffs()*=scalar; 

		return DualQuaternion(non_dual_part, dual_part);
	}

	/// Return a dual quaternion with no translation and no rotation
	static DualQuaternion identity()
	{
		return DualQuaternion(Eigen::Quaterniond(1.0, 0.0, 0.0, 0.0),
							Vec3(0.0, 0.0, 0.0) );
	}

	Eigen::Quaterniond get_dual_part() const { return _quat_e; }

	Eigen::Quaterniond get_non_dual_part() const { return _quat_0; }

	Eigen::Quaterniond translation() const { return _quat_e; }

	Eigen::Quaterniond rotation() const { return _quat_0; }

	void set_rotation( const Eigen::Quaterniond& q ){ _quat_0 = q; }



private:
	Eigen::Quaterniond _quat_0;
	Eigen::Quaterniond _quat_e;

};

} // namespace modeling

} // namespace cgogn

#endif // CGOGN_MODELING_TYPES_DUAL_QUATERNIONS_H_

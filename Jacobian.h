#pragma once
#include <array>
#include <vector>

#include <ManifestMath/Point3.h>
#include <ManifestMath/Quaternion.h>

using Manifest_Math::MFpoint3, Manifest_Math::MFvec3, Manifest_Math::MFquaternion, Manifest_Math::MFmat3;

namespace Manifest_Simulation
{
	//sorted with '>' to place friction contacts first
	enum class JacobianType
	{
		NORMAL,
		FIXED,
		FRICTION
	};
	MFbool operator>(const JacobianType& left,const JacobianType& right);		

	//represents the 1x12 row vector Ji, vi0 wi0 vi1 wi1 respectively
	using JacobianVector = std::array<MFvec3, 4>;		
	///computes Jconstraint to solve C*constraint = JV + b= 0 	
	JacobianVector ComputeNormalJacobian( const MFvec3 normal, const MFpoint3& r0, const MFpoint3& r1);
	JacobianVector ComputeFixedJacobian(const MFvec3& normal, const MFpoint3& localPoint);
	JacobianVector ComputeFrictionJacobian( const MFvec3& frictionVector, const MFpoint3& r0, const MFpoint3& r1);
}
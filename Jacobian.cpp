#include "Jacobian.h"

using namespace Manifest_Simulation;

MFbool Manifest_Simulation::operator>(const JacobianType& left, const JacobianType& right)
{
	return static_cast<MFu32>(left) > static_cast<MFu32>(right);
}

//returns J[-nT -(r0 x n)T nT (r1 x n)T] where,
//n = contact normal, rn = local contact point n
// Jacobian is stored as a row vector to implictily satisfy T
JacobianVector Manifest_Simulation::ComputeNormalJacobian( const MFvec3 normal, const MFpoint3& r0, const MFpoint3& r1)
{	
	const MFpoint3 torquePoint0{ Cross(normal, r0) };
	const MFpoint3 torquePoint1{ Cross(r1,normal) };

	return JacobianVector{-normal, torquePoint0, normal, torquePoint1 };
}

//returns J[-nT -(r0 x n)T 0 0] where,
//n = contact normal, r = local contact point 
//modified version of the normal constraint to 0 out the constant point being differentiated
JacobianVector Manifest_Simulation::ComputeFixedJacobian(const MFvec3& normal, const const MFpoint3& localPoint)
{
	const MFpoint3 torqueDirection{ Cross(normal,localPoint) };

	return JacobianVector{ -normal, torqueDirection,0,0 };
}

//returns J[-fT -(r0 x f)T fT (r1 x f)T] where,
//f = contact plane vector, rn = local contact point n
//modified version of the normal constraint to apply an impulse along the contact plane
JacobianVector Manifest_Simulation::ComputeFrictionJacobian(const MFvec3& frictionVector, const MFpoint3& r0, const MFpoint3& r1)
{
	const MFpoint3 slipDirection0{ Cross(frictionVector,r0) };
	const MFpoint3 slipDirection1{ Cross(r1,frictionVector) };
		
	return JacobianVector{ -frictionVector, slipDirection0,  frictionVector, slipDirection1 };
}
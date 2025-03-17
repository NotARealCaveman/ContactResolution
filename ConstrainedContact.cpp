#include "ConstrainedContact.h"

using namespace Manifest_Simulation;

MFvec3 Manifest_Simulation::ComputeRelativeContactVelocity(const ConstrainedContact& contact)
{
	const MFvec3 offset{ contact.normal * (contact.interpenetration * 0.5f) };

	const MFpoint3 offsetPoint0{ contact.contactPointsLocal[0] - offset };
	const MFpoint3 offsetPoint1{ contact.contactPointsLocal[1] + offset };

	const MFvec3 relativeVelocity0{ -(*contact.linearVelocityPtr[0]) - Cross(*contact.angularVelocityPtr[0],offsetPoint0)};
	const MFvec3 relativeVelocity1{ *contact.linearVelocityPtr[1] + Cross(*contact.angularVelocityPtr[1],offsetPoint1) };

	MFvec3 result{ relativeVelocity0 + relativeVelocity1 };

	for (MFu32 axis{ 0 }; axis < 3; ++axis)
		if (std::fabsf(result[axis]) < 0.0001f)
			result[axis] = 0;

	return result;
}

//computes default vectors based on the contacnt normal
//orthognalizes vectors with gram-schmit
void Manifest_Simulation::ComputeContactPlaneVectors(const MFvec3& relativeContactVelocity, MFvec3& normal, MFvec3& tangent, MFvec3& bitangent)
{ 		
	tangent =  relativeContactVelocity - normal * Dot(relativeContactVelocity, normal);
	MFfloat squaredMagnitude{ MagnitudeSquared(tangent) };	
	if (squaredMagnitude > 0.0001f)
	{
		tangent *= 1.0f / std::sqrtf(squaredMagnitude);
		bitangent = Normalize(Cross(tangent, normal));
	}
	else 
	{
		MFfloat x{ std::fabsf(normal.x) };
		MFfloat y{ std::fabsf(normal.y) };
		MFfloat z{ std::fabsf(normal.z) };

		if (x > y && x > z)
		{
			tangent = { 0,1,0 };
			bitangent = { 0,0,1 };
		}
		else if (y > z)
		{
			tangent = { 0,0,1 };
			bitangent = { 1,0,0 };
		}
		else
		{
			tangent = { 1,0,0 };
			bitangent = { 0,1,0 };
		}
	} 

	GramSchmidt(normal, tangent, bitangent);

}
#pragma once
#include <functional>

#include <ManifestPersistence/DatabaseTypes.h>
#include <ManifestMath/Vector3.h>
#include <ManifestMath/point3.h>
#include <ManifestMath/Plane.h>

#include <ManifestSimulation/PhysicsEngine/RigidBody.h>
#include "Jacobian.h"
#include "ContactManifold.h"//constraint types

using namespace Manifest_Math;
using Manifest_Persistence::PrimaryKey;

namespace Manifest_Simulation
{

	struct ConstrainedContact
	{	
		ConstraintType constraintType;		
		MFpoint3 contactPointsLocal[2];//local space contact point on body
		MFpoint3 contactObjectOrigin[2];//world space origin of contact objects
		MFvec3 normal;
		MFvec3 tangent;
		MFvec3 bitangent;
		MFvec3 relativeContactVelocity;
		InverseInertiaTensor* iTensorPtr[2];
		MFvec3* linearVelocityPtr[2];
		MFvec3* angularVelocityPtr[2];
		MFpoint3* positionPtr[2];
		MFquaternion* orientationPtr[2];
		MFfloat inverseObjectMass[2];//pure object mass
		MFfloat inverseContactMass[2];//divided object mass per contact point
		MFfloat interpenetration;//initial position error, C		
		MFfloat restitution;
		MFfloat friction;
		MFbool isDynamic[2];
	};
	//computes relative velocity of the contact bodies
	MFvec3 ComputeRelativeContactVelocity(const ConstrainedContact& contact);
	//computes orthonormal basis from the contact normal	
	//performs normalization of normal if magnitude is not satisfactory
	void ComputeContactPlaneVectors(const MFvec3& relativeContactVelocity, MFvec3& normal, MFvec3& tangent, MFvec3& bitangent);
}
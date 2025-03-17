#pragma once

#include <ManifestCore/Timer.h>

#include <ManifestSimulation/PhysicsEngine/PhysicsEngine.h>
#include <ManifestSimulation/EntityEngine/EntityEngine.h>
#include "ContactManifold.h"
#include "Constraint3D.h"

using namespace Manifest_Core;

namespace Manifest_Simulation
{
	//converts contact manifolds into ConstrainedContacts for resolution - returns the number of total contacts created
	void ConvertPhysicsContactManifolds(const PhysicsEngine& physicsEngine, const std::vector<ContactManifold>& contactManifolds, std::vector<ConstrainedContact>& contacts);
	void ConvertPhysicsManifoldData(const PhysicsEngine& physicsEngine, const PrimaryKey physicsID, const MFfloat massDivisor, MFtransform& inverseWorldSpace, MFmat3& inverseRotation, MFpoint3& position, InverseInertiaTensor*& iTensor, MFvec3*& linearVelocity, MFvec3*& angularVelocity, MFpoint3*& positionPtr, MFquaternion*& orientationPtr, MFfloat& inverseBodyMass, MFfloat& inverseContactMass, MFbool& isDynamic);
	//second phase - computes data unavailable during conversion		
	void PrepareContacts(std::vector<ConstrainedContact>& contacts);	
	void ResolveVelocityConstraints(const SparseJacobian& sparseJacobian, const InverseMassMatrix& inverseMassMatrix, const MFu32 MAX_ITERATIONS, const MFfloat dt, std::vector<Constraint>& constraintVector, std::vector<ConstrainedContact>& contacts);
	void ResolvePositionConstraints(const MFfloat dt, SparseJacobian& sparseJacobian, InverseMassMatrix& inverseMassMatrix, std::vector<Constraint>& constraintVector, std::vector<ConstrainedContact>& contacts);

	//visual debugging
	class ContactResolver
	{
	public:
		static std::atomic<MFvec3> frameNormal;
		static std::atomic<MFvec3> frameTangent;
		static std::atomic<MFvec3> frameBitangent;
		static std::atomic<MFpoint3> frameOrigin;
	};
}
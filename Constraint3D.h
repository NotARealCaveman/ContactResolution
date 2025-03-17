#pragma once
#include <ranges>
#include <unordered_set>

#include <ManifestSimulation/CollisionEngine/Colliders/Tensors.h>
#include "ConstrainedContact.h"

using namespace Manifest_Math;

namespace Manifest_Simulation
{
	struct JacobianEntry
	{	
		JacobianType jacobianType;
		ptrdiff_t contactIndex;
		JacobianVector jacobianVector;
	};
	using SparseJacobian = std::vector<JacobianEntry>;
	SparseJacobian ComputeSparseJacobian(const std::vector<ConstrainedContact>& contacts);

	struct InverseMassPair
	{
		MFfloat inverseContactMass;
		InverseInertiaTensor* iTensor;
	};
	using InverseMassVector = std::array<InverseMassPair, 2>;
	using InverseMassMatrix = std::vector<InverseMassVector>;
	InverseMassMatrix ComputeInverseMassMatrix(const std::vector<ConstrainedContact>& contacts);
		
	struct Constraint
	{			
		MFfloat positionalError;
		MFfloat velocityError;
		MFfloat effectiveMass;
		MFfloat totalImpulse;
		std::array<MFfloat, 2> impulseBounds;		
	};	
		
	std::vector<Constraint> ComputeConstraintVector(const SparseJacobian& sparseJacobian, const std::vector<ConstrainedContact>& contacts);
	//returns C*, where C* = JV
	MFfloat EvaluateGeneralVelocityConstraint(const ConstrainedContact& contact, const JacobianVector& jacobianVector);
	///specific position constraint evaluation functions
	MFfloat EvaluateNormalPositionConstraint(const ConstrainedContact& contact, const MFpoint3(&contactPointsLocal)[2], const MFvec3& normal, const MFfloat halfPositionError);
	MFfloat EvaluateFixedPositionConstraint(const ConstrainedContact& contact, const MFpoint3& localPoint, const MFpoint3& fixedPosition, const MFvec3& normal, const MFfloat positionError);
}
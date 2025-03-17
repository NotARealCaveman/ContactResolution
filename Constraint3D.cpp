#include "Constraint3D.h"

using namespace Manifest_Simulation;

//compute the sparse Jacobian of pairwise constrained bodies in the system
SparseJacobian Manifest_Simulation::ComputeSparseJacobian(const std::vector<ConstrainedContact>& contacts)
{
	SparseJacobian result;
	ConstrainedContact const* contactBegin{ contacts.data() };
	std::ranges::for_each(contacts, [&](const ConstrainedContact& contact)->void
		{
			const ptrdiff_t contactIndex{ std::distance(contactBegin,&contact) };
			switch (contact.constraintType)
			{			
				//NON ENTITY INVOLVED CONTACT CONSTRAINTS
			case ConstraintType::TERRAIN:
				result.emplace_back(JacobianEntry{ JacobianType::FRICTION, contactIndex, ComputeFrictionJacobian(contact.tangent,  contact.contactPointsLocal[0], 0) });
				result.emplace_back(JacobianEntry{ JacobianType::FRICTION, contactIndex, ComputeFrictionJacobian(contact.bitangent,  contact.contactPointsLocal[0], 0) });
				result.emplace_back(JacobianEntry{ JacobianType::FIXED, contactIndex,ComputeFixedJacobian(contact.normal, contact.contactPointsLocal[0]) });
				return;
			case ConstraintType::CONTACT:
				result.emplace_back(JacobianEntry{ JacobianType::FRICTION, contactIndex, ComputeFrictionJacobian(contact.tangent, contact.contactPointsLocal[0] , contact.contactPointsLocal[1]) });
				result.emplace_back(JacobianEntry{ JacobianType::FRICTION, contactIndex, ComputeFrictionJacobian(contact.bitangent, contact.contactPointsLocal[0] , contact.contactPointsLocal[1]) });
				result.emplace_back(JacobianEntry{ JacobianType::NORMAL, contactIndex, ComputeNormalJacobian(contact.normal, contact.contactPointsLocal[0], contact.contactPointsLocal[1]) });
				return;
			}
		});

	return result;
}
//compute the sparse Inverse Mass Matrix of pairwise constrained bodies in the system
InverseMassMatrix Manifest_Simulation::ComputeInverseMassMatrix(const std::vector<ConstrainedContact>& contacts)
{
	InverseMassMatrix result;
	std::ranges::transform(contacts, std::back_inserter(result),
		[&](const ConstrainedContact& contact)->InverseMassVector
		{
			InverseMassVector inverseMassVector;
			std::ranges::generate(inverseMassVector, [&, body = 0]() mutable ->InverseMassPair
				{
					const InverseMassPair inverseMassPair
					{
						 contact.inverseContactMass[body],
						 contact.iTensorPtr[body]
					};
					++body;
					return inverseMassPair;
				});
			return inverseMassVector;
		});

	return result;
}
 
 
//compute constraint vector C*=JV of each constrained body in the system
std::vector<Constraint>  Manifest_Simulation::ComputeConstraintVector(const SparseJacobian& sparseJacobian, const std::vector<ConstrainedContact>& contacts)
{
	const MFfloat TEMP_FRICTION{ 0.90f };

	const auto& ComputeBounds = [&](const JacobianType& jacobianType, const ptrdiff_t contactIndex)->std::array<MFfloat, 2>
	{
		switch (jacobianType)
		{
		case JacobianType::NORMAL: case JacobianType::FIXED:
			return { 0,std::numeric_limits<MFfloat>::infinity() };
		case JacobianType::FRICTION:
		{
			const MFfloat contactMass{ 1.0f/ (contacts[contactIndex].inverseContactMass[0] + contacts[contactIndex].inverseContactMass[1]) };
			const MFfloat gravitationFriction{ TEMP_FRICTION * 9.8f };
			const MFfloat applicableFriction{ contactMass*gravitationFriction };
			//sssssDLOG({ CONSOLE_ITALIC,CONSOLE_MAGENTA }, "Friction Force Bound:", applicableFriction);
			return { -applicableFriction,applicableFriction };
		}
		}

		return { 0.0f };
	};
	std::vector<Constraint> result(sparseJacobian.size());
	//evaluate each velocity constraint and compute potential violations
	JacobianEntry const* sparseBegin{ sparseJacobian.data() };
	std::ranges::transform(sparseJacobian, result.begin(), [&](const JacobianEntry& jacobianEntry)->Constraint
		{
			const ptrdiff_t contactIndex{ jacobianEntry.contactIndex };
			const JacobianVector& jacobianVector{ jacobianEntry.jacobianVector };

			//position error computed by contact detector
			const MFfloat positionalError{ contacts[contactIndex].interpenetration };
			//compute velocity error
			const MFfloat velocityError{ EvaluateGeneralVelocityConstraint(contacts[contactIndex],jacobianVector)};
			//compute inverse effective masses for each constraint using Meff = (JM-1JT)			


			const MFtransform mInverseJTranspose
			{//STRICTLY SPEAKING NOT A TRANSFORM BUT THE DATA LAYOUT WORKS
				jacobianVector[0] * contacts[contactIndex].inverseObjectMass[0],
				*contacts[contactIndex].iTensorPtr[0] * jacobianVector[1],
				jacobianVector[2] * contacts[contactIndex].inverseObjectMass[1],
				*contacts[contactIndex].iTensorPtr[1] * jacobianVector[3],
			};
			const MFfloat inverseEffectiveMass
			{
				Dot(jacobianVector[0],mInverseJTranspose[0]) +
				Dot(jacobianVector[1],mInverseJTranspose[1]) +
				Dot(jacobianVector[2],mInverseJTranspose[2]) +
				Dot(jacobianVector[3],mInverseJTranspose[3])
			};
			const MFfloat effectiveMass{ 1.0f / inverseEffectiveMass };		

			const std::array<MFfloat, 2> impulseBounds{ ComputeBounds(jacobianEntry.jacobianType,contactIndex) };

			return Constraint{ positionalError ,velocityError,effectiveMass, 0.0f,impulseBounds };
		});

	return result;
}

MFfloat Manifest_Simulation::EvaluateGeneralVelocityConstraint(const ConstrainedContact& contact, const JacobianVector& jacobianVector)
{
	//DLOG({ CONSOLE_RED }, "Body 0 Angular:", rigidBodies.angularVelocity[0]);

	//infinite masses/static bodies are set up in a way such that:
	//Jsp(i,n)V(bn) = 0, since Jsp(i,n) contains the [6x1]T block when multiplied by V(bn) which contains the [6x1] {0} block the result is 0
	const MFvec3 impulseDirection0{ jacobianVector[0] };
	const MFvec3 torqueDirection0{ jacobianVector[1] };
	const MFvec3 linearVelocity0{ *contact.linearVelocityPtr[0] };
	const MFvec3 angularVelocity0{ *contact.angularVelocityPtr[0] };
	MFfloat result = Dot(linearVelocity0, impulseDirection0) + Dot(angularVelocity0, torqueDirection0);
	//DLOG({ CONSOLE_DEFAULT }, "linearVelocity0", linearVelocity0, "angularVelocity0", angularVelocity0);
	const MFvec3 impulseDirection1{ jacobianVector[2] };
	const MFvec3 torqueDirection1{ jacobianVector[3] };
	const MFvec3 linearVelocity1{ *contact.linearVelocityPtr[1] };
	const MFvec3 angularVelocity1{ *contact.angularVelocityPtr[1] };
	result += Dot(linearVelocity1, impulseDirection1) + Dot(angularVelocity1, torqueDirection1);

	return result;
}
//returns the distance between two contacting bodies along the contact normal
//Cnrm = (x1+r1 - (x0+r0)) * n
MFfloat Manifest_Simulation::EvaluateNormalPositionConstraint(const ConstrainedContact& contact, const MFpoint3(&contactPointsLocal)[2], const MFvec3& normal, const MFfloat halfPositionError)
{
	const MFvec3 offset{ normal * halfPositionError };

	const MFpoint3 worldPoint0{ contact.contactObjectOrigin[0] + (contactPointsLocal[0]) };
	const MFpoint3 worldPoint1{ contact.contactObjectOrigin[1] + (contactPointsLocal[1]) };

	const MFpoint3 offsetPoint0{ worldPoint0 - offset };
	const MFpoint3 offsetPoint1{ worldPoint1 + offset };

	return Dot(offsetPoint1 - offsetPoint0, normal);
}

//modified version of the normal constraint to 0 out the constant point being differentiated
//returns the distance from a body to a fixed point along the contact normal
//Cnrm = (F - (x+r)) * n
MFfloat Manifest_Simulation::EvaluateFixedPositionConstraint(const ConstrainedContact& contact, const MFpoint3& localPoint, const MFpoint3& fixedPosition, const MFvec3& normal, const MFfloat positionError)
{
	const MFvec3 offset{ normal * positionError };

	const MFpoint3 worldPoint{ contact.contactObjectOrigin[0] + localPoint };
	const MFpoint3 bodyPoint{ worldPoint - offset };

	return Dot(fixedPosition - bodyPoint, normal);
}
#include "ContactResolver.h"

using namespace Manifest_Simulation;

//debug caputres
std::atomic<MFpoint3> ContactResolver::frameOrigin{};
std::atomic<MFvec3> ContactResolver::frameNormal{};
std::atomic<MFvec3> ContactResolver::frameTangent{};
std::atomic<MFvec3> ContactResolver::frameBitangent{};

void Manifest_Simulation::PrepareContacts(std::vector<ConstrainedContact>& contacts)
{
	std::ranges::for_each(contacts, [&](ConstrainedContact& contact)
		{
			contact.relativeContactVelocity = ComputeRelativeContactVelocity(contact);
			ComputeContactPlaneVectors(contact.relativeContactVelocity, contact.normal, contact.tangent, contact.bitangent);
		});
}

void Manifest_Simulation::ResolveVelocityConstraints(const SparseJacobian& sparseJacobian, const InverseMassMatrix& inverseMassMatrix, const MFu32 MAX_ITERATIONS, const MFfloat dt, std::vector<Constraint>& constraintVector, std::vector<ConstrainedContact>& contacts)
{
	constexpr MFfloat permissiblePenetration{ 0.005f};
	 
	Constraint* constraintBegin{ constraintVector.data() }; 
	//friction only first	
	for (MFu32 iteration{ 0 }; iteration < MAX_ITERATIONS; ++iteration)
		std::ranges::for_each(constraintVector, [&](Constraint& constraint)->void
			{
				const ptrdiff_t constraintIndex{ std::distance(constraintBegin,&constraint) };
				const ptrdiff_t contactIndex{ sparseJacobian[constraintIndex].contactIndex };

				//DLOG({ CONSOLE_RED }, "IDs:", contacts[contactIndex].physicsIDs[0], contacts[contactIndex].physicsIDs[1]);

				const JacobianVector& jacobianVector{ sparseJacobian[constraintIndex].jacobianVector };
				if (sparseJacobian[constraintIndex].jacobianType != JacobianType::FRICTION)
					return;
				const MFfloat velocityError{ EvaluateGeneralVelocityConstraint(contacts[contactIndex],jacobianVector)};
				//compute impulse λ to satisfy JV = 0, update total impulse applied 
				MFfloat& totalImpulse{ constraint.totalImpulse };
				MFfloat oldImpulse{ totalImpulse };
				//λ = -Meff * C* + b, C*<0 <=> JVT<0
				//DLOG({ CONSOLE_BLUE}, "b:", b);
				MFfloat impulse{ constraint.effectiveMass * -velocityError };
				totalImpulse = std::clamp(totalImpulse + impulse, constraint.impulseBounds[0], constraint.impulseBounds[1]);
				impulse = totalImpulse - oldImpulse;
				CONSOLE_CODE code{ &(CONSOLE_RESET)[iteration] };
				//DLOG({ code }, "pre impulse:",constraint.effectiveMass * (b -velocityError ), "total impulse:", totalImpulse, "new impulse:", impulse,"velocity error:",velocityError);
				//calculate ΔV by applying impulse λ to satisfy JV
				//ΔV = M-1JTλ, J(V + ΔV) = 0 :: J(V + M-1JTλ)= 0		
				//V0 + ΔV				
				const MFvec3 linearImpulse0{ jacobianVector[0] * impulse * inverseMassMatrix[contactIndex][0].inverseContactMass };
				const MFvec3 angularImpulse0{ *inverseMassMatrix[contactIndex][0].iTensor * jacobianVector[1] * impulse };
				*contacts[contactIndex].linearVelocityPtr[0] += linearImpulse0;
				*contacts[contactIndex].angularVelocityPtr[0] += angularImpulse0;
				//DLOG({ CONSOLE_BLUE }, "velocity contact:", contactIndex, "velocity error:", velocityError, "base impulse:", impulse, "totalImpulse:", totalImpulse, "linear0:", linearImpulse0); 
				//DLOG({ CONSOLE_GREEN }, "linearVelocity", rigidBodies.linearVelocity[bodyIndex0], "angular0:", angularImpulse0, "angularVelocity:", rigidBodies.angularVelocity[bodyIndex0]);
				//V1 + ΔV				
				const MFvec3 linearImpulse1{ jacobianVector[2] * impulse * inverseMassMatrix[contactIndex][1].inverseContactMass };
				const MFvec3 angularImpulse1{ *inverseMassMatrix[contactIndex][1].iTensor * jacobianVector[3] * impulse };
				*contacts[contactIndex].linearVelocityPtr[1] += linearImpulse1;
				*contacts[contactIndex].angularVelocityPtr[1] += angularImpulse1;
			});
	//then normal
	for (MFu32 iteration{ 0 }; iteration < MAX_ITERATIONS; ++iteration)
		std::ranges::for_each(constraintVector, [&](Constraint& constraint)->void
			{				
				const ptrdiff_t constraintIndex{ std::distance(constraintBegin,&constraint) };
				const ptrdiff_t contactIndex{ sparseJacobian[constraintIndex].contactIndex };		

				//DLOG({ CONSOLE_RED }, "IDs:", contacts[contactIndex].physicsIDs[0], contacts[contactIndex].physicsIDs[1]);

				const JacobianVector& jacobianVector{ sparseJacobian[constraintIndex].jacobianVector };
				//skip friction and entity based constrains
				if (sparseJacobian[constraintIndex].jacobianType == JacobianType::FRICTION )
					return;
				const MFfloat velocityError{ EvaluateGeneralVelocityConstraint(contacts[contactIndex],jacobianVector) };  				
				//compute impulse λ to satisfy JV = 0, update total impulse applied 
				MFfloat& totalImpulse{ constraint.totalImpulse };
				MFfloat oldImpulse{ totalImpulse };
				//λ = -Meff * C*, C*<0 <=> JVT<0											
				MFfloat impulse{ -constraint.effectiveMass * velocityError };
				totalImpulse = std::clamp(totalImpulse + impulse, constraint.impulseBounds[0], constraint.impulseBounds[1]);
				impulse = totalImpulse - oldImpulse;
				//DLOG({ CONSOLE_MAGENTA }, "pre impulse:",constraint.effectiveMass * (b -velocityError ), "total impulse:", totalImpulse, "new impulse:", impulse,"velocity error:",velocityError);
				//calculate ΔV by applying impulse λ to satisfy JV
				//ΔV = M-1JTλ, J(V + ΔV) = 0 :: J(V + M-1JTλ)= 0		
				//V0 + ΔV				
				const MFvec3 linearImpulse0{ jacobianVector[0] * impulse * inverseMassMatrix[contactIndex][0].inverseContactMass };
				const MFvec3 angularImpulse0{ *inverseMassMatrix[contactIndex][0].iTensor * jacobianVector[1] * impulse }; 
				//DLOG({ CONSOLE_DEFAULT },"Variables: idx:", constraintIndex, "iTensor:", *(inverseMassMatrix[contactIndex][0].iTensor), "J1:", jacobianVector[1], "l:", impulse);
				*contacts[contactIndex].linearVelocityPtr[0] += linearImpulse0;
				*contacts[contactIndex].angularVelocityPtr[0] += angularImpulse0;
				//DLOG({ CONSOLE_BLUE }, "velocity contact:", contactIndex, "velocity error:", velocityError, "base impulse:", impulse, "totalImpulse:", totalImpulse, "linear0:", linearImpulse0,"position error:",constraint.positionalError);
				//DLOG({ CONSOLE_GREEN }, "linearVelocity", *contacts[contactIndex].linearVelocityPtr[0], "angular0:", angularImpulse0, "angularVelocity:", *contacts[contactIndex].angularVelocityPtr[0]);
				//V1 + ΔV				
				const MFvec3 linearImpulse1{ jacobianVector[2] * impulse * inverseMassMatrix[contactIndex][1].inverseContactMass };
				const MFvec3 angularImpulse1{ *inverseMassMatrix[contactIndex][1].iTensor * jacobianVector[3] * impulse };
				*contacts[contactIndex].linearVelocityPtr[1] += linearImpulse1;
				*contacts[contactIndex].angularVelocityPtr[1] += angularImpulse1;
			});
}

void Manifest_Simulation::ResolvePositionConstraints(const MFfloat dt, SparseJacobian& sparseJacobian, InverseMassMatrix& inverseMassMatrix, std::vector<Constraint>& constraintVector, std::vector<ConstrainedContact>& contacts)
{
	constexpr MFfloat permissiblePenetration{ 0.005f };

	Constraint* constraintBegin{ constraintVector.data() };	
		std::ranges::for_each(constraintVector, [&](Constraint& constraint)->void
			{
				const ptrdiff_t constraintIndex{ std::distance(constraintBegin,&constraint) };
				const ptrdiff_t contactIndex{ sparseJacobian[constraintIndex].contactIndex };
				JacobianVector& jacobianVector{ sparseJacobian[constraintIndex].jacobianVector };
				if (sparseJacobian[constraintIndex].jacobianType == JacobianType::FRICTION)
					return;				
				constexpr MFfloat BETA{ 0.15f };
				MFfloat b = BETA * std::fminf(0.0f, (constraint.positionalError + permissiblePenetration));
				MFfloat impulse{ -constraint.effectiveMass * b };							
				const MFvec3 linearImpulse0{ jacobianVector[0] * impulse * inverseMassMatrix[contactIndex][0].inverseContactMass };
				const MFvec3 angularImpulse0{ *inverseMassMatrix[contactIndex][0].iTensor * jacobianVector[1] * impulse };								
				*contacts[contactIndex].positionPtr[0] += linearImpulse0;				
				MFquaternion& orientation0{ *contacts[contactIndex].orientationPtr[0] };
				orientation0 = Normalize(orientation0 + MFquaternion{ angularImpulse0 ,0.0 }*orientation0 * 0.5f * dt);

				const MFvec3 linearImpulse1{ jacobianVector[2] * impulse * inverseMassMatrix[contactIndex][1].inverseContactMass };
				const MFvec3 angularImpulse1{ *inverseMassMatrix[contactIndex][1].iTensor * jacobianVector[3] * impulse };
				*contacts[contactIndex].positionPtr[1] += linearImpulse1;
				MFquaternion& orientation1{ *contacts[contactIndex].orientationPtr[1] };
				orientation1 = Normalize(orientation1 + MFquaternion{ angularImpulse1 ,0.0 }*orientation1 * 0.5f * dt);

				//DLOG({ CONSOLE_BG_GREEN }, "linearImpulse0", linearImpulse0, "impulse",impulse,"constraintIndex", constraintIndex);

				ConstrainedContact& contact{ contacts[contactIndex] };				
				switch (sparseJacobian[constraintIndex].jacobianType)
				{
				case JacobianType::NORMAL:
					constraint.positionalError = EvaluateGeneralVelocityConstraint(contacts[contactIndex], jacobianVector);					
					break;
				case JacobianType::FIXED:
					constraint.positionalError = EvaluateGeneralVelocityConstraint(contacts[contactIndex], jacobianVector);					
					break;
				default:
					return;
				}				
				/*
				const MFtransform mInverseJTranspose
				{
					jacobianVector[0] * contacts[contactIndex].inverseContactMass[0],
					*contacts[contactIndex].iTensorPtr[0] * jacobianVector[1],
					jacobianVector[2] * contacts[contactIndex].inverseContactMass[1],
					*contacts[contactIndex].iTensorPtr[1] * jacobianVector[3],
				};
				const MFfloat inverseEffectiveMass
				{
					Dot(jacobianVector[0],mInverseJTranspose[0]) +
					Dot(jacobianVector[1],mInverseJTranspose[1]) +
					Dot(jacobianVector[2],mInverseJTranspose[2]) +
					Dot(jacobianVector[3],mInverseJTranspose[3])
				};
				constraint.effectiveMass = 1.0f / inverseEffectiveMass;
				*/
			});
}

void Manifest_Simulation::ConvertPhysicsContactManifolds(const PhysicsEngine& physicsEngine, const std::vector<ContactManifold>& contactManifolds, std::vector<ConstrainedContact>& contacts)
{
	std::ranges::for_each(contactManifolds, [&](const ContactManifold& contactManifold)
		{
			MFtransform inverseWorldSpace[2];
			MFmat3 inverseRotation[2];
			MFpoint3 position[2];						
			MFbool isDynamic[2];
			MFfloat inverseObjectMass[2];
			MFfloat inverseContactMass[2];
			InverseInertiaTensor* iTensor[2];
			MFvec3* linearVelocity[2];
			MFvec3* angularVelocity[2];
			MFpoint3* positionPtr[2];
			MFquaternion* orientationPtr[2];
			//C = a/b, C^-1 = b/a = b * a^-1 where C = contact mass, a = body mass, and b = n contact points
			const MFsize massDivisor{ contactManifold.contactPoints.size() };
			PrimaryKey ID{ contactManifold.collisionPair.a->physicsID };			
			for (MFu32 contactObject{ 0 }; contactObject < 2; ++contactObject)
			{								
				ConvertPhysicsManifoldData(physicsEngine, ID, massDivisor, inverseWorldSpace[contactObject], inverseRotation[contactObject], position[contactObject], iTensor[contactObject], linearVelocity[contactObject], angularVelocity[contactObject], positionPtr[contactObject],orientationPtr[contactObject], inverseObjectMass[contactObject], inverseContactMass[contactObject], isDynamic[contactObject]);				
				ID = contactManifold.collisionPair.b->physicsID;
			}
			

			std::ranges::transform(contactManifold.contactPoints, std::back_inserter(contacts), [&](const ContactPoint& manifoldContact)->ConstrainedContact
				{
					ConstrainedContact result;
					//copy tracking data					
					result.isDynamic[0] = isDynamic[0];
					result.isDynamic[1] = isDynamic[1];
					//copy physics data													
					result.normal = contactManifold.normal;
					result.linearVelocityPtr[0] = linearVelocity[0];
					result.linearVelocityPtr[1] = linearVelocity[1];
					result.angularVelocityPtr[0] = angularVelocity[0];
					result.angularVelocityPtr[1] = angularVelocity[1];
					result.iTensorPtr[0] = iTensor[0];
					result.iTensorPtr[1] = iTensor[1];
					result.positionPtr[0] = positionPtr[0];
					result.positionPtr[1] = positionPtr[1];	
					result.orientationPtr[0] = orientationPtr[0];
					result.orientationPtr[1] = orientationPtr[1];							
					result.constraintType = contactManifold.constraintType;
					result.interpenetration = manifoldContact.interpenetration;
					//convert contact points
					result.contactObjectOrigin[0] = position[0];
					result.contactPointsLocal[0] = (manifoldContact.collisionPointWorldSpace - result.contactObjectOrigin[0]);
					result.contactObjectOrigin[1] = position[1];
					result.contactPointsLocal[1] = (manifoldContact.collisionPointWorldSpace - result.contactObjectOrigin[1]);
					//divide mass across each contact point per body - [Catto05]
					result.inverseContactMass[0] = inverseContactMass[0];
					result.inverseContactMass[1] = inverseContactMass[1];
					result.inverseObjectMass[0] = inverseObjectMass[0];
					result.inverseObjectMass[1] = inverseObjectMass[1];

					//clean up detected values - clamp to avoid numerica 0l drift	 
					constexpr MFfloat CLAMP_EPSILON{ 1e-4 };
					Clamp_EPSILON(result.normal, CLAMP_EPSILON, 0.0f);
					Clamp_EPSILON(result.contactPointsLocal[0], CLAMP_EPSILON, 0.0f);
					Clamp_EPSILON(result.contactPointsLocal[1], CLAMP_EPSILON, 0.0f);

					return result;
				});
		});
}

void Manifest_Simulation::ConvertPhysicsManifoldData(const PhysicsEngine& physicsEngine, const PrimaryKey physicsID, const MFfloat massDivisor, MFtransform& inverseWorldSpace, MFmat3& inverseRotation, MFpoint3& position, InverseInertiaTensor*& iTensor, MFvec3*& linearVelocity, MFvec3*& angularVelocity, MFpoint3*& positionPtr, MFquaternion*& orientationPtr, MFfloat& inverseBodyMass, MFfloat& inverseContactMass, MFbool& isDynamic)
{
	isDynamic = physicsEngine.IsBodyDynamic(physicsID);
	const MFu32 bodyIndex = physicsEngine.GetBodyIndex(physicsID);

	const RigidBodyData& rigidBodies{ physicsEngine.rigidBodies };

	position = rigidBodies.GetPosition(bodyIndex);
	//store const compute data
	inverseWorldSpace = Inverse(rigidBodies.GetWorldSpace(bodyIndex));
	inverseRotation = Transpose(rigidBodies.GetOrientation(bodyIndex).GetRotation());
	//body might be static but given internal mass
	inverseBodyMass = rigidBodies.GetIMass(bodyIndex) * isDynamic;
	inverseContactMass = inverseBodyMass * massDivisor;	
	//store pointer update data	
	iTensor = &const_cast<InverseInertiaTensor&>(rigidBodies.GetWorldITensor(bodyIndex));
	linearVelocity = &const_cast<MFvec3&>(rigidBodies.GetLinearVelocity(bodyIndex));
	angularVelocity = &const_cast<MFvec3&>(rigidBodies.GetAngularVelocity(bodyIndex));
	positionPtr = &const_cast<MFpoint3&>(rigidBodies.GetPosition(bodyIndex));
	orientationPtr = &const_cast<MFquaternion&>(rigidBodies.GetOrientation(bodyIndex));	 
}

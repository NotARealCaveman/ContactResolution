#pragma once
#include "Collider.h"

using namespace Manifest_Math;

namespace Manifest_Simulation
{
	struct CollisionPair
	{
		Collider const* a;
		Collider const* b;				
		std::vector<MFtriangle*> triangleMesh;//used for processing mesh contacts
		MFbool operator==(const CollisionPair& other) const
		{
			if (other.a == a && other.b == b)
				return true;
			else if (other.b == a && other.a == b)
				return true;

			return false;
		}
	};
	//contact in world space
	struct ContactPoint
	{
		MFpoint3 collisionPointWorldSpace;
		MFfloat interpenetration;
	};

	//not a fan of this but it works for now
	enum class ConstraintType : MFu32
	{//terrain/plane probably to be renamed to cover cases like static meshes 
		TERRAIN,//physics vs terrain
		CONTACT,//physics vs physics				
	};
	//set of contacts in world space
	struct ContactManifold
	{		
		CollisionPair collisionPair;
		std::vector<ContactPoint> contactPoints;		
		MFvec3 normal;				
		ConstraintType constraintType;
	};	

	struct FrameCollisions
	{	
		std::vector<CollisionPair> collisionPairs;				
		std::vector<ContactManifold> contactManifolds;			
	};
}
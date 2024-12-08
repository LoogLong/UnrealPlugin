// Copyright LoogLong. All Rights Reserved.
#include "AnimNode_SpringMagic.h"
#include "AnimationRuntime.h"
#include "Animation/AnimInstanceProxy.h"
#include "Curves/CurveFloat.h"
#include "Runtime/Launch/Resources/Version.h"
#include "SceneInterface.h"
#include "PhysicsEngine/PhysicsAsset.h"


DECLARE_CYCLE_STAT(TEXT("SpringMagic_Eval"), STAT_SpringMagic_Eval, STATGROUP_Anim);

namespace SpringMagicIntersectionTest
{
	// take from https://gist.github.com/jdryg/ecde24d34aa0ce2d4d87

	/**
	 * return true if intersection happened
	 */
	bool IntersectRaySphere(const FRay& Ray, const FVector& SphereCenter, const float& SphereRadius, float& OutTMin, float& OutTMax)
	{
		FVector CO = Ray.Origin - SphereCenter;

		const float A = Ray.Direction.Dot(Ray.Direction);
		const float B = 2.0f * CO.Dot(Ray.Direction);
		const float C = CO.Dot(CO) - (SphereRadius * SphereRadius);

		float Discriminant = B * B - 4.0f * A * C;
		if (Discriminant < 0.0f)
			return false;

		OutTMin = (-B - FMath::Sqrt(Discriminant)) / (2.0f * A);
		OutTMax = (-B + FMath::Sqrt(Discriminant)) / (2.0f * A);
		if (OutTMin > OutTMax)
		{
			float temp = OutTMin;
			OutTMin       = OutTMax;
			OutTMax       = temp;
		}

		return true;
	}

	bool IntersectRayCapsule(const FRay& Ray, const FVector& CapsuleA, const FVector& CapsuleB, const float& CapsuleRadius, FVector& OutP1, FVector& OutP2, FVector& OutN1, FVector& OutN2, float& OutTMin, float& OutTMax)
	{
		// http://pastebin.com/2XrrNcxb

		// Substituting equ. (1) - (6) to equ. (I) and solving for t' gives:
		//
		// t' = (t * dot(AB, d) + dot(AB, AO)) / dot(AB, AB); (7) or
		// t' = t * m + n where 
		// m = dot(AB, d) / dot(AB, AB) and 
		// n = dot(AB, AO) / dot(AB, AB)
		//
		FVector AB = CapsuleB - CapsuleA;
		FVector AO = Ray.Origin - CapsuleA;

		float AB_dot_d  = AB.Dot(Ray.Direction);
		float AB_dot_AO = AB.Dot(AO);
		float AB_dot_AB = AB.Dot(AB);

		const float M = AB_dot_d / AB_dot_AB;
		const float N = AB_dot_AO / AB_dot_AB;

		// Substituting (7) into (II) and solving for t gives:
		//
		// dot(Q, Q)*t^2 + 2*dot(Q, R)*t + (dot(R, R) - r^2) = 0
		// where
		// Q = d - AB * m
		// R = AO - AB * n
		FVector Q = Ray.Direction - (AB * M);
		FVector R = AO - (AB * N);

		float A = Q.Dot(Q);
		float B = 2.0f * Q.Dot(R);
		float C = R.Dot(R) - (CapsuleRadius * CapsuleRadius);

		if (A == 0.0f)
		{
			// Special case: AB and ray direction are parallel. If there is an intersection it will be on the end spheres...
			// NOTE: Why is that?
			// Q = d - AB * m =>
			// Q = d - AB * (|AB|*|d|*cos(AB,d) / |AB|^2) => |d| == 1.0
			// Q = d - AB * (|AB|*cos(AB,d)/|AB|^2) =>
			// Q = d - AB * cos(AB, d) / |AB| =>
			// Q = d - unit(AB) * cos(AB, d)
			//
			// |Q| == 0 means Q = (0, 0, 0) or d = unit(AB) * cos(AB,d)
			// both d and unit(AB) are unit vectors, so cos(AB, d) = 1 => AB and d are parallel.
			// 

			float ATMin, ATMax, BTMin, BTMax;
			if (!IntersectRaySphere(Ray, CapsuleA, CapsuleRadius, ATMin, ATMax) ||
				!IntersectRaySphere(Ray, CapsuleB, CapsuleRadius, BTMin, BTMax))
			{
				// No intersection with one of the spheres means no intersection at all...
				return false;
			}

			if (ATMin < BTMin)
			{
				OutP1 = Ray.Origin + (Ray.Direction * ATMin);
				OutN1 = OutP1 - CapsuleA;
				OutTMin = ATMin;
				OutN1.Normalize();
			}
			else
			{
				OutP1 = Ray.Origin + (Ray.Direction * BTMin);
				OutN1 = OutP1 - CapsuleB;
				OutTMin = BTMin;
				OutN1.Normalize();
			}

			if (ATMax > BTMax)
			{
				OutP2 = Ray.Origin + (Ray.Direction * ATMax);
				OutN2 = OutP2 - CapsuleA;
				OutTMax = ATMax;
				OutN2.Normalize();
			}
			else
			{
				OutP2 = Ray.Origin + (Ray.Direction * BTMax);
				OutN2 = OutP2 - CapsuleB;
				OutTMax = BTMax;
				OutN2.Normalize();
			}

			return true;
		}

		float Discriminant = B * B - 4.0f * A * C;
		if (Discriminant < 0.0f)
		{
			// The ray doesn't hit the infinite cylinder defined by (A, B).
			// No intersection.
			return false;
		}

		float TMin = (-B - FMath::Sqrt(Discriminant)) / (2.0f * A);
		float TMax = (-B + FMath::Sqrt(Discriminant)) / (2.0f * A);
		if (TMin > TMax)
		{
			float Temp = TMin;
			TMin       = TMax;
			TMax       = Temp;
		}

		// Now check to see if K1 and K2 are inside the line segment defined by A,B
		float TK1 = TMin * M + N;
		if (TK1 < 0.0f)
		{
			// On sphere (A, r)...
			float STMin, STMax;
			if (IntersectRaySphere(Ray, CapsuleA, CapsuleRadius, STMin, STMax))
			{
				OutP1 = Ray.Origin + (Ray.Direction * STMin);
				OutN1 = OutP1 - CapsuleA;
				OutTMin = STMin;
				OutN1.Normalize();
			}
			else
				return false;
		}
		else if (TK1 > 1.0f)
		{
			// On sphere (B, r)...

			float STMin, STMax;
			if (IntersectRaySphere(Ray, CapsuleB, CapsuleRadius, STMin, STMax))
			{
				OutP1 = Ray.Origin + (Ray.Direction * STMin);
				OutN1 = OutP1 - CapsuleB;
				OutTMin = STMin;
				OutN1.Normalize();
			}
			else
				return false;
		}
		else
		{
			// On the cylinder...
			OutP1 = Ray.Origin + (Ray.Direction * TMin);
			OutTMin = TMin;
			FVector K1 = CapsuleA + AB * TK1;
			OutN1           = OutP1 - K1;
			OutN1.Normalize();
		}

		float TK2 = TMax * M + N;
		if (TK2 < 0.0f)
		{
			// On sphere (A, r)...

			float STMin, STMax;
			if (IntersectRaySphere(Ray, CapsuleA, CapsuleRadius, STMin, STMax))
			{
				OutP2 = Ray.Origin + (Ray.Direction * STMax);
				OutN2 = OutP2 - CapsuleA;
				OutTMax = STMax;
				OutN2.Normalize();
			}
			else
				return false;
		}
		else if (TK2 > 1.0f)
		{
			// On sphere (B, r)...

			float STMin, STMax;
			if (IntersectRaySphere(Ray, CapsuleB, CapsuleRadius, STMin, STMax))
			{
				OutP2 = Ray.Origin + (Ray.Direction * STMax);
				OutN2 = OutP2 - CapsuleB;
				OutTMax = STMax;
				OutN2.Normalize();
			}
			else
				return false;
		}
		else
		{
			OutP2 = Ray.Origin + (Ray.Direction * TMax);
			OutTMax = TMax;
			FVector k2 = CapsuleA + AB * TK2;
			OutN2           = OutP2 - k2;
			OutN2.Normalize();
		}

		return true;
	}
}

FAnimNode_SpringMagic::FAnimNode_SpringMagic()
{
}

void FAnimNode_SpringMagic::Initialize_AnyThread(const FAnimationInitializeContext& Context)
{
	Super::Initialize_AnyThread(Context);
	const FBoneContainer& RequiredBones = Context.AnimInstanceProxy->GetRequiredBones();
	InitializeBoneReferences(RequiredBones);
}

void FAnimNode_SpringMagic::CacheBones_AnyThread(const FAnimationCacheBonesContext& Context)
{
	Super::CacheBones_AnyThread(Context);
}

void FAnimNode_SpringMagic::ResetDynamics(ETeleportType InTeleportType)
{
	PendingTeleportType = InTeleportType;
}

void FAnimNode_SpringMagic::UpdateInternal(const FAnimationUpdateContext& Context)
{
	Super::UpdateInternal(Context);
}

void FAnimNode_SpringMagic::GatherDebugData(FNodeDebugData& DebugData)
{
	Super::GatherDebugData(DebugData);
}

FQuat SolveAimWithTwist(const FTransform& CurrentTransform, const FVector& TargetPosition, const FVector& AimVector, const FVector& TargetUpVector,const FVector& UpVector, float AimClampInDegree)
{
	if (!ensure(AimVector.IsNormalized()) || !ensure(TargetUpVector.IsNormalized()) || !ensure(UpVector.IsNormalized()))
	{
		return FQuat::Identity;
	}

	FVector ToTarget = TargetPosition - CurrentTransform.GetLocation();
	ToTarget.Normalize();

	if (AimClampInDegree > ZERO_ANIMWEIGHT_THRESH)
	{
		double AimClampInRadians = FMath::DegreesToRadians(FMath::Min(AimClampInDegree, 180.0));
		double DiffAngle = FMath::Acos(FVector::DotProduct(AimVector, ToTarget));

		if (DiffAngle > AimClampInRadians)
		{
			check(DiffAngle > 0.0);

			FVector DeltaTarget = ToTarget - AimVector;
			// clamp delta target to within the ratio
			DeltaTarget *= (AimClampInRadians / DiffAngle);
			// set new target
			ToTarget = AimVector + DeltaTarget;
			ToTarget.Normalize();
		}
	}
	FQuat Swing = FQuat::FindBetweenNormals(AimVector, ToTarget);

	FQuat NewRotation = Swing * CurrentTransform.GetRotation();
	FVector AimUpVector = NewRotation.RotateVector(UpVector);
	FVector AimTargetVector = FVector::VectorPlaneProject(TargetUpVector, ToTarget);
	AimTargetVector.Normalize();
	
	FQuat Twist = FQuat::FindBetweenNormals(AimUpVector, AimTargetVector);
	return Twist * Swing;
}

void FAnimNode_SpringMagic::EvaluateSkeletalControl_AnyThread(FComponentSpacePoseContext& Output,
                                                                TArray<FBoneTransform>& OutBoneTransforms)
{
	SCOPE_CYCLE_COUNTER(STAT_SpringMagic_Eval);
	const FTransform& ComponentTransform = Output.AnimInstanceProxy->GetComponentTransform();
	float DeltaTime = Output.AnimInstanceProxy->GetDeltaSeconds();
	if (PendingTeleportType == ETeleportType::TeleportPhysics)
	{
		
	}else if (PendingTeleportType == ETeleportType::ResetPhysics)
	{
		bNeedInitializeSimulation = true;
		PendingTeleportType = ETeleportType::None;
	}
	if (bNeedInitializeSimulation)
	{
		InitializeSimulation(Output, Output.AnimInstanceProxy->GetRequiredBones());
		if (bEnableCollision)
		{
			InitializeCollision(Output.AnimInstanceProxy->GetRequiredBones());
		}
		PreComponentTransform = ComponentTransform;
		bNeedInitializeSimulation = false;
	}

	const bool bNeedSimulate = DeltaTime > UE_KINDA_SMALL_NUMBER;
	if (bNeedSimulate)
	{
		{
			// update root node base on animation pose
			// Update Simulate Node Local Pose
			for (auto& Element : SimulationParticles)
			{
				if (Element.bRootBone)
				{
					Element.SimulationTransform = Output.Pose.GetComponentSpaceTransform(Element.BoneIndex);
				}
				else if (!Element.bDummyBone)
				{
					Element.LocalTransform = Output.Pose.GetLocalSpaceTransform(Element.BoneIndex);
				}
			}
			for (TArray<FSpringMagicJoint>& JointChain : SimulationJoints)
			{
				for (FSpringMagicJoint& Joint : JointChain)
				{
					Joint.ChildProxy = SimulationParticles[Joint.Child].LocalTransform * SimulationParticles[Joint.Parent].LocalTransform;
				}
			}

			if (bEnableCollision)
			{
				UpdateCollision(Output, Output.AnimInstanceProxy->GetRequiredBones());
			}
		}

		FVector ExternalForce = FVector::Zero();
		if (bUseWind)
		{
			SinCounter += DeltaTime * WindFrequency;
			SinCounter = FMath::Fmod(SinCounter, UE_PI);
			FVector WindDirSimSpace = ComponentTransform.InverseTransformVectorNoScale(WindDirection);
			ExternalForce += (FMath::Sin(SinCounter) * WindAmplitude + WindAmplitudeBias) * WindDirSimSpace;
		}

		FVector CompMoveVector = FVector::Zero();
		FQuat CompMoveRotation = FQuat::Identity;
		if (bUseComponentMovement)
		{
			CompMoveVector = ComponentTransform.InverseTransformPosition(PreComponentTransform.GetLocation()) * ComponentMovementTranslation;
			CompMoveRotation = ComponentTransform.GetRotation().Inverse() * PreComponentTransform.GetRotation();
		}
		// iterate from top to end
		for (int32 LayerIdx = 0; LayerIdx < MaxChainLength; ++LayerIdx)
		{
			for (TArray<FSpringMagicJoint>& JointChain : SimulationJoints)
			{
				if (LayerIdx < JointChain.Num())
				{
					auto& GrandParent = JointChain[LayerIdx].GrandParent;
					auto& Parent = JointChain[LayerIdx].Parent;
					auto& Child = JointChain[LayerIdx].Child;

					SimulationParticles[Parent].SimulationTransform = SimulationParticles[Parent].LocalTransform * SimulationParticles[GrandParent].SimulationTransform;
					FTransform ChildProxyTransform = JointChain[LayerIdx].ChildProxy * SimulationParticles[GrandParent].SimulationTransform;

					FVector NewChildPosition = ChildProxyTransform.GetTranslation();
					FVector CurrentChildPosition = SimulationParticles[Child].SimulationTransform.GetTranslation();

					//ComponentMovement
					if (bUseComponentMovement)
					{
						CurrentChildPosition += CompMoveVector + (CompMoveRotation.RotateVector(CurrentChildPosition) - CurrentChildPosition) * ComponentMovementRotation;
					}

					// Add Wind
					NewChildPosition += ExternalForce;

					// Add Inertia

					// Detected Collision
					if (bEnableCollision)
					{
						CollisionDetection(SimulationParticles[Child], NewChildPosition, CurrentChildPosition);
					}

					JointChain[LayerIdx].NewChildPosition = NewChildPosition;
					JointChain[LayerIdx].CurrentChildPosition = CurrentChildPosition;

					const FVector& DefaultAimVector = FVector::XAxisVector;
					const FVector& DefaultUpVector = FVector::YAxisVector;
					constexpr float AimClampInDegree = 90.f;

					// Twist
					FVector CurUPVector = ChildProxyTransform.GetRotation().RotateVector(DefaultUpVector);
					FVector PrevUPVector = SimulationParticles[Parent].SimulationTransform.GetRotation().RotateVector(DefaultUpVector);
					FVector UpVector = (PrevUPVector * TwistRatio) + (CurUPVector * (1 - TwistRatio));
					UpVector.Normalize();

					FVector AimVector = SimulationParticles[Parent].SimulationTransform.GetRotation().RotateVector(DefaultAimVector);
					FQuat Delta0 = SolveAimWithTwist(SimulationParticles[Parent].SimulationTransform, NewChildPosition, AimVector, UpVector, DefaultUpVector, AimClampInDegree);
					FQuat Delta1 = SolveAimWithTwist(SimulationParticles[Parent].SimulationTransform, CurrentChildPosition, AimVector, UpVector, DefaultUpVector, AimClampInDegree);

					FQuat AimRotation0 = Delta0 * SimulationParticles[Parent].SimulationTransform.GetRotation();
					FQuat AimRotation1 = Delta1 * SimulationParticles[Parent].SimulationTransform.GetRotation();
					FQuat BlendedRotation = FQuat::Slerp(AimRotation0, AimRotation1, SwingRatio);
					BlendedRotation.Normalize();

					SimulationParticles[Parent].SimulationTransform.SetRotation(BlendedRotation);

					SimulationParticles[Parent].LocalTransform.SetRotation(SimulationParticles[GrandParent].SimulationTransform.GetRotation().Inverse() * BlendedRotation);
					SimulationParticles[Parent].LocalTransform.NormalizeRotation();
				}
			}
		}

		// for debug
		for (TArray<FSpringMagicJoint>& JointChain : SimulationJoints)
		{
			SimulationParticles[JointChain.Last().Child].SimulationTransform = SimulationParticles[JointChain.Last().Child].LocalTransform * SimulationParticles[JointChain.Last().Parent].SimulationTransform;
		}
	}
	
	ApplySimulateResult(Output, OutBoneTransforms);
	PreComponentTransform = ComponentTransform;
}

bool FAnimNode_SpringMagic::IsValidToEvaluate(const USkeleton* Skeleton, const FBoneContainer& RequiredBones)
{
	for (auto& BoneChain : BoneChains)
	{
		if (BoneChain.RootBone.IsValidToEvaluate(RequiredBones) && BoneChain.EndBone.IsValidToEvaluate(RequiredBones))
		{
			return true;
		}
	}
	return false;
}

void FAnimNode_SpringMagic::InitializeBoneReferences(const FBoneContainer& RequiredBones)
{
	for (auto& BoneChain : BoneChains)
	{
		BoneChain.RootBone.Initialize(RequiredBones);
		BoneChain.EndBone.Initialize(RequiredBones);
	}
	bNeedInitializeSimulation = true;
}

void FAnimNode_SpringMagic::ApplySimulateResult(FComponentSpacePoseContext& Output, TArray<FBoneTransform>& OutBoneTransforms)
{
	for (auto& Particle : SimulationParticles)
	{
		if (Particle.bRootBone || Particle.bDummyBone)
		{
			continue;
		}
		if (Particle.BoneIndex == INDEX_NONE)
		{
			UE_LOG(LogAnimation, Error, TEXT("ApplySimulateResult: Bone Is INDEX_NONE!"));
			continue;
		}
		OutBoneTransforms.Add(FBoneTransform(Particle.BoneIndex, Particle.SimulationTransform));
	}

	// bones are already sorted in "Parents before Children" order.
	OutBoneTransforms.Sort(FCompareBoneTransformIndex());
}

void FAnimNode_SpringMagic::InitializeSimulation(FComponentSpacePoseContext& Output, const FBoneContainer& RequiredBones)
{
	SimulationParticles.Empty();
	MaxChainLength = 0;
	TMap<FCompactPoseBoneIndex, int32> BoneIndexToArrayIndex;
	SimulationJoints.Empty(BoneChains.Num());
	for (auto& BoneChain : BoneChains)
	{
		FCompactPoseBoneIndex RootBoneIndex = BoneChain.RootBone.GetCompactPoseIndex(RequiredBones);
		FCompactPoseBoneIndex EndBoneIndex = BoneChain.EndBone.GetCompactPoseIndex(RequiredBones);
		if (RootBoneIndex == INDEX_NONE || EndBoneIndex == INDEX_NONE)
		{
			continue;
		}
		if (RootBoneIndex == EndBoneIndex)
		{
			continue;
		}
		FCompactPoseBoneIndex BoneIndex = EndBoneIndex;
		TArray<FCompactPoseBoneIndex> BoneIndices;
		while (RootBoneIndex != BoneIndex)
		{
			BoneIndices.Add(BoneIndex);
			FCompactPoseBoneIndex ParentBoneIndex = RequiredBones.GetParentBoneIndex(BoneIndex);
			if (ParentBoneIndex == INDEX_NONE)
			{
				FString ErrorString = FString::Printf(TEXT("RootBone[%s] Is Not A Parent Of EndBone[%s] In Skeleton Hierarchy"), *BoneChain.RootBone.BoneName.ToString(), *BoneChain.EndBone.BoneName.ToString());
				FMessageDialog::Open(EAppMsgType::Ok, FText::FromString(ErrorString));
				return;
			}
			BoneIndex = ParentBoneIndex;
		}
		BoneIndices.Add(RootBoneIndex);
		int32 IndexCount = BoneIndices.Num();
		if (IndexCount > 1)
		{
			int32 Idx = IndexCount - 1;
			TArray<int32> Particles;
			for (int32 i = Idx; i >= 0; --i)
			{
				if (int32* ParticleIndex = BoneIndexToArrayIndex.Find(BoneIndices[i]))
				{
					Particles.Add(*ParticleIndex);
					continue;
				}
				int32 Index = SimulationParticles.Num();
				auto& Particle = SimulationParticles.AddZeroed_GetRef();
				Particle.BoneIndex = BoneIndices[i];
				Particle.SimulationTransform = Output.Pose.GetComponentSpaceTransform(BoneIndices[i]);
				Particle.Radius = 1.0f;
				Particle.bDummyBone = false;
				Particle.bRootBone = i == Idx;
				Particles.Add(Index);
			}
			{
				int32 Index = SimulationParticles.Num();
				auto& Particle = SimulationParticles.AddZeroed_GetRef();
				Particle.LocalTransform.SetRotation(FQuat::Identity);
				Particle.LocalTransform.SetTranslation(FVector::XAxisVector * BoneChain.DummyBoneLength);
				Particle.SimulationTransform = Particle.LocalTransform * Output.Pose.GetComponentSpaceTransform(BoneIndices[0]);
				Particle.Radius = 1.0f;
				Particle.BoneIndex = BoneIndices[0];
				Particle.bDummyBone = true;
				Particle.bRootBone = false;
				Particles.Add(Index);
			}

			auto& Chain = SimulationJoints.AddDefaulted_GetRef();
			for (int32 i = 1; i < IndexCount; ++i)
			{
				int32 Grandparent = Particles[i - 1];
				int32 Parent      = Particles[i];
				int32 Child       = Particles[i + 1];
				FSpringMagicJoint&    Joint       = Chain.AddZeroed_GetRef();
				Joint.GrandParent = Grandparent;
				Joint.Parent = Parent;
				Joint.Child = Child;
			}
			MaxChainLength = FMath::Max(MaxChainLength, Chain.Num());
		}
	}
}

void FAnimNode_SpringMagic::InitializeCollision(const FBoneContainer& RequiredBones)
{
	CapsuleColliders.Empty();
	SphereColliders.Empty();
	if (!PhysicsAsset)
	{
		return;
	}
	for (const TObjectPtr<USkeletalBodySetup>& SkeletalBodySetup : PhysicsAsset->SkeletalBodySetups)
	{
		if (SkeletalBodySetup->PhysicsType == EPhysicsType::PhysType_Kinematic)
		{
			for (const FKSphereElem& SphereElem : SkeletalBodySetup->AggGeom.SphereElems)
			{
				FSpringMagicSphere Sphere;
				Sphere.BindBone.BoneName = SkeletalBodySetup->BoneName;
				Sphere.BindBone.Initialize(RequiredBones);
				if (!Sphere.BindBone.IsValidToEvaluate(RequiredBones))
				{
					continue;
				}
				Sphere.Center = SphereElem.Center;
				Sphere.Radius = SphereElem.Radius;
				SphereColliders.Add(Sphere);
			}
			for (const auto& SphylElem : SkeletalBodySetup->AggGeom.SphylElems)
			{
				FSpringMagicCapsule Capsule;
				Capsule.BindBone.BoneName = SkeletalBodySetup->BoneName;
				Capsule.BindBone.Initialize(RequiredBones);
				if (!Capsule.BindBone.IsValidToEvaluate(RequiredBones))
				{
					continue;
				}
				FVector Dir = SphylElem.Rotation.RotateVector(FVector::ZAxisVector);
				Capsule.CenterA = SphylElem.Center + Dir * (SphylElem.Length * 0.5f);
				Capsule.CenterB = SphylElem.Center - Dir * (SphylElem.Length * 0.5f);
				Capsule.Radius = SphylElem.Radius;
				CapsuleColliders.Add(Capsule);
			}
		}
	}
}

void FAnimNode_SpringMagic::UpdateCollision(FComponentSpacePoseContext& Output, const FBoneContainer& RequiredBones)
{
	for (auto& SphereCollider : SphereColliders)
	{
		auto BoneIndex = SphereCollider.BindBone.GetCompactPoseIndex(RequiredBones);
		SphereCollider.CenterSimSpace = Output.Pose.GetComponentSpaceTransform(BoneIndex).TransformPosition(SphereCollider.Center);
	}
	for (auto& CapsuleCollider : CapsuleColliders)
	{
		auto BoneIndex = CapsuleCollider.BindBone.GetCompactPoseIndex(RequiredBones);
		CapsuleCollider.CenterASimSpace = Output.Pose.GetComponentSpaceTransform(BoneIndex).TransformPosition(CapsuleCollider.CenterA);
		CapsuleCollider.CenterBSimSpace = Output.Pose.GetComponentSpaceTransform(BoneIndex).TransformPosition(CapsuleCollider.CenterB);
	}
}

void FAnimNode_SpringMagic::CollisionDetection(const FSpringMagicParticle& Particle, FVector& InOutTargetChildPosition, FVector& InOutCurrentChildPosition)
{
	for (auto& SphereCollider : SphereColliders)
	{
		FRay Ray;
		Ray.Origin = InOutCurrentChildPosition;
		Ray.Direction = InOutTargetChildPosition - InOutCurrentChildPosition;
		float MinDistance = SphereCollider.Radius + Particle.Radius;
		if (!Ray.Direction.Normalize())
		{
			FVector Temp = InOutCurrentChildPosition - SphereCollider.CenterSimSpace;
			if (Temp.SquaredLength() < FMath::Square(MinDistance))
			{
				Temp.Normalize();
				InOutCurrentChildPosition = SphereCollider.CenterSimSpace + Temp * MinDistance;
				InOutTargetChildPosition = SphereCollider.CenterSimSpace + Temp * MinDistance;
			}
		}
		float TMin, TMax;
		if (SpringMagicIntersectionTest::IntersectRaySphere(Ray, SphereCollider.CenterSimSpace, SphereCollider.Radius, TMin, TMax))
		{
			FVector Temp = InOutTargetChildPosition - InOutCurrentChildPosition;
			float TargetT = Temp.Length();
			bool bCurrentInSphere = TMin < 0.f && TMax > 0.f;
			bool bTargetInSphere = TMin < TargetT && TMax > TargetT;
			if (bCurrentInSphere && !bTargetInSphere)
			{
				InOutCurrentChildPosition = Ray.Origin + Ray.Direction * TMax;
			}
			else if (!bCurrentInSphere && bTargetInSphere)
			{
				InOutTargetChildPosition = Ray.Origin + Ray.Direction * TMin;
			}
			else if (bCurrentInSphere && bTargetInSphere)
			{
				Temp = InOutTargetChildPosition - SphereCollider.CenterSimSpace;
				Temp.Normalize();
				InOutTargetChildPosition = Temp* MinDistance + SphereCollider.CenterSimSpace;

				Temp = InOutCurrentChildPosition - SphereCollider.CenterSimSpace;
				Temp.Normalize();
				InOutCurrentChildPosition = Temp * MinDistance + SphereCollider.CenterSimSpace;
			}
		}
	}
	for (auto& CapsuleCollider : CapsuleColliders)
	{
		FRay Ray;
		Ray.Origin = InOutCurrentChildPosition;
		Ray.Direction = InOutTargetChildPosition - InOutCurrentChildPosition;
		float MinDistance = CapsuleCollider.Radius + Particle.Radius;
		if (!Ray.Direction.Normalize())
		{
			FVector ClosestPointOnSegment = FMath::ClosestPointOnSegment(InOutCurrentChildPosition, CapsuleCollider.CenterASimSpace, CapsuleCollider.CenterBSimSpace);
			FVector Temp = InOutCurrentChildPosition - ClosestPointOnSegment;
			if (Temp.SquaredLength() < FMath::Square(MinDistance))
			{
				Temp.Normalize();
				InOutCurrentChildPosition = ClosestPointOnSegment + Temp * MinDistance;
				InOutTargetChildPosition = ClosestPointOnSegment + Temp * MinDistance;
			}
		}
		FVector P1, P2, N1, N2;
		float TMin, TMax;
		if (SpringMagicIntersectionTest::IntersectRayCapsule(Ray, CapsuleCollider.CenterASimSpace, CapsuleCollider.CenterBSimSpace, CapsuleCollider.Radius, P1, P2, N1, N2, TMin, TMax))
		{
			FVector Temp = InOutTargetChildPosition - InOutCurrentChildPosition;
			float TargetT = Temp.Length();
			bool bCurrentInSphere = TMin < 0.f && TMax > 0.f;
			bool bTargetInSphere = TMin < TargetT && TMax > TargetT;
			if (bCurrentInSphere && !bTargetInSphere)
			{
				InOutCurrentChildPosition = Ray.Origin + Ray.Direction * TMax;
			}
			else if (!bCurrentInSphere && bTargetInSphere)
			{
				InOutTargetChildPosition = Ray.Origin + Ray.Direction * TMin;
			}
			else if (bCurrentInSphere && bTargetInSphere)
			{
				float TMid = (TMin + TMax) * 0.5f;
				FVector MinPoint = Ray.Origin + Ray.Direction * TMid;
				FVector ClosestPointOnSegment = FMath::ClosestPointOnSegment(MinPoint, CapsuleCollider.CenterASimSpace, CapsuleCollider.CenterBSimSpace);

				Temp = InOutTargetChildPosition - ClosestPointOnSegment;
				Temp.Normalize();
				InOutTargetChildPosition = Temp * MinDistance + ClosestPointOnSegment;

				Temp = InOutCurrentChildPosition - ClosestPointOnSegment;
				Temp.Normalize();
				InOutCurrentChildPosition = Temp * MinDistance + ClosestPointOnSegment;
			}
		}
	}
}

// Copyright LoogLong. All Rights Reserved.
#include "AnimNode_LoogPhysics.h"
#include "AnimationRuntime.h"
#include "Animation/AnimInstanceProxy.h"
#include "Curves/CurveFloat.h"
#include "Runtime/Launch/Resources/Version.h"
#include "SceneInterface.h"
#include "PhysicsEngine/PhysicsAsset.h"


DECLARE_CYCLE_STAT(TEXT("LoogPhysics_Eval"), STAT_LoogPhysics_Eval, STATGROUP_Anim);

namespace LoogPhysicsIntersectionTest
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

FAnimNode_LoogPhysics::FAnimNode_LoogPhysics()
{
}

void FAnimNode_LoogPhysics::Initialize_AnyThread(const FAnimationInitializeContext& Context)
{
	Super::Initialize_AnyThread(Context);
	const FBoneContainer& RequiredBones = Context.AnimInstanceProxy->GetRequiredBones();
	InitializeBoneReferences(RequiredBones);
}

void FAnimNode_LoogPhysics::CacheBones_AnyThread(const FAnimationCacheBonesContext& Context)
{
	Super::CacheBones_AnyThread(Context);
}

void FAnimNode_LoogPhysics::ResetDynamics(ETeleportType InTeleportType)
{
	PendingTeleportType = InTeleportType;
}

void FAnimNode_LoogPhysics::UpdateInternal(const FAnimationUpdateContext& Context)
{
	Super::UpdateInternal(Context);
}

void FAnimNode_LoogPhysics::GatherDebugData(FNodeDebugData& DebugData)
{
	Super::GatherDebugData(DebugData);
}

void FAnimNode_LoogPhysics::EvaluateSkeletalControl_AnyThread(FComponentSpacePoseContext& Output,
                                                                TArray<FBoneTransform>& OutBoneTransforms)
{
	SCOPE_CYCLE_COUNTER(STAT_LoogPhysics_Eval);
	const auto& RequiredBones = Output.AnimInstanceProxy->GetRequiredBones();
	if (bNeedInitializeSimulation)
	{
		InitializeSimulation(Output, RequiredBones);
		InitializeCollision(RequiredBones);
		bNeedInitializeSimulation = false;
	}
	const float DeltaTimeSeconds = Output.AnimInstanceProxy->GetDeltaSeconds();
	const bool bPause = DeltaTimeSeconds < UE_KINDA_SMALL_NUMBER; // Now In Pie pause anim instance
	if (!bPause)
	{
		UpdateCollision(Output, RequiredBones);
		SimulatePhysics(Output, DeltaTimeSeconds);
	}
	ApplySimulateResult(Output, OutBoneTransforms);

}

bool FAnimNode_LoogPhysics::IsValidToEvaluate(const USkeleton* Skeleton, const FBoneContainer& RequiredBones)
{
	return bIsValidToSimulation;
}

void FAnimNode_LoogPhysics::InitializeBoneReferences(const FBoneContainer& RequiredBones)
{
	bIsValidToSimulation = true;
	for (auto& Particle : Particles)
	{
		if (!Particle.BindBone.Initialize(RequiredBones))
		{
			bIsValidToSimulation = false;
			break;
		}
	}
	bNeedInitializeSimulation = bIsValidToSimulation;
}



void FAnimNode_LoogPhysics::InitializeSimulation(FComponentSpacePoseContext& Output, const FBoneContainer& RequiredBones)
{
	RuntimeParticles.Empty(Particles.Num());
	RuntimeParticles.AddZeroed(Particles.Num());
	TMap<FCompactPoseBoneIndex, int32> BoneIndexToParticleIndex;
	for (int32 ParticleIndex = 0; ParticleIndex < Particles.Num(); ++ParticleIndex)
	{
		const auto& ConfigParticle = Particles[ParticleIndex];
		FCompactPoseBoneIndex BoneIndex = ConfigParticle.BindBone.GetCompactPoseIndex(RequiredBones);
		BoneIndexToParticleIndex.Add(BoneIndex, ParticleIndex);
		const auto& BoneTransform = Output.Pose.GetComponentSpaceTransform(BoneIndex);
		switch (ConfigParticle.ParticleType)
		{
		case ELoogPhysicsParticleType::Root:
			RuntimeParticles[ParticleIndex].BoneIndex = BoneIndex;
			RuntimeParticles[ParticleIndex].InvMass = 0.f;
			RuntimeParticles[ParticleIndex].PrevPosition = BoneTransform.GetTranslation();
			RuntimeParticles[ParticleIndex].Position = RuntimeParticles[ParticleIndex].PrevPosition;
			RuntimeParticles[ParticleIndex].Velocity = FVector::Zero();
			break;
		case ELoogPhysicsParticleType::Bone:
			RuntimeParticles[ParticleIndex].BoneIndex = BoneIndex;
			RuntimeParticles[ParticleIndex].InvMass = 1.0f / ConfigParticle.Mass;
			RuntimeParticles[ParticleIndex].PrevPosition = BoneTransform.GetTranslation();
			RuntimeParticles[ParticleIndex].Position = RuntimeParticles[ParticleIndex].PrevPosition;
			RuntimeParticles[ParticleIndex].Velocity = FVector::Zero();
			break;
		case ELoogPhysicsParticleType::VirtualBone:
			RuntimeParticles[ParticleIndex].BoneIndex = BoneIndex;
			RuntimeParticles[ParticleIndex].InvMass = 1.0f / ConfigParticle.Mass;
			RuntimeParticles[ParticleIndex].PrevPosition = BoneTransform.TransformPosition(ConfigParticle.VirtualBoneLocalPosition);
			RuntimeParticles[ParticleIndex].Position = RuntimeParticles[ParticleIndex].PrevPosition;
			RuntimeParticles[ParticleIndex].Velocity = FVector::Zero();
			break;
		case ELoogPhysicsParticleType::Max:
		default:
			break;
		}
	}
	for (int32 ParticleIndex = 0; ParticleIndex < Particles.Num(); ++ParticleIndex)
	{
		const auto& ConfigParticle = Particles[ParticleIndex];
		const auto& ParentBoneIndex = RequiredBones.GetParentBoneIndex(RuntimeParticles[ParticleIndex].BoneIndex);
		if (ParentBoneIndex == INDEX_NONE)
		{
			// impossible!
			check(false);
		}
		switch (ConfigParticle.ParticleType)
		{
		case ELoogPhysicsParticleType::Root:
			break;
		case ELoogPhysicsParticleType::Bone:
			RuntimeParticles[ParticleIndex].ParentParticleIndex = BoneIndexToParticleIndex[ParentBoneIndex];
			break;
		case ELoogPhysicsParticleType::VirtualBone:
		case ELoogPhysicsParticleType::Max:
		default:
			break;
		}
	}

	// first calculate constraints num
	int32 NumOfConstraints = 0;
	for (FLoogPhysicsClothSection& ClothSection : ClothSections)
	{
		for (auto& Chain : ClothSection.Chains)
		{
			NumOfConstraints += (Chain.ParticleIndices.Num() - 1) * 2;
		}
		NumOfConstraints += ClothSection.Constraints.Num();
	}
	
	RuntimeConstraints.Empty(NumOfConstraints);

	for (FLoogPhysicsClothSection& ClothSection : ClothSections)
	{
		int32 ChainLength = ClothSection.Chains[0].ParticleIndices.Num();
		for (int32 i = 1; i < ClothSection.Chains.Num(); ++i)
		{
			if (ClothSection.Chains[i].ParticleIndices.Num() != ChainLength)
			{
				FString ErrorString = FString::Printf(TEXT("Chain Length Must The Same!"));
				FMessageDialog::Open(EAppMsgType::Ok, FText::FromString(ErrorString));
				return;
			}
		}
	}
	for (FLoogPhysicsClothSection& ClothSection : ClothSections)
	{
		int32 ChainNum = ClothSection.Chains.Num();
		
		for (int32 ChainIndex = 0; ChainIndex < ChainNum; ++ChainIndex)
		{
			auto& Chain = ClothSection.Chains[ChainIndex];
			int32 ChainLength = Chain.ParticleIndices.Num();
			{ // global structure
				int32 RootIndex = Chain.ParticleIndices[Chain.RootParticleIndex];
				auto& RootRuntimeParticle = RuntimeParticles[RootIndex];
				auto& RootConfigParticle = Particles[RootIndex];
				for (int32 Idx = 1; Idx < ChainLength; ++Idx)
				{
					int32 ChildIndex = Chain.ParticleIndices[Idx];
					auto& ChildRuntimeParticle = RuntimeParticles[ChildIndex];
					auto& ChildConfigParticle = Particles[ChildIndex];
					auto& Constraint = RuntimeConstraints.AddZeroed_GetRef();
					Constraint.Particle0Index = RootIndex;
					Constraint.Particle1Index = ChildIndex;
					Constraint.RestLength = (RootRuntimeParticle.Position - ChildRuntimeParticle.Position).Length();
					Constraint.ShrinkCompliance = (RootConfigParticle.GlobalStructureCompliance + ChildConfigParticle.GlobalStructureCompliance) * 0.5f;
					Constraint.StretchCompliance = (RootConfigParticle.GlobalStructureCompliance + ChildConfigParticle.GlobalStructureCompliance) * 0.5f;
				}
			}
			{// local structure vertical
				for (int32 Idx = 1; Idx < ChainLength; ++Idx)
				{
					int32 AIndex = Chain.ParticleIndices[Idx - 1];
					auto& ARuntimeParticle = RuntimeParticles[AIndex];
					auto& AConfigParticle = Particles[AIndex];

					int32 BIndex = Chain.ParticleIndices[Idx];
					auto& BRuntimeParticle = RuntimeParticles[BIndex];
					auto& BConfigParticle = Particles[BIndex];
					auto& Constraint = RuntimeConstraints.AddZeroed_GetRef();
					Constraint.Particle0Index = AIndex;
					Constraint.Particle1Index = BIndex;
					Constraint.RestLength = (ARuntimeParticle.Position - BRuntimeParticle.Position).Length();
					Constraint.ShrinkCompliance = (AConfigParticle.LocalStructureCompliance + BConfigParticle.LocalStructureCompliance) * 0.5f;
					Constraint.StretchCompliance = (AConfigParticle.LocalStructureCompliance + BConfigParticle.LocalStructureCompliance) * 0.5f;
				}
			}

			{// local structure horizontal
				int32 NextChainIndex = ChainIndex + 1;
				if (NextChainIndex < ChainNum || ClothSection.bChainLoop)
				{
					if (NextChainIndex >= ChainNum)
					{
						NextChainIndex = 0;
					}
					auto& NextChain = ClothSection.Chains[NextChainIndex];
					for (int32 Idx = 1; Idx < ChainLength; ++Idx)
					{
						int32 AIndex = Chain.ParticleIndices[Idx];
						auto& ARuntimeParticle = RuntimeParticles[AIndex];
						auto& AConfigParticle = Particles[AIndex];

						int32 BIndex = NextChain.ParticleIndices[Idx];
						auto& BRuntimeParticle = RuntimeParticles[BIndex];
						auto& BConfigParticle = Particles[BIndex];
						auto& Constraint = RuntimeConstraints.AddZeroed_GetRef();
						Constraint.Particle0Index = AIndex;
						Constraint.Particle1Index = BIndex;
						Constraint.RestLength = (ARuntimeParticle.Position - BRuntimeParticle.Position).Length();
						Constraint.ShrinkCompliance = (AConfigParticle.LocalStructureCompliance + BConfigParticle.LocalStructureCompliance) * 0.5f;
						Constraint.StretchCompliance = (AConfigParticle.LocalStructureCompliance + BConfigParticle.LocalStructureCompliance) * 0.5f;
					}
				}
				
			}
		}

		for (FLoogPhysicsConstraint& Con : ClothSection.Constraints)
		{
			int32 AIndex = Con.ParticleAIndex;
			auto& ARuntimeParticle = RuntimeParticles[AIndex];

			int32 BIndex = Con.ParticleBIndex;
			auto& BRuntimeParticle = RuntimeParticles[BIndex];
			auto& Constraint = RuntimeConstraints.AddZeroed_GetRef();
			Constraint.Particle0Index = AIndex;
			Constraint.Particle1Index = BIndex;
			Constraint.RestLength = (ARuntimeParticle.Position - BRuntimeParticle.Position).Length();
			Constraint.ShrinkCompliance = Con.ShrinkCompliance;
			Constraint.StretchCompliance = Con.StretchCompliance;
		}
	}
}

void FAnimNode_LoogPhysics::InitializeCollision(const FBoneContainer& RequiredBones)
{
	RuntimeColliders.Empty();
	SetupColliders.Empty();
	if (!PhysicsAsset)
	{
		return;
	}

	for (const TObjectPtr<USkeletalBodySetup>& SkeletalBodySetup : PhysicsAsset->SkeletalBodySetups)
	{
		if (SkeletalBodySetup->PhysicsType == EPhysicsType::PhysType_Kinematic)
		{
			for (const auto& TaperedCapsuleElem : SkeletalBodySetup->AggGeom.TaperedCapsuleElems)
			{
				FLoogPhysicsCapsule Capsule;
				Capsule.BindBone.BoneName = SkeletalBodySetup->BoneName;
				Capsule.BindBone.Initialize(RequiredBones);
				if (!Capsule.BindBone.IsValidToEvaluate(RequiredBones))
				{
					continue;
				}
				FVector Dir = TaperedCapsuleElem.Rotation.RotateVector(FVector::ZAxisVector);
				Capsule.CenterA = TaperedCapsuleElem.Center + Dir * (TaperedCapsuleElem.Length * 0.5f);
				Capsule.CenterB = TaperedCapsuleElem.Center - Dir * (TaperedCapsuleElem.Length * 0.5f);
				RuntimeColliders.Add(Capsule);
				SetupColliders.Add(TaperedCapsuleElem);
			}
		}
	}
}



void FAnimNode_LoogPhysics::SimulatePhysics(FComponentSpacePoseContext& Output, const float& DeltaTimeSeconds)
{
	float SimulateDeltaTime = 1.0f / FrameRate;
	// compute how many times for simulation
	int32 SimTimes = static_cast<int32>(DeltaTimeSeconds * FrameRate) + 1;
	SimTimes = FMath::Min(SimTimes, MaxSimulationPerFrame);
	for (int32 i = 0; i < SimTimes; ++i)
	{
		SimulateOnce(Output, SimulateDeltaTime);
	}
}

void FAnimNode_LoogPhysics::SimulateOnce(FComponentSpacePoseContext& Output, const float& DeltaTimeSeconds)
{
	// first update particle dynamic
	const auto& ComponentTransform = Output.AnimInstanceProxy->GetComponentTransform();
	FVector GravityVelocity = ComponentTransform.InverseTransformVectorNoScale(GravityWorldSpace) * DeltaTimeSeconds;

	FVector WindVelocity = FVector::Zero();
	if (bUseWind)
	{
		WindSinCounter = WindSinCounter + DeltaTimeSeconds;
		WindSinCounter = FMath::Fmod(WindSinCounter, UE_PI);
		FVector WindDirSimSpace = ComponentTransform.InverseTransformVectorNoScale(WindDirectionWorldSpace);
		WindVelocity = (FMath::Sin(WindSinCounter) * WindAmplitude + WindAmplitudeBias) * WindDirSimSpace;
	}

	for (int32 ParticleIndex = 0; ParticleIndex < Particles.Num(); ++ParticleIndex)
	{
		const auto& ConfigParticle = Particles[ParticleIndex];
		auto& RuntimeParticle = RuntimeParticles[ParticleIndex];
		switch (ConfigParticle.ParticleType)
		{
		case ELoogPhysicsParticleType::Root:
			{
				const auto& BoneTransCS = Output.Pose.GetComponentSpaceTransform(RuntimeParticle.BoneIndex);
				RuntimeParticle.Position = BoneTransCS.GetTranslation();
			}
			break;
		case ELoogPhysicsParticleType::Bone:
		case ELoogPhysicsParticleType::VirtualBone:
			{
				RuntimeParticle.PrevPosition = RuntimeParticle.Position;
				FVector CurrentVelocity = GravityVelocity + WindVelocity + RuntimeParticle.Velocity * (1.0f - ConfigParticle.LinearDamping);
				FVector::FReal CurrentVelocityMagnitude = CurrentVelocity.Length();
				if (CurrentVelocityMagnitude > ConfigParticle.MaxVelocity)
				{
					CurrentVelocity = CurrentVelocity / CurrentVelocityMagnitude * ConfigParticle.MaxVelocity;
				}
				RuntimeParticle.Position += CurrentVelocity * DeltaTimeSeconds;
			}
			break;
		case ELoogPhysicsParticleType::Max:
		default:
			break;
		}
	}

	// second Update Constraint
	float InvDeltaTimeSquared = 1.0f / (DeltaTimeSeconds * DeltaTimeSeconds);
	for (const FLoogPhysicsRuntimeConstraint& RuntimeConstraint : RuntimeConstraints)
	{
		auto& Particle0 = RuntimeParticles[RuntimeConstraint.Particle0Index];
		auto& Particle1 = RuntimeParticles[RuntimeConstraint.Particle1Index];
		float W = Particle0.InvMass + Particle1.InvMass;
		if (W == 0.f)
		{
			continue;
		}
		FVector Delta = Particle0.Position - Particle1.Position;
		float CurrentLength = Delta.Length();
		if (CurrentLength <= 0.f)
		{
			continue;
		}
		FVector GradC = Delta / CurrentLength;
		float C = CurrentLength - RuntimeConstraint.RestLength;
		float Compliance = C > 0 ? RuntimeConstraint.ShrinkCompliance : RuntimeConstraint.StretchCompliance;
		float XPBDAlpha = Compliance * InvDeltaTimeSquared;
		float Lambda = -C / (W + XPBDAlpha);

		FVector Correction = GradC * Lambda;

		Particle0.Position += Correction * Particle0.InvMass;
		Particle1.Position -= Correction * Particle1.InvMass;
	}

	// handle collision
	if (bEnableCollision)
	{
		CollisionDetection();
	}

	// Update Velocity
	float InvDeltaTime = 1.0f / DeltaTimeSeconds;
	for (int32 ParticleIndex = 0; ParticleIndex < Particles.Num(); ++ParticleIndex)
	{
		auto& RuntimeParticle = RuntimeParticles[ParticleIndex];

		RuntimeParticle.Velocity = (RuntimeParticle.Position - RuntimeParticle.PrevPosition) * InvDeltaTime;
	}
}

void FAnimNode_LoogPhysics::UpdateCollision(FComponentSpacePoseContext& Output, const FBoneContainer& RequiredBones)
{
	int32 ColliderCount = RuntimeColliders.Num();
	for (int Index = 0; Index < ColliderCount; ++Index)
	{
		auto& SetupCollider = SetupColliders[Index];
		auto& RuntimeCollider = RuntimeColliders[Index];

		FTransform BoneTransform = Output.Pose.GetComponentSpaceTransform(RuntimeCollider.BindBone.GetCompactPoseIndex(RequiredBones));
		RuntimeCollider.TransformSimSpace = SetupCollider.GetTransform() * BoneTransform;
		RuntimeCollider.CenterASimSpace = BoneTransform.TransformPosition(RuntimeCollider.CenterA);
		RuntimeCollider.CenterBSimSpace = BoneTransform.TransformPosition(RuntimeCollider.CenterB);
	}
}

void FAnimNode_LoogPhysics::CollisionDetection()
{
	int32 ColliderCount = RuntimeColliders.Num();
	for (int Index = 0; Index < ColliderCount; ++Index)
	{
		auto& SetupCollider = SetupColliders[Index];
		auto& RuntimeCollider = RuntimeColliders[Index];

		const FVector& StartPoint = RuntimeCollider.CenterASimSpace;
		const FVector& EndPoint = RuntimeCollider.CenterBSimSpace;
		int32 ParticleCount = RuntimeParticles.Num();
		for (int ParticleIndex = 0; ParticleIndex < ParticleCount; ++ParticleIndex)
		{
			const auto& ConfigParticle = Particles[ParticleIndex];
			auto& Particle = RuntimeParticles[ParticleIndex];

			const FVector& Point = Particle.Position;

			const FVector Segment = EndPoint - StartPoint;
			const FVector VectToPoint = Point - StartPoint;

			FVector::FReal SphereRadius;
			FVector ClosestPoint;
			// See if closest point is before StartPoint
			const FVector::FReal Dot1 = VectToPoint | Segment;
			if (Dot1 <= 0)
			{
				ClosestPoint = StartPoint;
				SphereRadius = SetupCollider.Radius0;
			}
			else
			{
				// See if closest point is beyond EndPoint
				const FVector::FReal Dot2 = Segment | Segment;
				if (Dot2 <= Dot1)
				{
					ClosestPoint = EndPoint;
					SphereRadius = SetupCollider.Radius1;
				}
				else
				{
					// Closest Point is within segment
					ClosestPoint = StartPoint + Segment * (Dot1 / Dot2);
					SphereRadius = SetupCollider.Radius0 + (SetupCollider.Radius1 - SetupCollider.Radius0) * (Dot1 / Dot2);
				}
			}
			FVector DirVector = Point - ClosestPoint;
			FVector::FReal CurrentDistance = DirVector.Length();
			FVector::FReal RestMinDistance = ConfigParticle.Thickness + SphereRadius;
			if (CurrentDistance < RestMinDistance)
			{
				if (CurrentDistance < UE_SMALL_NUMBER)
				{
					UE_LOG(LogAnimation, Warning, TEXT("Zero Length In LoogPhysics!"));
				}
				else
				{
					Particle.Position = ClosestPoint + DirVector * RestMinDistance / CurrentDistance;
				}
			}
			
		}

		// capsule line
		for (auto& ClothSection : ClothSections)
		{
			for (FLoogPhysicsParticleCollider& Collider : ClothSection.Colliders)
			{
				auto& Particle0 = RuntimeParticles[Collider.Particle0Index];
				const auto& ConfigParticle0 = Particles[Collider.Particle0Index];
				auto& Particle1 = RuntimeParticles[Collider.Particle1Index];
				const auto& ConfigParticle1 = Particles[Collider.Particle1Index];

				FVector ClosestPoint1, ClosestPoint2;
				FMath::SegmentDistToSegmentSafe(RuntimeCollider.CenterASimSpace, RuntimeCollider.CenterBSimSpace, Particle0.Position, Particle1.Position, ClosestPoint1, ClosestPoint2);
				FVector CapsuleDir = RuntimeCollider.CenterBSimSpace - RuntimeCollider.CenterASimSpace;
				
				float S;
				if (CapsuleDir.Length() > 0)
				{
					S = (ClosestPoint1 - RuntimeCollider.CenterASimSpace).Dot(CapsuleDir / CapsuleDir.Length());
					S = S / CapsuleDir.Length();
				}
				else
				{
					S = 0.f;
				}

				FVector ParticleDir = Particle1.Position - Particle0.Position;
				float T;
				if (ParticleDir.Length() > 0)
				{
					T = (ClosestPoint2 - Particle0.Position).Dot(ParticleDir / ParticleDir.Length());
					T = T / ParticleDir.Length();
				}
				else
				{
					T = 0.f;
				}

				float Radius1 = FMath::Lerp(SetupCollider.Radius0, SetupCollider.Radius1, S);
				float Radius2 = FMath::Lerp(ConfigParticle0.Thickness, ConfigParticle1.Thickness, T);

				FVector ClosestDirection = ClosestPoint2 - ClosestPoint1;
				float Distance = ClosestDirection.Length();
				float RadiusSum = Radius1 + Radius2;

				if (Distance <= RadiusSum && Distance > 0.f)
				{
					FVector Normal = ClosestDirection / Distance;
					FVector Delta = Normal * (RadiusSum - Distance);
					Particle0.Position += Delta;
					Particle1.Position += Delta;
				}
			}
		}
	}

	
}

void FAnimNode_LoogPhysics::ApplySimulateResult(FComponentSpacePoseContext& Output, TArray<FBoneTransform>& OutBoneTransforms)
{
	const auto& RequiredBones = Output.AnimInstanceProxy->GetRequiredBones();
	
	for (int32 ParticleIndex = 0; ParticleIndex < Particles.Num(); ++ParticleIndex)
	{
		auto& RuntimeParticle = RuntimeParticles[ParticleIndex];
		auto& ParentParticle = RuntimeParticles[RuntimeParticle.ParentParticleIndex];
		FVector BoneSimDirection = RuntimeParticle.Position - ParentParticle.Position;
		if (!BoneSimDirection.Normalize())
		{
			continue;
		}
		// const FTransform& BoneRefPoseTransform = RequiredBones.GetRefPoseTransform(RuntimeParticle.BoneIndex);
		// float RefBoneLength = BoneRefPoseTransform.GetTranslation().Length();

		FTransform        BonePoseTransform = Output.Pose.GetComponentSpaceTransform(ParentParticle.BoneIndex);
		const FVector& DefaultForwardAxis = FVector::XAxisVector;
		FVector           BonePoseDirection = BonePoseTransform.GetRotation().RotateVector(DefaultForwardAxis);
		FQuat             AddRotation = FQuat::FindBetweenNormals(BonePoseDirection, BoneSimDirection);
		BonePoseTransform.SetRotation(AddRotation * BonePoseTransform.GetRotation());

		const auto& ConfigParticle = Particles[RuntimeParticle.ParentParticleIndex];
		switch (ConfigParticle.ParticleType)
		{
		case ELoogPhysicsParticleType::Root:
			{
				BonePoseTransform.SetRotation(AddRotation * BonePoseTransform.GetRotation());
				BonePoseTransform.SetTranslation(ParentParticle.Position);
				if (BonePoseTransform.ContainsNaN())
				{
					continue;
				}
				OutBoneTransforms.Add(FBoneTransform(ParentParticle.BoneIndex, BonePoseTransform));
			}
			break;
		case ELoogPhysicsParticleType::Bone:
			{
				
				BonePoseTransform.SetTranslation(ParentParticle.Position);
				if (BonePoseTransform.ContainsNaN())
				{
					continue;
				}
				OutBoneTransforms.Add(FBoneTransform(ParentParticle.BoneIndex, BonePoseTransform));
			}
			break;
		case ELoogPhysicsParticleType::VirtualBone:
		case ELoogPhysicsParticleType::Max:
		default:
			break;
		}
	}
	OutBoneTransforms.Sort(FCompareBoneTransformIndex());
}
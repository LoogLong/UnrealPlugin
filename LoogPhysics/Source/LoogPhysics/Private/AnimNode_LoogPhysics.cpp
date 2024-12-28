// Copyright LoogLong. All Rights Reserved.
#include "AnimNode_LoogPhysics.h"
#include "AnimationRuntime.h"
#include "Animation/AnimInstanceProxy.h"
#include "Curves/CurveFloat.h"
#include "Runtime/Launch/Resources/Version.h"
#include "SceneInterface.h"
#include "PhysicsEngine/PhysicsAsset.h"


DECLARE_CYCLE_STAT(TEXT("LoogPhysics_Eval"), STAT_LoogPhysics_Eval, STATGROUP_Anim);

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
		PreComponentTransform = Output.AnimInstanceProxy->GetComponentTransform();
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
	PreComponentTransform = Output.AnimInstanceProxy->GetComponentTransform();
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
			RuntimeParticles[ParticleIndex].BendingCompliance = ConfigParticle.BendingStructureCompliance;
			RuntimeParticles[ParticleIndex].PrevPosition = BoneTransform.GetTranslation();
			RuntimeParticles[ParticleIndex].Position = RuntimeParticles[ParticleIndex].PrevPosition;
			RuntimeParticles[ParticleIndex].Velocity = FVector::Zero();
			break;
		case ELoogPhysicsParticleType::Bone:
			RuntimeParticles[ParticleIndex].BoneIndex = BoneIndex;
			RuntimeParticles[ParticleIndex].InvMass = 1.0f / ConfigParticle.Mass;
			RuntimeParticles[ParticleIndex].BendingCompliance = ConfigParticle.BendingStructureCompliance;
			RuntimeParticles[ParticleIndex].PrevPosition = BoneTransform.GetTranslation();
			RuntimeParticles[ParticleIndex].Position = RuntimeParticles[ParticleIndex].PrevPosition;
			RuntimeParticles[ParticleIndex].Velocity = FVector::Zero();
			break;
		case ELoogPhysicsParticleType::VirtualBone:
			RuntimeParticles[ParticleIndex].BoneIndex = BoneIndex;
			RuntimeParticles[ParticleIndex].InvMass = 1.0f / ConfigParticle.Mass;
			RuntimeParticles[ParticleIndex].BendingCompliance = ConfigParticle.BendingStructureCompliance;
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
		check(ParentBoneIndex != INDEX_NONE);

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
	
	RuntimeDistanceConstraints.Empty(NumOfConstraints);

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
					auto& Constraint = RuntimeDistanceConstraints.AddZeroed_GetRef();
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
					auto& Constraint = RuntimeDistanceConstraints.AddZeroed_GetRef();
					Constraint.Particle0Index = AIndex;
					Constraint.Particle1Index = BIndex;
					Constraint.RestLength = (ARuntimeParticle.Position - BRuntimeParticle.Position).Length();
					Constraint.ShrinkCompliance = (AConfigParticle.LocalStructureCompliance + BConfigParticle.LocalStructureCompliance) * 0.5f;
					Constraint.StretchCompliance = (AConfigParticle.LocalStructureCompliance + BConfigParticle.LocalStructureCompliance) * 0.5f;
				}
			}
		}

		for (FLoogPhysicsConstraint& Con : ClothSection.Constraints)
		{
			int32 AIndex = Con.ParticleAIndex;
			auto& ARuntimeParticle = RuntimeParticles[AIndex];

			int32 BIndex = Con.ParticleBIndex;
			auto& BRuntimeParticle = RuntimeParticles[BIndex];
			auto& Constraint = RuntimeDistanceConstraints.AddZeroed_GetRef();
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
	// compute how many Steps for simulation
	int32 SimStepCount = static_cast<int32>(DeltaTimeSeconds * FrameRate) + 1;
	SimStepCount = FMath::Min(SimStepCount, MaxSimulationPerFrame);
	if (SimStepCount <= 0)
	{
		return;
	}
	const auto& ComponentTransform = Output.AnimInstanceProxy->GetComponentTransform();
	FVector PreComponentPosition = PreComponentTransform.GetTranslation();
	FVector CurrentComponentPosition = ComponentTransform.GetTranslation();
	FVector InvComponentTranslation = PreComponentPosition - CurrentComponentPosition;

	FQuat PreComponentRotation = PreComponentTransform.GetRotation();
	FQuat CurrentComponentRotation = ComponentTransform.GetRotation();
	FQuat InvComponentRotation = CurrentComponentRotation.Inverse() * PreComponentRotation;

	InvComponentTranslation = CurrentComponentRotation.Inverse().RotateVector(InvComponentTranslation);
	if (InvComponentTranslation.Length() > ComponentMovementMaxDistance)
	{
		InvComponentTranslation = InvComponentTranslation * (ComponentMovementMaxDistance / InvComponentTranslation.Length());
	}

	float InvRotationAngle = InvComponentRotation.GetAngle();
	if (InvRotationAngle > FMath::DegreesToRadians(ComponentMovementMaxAngle))
	{
		InvRotationAngle = FMath::DegreesToRadians(ComponentMovementMaxAngle);
	}
	FVector StepDeltaTranslation = InvComponentTranslation * ComponentMovementTranslation / SimStepCount;
	FQuat StepDeltaRotation = FQuat(InvComponentRotation.GetRotationAxis(), InvRotationAngle * ComponentMovementRotation / SimStepCount);

	for (int32 i = 0; i < SimStepCount; ++i)
	{
		SimulateOnce(Output, SimulateDeltaTime, StepDeltaTranslation, StepDeltaRotation);
	}
}

void FAnimNode_LoogPhysics::SimulateOnce(FComponentSpacePoseContext& Output, const float& DeltaTimeSeconds, const FVector& InStepDeltaTranslation, const FQuat& InStepDeltaRotation)
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
				RuntimeParticle.SimRotation = BoneTransCS.GetRotation();
			}
			break;
		case ELoogPhysicsParticleType::Bone:
		case ELoogPhysicsParticleType::VirtualBone:
			{
				// return to last step position in world space
				FVector LastStepPosition = InStepDeltaRotation.RotateVector(RuntimeParticle.Position);
				LastStepPosition = LastStepPosition + InStepDeltaTranslation ;
				RuntimeParticle.PrevPosition = LastStepPosition;

				FVector CurrentVelocity = GravityVelocity + WindVelocity + RuntimeParticle.Velocity * (1.0f - ConfigParticle.LinearDamping);
				FVector::FReal CurrentVelocityMagnitude = CurrentVelocity.Length();
				if (CurrentVelocityMagnitude > ConfigParticle.MaxVelocity)
				{
					CurrentVelocity = CurrentVelocity / CurrentVelocityMagnitude * ConfigParticle.MaxVelocity;
				}
				RuntimeParticle.Position = LastStepPosition + CurrentVelocity * DeltaTimeSeconds;
			}
			break;
		case ELoogPhysicsParticleType::Max:
		default:
			break;
		}
	}

	// second Update Constraint
	float InvDeltaTimeSquared = 1.0f / (DeltaTimeSeconds * DeltaTimeSeconds);
	for (const FLoogPhysicsDistanceConstraint& DistanceConstraint : RuntimeDistanceConstraints)
	{
		auto& Particle0 = RuntimeParticles[DistanceConstraint.Particle0Index];
		auto& Particle1 = RuntimeParticles[DistanceConstraint.Particle1Index];
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
		float C = CurrentLength - DistanceConstraint.RestLength;
		float Compliance = C > 0 ? DistanceConstraint.ShrinkCompliance : DistanceConstraint.StretchCompliance;
		float XPBDAlpha = Compliance * InvDeltaTimeSquared;
		float Lambda = -C / (W + XPBDAlpha);

		FVector Correction = GradC * Lambda;

		Particle0.Position += Correction * Particle0.InvMass;
		Particle1.Position -= Correction * Particle1.InvMass;
	}

	// bending constraint
	for (int32 ParticleIndex = 0; ParticleIndex < Particles.Num(); ++ParticleIndex)
	{
		auto& ChildConfigParticle = Particles[ParticleIndex];
		if (ChildConfigParticle.ParticleType == ELoogPhysicsParticleType::Root)
		{
			continue;
		}
		auto& ChildParticle = RuntimeParticles[ParticleIndex];
		auto& ParentParticle = RuntimeParticles[ChildParticle.ParentParticleIndex];
		float W = ChildParticle.InvMass + ParentParticle.InvMass;
		if (W == 0.f)
		{
			continue;
		}

		FTransform ChildPoseLocalTransform = Output.Pose.GetLocalSpaceTransform(ChildParticle.BoneIndex);
		
		FVector BonePoseSimDirection = ParentParticle.SimRotation.RotateVector(ChildPoseLocalTransform.GetTranslation());
		FVector BonePosePosition = BonePoseSimDirection + ParentParticle.Position;
		FVector Diff = ChildParticle.Position - BonePosePosition;
		float CurrentLength = Diff.Length();
		if (CurrentLength < UE_KINDA_SMALL_NUMBER)
		{
			continue;
		}
		FVector GradC      = Diff / CurrentLength;
		float   C          = CurrentLength;
		float   XPBDAlpha  = ParentParticle.BendingCompliance * InvDeltaTimeSquared;
		float   Lambda     = -C / (W + XPBDAlpha);
		FVector Correction = GradC * Lambda * ChildParticle.InvMass;
		ChildParticle.Position += Correction * ChildParticle.InvMass;

		// recompute parent rotation
		FVector BoneSimDirection = ChildParticle.Position - ParentParticle.Position;
		if (!BoneSimDirection.Normalize())
		{
			continue;
		}
		if (!BonePoseSimDirection.Normalize())
		{
			continue;
		}
		FQuat AddRotation = FQuat::FindBetweenNormals(BonePoseSimDirection, BoneSimDirection);
		ParentParticle.SimRotation = AddRotation * ParentParticle.SimRotation;
		ChildParticle.SimRotation = ParentParticle.SimRotation * ChildPoseLocalTransform.GetRotation();

		// FQuat AddRotation = FQuat::FindBetweenNormals(BonePoseDirection, BoneSimDirection);
		// FQuat::FReal DiffAngle = AddRotation.GetAngle();
		//
		// float C = DiffAngle;
		// float GradC = DiffAngle / FMath::Abs(DiffAngle);
		//
		// float XPBDAlpha = ParentParticle.BendingCompliance * InvDeltaTimeSquared;
		// float Lambda = -C / (W + XPBDAlpha);
		// FQuat Correction = FQuat(AddRotation.GetRotationAxis(), -GradC * Lambda * ChildParticle.InvMass);
		//  
		// ParentParticle.SimRotation = Correction * ParentParticle.SimRotation;
		// ChildParticle.SimRotation = ChildPoseLocalTransform.GetRotation() * ParentParticle.SimRotation;
		// ChildParticle.Position = ParentParticle.Position + ParentParticle.SimRotation.RotateVector(ChildPoseLocalTransform.GetTranslation());
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
	for (int32 ParticleIndex = 0; ParticleIndex < Particles.Num(); ++ParticleIndex)
	{
		auto& ChildParticle = RuntimeParticles[ParticleIndex];
		auto& ChildConfigParticle = Particles[ParticleIndex];
		if (ChildConfigParticle.ParticleType == ELoogPhysicsParticleType::Root)
		{
			continue;
		}
		auto& ParentParticle = RuntimeParticles[ChildParticle.ParentParticleIndex];
		FVector BoneSimDirection = ChildParticle.Position - ParentParticle.Position;
		if (!BoneSimDirection.Normalize())
		{
			continue;
		}
		const FTransform& BoneLocalPoseTransform = Output.Pose.GetLocalSpaceTransform(ChildParticle.BoneIndex);

		FTransform        BonePoseTransform = Output.Pose.GetComponentSpaceTransform(ParentParticle.BoneIndex);
		FVector ForwardAxis = BoneLocalPoseTransform.GetTranslation();
		FVector           BonePoseDirection = BonePoseTransform.GetRotation().RotateVector(ForwardAxis);
		if (!BonePoseDirection.Normalize())
		{
			continue;
		}
		FQuat             AddRotation = FQuat::FindBetweenNormals(BonePoseDirection, BoneSimDirection);
		BonePoseTransform.SetRotation(AddRotation * BonePoseTransform.GetRotation());

		const auto& ConfigParticle = Particles[ChildParticle.ParentParticleIndex];
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
		case ELoogPhysicsParticleType::VirtualBone:
			{
				BonePoseTransform.SetTranslation(ParentParticle.Position);
				if (BonePoseTransform.ContainsNaN())
				{
					continue;
				}
				OutBoneTransforms.Add(FBoneTransform(ParentParticle.BoneIndex, BonePoseTransform));
			}
			break;
		case ELoogPhysicsParticleType::Max:
		default:
			break;
		}
	}
	OutBoneTransforms.Sort(FCompareBoneTransformIndex());
}
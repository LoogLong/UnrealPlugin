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
	for (int32 ParticleIndex = 0; ParticleIndex < Particles.Num(); ++ParticleIndex)
	{
		auto& Particle = Particles[ParticleIndex];
		FCompactPoseBoneIndex BoneIndex = Particle.BindBone.GetCompactPoseIndex(RequiredBones);
		const auto& BoneTransform = Output.Pose.GetComponentSpaceTransform(BoneIndex);
		switch (Particle.ParticleType)
		{
		case ELoogPhysicsParticleType::Root:
			Particle.InvMass = 0.f;
			Particle.PrevPosition = BoneTransform.GetTranslation();
			Particle.Position = Particle.PrevPosition;
			Particle.Velocity = FVector::Zero();
			break;
		case ELoogPhysicsParticleType::Bone:
			Particle.InvMass = 1.0f / Particle.Mass;
			Particle.PrevPosition = BoneTransform.GetTranslation();
			Particle.Position = Particle.PrevPosition;
			Particle.Velocity = FVector::Zero();
			break;
		case ELoogPhysicsParticleType::VirtualBone:
			Particle.InvMass = 1.0f / Particle.Mass;
			Particle.PrevPosition = BoneTransform.TransformPosition(Particle.VirtualBoneLocalPosition);
			Particle.Position = Particle.PrevPosition;
			Particle.Velocity = FVector::Zero();
			break;
		case ELoogPhysicsParticleType::Max:
		default:
			break;
		}
	}
	
	for (FLoogPhysicsClothSection& ClothSection : ClothSections)
	{
		for (const auto& Chain : ClothSection.Chains)
		{
			for (int32 i = 1; i < Chain.ParticleIndices.Num(); ++i)
			{
				const int32& ParticleIndex = Chain.ParticleIndices[i];
				Particles[ParticleIndex].ParentParticleIndex = Chain.ParticleIndices[i - 1];
			}
		}

		for (FLoogPhysicsConstraint& Con : ClothSection.LocalVerticalConstraints)
		{
			int32 AIndex = Con.ParticleAIndex;
			auto& ARuntimeParticle = Particles[AIndex];

			int32 BIndex = Con.ParticleBIndex;
			auto& BRuntimeParticle = Particles[BIndex];
			Con.RestLength = (ARuntimeParticle.Position - BRuntimeParticle.Position).Length();
		}
		for (FLoogPhysicsConstraint& Con : ClothSection.GlobalStructureConstraints)
		{
			int32 AIndex = Con.ParticleAIndex;
			auto& ARuntimeParticle = Particles[AIndex];

			int32 BIndex = Con.ParticleBIndex;
			auto& BRuntimeParticle = Particles[BIndex];
			Con.RestLength = (ARuntimeParticle.Position - BRuntimeParticle.Position).Length();
		}

		for (FLoogPhysicsConstraint& Con : ClothSection.BendingStructureConstraints)
		{
			int32 AIndex = Con.ParticleAIndex;
			auto& ARuntimeParticle = Particles[AIndex];

			int32 BIndex = Con.ParticleBIndex;
			auto& BRuntimeParticle = Particles[BIndex];
			Con.RestLength = (ARuntimeParticle.Position - BRuntimeParticle.Position).Length();
		}
		for (FLoogPhysicsConstraint& Con : ClothSection.LocalHorizontalConstraints)
		{
			int32 AIndex = Con.ParticleAIndex;
			auto& ARuntimeParticle = Particles[AIndex];

			int32 BIndex = Con.ParticleBIndex;
			auto& BRuntimeParticle = Particles[BIndex];
			Con.RestLength = (ARuntimeParticle.Position - BRuntimeParticle.Position).Length();
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
	const auto& RequiredBones = Output.AnimInstanceProxy->GetRequiredBones();
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
		auto& RuntimeParticle = Particles[ParticleIndex];
		switch (ConfigParticle.ParticleType)
		{
		case ELoogPhysicsParticleType::Root:
			{
				const auto& BoneTransCS = Output.Pose.GetComponentSpaceTransform(RuntimeParticle.BindBone.GetCompactPoseIndex(RequiredBones));
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
	auto  UpdateConstraint    = [&](const FLoogPhysicsConstraint& DistanceConstraint)
	{
		auto& Particle0 = Particles[DistanceConstraint.ParticleAIndex];
		auto& Particle1 = Particles[DistanceConstraint.ParticleBIndex];
		float W         = Particle0.InvMass + Particle1.InvMass;
		if (W == 0.f)
		{
			return;
		}
		FVector Delta         = Particle0.Position - Particle1.Position;
		float   CurrentLength = Delta.Length();
		if (CurrentLength <= 0.f)
		{
			return;
		}
		FVector GradC      = Delta / CurrentLength;
		float   C          = CurrentLength - DistanceConstraint.RestLength;
		float   Compliance = C > 0 ? DistanceConstraint.ShrinkCompliance : DistanceConstraint.StretchCompliance;
		float   XPBDAlpha  = Compliance * InvDeltaTimeSquared;
		float   Lambda     = -C / (W + XPBDAlpha);

		FVector Correction = GradC * Lambda;

		Particle0.Position += Correction * Particle0.InvMass;
		Particle1.Position -= Correction * Particle1.InvMass;
	};

	// distance constraint
	for (auto& Section : ClothSections)
	{
		for (auto& Constraint : Section.LocalVerticalConstraints)
		{
			UpdateConstraint(Constraint);
		}
		for (auto& Constraint : Section.GlobalStructureConstraints)
		{
			UpdateConstraint(Constraint);
		}
		for (auto& Constraint : Section.LocalHorizontalConstraints)
		{
			UpdateConstraint(Constraint);
		}
		for (auto& Constraint : Section.BendingStructureConstraints)
		{
			UpdateConstraint(Constraint);
		}
	}

	// bending constraint
	if (bEnableChainBendingLocal)
	{
		for (FLoogPhysicsClothSection& ClothSection : ClothSections)
		{
			for (auto& Chain : ClothSection.Chains)
			{
				for (int32 i = 2; i < Chain.ParticleIndices.Num(); ++i)
				{
					int32 GrandparentIndex = Chain.ParticleIndices[i - 2];
					int32 ParentIndex = Chain.ParticleIndices[i - 1];
					int32 ChildIndex = Chain.ParticleIndices[i];
					auto& GrandparentParticle = Particles[GrandparentIndex];
					auto& ParentParticle = Particles[ParentIndex];
					auto& ChildParticle = Particles[ChildIndex];

					float W = GrandparentParticle.InvMass + ChildParticle.InvMass;
					if (W == 0.f)
					{
						continue;
					}
					float ThetaSim;
					float ThetaPose;
					FVector RotationAxis;
					{
						FVector Dir0 = ChildParticle.Position - ParentParticle.Position;
						FVector Dir1 = GrandparentParticle.Position - ParentParticle.Position;
						if (!Dir0.Normalize() || !Dir1.Normalize())
						{
							continue;
						}
						ThetaSim = FMath::Acos(FVector::DotProduct(Dir0, Dir1));
						RotationAxis = FVector::CrossProduct(Dir0, Dir1);
					}
					{
						FVector ChildPosePosition;
						if (ChildParticle.ParticleType == ELoogPhysicsParticleType::VirtualBone)
						{
							ChildPosePosition = Output.Pose.GetComponentSpaceTransform(ChildParticle.BindBone.GetCompactPoseIndex(RequiredBones)).TransformPosition(ChildParticle.VirtualBoneLocalPosition);
						}
						else
						{
							ChildPosePosition = Output.Pose.GetComponentSpaceTransform(ChildParticle.BindBone.GetCompactPoseIndex(RequiredBones)).GetTranslation();
						}
						FVector ParentPosePosition = Output.Pose.GetComponentSpaceTransform(ParentParticle.BindBone.GetCompactPoseIndex(RequiredBones)).GetTranslation();
						FVector GrandparentPosePosition = Output.Pose.GetComponentSpaceTransform(GrandparentParticle.BindBone.GetCompactPoseIndex(RequiredBones)).GetTranslation();
						FVector Dir0 = ChildPosePosition - ParentPosePosition;
						FVector Dir1 = GrandparentPosePosition - ParentPosePosition;
						if (!Dir0.Normalize() || !Dir1.Normalize())
						{
							continue;
						}
						ThetaPose = FMath::Acos(FVector::DotProduct(Dir0, Dir1));
					}

					FVector GradC     = RotationAxis;
					float   C         = ThetaSim - ThetaPose;
					float   XPBDAlpha = ParentParticle.LocalBendingCompliance * InvDeltaTimeSquared;
					float   Lambda    = -C / (W + XPBDAlpha);

					FQuat ChildDiffRot = FQuat(GradC, -Lambda * ChildParticle.InvMass);
					ChildParticle.Position = ChildDiffRot.RotateVector(ChildParticle.Position - ParentParticle.Position) + ParentParticle.Position;

					FQuat GrandparentDiffRot = FQuat(GradC, Lambda * GrandparentParticle.InvMass);
					GrandparentParticle.Position = GrandparentDiffRot.RotateVector(GrandparentParticle.Position - ParentParticle.Position) + ParentParticle.Position;
				}
			}

		}
	}
	if (bEnableChainBendingGlobal)
	{
		for (int32 ParticleIndex = 0; ParticleIndex < Particles.Num(); ++ParticleIndex)
		{
			auto& ChildParticle = Particles[ParticleIndex];
			if (ChildParticle.ParticleType == ELoogPhysicsParticleType::Root)
			{
				continue;
			}
			auto& ParentParticle = Particles[ChildParticle.ParentParticleIndex];
			float W = ChildParticle.InvMass + ParentParticle.InvMass;
			if (W == 0.f)
			{
				continue;
			}

			FTransform ChildPoseLocalTransform;
			if (ChildParticle.ParticleType == ELoogPhysicsParticleType::VirtualBone)
			{
				ChildPoseLocalTransform = FTransform(ChildParticle.VirtualBoneLocalPosition);
			}
			else
			{
				ChildPoseLocalTransform = Output.Pose.GetLocalSpaceTransform(ChildParticle.BindBone.GetCompactPoseIndex(RequiredBones));
			}

			FVector BonePoseSimDirection = ParentParticle.SimRotation.RotateVector(ChildPoseLocalTransform.GetTranslation());
			FVector BonePosePosition = BonePoseSimDirection + ParentParticle.Position;
			FVector Diff = ChildParticle.Position - BonePosePosition;
			float CurrentLength = Diff.Length();
			if (CurrentLength < UE_KINDA_SMALL_NUMBER)
			{
				continue;
			}
			FVector GradC = Diff / CurrentLength;
			float   C = CurrentLength;
			float   XPBDAlpha = ParentParticle.GlobalBendingCompliance * InvDeltaTimeSquared;
			float   Lambda = -C / (W + XPBDAlpha);
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
		auto& RuntimeParticle = Particles[ParticleIndex];

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
		int32 ParticleCount = Particles.Num();
		for (int ParticleIndex = 0; ParticleIndex < ParticleCount; ++ParticleIndex)
		{
			auto& Particle = Particles[ParticleIndex];

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
			FVector::FReal RestMinDistance = Particle.Thickness + SphereRadius;
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
				auto& Particle0 = Particles[Collider.Particle0Index];
				auto& Particle1 = Particles[Collider.Particle1Index];

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
				float Radius2 = FMath::Lerp(Particle0.Thickness, Particle1.Thickness, T);

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
		auto& ChildParticle = Particles[ParticleIndex];
		if (ChildParticle.ParticleType == ELoogPhysicsParticleType::Root)
		{
			continue;
		}
		auto& ParentParticle = Particles[ChildParticle.ParentParticleIndex];
		FVector BoneSimDirection = ChildParticle.Position - ParentParticle.Position;
		if (!BoneSimDirection.Normalize())
		{
			continue;
		}
		const FTransform& BoneLocalPoseTransform = Output.Pose.GetLocalSpaceTransform(ChildParticle.BindBone.GetCompactPoseIndex(RequiredBones));

		FTransform        BonePoseTransform = Output.Pose.GetComponentSpaceTransform(ParentParticle.BindBone.GetCompactPoseIndex(RequiredBones));
		FVector ForwardAxis = BoneLocalPoseTransform.GetTranslation();
		FVector           BonePoseDirection = BonePoseTransform.GetRotation().RotateVector(ForwardAxis);
		if (!BonePoseDirection.Normalize())
		{
			continue;
		}
		FQuat             AddRotation = FQuat::FindBetweenNormals(BonePoseDirection, BoneSimDirection);
		BonePoseTransform.SetRotation(AddRotation * BonePoseTransform.GetRotation());

		switch (ChildParticle.ParticleType)
		{
		case ELoogPhysicsParticleType::Root:
			{
				BonePoseTransform.SetRotation(AddRotation * BonePoseTransform.GetRotation());
				BonePoseTransform.SetTranslation(ParentParticle.Position);
				if (BonePoseTransform.ContainsNaN())
				{
					continue;
				}
				OutBoneTransforms.Add(FBoneTransform(ParentParticle.BindBone.GetCompactPoseIndex(RequiredBones), BonePoseTransform));
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
				OutBoneTransforms.Add(FBoneTransform(ParentParticle.BindBone.GetCompactPoseIndex(RequiredBones), BonePoseTransform));
			}
			break;
		case ELoogPhysicsParticleType::Max:
		default:
			break;
		}
	}
	OutBoneTransforms.Sort(FCompareBoneTransformIndex());
}
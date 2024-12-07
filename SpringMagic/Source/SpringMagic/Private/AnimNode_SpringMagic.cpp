// Copyright LoogLong. All Rights Reserved.
#include "AnimNode_SpringMagic.h"
#include "AnimationRuntime.h"
#include "Animation/AnimInstanceProxy.h"
#include "Curves/CurveFloat.h"
#include "Runtime/Launch/Resources/Version.h"
#include "SceneInterface.h"
#include "PhysicsEngine/PhysicsAsset.h"


DECLARE_CYCLE_STAT(TEXT("SpringMagic_Eval"), STAT_SpringMagic_Eval, STATGROUP_Anim);

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

	float DeltaTime = Context.GetDeltaTime();
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
	float   DeltaTime = Output.AnimInstanceProxy->GetDeltaSeconds();
	if (PendingTeleportType == ETeleportType::TeleportPhysics)
	{
		
	}else if (PendingTeleportType == ETeleportType::ResetPhysics)
	{
		bNeedInitializeSimulation = true;
		PendingTeleportType = ETeleportType::None;
	}
	if (bNeedInitializeSimulation)
	{
		InitializeSimulation(Output.AnimInstanceProxy->GetRequiredBones());
		PreComponentTransform = ComponentTransform;
		bNeedInitializeSimulation = false;
	}

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

				// TODO Add Inertia
				// TODO Detected Collision


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

void FAnimNode_SpringMagic::InitializeSimulation(const FBoneContainer& RequiredBones)
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
				Particle.bDummyBone = false;
				Particle.bRootBone = i == Idx;
				Particles.Add(Index);
			}
			{
				int32 Index = SimulationParticles.Num();
				auto& Particle = SimulationParticles.AddZeroed_GetRef();
				Particle.LocalTransform.SetRotation(FQuat::Identity);
				Particle.LocalTransform.SetTranslation(FVector::XAxisVector * BoneChain.DummyBoneLength);

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

				// FTransform GrandparentTransform = RequiredBones.GetRefPoseTransform(SimulationParticles[GrandParent].BoneIndex);
				FTransform& ParentTransform = SimulationParticles[Parent].LocalTransform;
				FTransform& ChildTransform = SimulationParticles[Child].LocalTransform;
				Joint.ChildProxy = ChildTransform * ParentTransform;
			}
			MaxChainLength = FMath::Max(MaxChainLength, Chain.Num());
		}
	}
}

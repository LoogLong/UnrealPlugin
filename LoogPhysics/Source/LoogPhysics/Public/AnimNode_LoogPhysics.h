﻿// Copyright LoogLong. All Rights Reserved.
#pragma once
#include "CoreMinimal.h"
#include "BoneContainer.h"
#include "BonePose.h"
#include "BoneControllers/AnimNode_SkeletalControlBase.h"
#include "PhysicsEngine/TaperedCapsuleElem.h"
#include "AnimNode_LoogPhysics.generated.h"

UENUM()
enum class ELoogPhysicsParticleType
{
	Root,
	Bone,
	VirtualBone,
	Max UMETA(Hidden)
};

UENUM()
enum class ELoogPhysicsSimulationMethod
{
	ForceBase,
	VelocityBase,
	PBD,
	XPBD,
	Max UMETA(Hidden)
};


USTRUCT()
struct FLoogPhysicsParticle
{
	GENERATED_USTRUCT_BODY()

	UPROPERTY(EditAnywhere)
	FBoneReference BindBone;

	UPROPERTY(EditAnywhere)
	ELoogPhysicsParticleType ParticleType = ELoogPhysicsParticleType::Root;

	UPROPERTY(EditAnywhere)
	FVector VirtualBoneLocalPosition = FVector::Zero();

	UPROPERTY(EditAnywhere)
	FVector VirtualBoneLocalNormal = FVector::Zero();
	
	UPROPERTY(EditAnywhere)
	float Thickness = 0.1f;
	
	UPROPERTY(EditAnywhere)
	float Mass = 1.0f;

	UPROPERTY(EditAnywhere)
	float LinearDamping = 0.01f;

	UPROPERTY(EditAnywhere)
	float AngularDamping = 0.01f;

	UPROPERTY(EditAnywhere)
	float GlobalStructureCompliance = 0.01f;

	UPROPERTY(EditAnywhere)
	float LocalStructureCompliance = 0.01f;

	UPROPERTY(EditAnywhere)
	float MaxVelocity = 1000.f;

	UPROPERTY(EditAnywhere)
	float MaxAngle = 180.f;
};

USTRUCT()
struct FLoogPhysicsConstraint
{
	GENERATED_USTRUCT_BODY()

	UPROPERTY(EditAnywhere)
	int32 ParticleAIndex;

	UPROPERTY(EditAnywhere)
	int32 ParticleBIndex;

	UPROPERTY(EditAnywhere)
	ELoogPhysicsSimulationMethod SimulationMethod = ELoogPhysicsSimulationMethod::XPBD;

	UPROPERTY(EditAnywhere)
	float ShrinkCompliance = 0.01f;

	UPROPERTY(EditAnywhere)
	float StretchCompliance = 0.01f;
};

USTRUCT()
struct FLoogPhysicsChain
{
	GENERATED_USTRUCT_BODY()

	UPROPERTY(EditAnywhere)
	TArray<int32> ParticleIndices;

	UPROPERTY(EditAnywhere)
	int32 RootParticleIndex = 0;

	UPROPERTY(EditAnywhere)
	float MaxRootVelocity = 1000.f;

	UPROPERTY(EditAnywhere)
	bool bHighLight = false;
};

USTRUCT()
struct FLoogPhysicsParticleCollider
{
	GENERATED_USTRUCT_BODY()

	UPROPERTY(EditAnywhere)
	int32 Particle0Index;

	UPROPERTY(EditAnywhere)
	int32 Particle1Index;

	// method to handle collision
};

USTRUCT()
struct FLoogPhysicsClothSection
{
	GENERATED_USTRUCT_BODY()

	UPROPERTY(EditAnywhere)
	FString SectionName;

	UPROPERTY(EditAnywhere)
	bool bChainLoop = true;

	UPROPERTY(EditAnywhere)
	TArray<FLoogPhysicsChain> Chains;

	UPROPERTY(EditAnywhere)
	TArray<FLoogPhysicsConstraint> Constraints;

	UPROPERTY(EditAnywhere)
	TArray<FLoogPhysicsParticleCollider> Colliders;
};

struct FLoogPhysicsRuntimeParticle
{
	FCompactPoseBoneIndex BoneIndex;
	float InvMass;
	int32 ParentParticleIndex;

	FVector PrevPosition;
	FVector Position;
	FVector Velocity;

};

struct FLoogPhysicsRuntimeConstraint
{
	int32 Particle0Index;
	int32 Particle1Index;

	float RestLength;
	float ShrinkCompliance;
	float StretchCompliance;
};

struct FLoogPhysicsCapsule
{
	FBoneReference BindBone;
	FVector CenterA;
	FVector CenterB;

	FTransform TransformSimSpace;
	FVector CenterASimSpace;
	FVector CenterBSimSpace;
};
USTRUCT(BlueprintType)
struct LOOGPHYSICS_API FAnimNode_LoogPhysics : public FAnimNode_SkeletalControlBase
{
	GENERATED_USTRUCT_BODY()

public:
	UPROPERTY(EditAnywhere, Category = "ClothSettings")
	TArray<FLoogPhysicsParticle> Particles;

	UPROPERTY(EditAnywhere, Category = "ClothSettings")
	TArray<FLoogPhysicsClothSection> ClothSections;

	UPROPERTY(EditAnywhere, Category = "SimulationSetting")
	float FrameRate = 30.f;

	UPROPERTY(EditAnywhere, Category = "SimulationSetting")
	int32 MaxSimulationPerFrame = 2;

	UPROPERTY(EditAnywhere, Category = "SimulationSetting")
	FVector GravityWorldSpace = FVector(0, 0, -980);

	UPROPERTY(EditAnywhere, Category = "Wind")
	bool bUseWind = true;

	UPROPERTY(EditAnywhere, Category = "Wind")
	float WindFrequency = 0.15f;

	UPROPERTY(EditAnywhere, Category = "Wind")
	float WindAmplitude = 1.0f;

	UPROPERTY(EditAnywhere, Category = "Wind")
	float WindAmplitudeBias = 0.5f;

	UPROPERTY(EditAnywhere, Category = "Wind", meta = (PinHiddenByDefault))
	FVector WindDirectionWorldSpace = FVector::Zero();

	UPROPERTY(EditAnywhere, Category = "ComponentMovenent")
	bool bUseComponentMovement = true;

	UPROPERTY(EditAnywhere, Category = "ComponentMovenent")
	float ComponentMovementRotation = 0.7f;

	UPROPERTY(EditAnywhere, Category = "ComponentMovenent")
	float ComponentMovementTranslation = 0.2;

	UPROPERTY(EditAnywhere, Category = "Collision")
	bool bEnableCollision = false;

	UPROPERTY(EditAnywhere, Category = "Collision")
	TObjectPtr<UPhysicsAsset> PhysicsAsset;

public:
	friend class FLoogPhysicsEditMode;
	FAnimNode_LoogPhysics();

	// FAnimNode_Base interface
	virtual void UpdateInternal(const FAnimationUpdateContext& Context) override;
	virtual void GatherDebugData(FNodeDebugData& DebugData) override;
	virtual void Initialize_AnyThread(const FAnimationInitializeContext& Context) override;
	virtual void CacheBones_AnyThread(const FAnimationCacheBonesContext& Context) override;
	virtual bool NeedsDynamicReset() const override { return true; }
	virtual void ResetDynamics(ETeleportType InTeleportType) override;
	// End of FAnimNode_Base interface

	// FAnimNode_SkeletalControlBase interface
	virtual void EvaluateSkeletalControl_AnyThread(FComponentSpacePoseContext& Output,
	                                               TArray<FBoneTransform>& OutBoneTransforms) override;
	virtual bool IsValidToEvaluate(const USkeleton* Skeleton, const FBoneContainer& RequiredBones) override;
	// virtual bool HasPreUpdate() const override;
	// virtual void PreUpdate(const UAnimInstance* InAnimInstance) override;
	// End of FAnimNode_SkeletalControlBase interface


protected:
	// FAnimNode_SkeletalControlBase interface
	virtual void InitializeBoneReferences(const FBoneContainer& RequiredBones) override;
	// End of FAnimNode_SkeletalControlBase interface

	

private:
	void InitializeSimulation(FComponentSpacePoseContext& Output, const FBoneContainer& RequiredBones);
	void InitializeCollision(const FBoneContainer& RequiredBones);
	void SimulatePhysics(FComponentSpacePoseContext& Output, const float& DeltaTimeSeconds);
	void SimulateOnce(FComponentSpacePoseContext& Output, const float& DeltaTimeSeconds);

	void UpdateCollision(FComponentSpacePoseContext& Output, const FBoneContainer& RequiredBones);
	void CollisionDetection();

	/**
	 * Applies the simulation results to the bone transforms.
	 *
	 * @param Output The pose context.
	 * @param OutBoneTransforms An array to store the resulting bone transforms.
	 */
	void ApplySimulateResult(FComponentSpacePoseContext& Output, TArray<FBoneTransform>& OutBoneTransforms);

	TArray<FLoogPhysicsRuntimeParticle> RuntimeParticles;
	TArray<FLoogPhysicsRuntimeConstraint> RuntimeConstraints;
	TArray<FKTaperedCapsuleElem> SetupColliders;
	TArray<FLoogPhysicsCapsule> RuntimeColliders;


	int32 MaxChainLength = 0;
	float WindSinCounter = 0.f;
	FTransform PreComponentTransform;
	ETeleportType PendingTeleportType = ETeleportType::None;
	bool bIsValidToSimulation = false;
	bool bNeedInitializeSimulation = true;
};
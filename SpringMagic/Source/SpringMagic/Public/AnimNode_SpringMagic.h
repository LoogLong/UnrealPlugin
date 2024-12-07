// Copyright LoogLong. All Rights Reserved.
#pragma once
#include "CoreMinimal.h"
#include "BoneContainer.h"
#include "BonePose.h"
#include "BoneControllers/AnimNode_SkeletalControlBase.h"
#include "AnimNode_SpringMagic.generated.h"

USTRUCT()
struct  FSpringMagicChain
{
	GENERATED_USTRUCT_BODY()

	UPROPERTY(EditAnywhere)
	FBoneReference RootBone;

	UPROPERTY(EditAnywhere)
	FBoneReference EndBone;

	UPROPERTY(EditAnywhere)
	float DummyBoneLength = 10.f;

};


struct FSpringMagicParticle
{
	// simulation space(component space)
	FTransform SimulationTransform;
	FTransform LocalTransform;
	FCompactPoseBoneIndex BoneIndex;


	bool bDummyBone;
	bool bRootBone;
};


struct FSpringMagicJoint
{
	int32 GrandParent;
	int32 Parent;
	int32 Child;

	FTransform ChildProxy;

	FVector NewChildPosition;
	FVector CurrentChildPosition;
};

/**
 * The implementation of FAnimNode_SpringMagic is highly inspired By Anim Bai's SpringMagic for Maya
 * Anim Bai:https://animbai.com/
 */
USTRUCT(BlueprintType)
struct SPRINGMAGIC_API FAnimNode_SpringMagic : public FAnimNode_SkeletalControlBase
{
	GENERATED_USTRUCT_BODY()

public:
	UPROPERTY(EditAnywhere, Category = "Settings")
	TArray<FSpringMagicChain> BoneChains;

	UPROPERTY(EditAnywhere)
	float SwingRatio = 0.7f;

	UPROPERTY(EditAnywhere)
	float TwistRatio = 0.7f;

	UPROPERTY(EditAnywhere, Category = "Wind")
	bool bUseWind = true;

	UPROPERTY(EditAnywhere, Category = "Wind")
	float WindFrequency = 0.15f;

	UPROPERTY(EditAnywhere, Category = "Wind")
	float WindAmplitude = 1.0f;

	UPROPERTY(EditAnywhere, Category = "Wind")
	float WindAmplitudeBias = 0.5f;

	UPROPERTY(EditAnywhere, Category = "Wind", meta = (PinHiddenByDefault))
	FVector WindDirection = FVector::Zero();

	UPROPERTY(EditAnywhere, Category = "ComponentMovenent")
	bool bUseComponentMovement = true;

	UPROPERTY(EditAnywhere, Category = "ComponentMovenent")
	float ComponentMovementRotation = 0.7f;

	UPROPERTY(EditAnywhere, Category = "ComponentMovenent")
	float ComponentMovementTranslation = 0.2;

public:
	friend class FSpringMagicEditMode;
	FAnimNode_SpringMagic();

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

	/**
	 * Applies the simulation results to the bone transforms.
	 *
	 * @param Output The pose context.
	 * @param OutBoneTransforms An array to store the resulting bone transforms.
	 */
	void ApplySimulateResult(FComponentSpacePoseContext& Output, TArray<FBoneTransform>& OutBoneTransforms);


private:
	void InitializeSimulation(const FBoneContainer& RequiredBones);

	TArray<FSpringMagicParticle> SimulationParticles;
	TArray<TArray<FSpringMagicJoint>> SimulationJoints;
	int32 MaxChainLength = 0;
	float SinCounter = 0.f;
	FTransform PreComponentTransform;

	ETeleportType PendingTeleportType = ETeleportType::None;
	bool bNeedInitializeSimulation = true;
};
